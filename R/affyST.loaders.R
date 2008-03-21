.idxToIdx <- function(pgfFrom, pgfTo, ids) {
  starts <- pgfFrom[ids]
  ends <- pgfFrom[ids+1] - 1
  ends[is.na(ends)] <- length(pgfTo)
  mapply(":", starts, ends, SIMPLIFY=FALSE)
}

.combineIdx <- function(psToAtomIdx, atomToProbeIdx) {
  probesPerAtom <- with(atomToProbeIdx,
                        sapply(split(atomIdx, atomIdx), length))
  cbind(probesetIdx=rep(psToAtomIdx$probesetIdx, probesPerAtom),
        atomToProbeIdx)
}

probesetIdxToAtomIdx <- function(pgf, probesetIdx) {
  atoms <- .idxToIdx(pgf[["probesetStartAtom"]],
                     pgf[["atomId"]], probesetIdx)
  data.frame(probesetIdx=rep(probesetIdx, sapply(atoms, length)),
             atomIdx=unlist(atoms))
}

atomIdxToProbeIdx <- function(pgf, atomIdx) {
  probes <- .idxToIdx(pgf[["atomStartProbe"]],
                      pgf[["probeId"]], atomIdx)
  data.frame(atomIdx=rep(atomIdx, sapply(probes, length)),
             probeIdx=unlist(probes))
}

probesetIdxToTripletIdx <- function(pgf, probesetIdx) {
  df1 <- probesetIdxToAtomIdx(pgf, probesetIdx)
  df2 <- atomIdxToProbeIdx(pgf, df1$atomIdx)
  .combineIdx(df1, df2)
}

loadUnits.affyST <- function(db, pgf, clf) {
    ## load featureSet table
    fset.table <- data.frame(fsetid = pgf[["probesetId"]],
                             man_fsetid = pgf[["probesetName"]],
                             type = pgf[["probesetType"]],
                             start_atom = pgf[["probesetStartAtom"]],
                             stringsAsFactors=FALSE)
    i <- fset.table[["man_fsetid"]] == ""
    fset.table[i, "man_fsetid"] <- fset.table[i, "fsetid"]
    rm(i)
    fset.table <- fset.table[order(fset.table[["man_fsetid"]]),]
    dbBeginTransaction(db)
    values <- "(:fsetid, :man_fsetid, :type, :start_atom)"
    sql <- "insert into featureSet (fsetid, man_fsetid, type, start_atom) values"
    dbGetPreparedQuery(db, paste(sql, values), bind.data=fset.table)
    dbCommit(db)

    ## load pmfeature / mmfeature / qc*feature
    ## (put together the "probe df")
    df <- probesetIdxToTripletIdx(pgf, 1:length(pgf[["probesetId"]]))
    fid <- pgf[["probeId"]][df[["probeIdx"]]]
    i <- match(fid, pgf[["probeId"]])
    ii <- match(fid, clf[["id"]])
    probes.table <- data.frame(fid=fid,
                               fsetid=pgf[["probesetId"]][df[["probesetIdx"]]],
                               pbase=as.integer(NA),
                               tbase=as.integer(NA),
                               atom=pgf[["atomId"]][df[["atomIdx"]]],
                               x=clf[["x"]][ii],
                               y=clf[["y"]][ii],
                               gc_count=pgf[["probeGcCount"]][i],
                               type=pgf[["probeType"]][i],
                               strand=SENSE,
                               pm=as.integer(NA),
                               stringsAsFactors=FALSE)
    rm(i, ii, df)
    i <- grep(":at$", tolower(probes.table[["type"]]))
    probes.table[i, "strand"] <- ANTISENSE
    i <- grep("^mm:", tolower(probes.table[["type"]]))
    pmonly <- length(i) == 0
    if (!pmonly)
      probes.table[i, "pm"] <- as.integer(0)
    i <- grep("^pm:", tolower(probes.table[["type"]]))
    probes.table[i, "pm"] <- as.integer(1)
    probes.table[["type"]] <- NULL
    lastCol <- match("pm", names(probes.table))
    probes.table <- probes.table[order(probes.table[["fid"]]),]

    ## loading pm
    sql <- paste("insert into pmfeature_tmp values",
                 "(:fid, :fsetid, :pbase, :tbase, :atom, :x, :y, :gc_count, :strand)")
    dbBeginTransaction(db)
    dbGetPreparedQuery(db, sql, bind.data=subset(probes.table, pm == 1)[,-lastCol])
    dbCommit(db)

    ## loading mm
    if (!pmonly){
      sql <- paste("insert into mmfeature_tmp values",
                   "(:fid, :fsetid, :pbase, :tbase, :atom, :x, :y, :gc_count, :strand)")
      dbBeginTransaction(db)
      dbSendPreparedQuery(db, sql, subset(probes.table, pm == 0)[,-lastCol])
      dbCommit(db)
    }
    
    #### FIXME!!!
    #### ADD PM-MM LINK FOR ARRAYS THAT CONTAIN MM PROBES
    ## Insert pm <--> mm link
##     values <- "(:pmi, :mmi)"
##     sql <- paste("insert into", pmmm, "values", values)
##     dbBeginTransaction(db)
##     rset <- dbSendPreparedQuery(db, sql,
##                                 data.frame(pmi=batchMat[isPm, "indices"],
##                                            mmi=batchMat[!isPm, "indices"]))
##     dbClearResult(rset)
##     dbCommit(db)
}

readProbeFile <- function(filename){
  header <- as.character(read.table(filename, nrow=1, stringsAsFactors=FALSE, header=FALSE, sep="\t"))
  header <- tolower(gsub(" ", "_", header))
  expected <- c("probe_id", "transcript_cluster_id", "probe_x",
                "probe_y", "assembly", "seqname", "start", "stop",
                "strand", "probe_sequence", "target_strandedness",
                "category")
  stopifnot(identical(sort(header), sort(expected)))
  reorder <- match(expected, header)
  classes <- rep(rep(c("integer", "character"), 2), c(4, 2, 2, 4))
  probes <- read.table(filename, skip=1, header=FALSE, colClasses=classes, sep="\t", na.string="---")
  names(probes) <- header
  probes <- probes[order(probes[["probe_id"]]),]
  return(probes[,reorder])
}


loadAffySeqCsv.affyST <- function(db, probeFile, pgf, batch_size=5000) {
  probeSeq <- readProbeFile(probeFile)
  i <- which(probeSeq[["probe_id"]] %in% pgf[["probeId"]])
  probeSeq[["interrogation_position"]] <- NA
  probeSeq[i, "interrogation_position"] <- pgf[["probeInterrogationPosition"]][match(probeSeq[["probe_id"]][i], pgf[["probeId"]])]
  probeSet <- probeSeq[, c("probe_id", "target_strandedness",
                           "interrogation_position", "transcript_cluster_id",
                           "seqname", "start", "stop", "probe_sequence",
                           "category")]
  ts <- ifelse(probeSet[["target_strandedness"]] == "Sense", SENSE, ANTISENSE)
  ts <- as.integer(ts)
  probeSet[["target_strandedness"]] <- ts
  rm(ts,i, probeSeq)
  probeSet <- probeSet[order(probeSet[["probe_id"]]),]
  
  sql <- paste("INSERT INTO sequence VALUES",
               "(:probe_id, :target_strandedness,",
               ":interrogation_position, :transcript_cluster_id,",
               ":seqname, :start, :stop, :probe_sequence, :category)")
  
  ## anything specific for MM probes?
  dbBeginTransaction(db)
  dbGetPreparedQuery(db, sql, bind.data=probeSet)
  dbCommit(db)

}

readTranscriptFile <- function(filename){
  transcript <- read.csv(filename, comment.char="#", stringsAsFactors=FALSE, na.string="---")
  expected <- c("transcript_cluster_id", "probeset_id", "seqname",
                "strand", "start", "stop", "total_probes",
                "gene_assignment", "mrna_assignment", "swissprot",
                "unigene", "GO_biological_process",
                "GO_cellular_component", "GO_molecular_function",
                "pathway", "protein_domains", "crosshyb_type",
                "category")
  stopifnot(identical(sort(names(transcript)), sort(expected)))
  transcript[["strand"]] <- as.integer(ifelse(transcript[["strand"]] == "+", SENSE, ANTISENSE))
  reorder <- match(expected, names(transcript))
  transcript <- transcript[order(transcript[["probeset_id"]]),]
  return(transcript[, reorder])
}

loadAffyCsv.affyST <- function(db, transFile, batch_size=5000) {
  transcript <- readTranscriptFile(transFile)

  db_cols <- c("transcript_cluster_id", "fsetid", "seqname",
               "strand", "start", "stop", "total_probes",
               "gene_assignment", "mrna_assignment", "swissprot",
               "unigene", "GO_biological_process",
               "GO_cellular_component", "GO_molecular_function",
               "pathway", "protein_domains", "crosshyb_type",
               "category")

  val_holders <- paste(":", db_cols, sep="")
  val_holders <- gsub("fsetid", "probeset_id", val_holders)
  exprs <- paste(db_cols, " = ", val_holders, sep="", collapse=", ")
  sql <- paste("UPDATE featureSet set", exprs, "WHERE fsetid = :probeset_id")
  
  dbBeginTransaction(db)
  dbGetPreparedQuery(db, sql, bind.data=transcript)
  dbCommit(db)
}

buildPdInfoDb.affyST <- function(pgfFile, clfFile, probeFile, transFile, dbFile, matFile,
                          batch_size=10000, verbose=FALSE) {

    ST <- system.time
    printTime <- function(msg, t) {
        if (verbose) {
            m <- paste(msg, "took %.2f sec\n")
            cat(sprintf(m, t))
        }
    }
    clf <- readClf(clfFile)
    pgf <- readPgf(pgfFile)

    db <- initDb.affyST(dbFile)

    t <- ST(loadUnits.affyST(db, pgf, clf))
    printTime("loadUnitsByBatch", t[3])
    t <- ST(loadAffyCsv.affyST(db, transFile, batch_size=batch_size))
    printTime("loadAffyCsv", t[3])
    t <- ST(loadAffySeqCsv.affyST(db, probeFile, pgf, batch_size=batch_size))
    printTime("loadAffySeqCsv", t[3])
    t <- ST({
        sortFeatureTables.affyST(db)
        createIndicesDb.affyST(db)
        createTableInfoTable.affyST(db)
        createFeatureTableInfo.affyST(db)
    })
    printTime("DB sort, index creation", t[3])

##     t <- ST({
##         seqMat <- createSeqMat(db)
##         save(seqMat, file=matFile, compress=TRUE)
##     })
##     printTime("sequence matrix", t[3])
    closeDb(db)
}
