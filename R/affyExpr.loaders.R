loadUnitNames.affyExpr <- function(db, unames) {
    dbBeginTransaction(db)
    ## To use an auto-incrementing field via RSQLite, you need
    ## to be careful to pass integer NA's
    df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
    values <- "(:id, :name)"
    sql <- "insert into featureSet (fsetid, man_fsetid) values"
    dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
    dbCommit(db)
}

loadUnits.affyExpr <- function(db, batch, isQc=FALSE) {
    pmfeature <- "pmfeature_tmp"
    mmfeature <- "mmfeature_tmp"
    pmmm <- "pm_mm"
    if (isQc) {
        pmfeature <- "qcpmfeature"
        mmfeature <- "qcmmfeature"
        pmmm <- "qcpm_qcmm"
    }
    batchMat <- do.call(rbind, lapply(batch, readCdfUnitToMat.affyExpr))

    loadUnitNames.affyExpr(db, names(batch))

    ## Find internal featureSet IDs for these features
    batchLens <- sapply(batch, function(x)
                        sum(sapply(x$groups, function(y)
                                   length(y$indices))))
    sql <- paste("select fsetid from featureSet where man_fsetid in (",
                 paste('"', names(batch), '"', sep="", collapse=","),
                 ") order by fsetid")
    batchIds <- dbGetQuery(db, sql)[[1]]
    batchIds <- rep(batchIds, batchLens)
    batchMat <- cbind(batchMat, fsetid=batchIds)

    ## Insert pm and mm into respective tables
    isPm <- as.logical(batchMat[, "ispm"])
    values <- "(:indices, :fsetid, :pbase, :tbase,  :atom, :x, :y)"
    sql <- paste("insert into", pmfeature, "values", values)
    dbBeginTransaction(db)
    rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[isPm, ]))
    dbClearResult(rset)

    sql <- paste("insert into", mmfeature, "values", values)
    rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[!isPm, ]))
    dbClearResult(rset)
    dbCommit(db)

    ## Insert pm <--> mm link
    values <- "(:pmi, :mmi)"
    sql <- paste("insert into", pmmm, "values", values)
    dbBeginTransaction(db)
    rset <- dbSendPreparedQuery(db, sql,
                                data.frame(pmi=batchMat[isPm, "indices"],
                                           mmi=batchMat[!isPm, "indices"]))
    dbClearResult(rset)
    dbCommit(db)
}


loadUnitsByBatch.affyExpr <- function(db, cdfFile, batch_size=10000,
                             max_units=NULL, verbose=FALSE) {
    unames <- readCdfUnitNames(cdfFile)
    if (is.null(max_units))
      max_units <- length(unames)
    offset <- 1
    whQc <- grep("^AFFX", unames)
    if (length(whQc)) {                 # load all QC at once
        qcunits <- readCdf(cdfFile, units=whQc, readGroupDirection=TRUE,
                           readIndices=TRUE, readIsPm=TRUE)
        loadUnits.affyExpr(db, qcunits, isQc=TRUE)
        offset <- max(whQc) + 1
    }
    extra <- (length(unames)-length(whQc)) %% batch_size
    num_batches <- (length(unames)-length(whQc)) %/% batch_size
    if (extra != 0)
      num_batches <- num_batches + 1
    done <- FALSE
    while (!done) {
        end <- min(offset + batch_size, max_units)
        if (end == max_units)
          done <- TRUE
        wanted <- seq.int(offset, end)
        offset <- offset + batch_size + 1
        vvunits <- readCdf(cdfFile, units=wanted, readGroupDirection=TRUE,
                           readIndices=TRUE, readIsPm=TRUE)
        loadUnits.affyExpr(db, vvunits)
    }
}


loadAffySeqCsv.affyExpr <- function(db, csvFile, cdfFile, batch_size=5000) {
    cdfHeader <- readCdfHeader(cdfFile)
    ncol <- cdfHeader$ncol
    xy2i <- function(x, y) x + 1 + y * ncol

    complementBase <- function(x, special=FALSE){
      bases <- c("A", "C", "G", "T")
      if (!special)
        comp <- c("T", "G", "C", "A")
      else
        comp <- c("G", "T", "A", "C")
      comp[match(x, bases)]
    }

    con <- file(csvFile, open="r")
    on.exit(close(con))
    header <- c("fset.name", "x", "y", "interrogation_position", "seq", "tstrand")
    done <- FALSE
    tmp <- read.table(con, sep="\t", stringsAsFactors=FALSE,
                      nrows=1, na.strings="---", header=FALSE)
    while (!done) {
        pmdf <- read.table(con, sep="\t", stringsAsFactors=FALSE,
                           nrows=batch_size, na.strings="---",
                           header=FALSE)
        if (nrow(pmdf) < batch_size) {
            done <- TRUE
            if (nrow(pmdf) == 0)
              break
        }
        names(pmdf) <- header
        pmdf[["fid"]] <- xy2i(pmdf[["x"]], pmdf[["y"]])
        pmdf[["tstrand"]] <- ifelse(tolower(pmdf[["tstrand"]]) == "sense", SENSE, ANTISENSE)
        N <- nrow(pmdf)

        ## FIXME: Are these files PM-only???
 
        values <- "(:fid, :tstrand, :interrogation_position, :seq)"
        sql <- paste("insert into sequence values", values)
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, sql, bind.data=pmdf)
##        dbGetPreparedQuery(db, sql, bind.data=mmdf)
        dbCommit(db)
    }

}


buildPdInfoDb.affyExpr <- function(cdfFile, csvFile, csvSeqFile, dbFile, matFile,
                          batch_size=10000, verbose=FALSE) {

    ST <- system.time
    printTime <- function(msg, t) {
        if (verbose) {
            m <- paste(msg, "took %.2f sec\n")
            cat(sprintf(m, t))
        }
    }

    db <- initDb.affyExpr(dbFile)

    t <- ST(loadUnitsByBatch.affyExpr(db, cdfFile, batch_size=batch_size))
    printTime("loadUnitsByBatch", t[3])
    t <- ST(loadAffyCsv.affyExpr(db, csvFile, batch_size=batch_size))
    printTime("loadAffyCsv", t[3])
    t <- ST(loadAffySeqCsv.affyExpr(db, csvSeqFile, cdfFile, batch_size=batch_size))
    printTime("loadAffySeqCsv", t[3])
    t <- ST({
        sortFeatureTables.affyExpr(db)
        createIndicesDb.affyExpr(db)
        createTableInfoTable.affyExpr(db)
        createFeatureTableInfo.affyExpr(db)
    })
    printTime("DB sort, index creation", t[3])

##     t <- ST({
##         seqMat <- createSeqMat(db)
##         save(seqMat, file=matFile, compress=TRUE)
##     })
##     printTime("sequence matrix", t[3])
    closeDb(db)
}

loadAffyCsv.affyExpr <- function(db, csvFile, batch_size=5000) {
    con <- file(csvFile, open="r")
    on.exit(close(con))

    df <- read.table(con, sep=",", stringsAsFactors=FALSE, nrows=10,
                     na.strings="---", header=TRUE)
    wantedCols <- match(c("Probe.Set.ID", "Alignments", "Gene.Symbol", "Chromosomal.Location", "Ensembl"), names(df))
    df <- df[, wantedCols]
    header <- gsub(".", "_", names(df), fixed=TRUE)
    names(df) <- header

    db_cols <- c("alignment", "gene_symbol",  "chrom", "ensembl")
    val_holders <- c(":Alignments", ":Gene_Symbol", ":Chromosomal_Location", ":Ensembl")
    exprs <- paste(db_cols, " = ", val_holders, sep="", collapse=", ")
    sql <- paste("update featureSet set ", exprs,
                 "where man_fsetid = :Probe_Set_ID")

    dbBeginTransaction(db)
    dbGetPreparedQuery(db, sql, bind.data=df)
    dbCommit(db)

    ## Now do the rest in batches
    done <- FALSE
    while (!done) {
        df <- read.table(con, sep=",", stringsAsFactors=FALSE,
                         nrows=batch_size, na.strings="---",
                         header=FALSE)[, wantedCols]
        if (nrow(df) < batch_size) {
            done <- TRUE
            if (nrow(df) == 0)
              break
        }
        names(df) <- header
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, sql, bind.data=df)
        dbCommit(db)
    }
}

