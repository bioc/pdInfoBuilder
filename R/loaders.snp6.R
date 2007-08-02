snp6.loadUnitNames <- function(db, unames) {
    dbBeginTransaction(db)
    ## To use an auto-incrementing field via RSQLite, you need
    ## to be careful to pass integer NA's
    df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
    values <- "(:id, :name)"
    sql <- "insert into featureSet (fsetid, man_fsetid) values"
    dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
    dbCommit(db)
}

snp6.loadUnitNames.cnv <- function(db, unames) {
    dbBeginTransaction(db)
    ## To use an auto-incrementing field via RSQLite, you need
    ## to be careful to pass integer NA's
    df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
    values <- "(:id, :name)"
    sql <- "insert into featureSetCNV (fsetid, man_fsetid) values"
    dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
    dbCommit(db)
}

snp6.loadUnits.snp <- function(db, batch, isQc=FALSE) {
  pmfeature <- "pmfeature_tmp"

  ## Don't check PM/MM matching.
  ## SNP 5/6 is PM-only
  batchMat <- do.call(rbind, lapply(batch, readCdfUnitToMat, verify.pmmm=FALSE))

  loadUnitNames(db, names(batch))

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
  
  ## Insert pm 
  isPm <- as.logical(batchMat[, "ispm"])
  values <- "(:indices, :strand, :allele, :fsetid, :indexpos, :x, :y)"
  sql <- paste("insert into", pmfeature, "values", values)
  dbBeginTransaction(db)
  rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[isPm, ]))
  dbClearResult(rset)
  dbCommit(db)
}

readCdfUnitToMat.cnv <- function(u){
  mat <- t(sapply(u$groups,
                  function(vv){
                    theDir <- groupDirectionToInt(vv$groupdirection)
                    c(indices=vv$indice, strand=theDir, x=vv$x, y=vv$y, ispm=vv$ispm)
                  }))
}

snp6.loadUnits.cnv <- function(db, batch, isQc=FALSE) {
  pmfeature <- "pmfeatureCNV_tmp"

  ## Don't check PM/MM matching.
  ## SNP 5/6 is PM-only
  batchMat <- do.call(rbind, lapply(batch, readCdfUnitToMat.cnv))

  snp6.loadUnitNames.cnv(db, names(batch))

  ## Find internal featureSet IDs for these features
  batchLens <- sapply(batch, function(x)
                      sum(sapply(x$groups, function(y)
                                 length(y$indices))))
  sql <- paste("select fsetid from featureSetCNV where man_fsetid in (",
               paste('"', names(batch), '"', sep="", collapse=","),
               ") order by fsetid")
  batchIds <- dbGetQuery(db, sql)[[1]]
  batchIds <- rep(batchIds, batchLens)
  batchMat <- cbind(batchMat, fsetid=batchIds)
  
  ## Insert pm 
  isPm <- as.logical(batchMat[, "ispm"])
  values <- "(:indices, :strand, :fsetid, :x, :y)"
  sql <- paste("insert into", pmfeature, "values", values)
  dbBeginTransaction(db)
  rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[isPm, ]))
  dbClearResult(rset)
  dbCommit(db)
}

snp6.loadUnitsByBatch <- function(db, cdfFile, batch_size=10000,
                                  max_units=NULL, verbose=FALSE) {
  unames <- readCdfUnitNames(cdfFile)
  if (is.null(max_units))
    max_units <- length(unames)

  snp.probes <- grep("^SNP", unames)
  cnv.probes <- grep("^CN", unames)
  whQc <- (1:length(unames))[-c(snp.probes, cnv.probes)]

  ## SNP probes
  done <- FALSE
  while(!done){
    if (length(snp.probes) >= batch_size){
      wanted <- snp.probes[1:batch_size]
    }else{
      wanted <- snp.probes
      done <- TRUE
    }
    vvunits <- readCdf(cdfFile, units=wanted, readGroupDirection=TRUE,
                       readIndices=TRUE, readIsPm=TRUE)
    snp6.loadUnits.snp(db, vvunits)
    if (!done)
      snp.probes <- snp.probes[-(1:batch_size)]
  }
  
  ## CNV probes
  done <- FALSE
  while(!done){
    if (length(snp.probes) >= batch_size){
      wanted <- cnv.probes[1:batch_size]
    }else{
      wanted <- cnv.probes
      done <- TRUE
    }
    vvunits <- readCdf(cdfFile, units=wanted, readGroupDirection=TRUE,
                       readIndices=TRUE, readIsPm=TRUE)
    snp6.loadUnits.cnv(db, vvunits)
    if (!done)
      cnv.probes <- cnv.probes[-(1:batch_size)]
  }
}

snp6.loadAffySeqCsv <- function(db, csvFile, cdfFile, batch_size=5000) {
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
    header <- c("fset.name", "x", "y", "offset", "seq", "tstrand", "type",
                "tallele")
    done <- FALSE
    pmdf <- read.table(con, sep="\t", stringsAsFactors=FALSE, nrows=1, header=FALSE)
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
        N <- nrow(pmdf)

        ## process MMs
##         mmSql <- paste("select mm_fid, pm_fid from pm_mm where pm_fid in (",
##                        paste(pmdf[["fid"]], collapse=","), ")")
##         pairedIds <- dbGetQuery(db, mmSql)
##         foundIdIdx <- match(pmdf[["fid"]], pairedIds[["pm_fid"]], 0)
##         mmdf <- pmdf[foundIdIdx, ]
##         mmdf[["fid"]] <-  pairedIds[["mm_fid"]]
## 
##         ## Assuming 25mers
##         midbase <- substr(mmdf$seq, 13, 13)
##         types <- apply(table(mmdf$fset.name, mmdf$tallele)>0, 1,
##                        function(v) paste(c("A", "C", "G", "T")[v], collapse=""))
##         types <- rep(types, as.integer(table(mmdf$fset.name)))
##         isSpecial <- (types == "AT" | types == "AG") & mmdf$offset == 0
##         rm(types)
##         midbase[isSpecial] <- complementBase(midbase[isSpecial], T)
##         midbase[!isSpecial] <- complementBase(midbase[!isSpecial])
##         rm(isSpecial)
##         mmdf$seq <- paste(substr(mmdf$seq, 1, 12), midbase,
##                           substr(mmdf$seq, 14, 25), sep="")
##         rm(midbase)
##         ## end MM seq

        values <- "(:fid, :offset, :tstrand, :tallele, :seq)"
        sql <- paste("insert into sequence values", values)
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, sql, bind.data=pmdf)
##        dbGetPreparedQuery(db, sql, bind.data=mmdf)
        dbCommit(db)
    }

}


snp6.buildPdInfoDb <- function(cdfFile, csvFile, csvSeqFile, csvFileCnv,
                               csvSeqFileCnv, dbFile, matFile, matFileCnv,
                               batch_size=10000, verbose=FALSE) {
  ST <- system.time
  printTime <- function(msg, t) {
    if (verbose) {
      m <- paste(msg, "took %.2f sec\n")
      cat(sprintf(m, t))
    }
  }

  db <- snp6.initDb(dbFile)
  t <- ST(snp6.loadUnitsByBatch(db, cdfFile, batch_size=batch_size))
  printTime("loadUnitsByBatch", t[3])
  t <- ST(snp6.loadAffyCsv(db, csvFile, batch_size=batch_size))
  printTime("loadAffyCsv", t[3])
  t <- ST(snp6.loadAffySeqCsv(db, csvSeqFile, cdfFile, batch_size=batch_size))
  printTime("loadAffySeqCsv", t[3])
  t <- ST(snp6.loadAffyCsv.cnv(db, csvFileCnv, batch_size=batch_size))
  printTime("loadAffyCsv-CNV", t[3])
  t <- ST(snp6.loadAffySeqCsv.cnv(db, csvSeqFileCnv, cdfFile, batch_size=batch_size))
  printTime("loadAffySeqCsv-CNV", t[3])
  t <- ST({
    snp6.sortFeatureTables(db)
    snp6.createIndicesDb(db)
    snp6.createTableInfoTable(db)
    snp6.createFeatureTableInfo(db)
  })
  printTime("DB sort, index creation", t[3])
  
  t <- ST({
    seqMat <- createSeqMat(db)
    save(seqMat, file=matFile, compress=TRUE)
    seqMatCnv <- createSeqMat(db)
    save(seqMatCnv, file=matFileCnv, compress=TRUE)
  })
  printTime("sequence matrix", t[3])
  closeDb(db)
}

# hacked by VC -- original code up above in loadAffyCsvNOCYTOBAND
#
snp6.loadAffyCsv <- function(db, csvFile, batch_size=5000) {
  con <- file(csvFile, open="r")
  on.exit(close(con))

  getFragLength <- function(v){
    tmp <- sapply(strsplit(v, " // "), function(obj) obj[[1]])
    tmp[tmp == "---"] <- NA
    as.integer(tmp)
  }
    
##    wantedCols <- c(1,2,3,4,7,8,10,12,13,14,17) # added 10/14

  wantedCols <- c(1, 2, 3, 4, 5, 6, 8, 10, 11, 12, 15)
  df <- read.table(con, sep=",", stringsAsFactors=FALSE, nrows=10,
                   na.strings="---", header=TRUE)[, wantedCols]
  header <- gsub(".", "_", names(df), fixed=TRUE)
  names(df) <- header

  FRAG_COL <- "Fragment_Length_Start_Stop"
  df[ , FRAG_COL] <- getFragLength(df[ , FRAG_COL])
  
  db_cols <- c("affy_snp_id", "dbsnp_rs_id", "chrom",
               "physical_pos", "strand", "cytoband", "allele_a",
               "allele_b", "gene_assoc", "fragment_length")

  val_holders <- c(":Affy_SNP_ID", ":dbSNP_RS_ID", ":Chromosome",
                   ":Physical_Position", ":Strand", ":Cytoband", ":Allele_A",
                   ":Allele_B", ":Associated_Gene", ":Fragment_Length_Start_Stop")

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
    df[ , FRAG_COL] <- getFragLength(df[ , FRAG_COL])
    dbBeginTransaction(db)
    dbGetPreparedQuery(db, sql, bind.data=df)
    dbCommit(db)
  }
}

snp6.loadAffyCsv.cnv <- function(db, csvFile, batch_size=5000) {
  getFragLength <- function(v){
    tmp <- sapply(strsplit(v, " // "), function(obj) obj[[1]])
    tmp[tmp == "---"] <- NA
    as.integer(tmp)
  }

  getPAR <- function(theDF){
    ## PAR: Pseudo-Autosomal Region
    ##   0: No / 1: PAR1 / 2: PAR2
    PAR <- rep(0, nrow(theDF))
    PAR[theDF[["ChrX_pseudo_autosomal_region_1"]] == 1] <- 1
    PAR[theDF[["ChrX_pseudo_autosomal_region_2"]] == 2] <- 2
    theDF[["XPAR"]] <- PAR
    theDF[["ChrX_pseudo_autosomal_region_1"]] <- NULL
    theDF[["ChrX_pseudo_autosomal_region_2"]] <- NULL
    theDF
  }

  ## CNV probes
  con <- file(csvFile, open="r")
  on.exit(close(con))
  
  wantedCols <- c(1, 2, 3, 4, 5, 6, 7, 8, 10, 13)
  df <- read.table(con, sep=",", stringsAsFactors=FALSE, nrows=10,
                   na.strings="---", header=TRUE)[, wantedCols]
  header <- gsub(".", "_", names(df), fixed=TRUE)
  names(df) <- header

  FRAG_COL <- "Fragment_Length_Start_Stop"
  df[ , FRAG_COL] <- getFragLength(df[ , FRAG_COL])

  df <- getPAR(df)
  df[["Strand"]] <- ifelse(df[["Strand"]] == "+", SENSE, ANTISENSE)
  db_cols <- c("chrom", "chrom_start", "chrom_stop", "strand",
               "cytoband", "gene_assoc", "fragment_length", "xpar")
  
  val_holders <- c(":Chromosome", ":Chromosome_Start", ":Chromosome_Stop",
                   ":Strand", ":Cytoband", ":Associated_Gene",
                   ":Fragment_Length_Start_Stop", ":XPAR")
  
  exprs <- paste(db_cols, " = ", val_holders, sep="", collapse=", ")
  sql <- paste("update featureSetCNV set ", exprs,
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
    df[ , FRAG_COL] <- getFragLength(df[ , FRAG_COL])
    df <- getPAR(df)
    df[["Strand"]] <- ifelse(df[["Strand"]] == "+", SENSE, ANTISENSE)
    dbBeginTransaction(db)
    dbGetPreparedQuery(db, sql, bind.data=df)
    dbCommit(db)
  }


}

snp6.loadAffySeqCsv.cnv <- function(db, csvFile, cdfFile, batch_size=5000) {
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
    header <- c("fset.name", "x", "y", "offset", "seq", "tstrand",
                "type", "chromosome", "strand", "probe_start_pos",
                "in_xpar1", "in_xpar2", "nsp_frag_start", "nsp_frag_end",
                "sty_frag_start", "sty_frag_end")
    wanted <- 1:7
    header <- header[wanted]
    done <- FALSE
    pmdf <- read.table(con, sep="\t", stringsAsFactors=FALSE, nrows=1, header=FALSE)[, wanted]
    while (!done) {
        pmdf <- read.table(con, sep="\t", stringsAsFactors=FALSE,
                           nrows=batch_size, na.strings="---",
                           header=FALSE)[, wanted]
        if (nrow(pmdf) < batch_size) {
            done <- TRUE
            if (nrow(pmdf) == 0)
              break
        }
        names(pmdf) <- header
        pmdf[["fid"]] <- xy2i(pmdf[["x"]], pmdf[["y"]])
        N <- nrow(pmdf)
        values <- "(:fid, :offset, :tstrand, :seq)"
        sql <- paste("insert into sequenceCNV values", values)
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, sql, bind.data=pmdf)
        dbCommit(db)
    }

}

