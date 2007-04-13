loadUnitNames <- function(db, unames) {
    dbBeginTransaction(db)
    ## To use an auto-incrementing field via RSQLite, you need
    ## to be careful to pass integer NA's
    df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
    values <- "(:id, :name)"
    sql <- "insert into featureSet (fsetid, man_fsetid) values"
    dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
    dbCommit(db)
}

loadUnits <- function(db, batch, isQc=FALSE) {
    pmfeature <- "pmfeature_tmp"
    mmfeature <- "mmfeature_tmp"
    pmmm <- "pm_mm"
    if (isQc) {
        pmfeature <- "qcpmfeature"
        mmfeature <- "qcmmfeature"
        pmmm <- "qcpm_qcmm"
    }
    batchMat <- do.call(rbind, lapply(batch, readCdfUnitToMat))

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

    ## Insert pm and mm into respective tables
    isPm <- as.logical(batchMat[, "ispm"])
    values <- "(:indices, :strand, :allele, :fsetid, :indexpos, :x, :y)"
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


loadUnitsByBatch <- function(db, cdfFile, batch_size=10000,
                             max_units=NULL, verbose=FALSE) {
    unames <- readCdfUnitNames(cdfFile)
    if (is.null(max_units))
      max_units <- length(unames)
    offset <- 1
    whQc <- grep("^AFFX", unames)
    if (length(whQc)) {                 # load all QC at once
        qcunits <- readCdf(cdfFile, units=whQc, readGroupDirection=TRUE,
                           readIndices=TRUE, readIsPm=TRUE)
        loadUnits(db, qcunits, isQc=TRUE)
##        unames <- unames[-whQc]
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
        loadUnits(db, vvunits)
    }
}


loadAffyCsvNOCYTOBAND <- function(db, csvFile, batch_size=5000) {
    con <- file(csvFile, open="r")
    on.exit(close(con))

    getFragLength <- function(v){
      tmp <- sapply(strsplit(v, " // "), function(obj) obj[[1]])
      tmp[tmp == "---"] <- NA
      as.integer(tmp)
    }
    
    wantedCols <- c(1,2,3,4,7,8,12,13,17)
    df <- read.table(con, sep=",", stringsAsFactors=FALSE, nrows=10,
                     na.strings="---", header=TRUE)[, wantedCols]
    header <- gsub(".", "_", names(df), fixed=TRUE)
    names(df) <- header

    FRAG_COL <- "Fragment_Length_Start_Stop"
    df[ , FRAG_COL] <- getFragLength(df[ , FRAG_COL])

    db_cols <- c("affy_snp_id", "dbsnp_rs_id", "chrom",
                 "physical_pos", "strand", "allele_a",
                 "allele_b", "fragment_length")

    val_holders <- c(":Affy_SNP_ID", ":dbSNP_RS_ID", ":Chromosome",
                     ":Physical_Position", ":Strand", ":Allele_A",
                     ":Allele_B", ":Fragment_Length_Start_Stop")

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


loadAffySeqCsv <- function(db, csvFile, cdfFile, batch_size=5000) {
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
        mmSql <- paste("select mm_fid, pm_fid from pm_mm where pm_fid in (",
                       paste(pmdf[["fid"]], collapse=","), ")")
        pairedIds <- dbGetQuery(db, mmSql)
        foundIdIdx <- match(pmdf[["fid"]], pairedIds[["pm_fid"]], 0)
        mmdf <- pmdf[foundIdIdx, ]
        mmdf[["fid"]] <-  pairedIds[["mm_fid"]]

        ## Assuming 25mers
        midbase <- substr(mmdf$seq, 13, 13)
        types <- apply(table(mmdf$fset.name, mmdf$tallele)>0, 1,
                       function(v) paste(c("A", "C", "G", "T")[v], collapse=""))
        types <- rep(types, as.integer(table(mmdf$fset.name)))
        isSpecial <- (types == "AT" | types == "AG") & mmdf$offset == 0
        rm(types)
        midbase[isSpecial] <- complementBase(midbase[isSpecial], T)
        midbase[!isSpecial] <- complementBase(midbase[!isSpecial])
        rm(isSpecial)
        mmdf$seq <- paste(substr(mmdf$seq, 1, 12), midbase,
                          substr(mmdf$seq, 14, 25), sep="")
        rm(midbase)
        ## end MM seq

        values <- "(:fid, :offset, :tstrand, :tallele, :seq)"
        sql <- paste("insert into sequence values", values)
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, sql, bind.data=pmdf)
        dbGetPreparedQuery(db, sql, bind.data=mmdf)
        dbCommit(db)
    }

}


buildPdInfoDb <- function(cdfFile, csvFile, csvSeqFile, dbFile, matFile,
                          batch_size=10000, verbose=FALSE) {

    ST <- system.time
    printTime <- function(msg, t) {
        if (verbose) {
            m <- paste(msg, "took %.2f sec\n")
            cat(sprintf(m, t))
        }
    }

    db <- initDb(dbFile)

    t <- ST(loadUnitsByBatch(db, cdfFile, batch_size=batch_size))
    printTime("loadUnitsByBatch", t[3])
    t <- ST(loadAffyCsv(db, csvFile, batch_size=batch_size))
    printTime("loadAffyCsv", t[3])
    t <- ST(loadAffySeqCsv(db, csvSeqFile, cdfFile, batch_size=batch_size))
    printTime("loadAffySeqCsv", t[3])
    t <- ST({
        sortFeatureTables(db)
        createIndicesDb(db)
        createTableInfoTable(db)
        createFeatureTableInfo(db)
    })
    printTime("DB sort, index creation", t[3])

    t <- ST({
        seqMat <- createSeqMat(db)
        save(seqMat, file=matFile, compress=TRUE)
    })
    printTime("sequence matrix", t[3])
    closeDb(db)
}

# hacked by VC -- original code up above in loadAffyCsvNOCYTOBAND
#
loadAffyCsv <- function(db, csvFile, batch_size=5000) {
    con <- file(csvFile, open="r")
    on.exit(close(con))

    getFragLength <- function(v){
      tmp <- sapply(strsplit(v, " // "), function(obj) obj[[1]])
      tmp[tmp == "---"] <- NA
      as.integer(tmp)
    }
    
    wantedCols <- c(1,2,3,4,7,8,10,12,13,17) # added 10
    df <- read.table(con, sep=",", stringsAsFactors=FALSE, nrows=10,
                     na.strings="---", header=TRUE)[, wantedCols]
    header <- gsub(".", "_", names(df), fixed=TRUE)
    names(df) <- header

    FRAG_COL <- "Fragment_Length_Start_Stop"
    df[ , FRAG_COL] <- getFragLength(df[ , FRAG_COL])

    db_cols <- c("affy_snp_id", "dbsnp_rs_id", "chrom",
                 "physical_pos", "strand", "cytoband", "allele_a",
                 "allele_b", "fragment_length")

    val_holders <- c(":Affy_SNP_ID", ":dbSNP_RS_ID", ":Chromosome",
                     ":Physical_Position", ":Strand", ":Cytoband", ":Allele_A",
                     ":Allele_B", ":Fragment_Length_Start_Stop")

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

