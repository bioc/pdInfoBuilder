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
                 ")")
    batchIds <- dbGetQuery(db, sql)[[1]]
    batchIds <- rep(batchIds, batchLens)
    batchMat <- cbind(batchMat, fsetid=batchIds, offset=NA)

    ## Insert pm and mm into respective tables
    isPm <- as.logical(batchMat[, "ispm"])
    values <- "(:indices, :strand, :allele, :fsetid, :indexpos, :x, :y, :offset)"
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
    offset <- 1
    whQc <- grep("^AFFX", unames)
    if (length(whQc)) {                 # load all QC at once
        qcunits <- readCdf(cdfFile, units=whQc, readGroupDirection=TRUE,
                           readIndices=TRUE, readIsPm=TRUE)
        loadUnits(db, qcunits, isQc=TRUE)
        unames <- unames[-whQc]
        offset <- max(whQc) + 1
    }
    extra <- length(unames) %% batch_size
    num_batches <- length(unames) %/% batch_size
    if (extra != 0)
      num_batches <- num_batches + 1
    done <- FALSE
    if (is.null(max_units))
      max_units <- length(unames)
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


loadAffyCsv <- function(db, csvFile, batch_size=5000) {
    con <- file(csvFile, open="r")
    on.exit(close(con))

    getLength <- function(v){
      tmp <- sapply(strsplit(v, "//"), function(obj) obj[[1]])
      tmp[tmp == "---"] <- NA
      as.integer(tmp)
    }
    
    ## BC: Added column 17 (need fragment length)
    wantedCols <- c(1,2,3,4,7,8,12,13,17)
    df <- read.table(con, sep=",", stringsAsFactors=FALSE, nrows=10,
                     na.strings="---", header=TRUE)[, wantedCols]
    df[,9] <- getLength(df[,9])
    header <- gsub(".", "_", names(df), fixed=TRUE)
    names(df) <- header

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
        df[,9] <- getLength(df[,9])
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


loadAffySeqCsv <- function(db, csvFile, cdfFile, batch_size=5000) {
    cdfHeader <- readCdfHeader(cdfFile)
    ncol <- cdfHeader$ncol
    xy2i <- function(x, y) x + 1 + y * ncol

    complementBase <- function(x, special=FALSE){
      bases <- c("A", "C", "G", "T")
      if (!special) comp <- c("T", "G", "C", "A")
      else comp <- c("G", "T", "A", "C")
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

        ## BC: adding offset to pmfeature
        offsetSql <- paste("update pmfeature_tmp set offset = :offset where pmfeature_tmp.fid = :fid")
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, offsetSql, bind.data=pmdf)
        dbCommit(db)
        
        ## process MMs
        mmSql <- paste("select mm_fid, pm_fid from pm_mm where pm_fid in (",
                       paste(pmdf[["fid"]], collapse=","), ")")
        pairedIds <- dbGetQuery(db, mmSql)
        foundIdIdx <- match(pmdf[["fid"]], pairedIds[["pm_fid"]], 0)
        mmdf <- pmdf[foundIdIdx, ]
        mmdf[["fid"]] <-  pairedIds[["mm_fid"]]

        ## BC: adding offset to mmfeature
        offsetSql <- paste("update mmfeature_tmp set offset = :offset where mmfeature_tmp.fid = :fid")
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, offsetSql, bind.data=mmdf)
        dbCommit(db)

        
        ## Assuming 25mers
        midbase <- substr(mmdf$seq, 13, 13)
        types <- aggregate(mmdf$tallele, by=list(mmdf$fset.name),
                             FUN=function(v) paste(sort(unique(v)), collapse=""))[,2]
        types <- rep(types, as.integer(table(mmdf$fset.name)))
        isSpecial <- (types == "AT" | types == "AG") & mmdf$offset == 0
        rm(types)
        midbase[isSpecial] <- complementBase(midbase[isSpecial], T)
        midbase[!isSpecial] <- complementBase(midbase[!isSpecial])
        rm(isSpecial)
        mmdf$seq <- paste(substr(mmdf$seq, 1, 12), midbase, substr(mmdf$seq, 14, 25), sep="")
        rm(midbase)
        ## end MM seq

        values <- "(:fid, :offset, :tstrand, :tallele, :type, :seq)"
        sql <- paste("insert into sequence values", values)
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, sql, bind.data=pmdf)
        dbGetPreparedQuery(db, sql, bind.data=mmdf)
        dbCommit(db)
    }

}


buildPdInfoDb <- function(cdfFile, csvFile, csvSeqFile, dbFile, matFile) {
    db <- initDb(dbFile)

    loadUnitsByBatch(db, cdfFile, batch_size=600000)
    loadAffyCsv(db, csvFile, batch_size=300000)
    loadAffySeqCsv(db, csvSeqFile, cdfFile, batch_size=300000)

    sortFeatureTables(db)
    createIndicesDb(db)
    createTableInfoTable(db)
    createFeatureTableInfo(db)

    seqMat <- createSeqMat(db)
    save(seqMat, file=matFile, compress=TRUE)
    closeDb(db)
}
