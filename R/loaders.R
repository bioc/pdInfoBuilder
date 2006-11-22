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

loadUnits <- function(db, batch) {
    batchMat <- do.call(rbind, lapply(batch, pmmmBlockToMat))

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
    batchMat <- cbind(batchMat, fsetid=batchIds)

    ## Insert pm and mm into respective tables
    isPm <- as.logic(batchmat[, "ispm"])
    values <- "(:indices, :strand, :allele, :fsetid, :expos, :x, :y)"
    sql <- paste("insert into pmfeature_tmp values", values)
    dbBeginTransaction(db)
    rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[isPm, ]))
    dbClearResult(rset)

    sql <- paste("insert into mmfeature_tmp values", values)
    rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[!isPm, ]))
    dbClearResult(rset)
    dbCommit(db)

    ## Insert pm <--> mm link
    values <- "(:pmi, :mmi)"
    sql <- paste("insert into pm_mm values", values)
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

    wantedCols <- c(1,2,3,4,7,8,12,13)
    df <- read.table(con, sep=",", stringsAsFactors=FALSE, nrows=10,
                     na.strings="---", header=TRUE)[, wantedCols]
    header <- gsub(".", "_", names(df), fixed=TRUE)
    names(df) <- header
    
    db_cols <- c("affy_snp_id", "dbsnp_rs_id", "chrom",
                 "physical_pos", "strand", "allele_a", "allele_b")
    
    val_holders <- c(":Affy_SNP_ID", ":dbSNP_RS_ID", ":Chromosome",
                     ":Physical_Position", ":Strand", ":Allele_A",
                     ":Allele_B")

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


buildPdInfoDb <- function(cdfFile, csvFile, dbFile) {
    db <- initDb(dbFile)
    loadUnitsByBatch(db, cdfFile)
    loadAffyCsv(db, csvFile)
    createIndicesDb(db)
    closeDb(db)
}
