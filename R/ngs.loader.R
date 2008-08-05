"loadUnitNames.ngs" <- 
function(db, unames) {
    dbBeginTransaction(db)
    ## To use an auto-incrementing field via RSQLite, you need
    ## to be careful to pass integer NA's
    df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
    values <- "(:id, :name)"
    sql <- "insert into featureSet (fsetid, man_fsetid) values"
    dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
    dbCommit(db)
}

"loadUnits.ngs" <- 
function(db, batch, isQc=FALSE) {
## debug
#batch = ndfdata[-controls,]
#isQc=FALSE
## create a couple of MISMATCH Probes for testing
#batch[c(4,5,6),"MISMATCH"] <- 1
#batch[1:6,"MATCH_INDEX"] <- rep(1:3,2)
## end debug

    pmfeature <- "pmfeature_tmp"
    mmfeature <- "mmfeature_tmp"
    pmmm      <- "pm_mm"
    if (isQc) {
        pmfeature <- "qcpmfeature"
        mmfeature <- "qcmmfeature"
        pmmm      <- "qcpm_qcmm"
    }

    ## fill featureSet object with SEQ_ID
    loadUnitNames.ngs(db, unique(batch[["SEQ_ID"]]))
   
    batch <- batch[order(batch$CONTAINER,batch$SEQ_ID,batch$POSITION),]
    pcounts <- table(paste(batch$CONTAINER,batch$SEQ_ID,sep="."))
    batch$ATOM <- unlist(sapply(pcounts,function(x) 1:x))
    batch <- batch[order(batch$fid),]

    batchMat <- batch[, c("fid", "CONTAINER", "FEATURE_ID", "MISMATCH", "MATCH_INDEX", "ATOM", "PROBE_ID", "X", "Y")]
    ## Find internal featureSet IDs for these feature
    fset <- dbGetQuery(db, "SELECT fsetid, man_fsetid FROM featureSet")
    batchIds <- fset[match(batch[["SEQ_ID"]], fset[["man_fsetid"]]), "fsetid"]
    batchMat <- cbind(batchMat, fsetid=batchIds)
    
    ## Insert pm and mm into respective tables
    isPm <- batchMat[["MISMATCH"]] == 0
    isMm <- batchMat[["MISMATCH"]] > 0 & batchMat[["MISMATCH"]] < 10000 ## NGS can have mismatch values > 10000; May be overridden by users if MISMATCH in probe file >= 10,000.
    
    batchMat[["MISMATCH"]] <- NULL
    names(batchMat) <- c("fid", "container", "unit_id", "match_index", "atom", "probe_id", "x", "y", "fsetid")

    ## insert Pm value
    values <- "(:fid, :container, :unit_id, :match_index, :atom, :probe_id, :x, :y, :fsetid)"
    sql <- paste("insert into", pmfeature, "values", values)
    dbBeginTransaction(db)
    rset <- dbSendPreparedQuery(db, sql, batchMat[isPm, ])
    dbClearResult(rset)
    
    if (sum(isMm) > 0){
        sql <- paste("insert into", mmfeature, "values", values)
        rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[isMm, ]))
        dbClearResult(rset)
    }
    dbCommit(db)
    
    ## Insert pm <--> mm link
    ## not necessarily exists for NGS (or any other PM-only arrays) 
    ## NGS can have MM, but usually not
    ## MAKE ME MORE ROBUST!!!!
    int <- intersect(batchMat[isPm, "match_index"], batchMat[isMm, "match_index"])
    if (length(int) > 0){
        values <- "(:pmi, :mmi)"
        sql <- paste("insert into", pmmm, "values", values)
        df1 <- subset(batchMat[isPm,],  batchMat[isPm,"match_index"] %in% int)
        df2 <- subset(batchMat[isMm,],  batchMat[isPm,"match_index"] %in% int)
        dbBeginTransaction(db)
        rset <- dbSendPreparedQuery(db, sql,
            data.frame(pmi=df1[["fid"]],mmi=df2[["fid"]]))
        dbClearResult(rset)
        dbCommit(db)
    }
}


## loadUnitsByBatch.ngs <- function(db, ndfFile, batch_size=10000,
##                              max_units=NULL, verbose=FALSE) {
"loadUnitsByBatch.ngs" <- 
function(db, ndfdata, batch_size=10000,
                             max_units=NULL, verbose=FALSE) {
## debug
#max_units=NULL
#verbose=FALSE
## end debug

    ctrl.containers <- c("CODE", "RANDOM", "CONTROL")
    controls <- sort(as.integer(unlist(sapply(ctrl.containers, grep, toupper(ndfdata[["CONTAINER"]])))))
    ## fill qcpmfeature, qcmmfeature, qcpm_qcmm tables    
    loadUnits.ngs(db, ndfdata[controls,], isQc=TRUE)
    ## fill pmfeature, mmfeature, pm_qcmm tables    
    loadUnits.ngs(db, ndfdata[-controls,], isQc=FALSE)
    ## fill sequences table    
    seq.df <- data.frame(fid=ndfdata$fid, seq=ndfdata$PROBE_SEQUENCE, interrogation_position=ndfdata$POSITION)
    sql <- "INSERT INTO sequence (fid, seq, interrogation_position) VALUES (:fid, :seq, :interrogation_position)"
    dbBeginTransaction(db)
    dbGetPreparedQuery(db, sql, bind.data=seq.df)
    dbCommit(db)
}

"loadNgd.ngs" <- 
function(db, ngdfile){ ## I had to fix ngdfile, there as a bad character somewhere and the whole table wasn't being loaded, specified comment column as character
    ngddata <- read.delim(ngdfile, sep="\t", as.is=TRUE, header=TRUE)
    ngddata <- as.data.frame(ngddata)
    colnames(ngddata) <- c("man_fsetid","comment")
 
    fset <- dbGetQuery(db, "SELECT fsetid, man_fsetid FROM featureSet")
    ngdIds <- match(fset[["man_fsetid"]],ngddata[["man_fsetid"]])
    ngddata <- cbind(fset,comment=ngddata[ngdIds,"comment"])

    db_cols <- c("comment")
    val_holders <- c(":comment")
    exprs <- paste(db_cols, " = ", val_holders, sep="", collapse=", ")
    sql <- paste("update featureSet set ", exprs,
                 "where man_fsetid = :man_fsetid")

    dbBeginTransaction(db)
    dbGetPreparedQuery(db, sql, bind.data=ngddata)
    dbCommit(db)
}
## TODO: fix for tiling arrays
"loadPos.ngs" <- 
function(db, posfile){
    posdata <- read.delim(posfile, as.is=TRUE, header=TRUE)
    posdata <- posdata[, c("SEQ_ID", "PROBE_ID", "POSITION")]
    posdata[["id"]] <- paste(posdata[["SEQ_ID"]], posdata[["PROBE_ID"]], sep=":::")
    posdata <- posdata[order(posdata[["id"]]),]
    df1 <- dbGetQuery(db, "SELECT fid, man_fsetid, probe_id FROM pmfeature_tmp, featureSet WHERE pmfeature_tmp.fsetid = featureSet.fsetid")
    df1[["id"]] <- paste(df1[["man_fsetid"]], df1[["probe_id"]], sep=":::")
    idx <- match(posdata[["id"]], df1[["id"]])
    df1 <- df1[idx,]
    tmp1 <- data.frame(fid=as.integer(df1$fid), position=as.integer(posdata$POSITION))
    sql <- "UPDATE pmfeature_tmp SET position = :position WHERE fid = :fid"
    dbBeginTransaction(db)
    dbGetPreparedQuery(db, sql, bind.data=tmp1)
    dbCommit(db)
    
    df1 <- dbGetQuery(db, "SELECT fid, man_fsetid, probe_id FROM mmfeature_tmp, featureSet WHERE mmfeature_tmp.fsetid = featureSet.fsetid")
    if (nrow(df1) > 0){
        df1[["id"]] <- paste(df1[["man_fsetid"]], df1[["probe_id"]], sep=":::")
        idx <- match(posdata[["id"]], df1[["id"]])
        df1 <- df1[idx,]
        tmp1 <- data.frame(fid=as.integer(df1$fid), position=as.integer(posdata$POSITION))
        sql <- "UPDATE mmfeature_tmp SET position = :position WHERE fid = :fid"
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, sql, bind.data=tmp1)
        dbCommit(db)
    }
}


## buildPdInfoDb.ngs <- function(ndfFile, posFile, csvSeqFile, dbFile, matFile,
##                           batch_size=10000, verbose=FALSE) {
"buildPdInfoDb.ngsExprs" <- 
function(ndfdata, ngdFile, dbFile,
                          batch_size=10000, verbose=FALSE) {
## debug
#ndfdata
#ngdFile = object@ngdFile
#dbFile = dbFilePath 
#batch_size=batch_size
#verbose=!quiet
## end debug
    ST <- system.time
    printTime <- function(msg, t) {
        if (verbose) {
            m <- paste(msg, "took %.2f sec\n")
            cat(sprintf(m, t))
        }
    }
  
    # Initialize DB, create tables
    db <- initDb.ngs(dbFile)
    # Load primary feature information and sequence data    
    t <- ST(loadUnitsByBatch.ngs(db, ndfdata, batch_size=batch_size))
    printTime("loadUnitsByBatch", t[3])
    # load 
    if (file.exists(ngdFile)){
        t <- ST(loadNgd.ngs(db, ngdFile))
        printTime("loadNgd.ngs", t[3])
    }
    t <- ST({
        sortFeatureTables.ngs(db)
        createIndicesDb.ngs(db)
        createTableInfoTable.ngs(db)
        createFeatureTableInfo.ngs(db)
    })
    printTime("DB sort, index creation", t[3])
    closeDb(db)
}

## TODO: Fix for tiling arrays
"buildPdInfoDb.ngsTiling" <- 
function(ndfdata, posFile, dbFile,
                          batch_size=10000, verbose=FALSE) {
    ST <- system.time
    printTime <- function(msg, t) {
      if (verbose) {
        m <- paste(msg, "took %.2f sec\n")
        cat(sprintf(m, t))
      }
    }
   
    db <- initDb.ngs(dbFile)
    t <- ST(loadUnitsByBatch.ngs(db, ndfdata, batch_size=batch_size))
    printTime("loadUnitsByBatch", t[3])
    t <- ST(loadPos.ngs(db, posFile))
    printTime("loadAffySeqCsv", t[3])
    t <- ST({
        sortFeatureTables.ngs(db)
        createIndicesDb.ngs(db)
        createTableInfoTable.ngs(db)
        createFeatureTableInfo.ngs(db)
    })
    printTime("DB sort, index creation", t[3])
    closeDb(db)
}
