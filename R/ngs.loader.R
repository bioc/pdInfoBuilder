loadUnitNames.ngs <- function(db, unames) {
    dbBeginTransaction(db)
    ## To use an auto-incrementing field via RSQLite, you need
    ## to be careful to pass integer NA's
    df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
    values <- "(:id, :name)"
    sql <- "insert into featureSet (fsetid, man_fsetid) values"
    dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
    dbCommit(db)
}

loadUnits.ngs <- function(db, batch, isQc=FALSE, ncol) {
    pmfeature <- "pmfeature_tmp"
    mmfeature <- "mmfeature_tmp"
    pmmm <- "pm_mm"
    if (isQc) {
        pmfeature <- "qcpmfeature"
        mmfeature <- "qcmmfeature"
        pmmm <- "qcpm_qcmm"
    }
    batchMat <- batch[, c("fid", "CONTAINER", "FEATURE_ID", "MISMATCH", "MATCH_INDEX", "PROBE_ID", "POSITION", "X", "Y")]

    loadUnitNames.ngs(db, unique(batch[["SEQ_ID"]]))

    ## Find internal featureSet IDs for these feature
    fset <- dbGetQuery(db, "SELECT fsetid, man_fsetid FROM featureSet")
    batchIds <- fset[match(batch[["SEQ_ID"]], fset[["man_fsetid"]]), "fsetid"]
    batchMat <- cbind(batchMat, fsetid=batchIds)

    ## Insert pm and mm into respective tables
    isPm <- batchMat[["MISMATCH"]] == 0
    batchMat[["MISMATCH"]] <- NULL
##     xy2i <- function(x, y) x  + (y-1) * ncol
##     batchMat <- cbind(fid=as.integer(xy2i(batchMat$X, batchMat$Y)), batchMat)
##     batchMat <- batchMat[order(batchMat[,1]),]
##     batchMat[,1] <- 1:nrow(batchMat)
    names(batchMat) <- c("fid", "container", "unit_id", "match_index", "probe_id", "position", "x", "y", "fsetid")

    values <- "(:fid, :container, :unit_id, :match_index, :probe_id, :position, :x, :y, :fsetid)"
    sql <- paste("insert into", pmfeature, "values", values)

    dbBeginTransaction(db)
    rset <- dbSendPreparedQuery(db, sql, batchMat[isPm, ])
    dbClearResult(rset)

    if (sum(!isPm) > 0){
      sql <- paste("insert into", mmfeature, "values", values)
      rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[!isPm, ]))
      dbClearResult(rset)
    }
    dbCommit(db)

    ## Insert pm <--> mm link
    ## not necessarily exists for NGS (or any other PM-only arrays)
    ## MAKE ME MORE ROBUST!!!!

    int <- intersect(batchMat[isPm, "match_index"], batch[!isPm, "match_index"])
    if (length(int) > 0){
      values <- "(:pmi, :mmi)"
      sql <- paste("insert into", pmmm, "values", values)
      df1 <- subset(batchMat[isPm,], match_index %in% int)
      df2 <- subset(batchMat[!isPm,], match_index %in% int)
      dbBeginTransaction(db)
      rset <- dbSendPreparedQuery(db, sql,
                                  data.frame(pmi=df1[["fid"]],
                                             mmi=df2[["fid"]]))
      dbClearResult(rset)
      dbCommit(db)
    }
}


## loadUnitsByBatch.ngs <- function(db, ndfFile, batch_size=10000,
##                              max_units=NULL, verbose=FALSE) {

loadUnitsByBatch.ngs <- function(db, ndfdata, batch_size=10000,
                             max_units=NULL, verbose=FALSE) {

  ## not running on batches yet
  ## i have the feeling we can handle all in memmory

##  ndfdata <- read.delim(ndfFile, as.is=TRUE, header=TRUE)
##  ctrl.containers <- c("NGS_CONTROLS", "V_CODE", "H_CODE", "RANDOM", "OLIGO_CONTROL")
  ctrl.containers <- c("CODE", "RANDOM", "CONTROL")
  nx <- max(ndfdata$X)
##  controls <- which(toupper(ndfdata[["CONTAINER"]]) %in% ctrl.containers)
  controls <- sort(as.integer(unlist(sapply(ctrl.containers, grep, toupper(ndfdata[["CONTAINER"]])))))
  loadUnits.ngs(db, ndfdata[controls,], isQc=TRUE, ncol=nx)
  loadUnits.ngs(db, ndfdata[-controls,], isQc=FALSE, ncol=nx)
  seq.df <- data.frame(fid=ndfdata$fid, seq=ndfdata$PROBE_SEQUENCE)
  sql <- "INSERT INTO sequence (fid, seq) VALUES (:fid, :seq)"
  dbBeginTransaction(db)
  dbGetPreparedQuery(db, sql, bind.data=seq.df)
  dbCommit(db)
}


loadPos.ngs <- function(db, posfile){
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

buildPdInfoDb.ngs <- function(ndfdata, posFile, dbFile,
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

##     t <- ST({
##         seqMat <- createSeqMat(db)
##         save(seqMat, file=matFile, compress=TRUE)
##     })
##     printTime("sequence matrix", t[3])
    closeDb(db)
}
