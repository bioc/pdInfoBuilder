loadUnitNames.affyTiling <- function(db, unames) {
    dbBeginTransaction(db)
    ## To use an auto-incrementing field via RSQLite, you need
    ## to be careful to pass integer NA's
    df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
    values <- "(:id, :name)"
    sql <- "insert into featureSet (fsetid, man_fsetid) values"
    dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
    dbCommit(db)
}

loadUnits.affyTiling <- function(db, batch, nx, isQc=FALSE) {
    pmfeature <- "pmfeature_tmp"
    mmfeature <- "mmfeature_tmp"
    pmmm <- "pm_mm"
    if (isQc) {
        pmfeature <- "qcpmfeature"
        mmfeature <- "qcmmfeature"
        pmmm <- "qcpm_qcmm"
    }

    batchMat <- do.call(rbind, lapply(batch, readCdfUnitToMat.affyTiling, nx=nx))

    sets <- as.character(unlist(sapply(batch, function(x) paste(x$seqInfo$fullname, x$startpos, sep=":"))))

    tmp.df <- function(x){
      tmp <- data.frame(id=as.integer(NA),
                        man_fsetid=paste(x$seqInfo$fullname, x$startpos, sep=":"),
                        groupname=x$seqInfo$groupname,
                        version=x$seqInfo$version,
                        fullname=x$seqInfo$fullname,
                        name=x$seqInfo$name,
                        stringsAsFactors=FALSE)
      tmp[!duplicated(tmp[["man_fsetid"]]), ]
    }

    super.set <- do.call("rbind", lapply(batch, tmp.df))
    values <- "(:id, :man_fsetid, :groupname, :version, :fullname, :name)"
    sql <- "INSERT INTO featureSet VALUES"
    dbBeginTransaction(db)
    dbGetPreparedQuery(db, paste(sql, values), bind.data=super.set)
    dbCommit(db)
    rm(super.set)
    
##    loadUnitNames.affyTiling(db, unique(sets))

    ## Find internal featureSet IDs for these features
    ## Using subsets of 1000 units
    tmp <- NULL
    us <- unique(sets)
    subsets <- rep(1:length(us), length.out=length(us), each=1000)
    subgrps <- split(us, subsets)
    for (i in 1:length(subgrps)){
      sql <- paste("select man_fsetid, fsetid from featureSet where man_fsetid in (",
                   paste('"', subgrps[[i]], '"', sep="", collapse=","),
                 ") order by fsetid")
      tmp <- rbind(tmp, dbGetQuery(db, sql))
    }
    
##    batchIds <- dbGetQuery(db, sql)
    batchIds <- tmp
    idx <- match(sets, batchIds$man_fsetid)
    ids <- batchIds$fsetid[idx]
#    if (length(ids) %in% c(622, 417497)) browser()
    pmmmtest <- any(batchMat[, "ispm"] == 0)
    if (pmmmtest){
      batchMat <- rbind(cbind(batchMat[batchMat[,"ispm"]==1,], fsetid=ids),
                        cbind(batchMat[batchMat[,"ispm"]==0,], fsetid=ids))
    }else{
      batchMat <- cbind(batchMat[batchMat[,"ispm"]==1,], fsetid=ids)
    }
    theDups <- duplicated(batchMat[["indices"]])
    batchMat <- batchMat[!theDups,]
    rm(theDups)
    
    ## Insert pm and mm into respective tables
    isPm <- as.logical(batchMat[, "ispm"]) 
    values <- "(:indices, :fsetid, :strand, :startpos,  :atom, :x, :y)"
    sql <- paste("insert into", pmfeature, "values", values)
    dbBeginTransaction(db)
    rset <- dbSendPreparedQuery(db, sql, batchMat[isPm, ])
    dbClearResult(rset)

    if (pmmmtest){
      sql <- paste("insert into", mmfeature, "values", values)
      rset <- dbSendPreparedQuery(db, sql, batchMat[!isPm, ])
      dbClearResult(rset)
    }
    dbCommit(db)

    ## Insert pm <--> mm link
    if (pmmmtest){
      pm.fsetid.atom.set = subset(batchMat, ispm == 1)[, c("fsetid", "atom", "indices")]
      pm.fsetid.atom = paste(pm.fsetid.atom.set[["fsetid"]], pm.fsetid.atom.set[["atom"]], sep=":")
      mm.fsetid.atom.set = subset(batchMat, ispm == 0)[, c("fsetid", "atom", "indices")]
      mm.fsetid.atom = paste(mm.fsetid.atom.set[["fsetid"]], mm.fsetid.atom.set[["atom"]], sep=":")
      common <- sort(intersect(pm.fsetid.atom, mm.fsetid.atom))

      pmset <- match(common, pm.fsetid.atom)
      mmset <- match(common, mm.fsetid.atom)
    
      values <- "(:pmi, :mmi)"
      sql <- paste("insert into", pmmm, "values", values)
      dbBeginTransaction(db)
    
##     rset <- dbSendPreparedQuery(db, sql,
##                                 data.frame(pmi=batchMat[isPm, "indices"],
##                                            mmi=batchMat[!isPm, "indices"]))
    
      rset <- dbSendPreparedQuery(db, sql,
                                  data.frame(pmi=pm.fsetid.atom.set[pmset, "indices"],
                                             mmi=mm.fsetid.atom.set[mmset, "indices"]))
      dbClearResult(rset)
      dbCommit(db)
    }

    ## Insert sequence info
    values <- "(:indices, :probeseq)"
    sql <- paste("INSERT INTO sequence VALUES", values)
    dbBeginTransaction(db)
    rset <- dbSendPreparedQuery(db, sql, batchMat)
    dbClearResult(rset)
    dbCommit(db)
}


loadUnitsByBatch.affyTiling <- function(db, bpmapFile, batch_size=1, nx,
                             max_units=NULL, verbose=FALSE) {
    seqs <- readBpmapSeqinfo(bpmapFile)
    grps <- sapply(seqs, "[[", "parameters")
    whQc <- which(tolower(grps) != "tiling")

    offset <- 1
    if (length(whQc)) {                 # load all QC at once
      qcunits <- readBpmap(bpmapFile, whQc)
      loadUnits.affyTiling(db, qcunits, nx=nx, isQc=TRUE)
      offset <- max(whQc) + 1
    }

    for (wanted in offset:length(grps)){
      message("Processing unit ", wanted, " out of ", length(grps), ".")
      vvunits <- readBpmap(bpmapFile, wanted)
      loadUnits.affyTiling(db, vvunits, nx=nx)
    }
    
##     if (is.null(max_units))
##       max_units <- length(grps)
## 
##     extra <- (max_units-length(whQc)) %% batch_size
##     num_batches <- (max_units-length(whQc)) %/% batch_size
##     if (extra != 0)
##       num_batches <- num_batches + 1
##     done <- FALSE
##     while (!done) {
##         end <- min(offset + batch_size, max_units)
##         if (end == max_units)
##           done <- TRUE
##         wanted <- seq.int(offset, end)
##         offset <- offset + batch_size + 1
##         vvunits <- readBpmap(bpmapFile, wanted)
##         loadUnits.affyTiling(db, vvunits, nx=nx)
##     }
}


buildPdInfoDb.affyTiling <- function(bpmapFile,  cifFile, dbFile, matFile,
                          batch_size=10000, verbose=FALSE) {

    ST <- system.time
    printTime <- function(msg, t) {
        if (verbose) {
            m <- paste(msg, "took %.2f sec\n")
            cat(sprintf(m, t))
        }
    }

    db <- initDb.affyTiling(dbFile)

    cif <- scan(cifFile, what="c", quiet=TRUE)
    nx <- as.integer(unlist(strsplit(cif[grep("Cols", cif)], "="))[2])
    rm(cif)
    
    t <- ST(loadUnitsByBatch.affyTiling(db, bpmapFile, batch_size=batch_size, nx=nx))
    printTime("loadUnitsByBatch", t[3])
    t <- ST({
        sortFeatureTables.affyTiling(db)
        createIndicesDb.affyTiling(db)
        createTableInfoTable.affyTiling(db)
        createFeatureTableInfo.affyTiling(db)
    })
    printTime("DB sort, index creation", t[3])

##     t <- ST({
##         seqMat <- createSeqMat(db)
##         save(seqMat, file=matFile, compress=TRUE)
##     })
##     printTime("sequence matrix", t[3])
    closeDb(db)
}

