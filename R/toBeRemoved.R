######################################################
######################################################
#### EVERYTHING IS COMMENTED WITH 5 #'S
######################################################
######################################################

##### ################################################
##### ## affyExpr.initDb.R
##### ################################################
##### initDb.affyExpr <- function(dbname) {
#####     db <- dbConnect(dbDriver("SQLite"), dbname)
##### 
#####     ## Set page size
#####     dbGetQuery(db, setPageSizeSql)
##### 
#####     ## Create tables
#####     ## BC: Soon we need to add a table for the control probes
#####     dbGetQuery(db, createAffyExprFeatureSetSql)
##### 
#####     dbGetQuery(db, sprintf(createAffyExprFeatureSql, "pmfeature_tmp"))
##### 
#####     dbGetQuery(db, sprintf(createAffyExprFeatureSql, "mmfeature_tmp"))
##### 
#####     dbGetQuery(db, sprintf(createAffyExprFeatureSql, "qcpmfeature"))
#####     dbGetQuery(db, sprintf(createAffyExprFeatureSql, "qcmmfeature"))
##### 
#####     dbGetQuery(db, sprintf(createAffyExprPm_MmSql, "pm_mm"))
#####     dbGetQuery(db, sprintf(createAffyExprPm_MmSql, "qcpm_qcmm"))
##### 
#####     dbGetQuery(db, createAffyExprSequenceSql)
##### 
#####     ## Create index
#####     ## NOTE: might be more efficient to create this after loading,
#####     ## but current perf is ok.
#####     sql <- 'create index man_fsetid_idx on featureSet ("man_fsetid")'
#####     dbGetQuery(db, sql)
#####     db
##### }
##### 
##### sortFeatureTables.affyExpr <- function(db) {
#####     dbGetQuery(db, sprintf(createAffyExprFeatureSql, "pmfeature"))
#####     dbGetQuery(db, sprintf(createAffyExprFeatureSql, "mmfeature"))
##### 
#####     ## Reorder XXfeature tables
#####     fillSql <- paste("insert into %s select * from %s order by",
#####                      "fsetid, atom")
#####     dbBeginTransaction(db)
#####     dbGetQuery(db, sprintf(fillSql, "pmfeature", "pmfeature_tmp"))
#####     dbGetQuery(db, sprintf(fillSql, "mmfeature", "mmfeature_tmp"))
#####     dbCommit(db)
#####     ## drop temp tables
#####     dbBeginTransaction(db)
#####     dbGetQuery(db, "drop table pmfeature_tmp")
#####     dbGetQuery(db, "drop table mmfeature_tmp")
#####     dbCommit(db)
##### }
##### 
##### 
##### createIndicesDb.affyExpr <- function(db) {
#####     makeIndex <- function(name, t, cols) {
#####         sql <- paste("create index", name, "on", t,
#####                      paste("(", paste(cols, collapse=","), ")"))
#####         dbGetQuery(db, sql)
#####     }
##### 
#####     ## Create DB indices and fix ordering
#####     makeIndex("pmf_idx_fsetid", "pmfeature", "fsetid")
#####     makeIndex("mmf_idx_fsetid", "mmfeature", "fsetid")
##### 
##### ##    makeIndex("fset_idx_chrom", "featureSet", "chrom")
#####     makeIndex("fset_idx_fsetid", "featureSet", "fsetid")
##### 
#####     ## finally, run analyze (SQLite specific?)
#####     dbGetQuery(db, "analyze")
##### }
##### 
##### 
##### createTableInfoTable.affyExpr <- function(db, verbose=FALSE) {
#####     tables <- dbListTables(db)
#####     counts <- integer(length(tables))
#####     sql <- "select count(*) from %s"
#####     for (i in seq(along=counts)) {
#####         if (verbose)
#####           cat("counting rows in ", tables[i], "\n")
#####         counts[i] <- dbGetQuery(db, sprintf(sql, tables[i]))[[1]][1]
#####     }
##### 
#####     df <- data.frame(tbl=tables, row_count=counts,
#####                      stringsAsFactors=FALSE)
#####     dbWriteTable(db, "table_info", df, row.names=FALSE)
##### }
##### 
##### 
##### createFeatureTableInfo.affyExpr <- function(db, tname) {
#####     return(FALSE)
#####     ## FIXME: add code to determine offsets of sorted
#####     ## strand and allele
##### }
##### 
##### 
##### ################################################
##### ## affyExpr.loaders.R
##### ################################################
##### loadUnitNames.affyExpr <- function(db, unames) {
#####     dbBeginTransaction(db)
#####     ## To use an auto-incrementing field via RSQLite, you need
#####     ## to be careful to pass integer NA's
#####     df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
#####     values <- "(:id, :name)"
#####     sql <- "insert into featureSet (fsetid, man_fsetid) values"
#####     dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
#####     dbCommit(db)
##### }
##### 
##### loadUnits.affyExpr <- function(db, batch, isQc=FALSE) {
#####     pmfeature <- "pmfeature_tmp"
#####     mmfeature <- "mmfeature_tmp"
#####     pmmm <- "pm_mm"
#####     if (isQc) {
#####         pmfeature <- "qcpmfeature"
#####         mmfeature <- "qcmmfeature"
#####         pmmm <- "qcpm_qcmm"
#####     }
#####     batchMat <- do.call(rbind, lapply(batch, readCdfUnitToMat.affyExpr))
##### 
#####     loadUnitNames.affyExpr(db, names(batch))
##### 
#####     ## Find internal featureSet IDs for these features
#####     batchLens <- sapply(batch, function(x)
#####                         sum(sapply(x$groups, function(y)
#####                                    length(y$indices))))
#####     sql <- paste("select fsetid from featureSet where man_fsetid in (",
#####                  paste('"', names(batch), '"', sep="", collapse=","),
#####                  ") order by fsetid")
#####     batchIds <- dbGetQuery(db, sql)[[1]]
#####     batchIds <- rep(batchIds, batchLens)
#####     batchMat <- cbind(batchMat, fsetid=batchIds)
##### 
#####     ## Insert pm and mm into respective tables
#####     isPm <- as.logical(batchMat[, "ispm"])
#####     values <- "(:indices, :fsetid, :pbase, :tbase,  :atom, :x, :y)"
#####     sql <- paste("insert into", pmfeature, "values", values)
#####     dbBeginTransaction(db)
#####     rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[isPm, ]))
#####     dbClearResult(rset)
##### 
#####     sql <- paste("insert into", mmfeature, "values", values)
#####     rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[!isPm, ]))
#####     dbClearResult(rset)
#####     dbCommit(db)
##### 
#####     ## Insert pm <--> mm link
#####     values <- "(:pmi, :mmi)"
#####     sql <- paste("insert into", pmmm, "values", values)
#####     dbBeginTransaction(db)
#####     rset <- dbSendPreparedQuery(db, sql,
#####                                 data.frame(pmi=batchMat[isPm, "indices"],
#####                                            mmi=batchMat[!isPm, "indices"]))
#####     dbClearResult(rset)
#####     dbCommit(db)
##### }
##### 
##### 
##### loadUnitsByBatch.affyExpr <- function(db, cdfFile, batch_size=10000,
#####                              max_units=NULL, verbose=FALSE) {
#####     unames <- readCdfUnitNames(cdfFile)
#####     if (is.null(max_units))
#####       max_units <- length(unames)
#####     offset <- 1
#####     whQc <- grep("^AFFX", unames)
#####     if (length(whQc)) {                 # load all QC at once
#####         qcunits <- readCdf(cdfFile, units=whQc, readGroupDirection=TRUE,
#####                            readIndices=TRUE, readIsPm=TRUE)
#####         loadUnits.affyExpr(db, qcunits, isQc=TRUE)
#####         offset <- max(whQc) + 1
#####     }
#####     extra <- (length(unames)-length(whQc)) %% batch_size
#####     num_batches <- (length(unames)-length(whQc)) %/% batch_size
#####     if (extra != 0)
#####       num_batches <- num_batches + 1
#####     done <- FALSE
#####     while (!done) {
#####         end <- min(offset + batch_size, max_units)
#####         if (end == max_units)
#####           done <- TRUE
#####         wanted <- seq.int(offset, end)
#####         offset <- offset + batch_size + 1
#####         vvunits <- readCdf(cdfFile, units=wanted, readGroupDirection=TRUE,
#####                            readIndices=TRUE, readIsPm=TRUE)
#####         loadUnits.affyExpr(db, vvunits)
#####     }
##### }
##### 
##### 
##### loadAffySeqCsv.affyExpr <- function(db, csvFile, cdfFile, batch_size=5000) {
#####     cdfHeader <- readCdfHeader(cdfFile)
#####     ncol <- cdfHeader$ncol
#####     xy2i <- function(x, y) x + 1 + y * ncol
##### 
#####     complementBase <- function(x, special=FALSE){
#####       bases <- c("A", "C", "G", "T")
#####       if (!special)
#####         comp <- c("T", "G", "C", "A")
#####       else
#####         comp <- c("G", "T", "A", "C")
#####       comp[match(x, bases)]
#####     }
##### 
#####     con <- file(csvFile, open="r")
#####     on.exit(close(con))
#####     header <- c("fset.name", "x", "y", "interrogation_position", "seq", "tstrand")
#####     done <- FALSE
#####     tmp <- read.table(con, sep="\t", stringsAsFactors=FALSE,
#####                       nrows=1, na.strings="---", header=FALSE)
#####     while (!done) {
#####         pmdf <- read.table(con, sep="\t", stringsAsFactors=FALSE,
#####                            nrows=batch_size, na.strings="---",
#####                            header=FALSE)
#####         if (nrow(pmdf) < batch_size) {
#####             done <- TRUE
#####             if (nrow(pmdf) == 0)
#####               break
#####         }
#####         names(pmdf) <- header
#####         pmdf[["fid"]] <- xy2i(pmdf[["x"]], pmdf[["y"]])
#####         pmdf[["tstrand"]] <- ifelse(tolower(pmdf[["tstrand"]]) == "sense", SENSE, ANTISENSE)
#####         N <- nrow(pmdf)
##### 
#####         ## FIXME: Are these files PM-only???
#####  
#####         values <- "(:fid, :tstrand, :interrogation_position, :seq)"
#####         sql <- paste("insert into sequence values", values)
#####         dbBeginTransaction(db)
#####         dbGetPreparedQuery(db, sql, bind.data=pmdf)
##### ##        dbGetPreparedQuery(db, sql, bind.data=mmdf)
#####         dbCommit(db)
#####     }
##### 
##### }
##### 
##### 
##### buildPdInfoDb.affyExpr <- function(cdfFile, csvFile, csvSeqFile, dbFile, matFile,
#####                           batch_size=10000, verbose=FALSE) {
##### 
#####     ST <- system.time
#####     printTime <- function(msg, t) {
#####         if (verbose) {
#####             m <- paste(msg, "took %.2f sec\n")
#####             cat(sprintf(m, t))
#####         }
#####     }
##### 
#####     db <- initDb.affyExpr(dbFile)
##### 
#####     t <- ST(loadUnitsByBatch.affyExpr(db, cdfFile, batch_size=batch_size))
#####     printTime("loadUnitsByBatch", t[3])
#####     t <- ST(loadAffyCsv.affyExpr(db, csvFile, batch_size=batch_size))
#####     printTime("loadAffyCsv", t[3])
#####     t <- ST(loadAffySeqCsv.affyExpr(db, csvSeqFile, cdfFile, batch_size=batch_size))
#####     printTime("loadAffySeqCsv", t[3])
#####     t <- ST({
#####         sortFeatureTables.affyExpr(db)
#####         createIndicesDb.affyExpr(db)
#####         createTableInfoTable.affyExpr(db)
#####         createFeatureTableInfo.affyExpr(db)
#####     })
#####     printTime("DB sort, index creation", t[3])
##### 
##### ##     t <- ST({
##### ##         seqMat <- createSeqMat(db)
##### ##         save(seqMat, file=matFile, compress=TRUE)
##### ##     })
##### ##     printTime("sequence matrix", t[3])
#####     closeDb(db)
##### }
##### 
##### loadAffyCsv.affyExpr <- function(db, csvFile, batch_size=5000) {
#####     con <- file(csvFile, open="r")
#####     on.exit(close(con))
##### 
#####     df <- read.table(con, sep=",", stringsAsFactors=FALSE, nrows=10,
#####                      na.strings="---", header=TRUE)
#####     wantedCols <- match(c("Probe.Set.ID", "Alignments", "Gene.Symbol", "Chromosomal.Location", "Ensembl"), names(df))
#####     df <- df[, wantedCols]
#####     header <- gsub(".", "_", names(df), fixed=TRUE)
#####     names(df) <- header
##### 
#####     db_cols <- c("alignment", "gene_symbol",  "chrom", "ensembl")
#####     val_holders <- c(":Alignments", ":Gene_Symbol", ":Chromosomal_Location", ":Ensembl")
#####     exprs <- paste(db_cols, " = ", val_holders, sep="", collapse=", ")
#####     sql <- paste("update featureSet set ", exprs,
#####                  "where man_fsetid = :Probe_Set_ID")
##### 
#####     dbBeginTransaction(db)
#####     dbGetPreparedQuery(db, sql, bind.data=df)
#####     dbCommit(db)
##### 
#####     ## Now do the rest in batches
#####     done <- FALSE
#####     while (!done) {
#####         df <- read.table(con, sep=",", stringsAsFactors=FALSE,
#####                          nrows=batch_size, na.strings="---",
#####                          header=FALSE)[, wantedCols]
#####         if (nrow(df) < batch_size) {
#####             done <- TRUE
#####             if (nrow(df) == 0)
#####               break
#####         }
#####         names(df) <- header
#####         dbBeginTransaction(db)
#####         dbGetPreparedQuery(db, sql, bind.data=df)
#####         dbCommit(db)
#####     }
##### }
##### 
##### 
##### ################################################
##### ## affyExpr.schema.R
##### ################################################
##### 
##### BASE_A <- 1
##### BASE_C <- 2
##### BASE_G <- 3
##### BASE_T <- 4
##### 
##### createAffyExprFeatureSetSql <- ('
##### create table featureSet (
#####     fsetid integer primary key,
#####     man_fsetid text,
#####     alignment text,
#####     gene_symbol text,
#####     chrom text,
#####     ensembl text,
#####     strand integer)
##### ')
##### 
##### createAffyExprFeatureSql <- ('
##### create table %s (
#####     fid integer primary key,
#####     fsetid integer not null references "featureSet" ("fsetid"),
#####     pbase integer,
#####     tbase integer,
#####     atom integer,
#####     x integer,
#####     y integer)
##### ')
##### 
##### createAffyExprPm_MmSql <- ('
##### create table %s (
#####     pm_fid integer primary key references "pmfeature" ("fid"),
#####     mm_fid integer references "mmfeature" ("fid"))
##### ')
##### 
##### createAffyExprSequenceSql <- ('
##### create table sequence (
#####     fid integer primary key,
#####         tstrand integer,
#####     interrogation_position integer,
#####     seq text)
##### ')
##### 
##### 
##### ################################################
##### ## affyST.initDb.R
##### ################################################
##### 
##### initDb.affyST <- function(dbname) {
#####     db <- dbConnect(dbDriver("SQLite"), dbname)
##### 
#####     ## Set page size
#####     dbGetQuery(db, setPageSizeSql)
##### 
#####     ## Create tables
##### 
#####     dbGetQuery(db, createAffySTFeatureSetSql)
#####     dbGetQuery(db, sprintf(createAffySTFeatureSql, "pmfeature_tmp"))
#####     dbGetQuery(db, sprintf(createAffySTFeatureSql, "mmfeature_tmp"))
#####     dbGetQuery(db, sprintf(createAffySTPm_MmSql, "pm_mm"))
#####     
#####     
#####     dbGetQuery(db, sprintf(createAffySTFeatureSql, "qcpmfeature")) ## never used
#####     dbGetQuery(db, sprintf(createAffySTFeatureSql, "qcmmfeature")) ## never used
#####     dbGetQuery(db, sprintf(createAffySTPm_MmSql, "qcpm_qcmm")) ## never used
#####     dbGetQuery(db, createAffySTSequenceSql)
##### 
#####     ## Create index
#####     ## NOTE: might be more efficient to create this after loading,
#####     ## but current perf is ok.
##### 
#####     sql <- 'create index man_fsetid_idx on featureSet ("man_fsetid")'
#####     dbGetQuery(db, sql)
#####     db
##### }
##### 
##### sortFeatureTables.affyST <- function(db) {
#####     dbGetQuery(db, sprintf(createAffySTFeatureSql, "pmfeature"))
#####     dbGetQuery(db, sprintf(createAffySTFeatureSql, "mmfeature"))
##### 
#####     ## Reorder XXfeature tables
##### 
#####     ## on other arrays, fid is primary key (therefore, tbls ordered by
#####     ## fid, no matter what)
#####     
#####     ## fillSql <- paste("insert into %s select * from %s order by",
#####     ##                  "fsetid, atom")
#####     fillSql <- paste("insert into %s select * from %s order by fid")
#####     
#####     dbBeginTransaction(db)
#####     dbGetQuery(db, sprintf(fillSql, "pmfeature", "pmfeature_tmp"))
#####     dbGetQuery(db, sprintf(fillSql, "mmfeature", "mmfeature_tmp"))
#####     dbCommit(db)
#####     ## drop temp tables
#####     dbBeginTransaction(db)
#####     dbGetQuery(db, "drop table pmfeature_tmp")
#####     dbGetQuery(db, "drop table mmfeature_tmp")
#####     dbCommit(db)
##### }
##### 
##### 
##### createIndicesDb.affyST <- function(db) {
#####     makeIndex <- function(name, t, cols) {
#####         sql <- paste("create index", name, "on", t,
#####                      paste("(", paste(cols, collapse=","), ")"))
#####         dbGetQuery(db, sql)
#####     }
##### 
#####     ## Create DB indices and fix ordering
#####     makeIndex("pmf_idx_fsetid", "pmfeature", "fsetid")
#####     makeIndex("mmf_idx_fsetid", "mmfeature", "fsetid")
##### 
##### ##    makeIndex("fset_idx_chrom", "featureSet", "chrom")
#####     makeIndex("fset_idx_fsetid", "featureSet", "fsetid")
##### 
#####     ## finally, run analyze (SQLite specific?)
#####     dbGetQuery(db, "analyze")
##### }
##### 
##### 
##### createTableInfoTable.affyST <- function(db, verbose=FALSE) {
#####     tables <- dbListTables(db)
#####     counts <- integer(length(tables))
#####     sql <- "select count(*) from %s"
#####     for (i in seq(along=counts)) {
#####         if (verbose)
#####           cat("counting rows in ", tables[i], "\n")
#####         counts[i] <- dbGetQuery(db, sprintf(sql, tables[i]))[[1]][1]
#####     }
##### 
#####     df <- data.frame(tbl=tables, row_count=counts,
#####                      stringsAsFactors=FALSE)
#####     dbWriteTable(db, "table_info", df, row.names=FALSE)
##### }
##### 
##### 
##### createFeatureTableInfo.affyST <- function(db, tname) {
#####     return(FALSE)
#####     ## FIXME: add code to determine offsets of sorted
#####     ## strand and allele
##### }
##### 
##### 
##### ################################################
##### ## affyST.loaders.R
##### ################################################
##### 
##### .idxToIdx <- function(pgfFrom, pgfTo, ids) {
#####   starts <- pgfFrom[ids]
#####   ends <- pgfFrom[ids+1] - 1
#####   ends[is.na(ends)] <- length(pgfTo)
#####   mapply(":", starts, ends, SIMPLIFY=FALSE)
##### }
##### 
##### .combineIdx <- function(psToAtomIdx, atomToProbeIdx) {
#####   probesPerAtom <- with(atomToProbeIdx,
#####                         sapply(split(atomIdx, atomIdx), length))
#####   cbind(probesetIdx=rep(psToAtomIdx$probesetIdx, probesPerAtom),
#####         atomToProbeIdx)
##### }
##### 
##### probesetIdxToAtomIdx <- function(pgf, probesetIdx) {
#####   atoms <- .idxToIdx(pgf[["probesetStartAtom"]],
#####                      pgf[["atomId"]], probesetIdx)
#####   data.frame(probesetIdx=rep(probesetIdx, sapply(atoms, length)),
#####              atomIdx=unlist(atoms))
##### }
##### 
##### atomIdxToProbeIdx <- function(pgf, atomIdx) {
#####   probes <- .idxToIdx(pgf[["atomStartProbe"]],
#####                       pgf[["probeId"]], atomIdx)
#####   data.frame(atomIdx=rep(atomIdx, sapply(probes, length)),
#####              probeIdx=unlist(probes))
##### }
##### 
##### probesetIdxToTripletIdx <- function(pgf, probesetIdx) {
#####   df1 <- probesetIdxToAtomIdx(pgf, probesetIdx)
#####   df2 <- atomIdxToProbeIdx(pgf, df1$atomIdx)
#####   .combineIdx(df1, df2)
##### }
##### 
##### loadUnits.affyST <- function(db, pgf, clf) {
##### 
#####     ## create tables table
#####     sqliteQuickSQL(db, createAffySTFeatureSetSql)
#####     sqliteQuickSQL(db, sprintf(createAffySTFeatureSql, "pmfeature_tmp"))
#####     sqliteQuickSQL(db, sprintf(createAffySTFeatureSql, "mmfeature_tmp"))
#####     sqliteQuickSQL(db, sprintf(createAffySTPm_MmSql, "pm_mm"))
#####     
#####     ## load featureSet table
#####     fset.table <- data.frame(fsetid = pgf[["probesetId"]],
#####                              man_fsetid = pgf[["probesetName"]],
#####                              type = pgf[["probesetType"]],
#####                              start_atom = pgf[["probesetStartAtom"]],
#####                              stringsAsFactors=FALSE)
#####     i <- fset.table[["man_fsetid"]] == ""
#####     fset.table[i, "man_fsetid"] <- fset.table[i, "fsetid"]
#####     rm(i)
#####     fset.table <- fset.table[order(fset.table[["man_fsetid"]]),]
#####     
#####     
#####     dbBeginTransaction(db)
#####     values <- "(:fsetid, :man_fsetid, :type, :start_atom)"
#####     sql <- "insert into featureSet (fsetid, man_fsetid, type, start_atom) values"
#####     dbGetPreparedQuery(db, paste(sql, values), bind.data=fset.table)
#####     dbCommit(db)
##### 
#####     ## load pmfeature / mmfeature / qc*feature
#####     ## (put together the "probe df")
#####     df <- probesetIdxToTripletIdx(pgf, 1:length(pgf[["probesetId"]]))
#####     fid <- pgf[["probeId"]][df[["probeIdx"]]]
#####     i <- match(fid, pgf[["probeId"]])
#####     ii <- match(fid, clf[["id"]])
#####     probes.table <- data.frame(fid=fid,
#####                                fsetid=pgf[["probesetId"]][df[["probesetIdx"]]],
#####                                pbase=as.integer(NA),
#####                                tbase=as.integer(NA),
#####                                atom=pgf[["atomId"]][df[["atomIdx"]]],
#####                                x=clf[["x"]][ii],
#####                                y=clf[["y"]][ii],
#####                                gc_count=pgf[["probeGcCount"]][i],
#####                                type=pgf[["probeType"]][i],
#####                                strand=SENSE,
#####                                pm=as.integer(NA),
#####                                stringsAsFactors=FALSE)
#####     rm(i, ii, df)
#####     i <- grep(":at$", tolower(probes.table[["type"]]))
#####     probes.table[i, "strand"] <- ANTISENSE
#####     i <- grep("^mm:", tolower(probes.table[["type"]]))
#####     pmonly <- length(i) == 0
#####     if (!pmonly)
#####       probes.table[i, "pm"] <- as.integer(0)
#####     i <- grep("^pm:", tolower(probes.table[["type"]]))
#####     probes.table[i, "pm"] <- as.integer(1)
#####     probes.table[["type"]] <- NULL
#####     lastCol <- match("pm", names(probes.table))
#####     probes.table <- probes.table[order(probes.table[["fid"]]),]
##### 
#####     ## loading pm
#####     sql <- paste("insert into pmfeature_tmp values",
#####                  "(:fid, :fsetid, :pbase, :tbase, :atom, :x, :y, :gc_count, :strand)")
#####     dbBeginTransaction(db)
#####     dbGetPreparedQuery(db, sql, bind.data=subset(probes.table, pm == 1)[,-lastCol])
#####     dbCommit(db)
##### 
#####     ## loading mm
#####     if (!pmonly){
#####       sql <- paste("insert into mmfeature_tmp values",
#####                    "(:fid, :fsetid, :pbase, :tbase, :atom, :x, :y, :gc_count, :strand)")
#####       dbBeginTransaction(db)
#####       dbSendPreparedQuery(db, sql, subset(probes.table, pm == 0)[,-lastCol])
#####       dbCommit(db)
##### 
#####       pm.info <- probes.table[probes.table[["pm"]]==1, "fsetid"]
#####       mm.info <- probes.table[probes.table[["pm"]]==0, "fsetid"]
#####       set <- intersect(pm.info, mm.info)
#####       f.split <- function(i){
#####         tmp <- subset(probes.table[, c("fid", "atom", "pm")], probes.table$fsetid == i)
#####         pm.fid <- subset(tmp, pm==1)[,-3]
#####         mm.fid <- subset(tmp, pm==0)[,-3]
#####         common <- sort(intersect(pm.fid[["atom"]], mm.fid[["atom"]]))
#####         pm.i <- match(common, pm.fid[["atom"]])
#####         mm.i <- match(common, mm.fid[["atom"]])
#####         data.frame(pmi=pm.fid[["fid"]][pm.i], mmi=mm.fid[["fid"]][mm.i])
#####       }
#####       link <- do.call(rbind, lapply(set, f.split))
#####       
#####       ## Insert pm <--> mm link
#####       values <- "(:pmi, :mmi)"
#####       sql <- paste("insert into pm_mm values", values)
#####       dbBeginTransaction(db)
#####       rset <- dbSendPreparedQuery(db, sql, link)
#####       dbClearResult(rset)
#####       dbCommit(db)
#####     }
##### }
##### 
##### readProbeFile <- function(filename){
#####   probes <- read.table(filename, stringsAsFactors=FALSE, header=TRUE, sep="\t")
#####   header <- tolower(gsub("\\.", "_", names(probes)))
#####   names(probes) <- header
#####   expected <- c("probe_id", "transcript_cluster_id", "probe_x",
#####                 "probe_y", "assembly", "seqname", "start", "stop",
#####                 "strand", "probe_sequence", "target_strandedness",
#####                 "category")
#####   missing.fields <- expected[!(expected %in% header)]
#####   if (length(missing.fields)>0) for (i in missing.fields) probes[[i]] <- NA
#####   not.expected <- names(probes)[!names(probes) %in% expected]
#####   if (length(not.expected)>0) for (i in  not.expected) probes[[i]] <- NULL
#####   reorder <- match(expected, names(probes))
#####   probes <- probes[order(probes[["probe_id"]]),]
#####   probes[["transcript_cluster_id"]] <- as.integer(probes[["transcript_cluster_id"]])
#####   return(probes[,reorder])
##### }
##### 
##### 
##### loadAffySeqCsv.affyST <- function(db, probeFile, pgf, batch_size=5000) {
#####   
#####   ## create tables table
#####   sqliteQuickSQL(db, createAffySTSequenceSql)
#####     
#####   probeSeq <- readProbeFile(probeFile)
#####   i <- which(probeSeq[["probe_id"]] %in% pgf[["probeId"]])
#####   probeSeq[["interrogation_position"]] <- NA
#####   probeSeq[i, "interrogation_position"] <- pgf[["probeInterrogationPosition"]][match(probeSeq[["probe_id"]][i], pgf[["probeId"]])]
#####   probeSet <- probeSeq[, c("probe_id", "target_strandedness",
#####                            "interrogation_position", "transcript_cluster_id",
#####                            "seqname", "start", "stop", "probe_sequence",
#####                            "category")]
#####   ts <- ifelse(probeSet[["target_strandedness"]] == "Sense", SENSE, ANTISENSE)
#####   ts <- as.integer(ts)
#####   probeSet[["target_strandedness"]] <- ts
#####   rm(ts,i, probeSeq)
#####   probeSet <- probeSet[order(probeSet[["probe_id"]]),]
#####   
#####   sql <- paste("insert into sequence VALUES",
#####                "(:probe_id, :target_strandedness,",
#####                ":interrogation_position, :transcript_cluster_id,",
#####                ":seqname, :start, :stop, :probe_sequence, :category)")
#####   
#####   ## anything specific for MM probes?
#####   dbBeginTransaction(db)
#####   dbGetPreparedQuery(db, sql, bind.data=probeSet)
#####   dbCommit(db)
##### 
##### }
##### 
##### readTranscriptFile <- function(filename){
#####   transcript <- read.csv(filename, comment.char="#", stringsAsFactors=FALSE, na.string="---")
#####   expected <- c("transcript_cluster_id", "probeset_id", "seqname",
#####                 "strand", "start", "stop", "total_probes",
#####                 "gene_assignment", "mrna_assignment", "swissprot",
#####                 "unigene", "GO_biological_process",
#####                 "GO_cellular_component", "GO_molecular_function",
#####                 "pathway", "protein_domains", "crosshyb_type",
#####                 "category")
#####   transcript[["strand"]] <- as.integer(ifelse(transcript[["strand"]] == "+", SENSE, ANTISENSE))
#####   missing.fields <- expected[!(expected %in% names(transcript))]
#####   if (length(missing.fields)>0) for (i in missing.fields) transcript[[i]] <- NA
#####   reorder <- match(expected, names(transcript))
#####   transcript <- transcript[order(transcript[["probeset_id"]]),]
#####   transcript[["crosshyb_type"]] <- as.integer(transcript[["crosshyb_type"]])
#####   return(transcript[, reorder])
##### }
##### 
##### loadAffyCsv.affyST <- function(db, transFile, batch_size=5000) {
#####   transcript <- readTranscriptFile(transFile)
##### 
#####   db_cols <- c("transcript_cluster_id", "fsetid", "seqname",
#####                "strand", "start", "stop", "total_probes",
#####                "gene_assignment", "mrna_assignment", "swissprot",
#####                "unigene", "GO_biological_process",
#####                "GO_cellular_component", "GO_molecular_function",
#####                "pathway", "protein_domains", "crosshyb_type",
#####                "category")
##### 
#####   val_holders <- paste(":", db_cols, sep="")
#####   val_holders <- gsub("fsetid", "probeset_id", val_holders)
#####   exprs <- paste(db_cols, " = ", val_holders, sep="", collapse=", ")
#####   sql <- paste("UPDATE featureSet set", exprs, "WHERE fsetid = :probeset_id")
#####   
#####   dbBeginTransaction(db)
#####   dbGetPreparedQuery(db, sql, bind.data=transcript)
#####   dbCommit(db)
##### }
##### 
##### buildPdInfoDb.affyST <- function(pgfFile, clfFile, probeFile, transFile, dbFile, matFile,
#####                           batch_size=10000, verbose=FALSE) {
##### 
#####     cat("starting Build process\n")
#####     ST <- system.time
#####     printTime <- function(msg, t) {
#####         if (verbose) {
#####             m <- paste(msg, "took %.2f sec\n")
#####             cat(sprintf(m, t))
#####         }
#####     }
#####     cat("ReadClfEnv\n")
#####     clf <- readClfEnv(clfFile)
#####     cat("ReadpgfEnv\n")
#####     pgf <- readPgfEnv(pgfFile)
##### 
#####     cat("initDB\n")
#####     db <- initDb.affyST(dbFile)
#####     cat("loadUnits\n") 
#####     t <- ST(loadUnits.affyST(db, pgf, clf))
#####     printTime("loadUnitsByBatch", t[3])
#####     cat("lodAffyCsv\n")
#####     t <- ST(loadAffyCsv.affyST(db, transFile, batch_size=batch_size))
#####     printTime("loadAffyCsv", t[3])
#####     cat("loadAffySeq\n")
#####     t <- ST(loadAffySeqCsv.affyST(db, probeFile, pgf, batch_size=batch_size))
#####     printTime("loadAffySeqCsv", t[3])
#####     cat("Sort\n")
#####     t <- ST({
#####         sortFeatureTables.affyST(db)
#####         createIndicesDb.affyST(db)
#####         createTableInfoTable.affyST(db)
#####         createFeatureTableInfo.affyST(db)
#####     })
#####     printTime("DB sort, index creation", t[3])
##### 
##### ##     t <- ST({
##### ##         seqMat <- createSeqMat(db)
##### ##         save(seqMat, file=matFile, compress=TRUE)
##### ##     })
##### ##     printTime("sequence matrix", t[3])
#####     closeDb(db)
##### }
##### 
##### 
##### ################################################
##### ## affyST.schema.R
##### ################################################
##### 
##### ## This is the schema file
##### ## for the Affy Gene ST arrays
##### ## Started: March/08 - Benilton Carvalho
##### 
##### ## Strategy:
##### ## featureSet table contains only *required* fields from PGF
##### ## sequence table contains info from probe.tab
##### ## likely to have only pmfeature table.
##### 
##### createAffySTFeatureSetSql <- paste('
##### create table featureSet (
#####     fsetid integer primary key,
#####     man_fsetid text,
#####     type text,
#####     start_atom integer,
#####     transcript_cluster_id integer,
#####     seqname text,
#####     strand integer,
#####     start integer,
#####     stop integer,
#####     total_probes integer,
#####     gene_assignment text,
#####     mrna_assignment text,
#####     swissprot text,
#####     unigene text,
#####     GO_biological_process text,
#####     GO_cellular_component text,
#####     GO_molecular_function text,
#####     pathway text,
#####     protein_domains text,
#####     crosshyb_type integer,
#####     category text
##### )
##### ')
##### 
##### createAffySTFeatureSql <- paste('
##### create table %s (
#####     fid integer,
#####     fsetid integer not null references "featureSet" ("fsetid"),
#####     pbase integer,
#####     tbase integer,
#####     atom integer,
#####     x integer,
#####     y integer,
#####     gc_count integer,
#####     strand integer,
#####     PRIMARY KEY (atom, fid))
##### ')
##### 
##### createAffySTPm_MmSql <- paste('
##### create table %s (
#####     pm_fid integer primary key references "pmfeature" ("fid"),
#####     mm_fid integer references "mmfeature" ("fid"))
##### ')
##### 
##### 
##### ##  info here could actually on on the *feature tables but on ST arrays,
##### ##  probe ID (fid) are not necessarily unique to a probeset (fsetid) and
##### ##  many of the things would be duplicated (but it is still duplicated
##### ##  here with transcript_cluster and trancript+probeid do not form a key
##### ##  =(
##### 
##### createAffySTSequenceSql <- paste('
##### create table sequence (
#####     fid integer,
#####     tstrand integer,
#####     interrogation_position integer,
#####     transcript_cluster integer,
#####     seqname text,
#####     start integer,
#####     stop integer,
#####     seq text,
#####     category text)
##### ')
##### 
##### 
##### ################################################
##### ## affyTiling.initDb.R
##### ################################################
##### 
##### initDb.affyTiling <- function(dbname) {
#####     db <- dbConnect(dbDriver("SQLite"), dbname)
##### 
#####     ## Set page size
#####     dbGetQuery(db, setPageSizeSql)
##### 
#####     ## Create tables
#####     ## BC: Soon we need to add a table for the control probes
#####     dbGetQuery(db, createAffyTilingFeatureSetSql)
##### 
#####     dbGetQuery(db, sprintf(createAffyTilingFeatureSql, "pmfeature_tmp"))
##### 
#####     dbGetQuery(db, sprintf(createAffyTilingFeatureSql, "mmfeature_tmp"))
##### 
#####     dbGetQuery(db, sprintf(createAffyTilingFeatureSql, "qcpmfeature"))
#####     dbGetQuery(db, sprintf(createAffyTilingFeatureSql, "qcmmfeature"))
##### 
#####     dbGetQuery(db, sprintf(createAffyTilingPm_MmSql, "pm_mm"))
#####     dbGetQuery(db, sprintf(createAffyTilingPm_MmSql, "qcpm_qcmm"))
##### 
#####     dbGetQuery(db, createAffyTilingSequenceSql)
##### 
#####     ## Create index
#####     ## NOTE: might be more efficient to create this after loading,
#####     ## but current perf is ok.
#####     sql <- 'create index man_fsetid_idx on featureSet ("man_fsetid")'
#####     dbGetQuery(db, sql)
#####     db
##### }
##### 
##### sortFeatureTables.affyTiling <- function(db) {
#####     dbGetQuery(db, sprintf(createAffyTilingFeatureSql, "pmfeature"))
#####     dbGetQuery(db, sprintf(createAffyTilingFeatureSql, "mmfeature"))
##### 
#####     ## Reorder XXfeature tables
#####     fillSql <- paste("insert into %s select * from %s order by",
#####                      "fsetid, startpos")
#####     dbBeginTransaction(db)
#####     dbGetQuery(db, sprintf(fillSql, "pmfeature", "pmfeature_tmp"))
#####     dbGetQuery(db, sprintf(fillSql, "mmfeature", "mmfeature_tmp"))
#####     dbCommit(db)
#####     ## drop temp tables
#####     dbBeginTransaction(db)
#####     dbGetQuery(db, "drop table pmfeature_tmp")
#####     dbGetQuery(db, "drop table mmfeature_tmp")
#####     dbCommit(db)
##### }
##### 
##### 
##### createIndicesDb.affyTiling <- function(db) {
#####     makeIndex <- function(name, t, cols) {
#####         sql <- paste("create index", name, "on", t,
#####                      paste("(", paste(cols, collapse=","), ")"))
#####         dbGetQuery(db, sql)
#####     }
##### 
#####     ## Create DB indices and fix ordering
#####     makeIndex("pmf_idx_fsetid", "pmfeature", "fsetid")
#####     makeIndex("mmf_idx_fsetid", "mmfeature", "fsetid")
##### 
##### ##    makeIndex("fset_idx_chrom", "featureSet", "chrom")
#####     makeIndex("fset_idx_fsetid", "featureSet", "fsetid")
##### 
#####     ## finally, run analyze (SQLite specific?)
#####     dbGetQuery(db, "analyze")
##### }
##### 
##### 
##### createTableInfoTable.affyTiling <- function(db, verbose=FALSE) {
#####     tables <- dbListTables(db)
#####     counts <- integer(length(tables))
#####     sql <- "select count(*) from %s"
#####     for (i in seq(along=counts)) {
#####         if (verbose)
#####           cat("counting rows in ", tables[i], "\n")
#####         counts[i] <- dbGetQuery(db, sprintf(sql, tables[i]))[[1]][1]
#####     }
##### 
#####     df <- data.frame(tbl=tables, row_count=counts,
#####                      stringsAsFactors=FALSE)
#####     dbWriteTable(db, "table_info", df, row.names=FALSE)
##### }
##### 
##### 
##### createFeatureTableInfo.affyTiling <- function(db, tname) {
#####     return(FALSE)
#####     ## FIXME: add code to determine offsets of sorted
#####     ## strand and allele
##### }
##### 
##### 
##### 
##### ################################################
##### ## affyTiling.loaders.R
##### ################################################
##### 
##### loadUnitNames.affyTiling <- function(db, unames) {
#####     dbBeginTransaction(db)
#####     ## To use an auto-incrementing field via RSQLite, you need
#####     ## to be careful to pass integer NA's
#####     df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
#####     values <- "(:id, :name)"
#####     sql <- "insert into featureSet (fsetid, man_fsetid) values"
#####     dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
#####     dbCommit(db)
##### }
##### 
##### loadUnits.affyTiling <- function(db, batch, nx, isQc=FALSE) {
#####     pmfeature <- "pmfeature_tmp"
#####     mmfeature <- "mmfeature_tmp"
#####     pmmm <- "pm_mm"
#####     if (isQc) {
#####         pmfeature <- "qcpmfeature"
#####         mmfeature <- "qcmmfeature"
#####         pmmm <- "qcpm_qcmm"
#####     }
##### 
#####     batchMat <- do.call(rbind, lapply(batch, readCdfUnitToMat.affyTiling, nx=nx))
##### 
#####     sets <- as.character(unlist(sapply(batch, function(x) paste(x$seqInfo$fullname, x$startpos, sep=":"))))
##### 
#####     tmp.df <- function(x){
#####       tmp <- data.frame(id=as.integer(NA),
#####                         man_fsetid=paste(x$seqInfo$fullname, x$startpos, sep=":"),
#####                         groupname=x$seqInfo$groupname,
#####                         version=x$seqInfo$version,
#####                         fullname=x$seqInfo$fullname,
#####                         name=x$seqInfo$name,
#####                         stringsAsFactors=FALSE)
#####       tmp[!duplicated(tmp[["man_fsetid"]]), ]
#####     }
##### 
#####     super.set <- do.call(rbind, lapply(batch, tmp.df))
#####     values <- "(:id, :man_fsetid, :groupname, :version, :fullname, :name)"
#####     sql <- "INSERT INTO featureSet VALUES"
#####     dbBeginTransaction(db)
#####     dbGetPreparedQuery(db, paste(sql, values), bind.data=super.set)
#####     dbCommit(db)
#####     rm(super.set)
##### 
#####     ## Find internal featureSet IDs for these features
#####     ## Using subsets of 1000 units
#####     tmp <- NULL
#####     us <- unique(sets)
#####     subsets <- rep(1:length(us), length.out=length(us), each=1000)
#####     subgrps <- split(us, subsets)
#####     for (i in 1:length(subgrps)){
#####       sql <- paste("select man_fsetid, fsetid from featureSet where man_fsetid in (",
#####                    paste('"', subgrps[[i]], '"', sep="", collapse=","),
#####                  ") order by fsetid")
#####       tmp <- rbind(tmp, dbGetQuery(db, sql))
#####     }
#####     
#####     batchIds <- tmp
#####     idx <- match(sets, batchIds$man_fsetid)
#####     ids <- batchIds$fsetid[idx]
#####     pmmmtest <- any(batchMat[, "ispm"] == 0)
#####     if (pmmmtest){
#####       batchMat <- rbind(cbind(batchMat[batchMat[,"ispm"]==1,], fsetid=ids),
#####                         cbind(batchMat[batchMat[,"ispm"]==0,], fsetid=ids))
#####     }else{
#####       batchMat <- cbind(batchMat[batchMat[,"ispm"]==1,], fsetid=ids)
#####     }
#####     theDups <- duplicated(batchMat[["indices"]])
#####     batchMat <- batchMat[!theDups,]
#####     rm(theDups)
#####     
#####     ## Insert pm and mm into respective tables
#####     isPm <- as.logical(batchMat[, "ispm"]) 
#####     values <- "(:indices, :fsetid, :strand, :startpos,  :atom, :x, :y)"
#####     sql <- paste("insert into", pmfeature, "values", values)
#####     dbBeginTransaction(db)
#####     rset <- dbSendPreparedQuery(db, sql, batchMat[isPm, ])
#####     dbClearResult(rset)
##### 
#####     if (pmmmtest){
#####       sql <- paste("insert into", mmfeature, "values", values)
#####       rset <- dbSendPreparedQuery(db, sql, batchMat[!isPm, ])
#####       dbClearResult(rset)
#####     }
#####     dbCommit(db)
##### 
#####     ## Insert pm <--> mm link
#####     if (pmmmtest){
#####       pm.fsetid.atom.set = subset(batchMat, ispm == 1)[, c("fsetid", "atom", "indices")]
#####       pm.fsetid.atom = paste(pm.fsetid.atom.set[["fsetid"]], pm.fsetid.atom.set[["atom"]], sep=":")
#####       mm.fsetid.atom.set = subset(batchMat, ispm == 0)[, c("fsetid", "atom", "indices")]
#####       mm.fsetid.atom = paste(mm.fsetid.atom.set[["fsetid"]], mm.fsetid.atom.set[["atom"]], sep=":")
#####       common <- sort(intersect(pm.fsetid.atom, mm.fsetid.atom))
##### 
#####       pmset <- match(common, pm.fsetid.atom)
#####       mmset <- match(common, mm.fsetid.atom)
#####     
#####       values <- "(:pmi, :mmi)"
#####       sql <- paste("insert into", pmmm, "values", values)
#####       dbBeginTransaction(db)
#####     
#####       rset <- dbSendPreparedQuery(db, sql,
#####                                   data.frame(pmi=pm.fsetid.atom.set[pmset, "indices"],
#####                                              mmi=mm.fsetid.atom.set[mmset, "indices"]))
#####       dbClearResult(rset)
#####       dbCommit(db)
#####     }
##### 
#####     ## Insert sequence info
#####     values <- "(:indices, :probeseq)"
#####     sql <- paste("INSERT INTO sequence VALUES", values)
#####     dbBeginTransaction(db)
#####     rset <- dbSendPreparedQuery(db, sql, batchMat)
#####     dbClearResult(rset)
#####     dbCommit(db)
##### }
##### 
##### 
##### loadUnitsByBatch.affyTiling <- function(db, bpmapFile, batch_size=1, nx,
#####                              max_units=NULL, verbose=FALSE) {
#####     seqs <- readBpmapSeqinfo(bpmapFile)
#####     grps <- sapply(seqs, "[[", "parameters")
#####     whQc <- which(tolower(grps) != "tiling")
##### 
#####     offset <- 1
#####     if (length(whQc)) {                 # load all QC at once
#####       qcunits <- readBpmap(bpmapFile, whQc)
#####       loadUnits.affyTiling(db, qcunits, nx=nx, isQc=TRUE)
#####       offset <- max(whQc) + 1
#####     }
##### 
#####     for (wanted in offset:length(grps)){
#####       message("Processing unit ", wanted, " out of ", length(grps), ".")
#####       vvunits <- readBpmap(bpmapFile, wanted)
#####       loadUnits.affyTiling(db, vvunits, nx=nx)
#####     }
#####     
##### ##     if (is.null(max_units))
##### ##       max_units <- length(grps)
##### ## 
##### ##     extra <- (max_units-length(whQc)) %% batch_size
##### ##     num_batches <- (max_units-length(whQc)) %/% batch_size
##### ##     if (extra != 0)
##### ##       num_batches <- num_batches + 1
##### ##     done <- FALSE
##### ##     while (!done) {
##### ##         end <- min(offset + batch_size, max_units)
##### ##         if (end == max_units)
##### ##           done <- TRUE
##### ##         wanted <- seq.int(offset, end)
##### ##         offset <- offset + batch_size + 1
##### ##         vvunits <- readBpmap(bpmapFile, wanted)
##### ##         loadUnits.affyTiling(db, vvunits, nx=nx)
##### ##     }
##### }
##### 
##### 
##### buildPdInfoDb.affyTiling <- function(bpmapFile,  celFile, dbFile, matFile,
#####                           batch_size=10000, verbose=FALSE) {
##### 
#####     ST <- system.time
#####     printTime <- function(msg, t) {
#####         if (verbose) {
#####             m <- paste(msg, "took %.2f sec\n")
#####             cat(sprintf(m, t))
#####         }
#####     }
##### 
#####     db <- initDb.affyTiling(dbFile)
##### 
#####     celh <- readCelHeader(celFile)
#####     nx <- as.integer(celh$cols)
#####     rm(celh)
#####     
#####     t <- ST(loadUnitsByBatch.affyTiling(db, bpmapFile, batch_size=batch_size, nx=nx))
#####     printTime("loadUnitsByBatch", t[3])
#####     t <- ST({
#####         sortFeatureTables.affyTiling(db)
#####         createIndicesDb.affyTiling(db)
#####         createTableInfoTable.affyTiling(db)
#####         createFeatureTableInfo.affyTiling(db)
#####     })
#####     printTime("DB sort, index creation", t[3])
##### 
##### ##     t <- ST({
##### ##         seqMat <- createSeqMat(db)
##### ##         save(seqMat, file=matFile, compress=TRUE)
##### ##     })
##### ##     printTime("sequence matrix", t[3])
#####     closeDb(db)
##### }
##### 
##### 
##### 
##### ################################################
##### ## affyTiling.schema.R
##### ################################################
##### 
##### createAffyTilingFeatureSetSql <- ('
##### create table featureSet (
#####     fsetid integer primary key,
#####     man_fsetid text,
#####     groupname text,
#####     version text,
#####     fullname text,
#####     name text)
##### ')
##### 
##### createAffyTilingFeatureSql <- ('
##### create table %s (
#####     fid integer primary key,
#####     fsetid integer not null references "featureSet" ("fsetid"),
#####     strand integer,
#####     startpos integer,
#####     atom integer,
#####     x integer,
#####     y integer)
##### ')
##### 
##### createAffyTilingPm_MmSql <- ('
##### create table %s (
#####     pm_fid integer primary key references "pmfeature" ("fid"),
#####     mm_fid integer references "mmfeature" ("fid"))
##### ')
##### 
##### createAffyTilingSequenceSql <- ('
##### create table sequence (
#####     fid integer primary key,
#####     seq text)
##### ')
##### 
##### 
##### ################################################
##### ## ngs.InitDb.R
##### ################################################
##### 
##### "initDb.ngs" <- 
##### function(dbname) {
#####     db <- dbConnect(dbDriver("SQLite"), dbname)
#####     
#####     ## Set page size
#####     dbGetQuery(db, setPageSizeSql)
#####     
#####     ## Create tables
#####     dbGetQuery(db, createNgsFeatureSetSql)
##### 
#####     dbGetQuery(db, sprintf(createNgsFeatureSql, "pmfeature_tmp"))
#####     dbGetQuery(db, sprintf(createNgsFeatureSql, "mmfeature_tmp"))
#####     dbGetQuery(db, sprintf(createNgsFeatureSql, "qcpmfeature"))
#####     dbGetQuery(db, sprintf(createNgsFeatureSql, "qcmmfeature"))
##### 
#####     dbGetQuery(db, sprintf(createNgsPm_MmSql, "pm_mm"))
#####     dbGetQuery(db, sprintf(createNgsPm_MmSql, "qcpm_qcmm"))
##### 
#####     dbGetQuery(db, createNgsSequenceSql)
#####     
#####     ## Create index
#####     ## NOTE: might be more efficient to create this after loading,
#####     ## but current perf is ok.
#####     sql <- 'create index man_fsetid_idx on featureSet ("man_fsetid")'
#####     dbGetQuery(db, sql)
#####     db
##### }
##### 
##### "sortFeatureTables.ngs" <- 
##### function(db) {
#####     dbGetQuery(db, sprintf(createNgsFeatureSql, "pmfeature"))
#####     dbGetQuery(db, sprintf(createNgsFeatureSql, "mmfeature"))
#####     
#####     ## Reorder XXfeature tables
#####     fillSql <- paste("insert into %s select * from %s order by",
#####         "fsetid, atom")
#####     dbBeginTransaction(db)
#####     dbGetQuery(db, sprintf(fillSql, "pmfeature", "pmfeature_tmp"))
#####     dbGetQuery(db, sprintf(fillSql, "mmfeature", "mmfeature_tmp"))
#####     dbCommit(db)
#####     ## drop temp tables
#####     dbBeginTransaction(db)
#####     dbGetQuery(db, "drop table pmfeature_tmp")
#####     dbGetQuery(db, "drop table mmfeature_tmp")
#####     dbCommit(db)
##### }
##### 
##### 
##### "createIndicesDb.ngs" <- 
##### function(db) {
#####     makeIndex.ngs <- function(name, t, cols) {
#####         sql <- paste("create index", name, "on", t,
#####             paste("(", paste(cols, collapse=","), ")"))
#####         dbGetQuery(db, sql)
#####     }
#####     
#####     ## Create DB indices and fix ordering
#####     makeIndex.ngs("pmf_idx_fsetid", "pmfeature", "fsetid")
#####     makeIndex.ngs("mmf_idx_fsetid", "mmfeature", "fsetid")
#####     
#####     ## makeIndex("fset_idx_chrom", "featureSet", "chrom")
#####     makeIndex.ngs("fset_idx_fsetid", "featureSet", "fsetid")
#####     
#####     ## finally, run analyze (SQLite specific?)
#####     dbGetQuery(db, "analyze")
##### }
##### 
##### 
##### "createTableInfoTable.ngs" <- 
##### function(db, verbose=FALSE) {
#####     tables <- dbListTables(db)
#####     counts <- integer(length(tables))
#####     sql <- "select count(*) from %s"
#####     for (i in seq(along=counts)) {
#####         if (verbose)
#####             cat("counting rows in ", tables[i], "\n")
#####         counts[i] <- dbGetQuery(db, sprintf(sql, tables[i]))[[1]][1]
#####     }
#####     
#####     df <- data.frame(tbl=tables, row_count=counts,
#####     stringsAsFactors=FALSE)
#####     dbWriteTable(db, "table_info", df, row.names=FALSE)
##### }
##### 
##### 
##### "createFeatureTableInfo.ngs" <- 
##### function(db, tname) {
#####     return(FALSE)
#####     ## FIXME: add code to determine offsets of sorted
#####     ## strand and allele
##### }
##### 
##### 
##### ################################################
##### ## ngs.loader.R
##### ################################################
##### 
##### loadUnitNames.ngs <- function(db, unames) {
#####   dbBeginTransaction(db)
#####   ## To use an auto-incrementing field via RSQLite, you need
#####   ## to be careful to pass integer NA's
#####   df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
#####   values <- "(:id, :name)"
#####   sql <- "insert into featureSet (fsetid, man_fsetid) values"
#####   dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
#####   dbCommit(db)
##### }
##### 
##### 
##### ## URGENT
##### loadUnits.ngs <- function(db, batch, isQc=FALSE) {
##### ## debug
##### #batch = ndfdata[-controls,]
##### #isQc=FALSE
##### ## create a couple of MISMATCH Probes for testing
##### #batch[c(4,5,6),"MISMATCH"] <- 1
##### #batch[1:6,"MATCH_INDEX"] <- rep(1:3,2)
##### ## end debug
##### 
#####     pmfeature <- "pmfeature_tmp"
#####     mmfeature <- "mmfeature_tmp"
#####     pmmm      <- "pm_mm"
#####     if (isQc) {
#####         pmfeature <- "qcpmfeature"
#####         mmfeature <- "qcmmfeature"
#####         pmmm      <- "qcpm_qcmm"
#####     }
##### 
#####     ## fill featureSet object with SEQ_ID
#####     loadUnitNames.ngs(db, unique(batch[["SEQ_ID"]]))
#####    
#####     batch <- batch[order(batch$CONTAINER, batch$SEQ_ID, batch$POSITION),]
#####     pcounts <- table(paste(batch$CONTAINER, batch$SEQ_ID, sep="."))
#####     batch$ATOM <- unlist(sapply(pcounts, function(x) 1:x))
#####     batch <- batch[order(batch$fid),]
##### 
#####     batchMat <- batch[, c("fid", "CONTAINER", "FEATURE_ID", "MISMATCH", "MATCH_INDEX", "ATOM", "PROBE_ID", "X", "Y", "POSITION")]
#####     ## Find internal featureSet IDs for these feature
#####     fset <- dbGetQuery(db, "SELECT fsetid, man_fsetid FROM featureSet")
#####     batchIds <- fset[match(batch[["SEQ_ID"]], fset[["man_fsetid"]]), "fsetid"]
#####     batchMat <- cbind(batchMat, fsetid=batchIds)
#####     
#####     ## Insert pm and mm into respective tables
#####     isPm <- batchMat[["MISMATCH"]] == 0
#####     isMm <- batchMat[["MISMATCH"]] > 0 & batchMat[["MISMATCH"]] < 10000 ## NGS can have mismatch values > 10000; May be overridden by users if MISMATCH in probe file >= 10,000.
#####     
#####     batchMat[["MISMATCH"]] <- NULL
#####     names(batchMat) <- c("fid", "container", "unit_id", "match_index", "atom", "probe_id", "x", "y", "position", "fsetid")
##### 
#####     ## insert Pm value
#####     values <- "(:fid, :container, :unit_id, :match_index, :atom, :probe_id, :x, :y, :position, :fsetid)"
#####     sql <- paste("insert into", pmfeature, "values", values)
#####     dbBeginTransaction(db)
#####     rset <- dbSendPreparedQuery(db, sql, batchMat[isPm, ])
#####     dbClearResult(rset)
#####     
#####     if (sum(isMm) > 0){
#####         sql <- paste("insert into", mmfeature, "values", values)
#####         rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[isMm, ]))
#####         dbClearResult(rset)
#####     }
#####     dbCommit(db)
#####     
#####     ## Insert pm <--> mm link
#####     ## not necessarily exists for NGS (or any other PM-only arrays) 
#####     ## NGS can have MM, but usually not
#####     ## MAKE ME MORE ROBUST!!!!
#####     int <- intersect(batchMat[isPm, "match_index"], batchMat[isMm, "match_index"])
#####     if (length(int) > 0){
#####         values <- "(:pmi, :mmi)"
#####         sql <- paste("insert into", pmmm, "values", values)
#####         df1 <- subset(batchMat[isPm,],  batchMat[isPm,"match_index"] %in% int)
#####         df2 <- subset(batchMat[isMm,],  batchMat[isPm,"match_index"] %in% int)
#####         dbBeginTransaction(db)
#####         rset <- dbSendPreparedQuery(db, sql,
#####             data.frame(pmi=df1[["fid"]],mmi=df2[["fid"]]))
#####         dbClearResult(rset)
#####         dbCommit(db)
#####     }
##### }
##### 
##### loadUnitsByBatch.ngs <- function(db, ndfdata, batch_size=10000,
#####                              max_units=NULL, verbose=FALSE) {
##### ## debug
##### #max_units=NULL
##### #verbose=FALSE
##### ## end debug
##### 
#####     ctrl.containers <- c("CODE", "RANDOM", "CONTROL")
#####     controls <- sort(as.integer(unlist(sapply(ctrl.containers, grep, toupper(ndfdata[["CONTAINER"]])))))
#####     ## fill qcpmfeature, qcmmfeature, qcpm_qcmm tables    
#####     loadUnits.ngs(db, ndfdata[controls,], isQc=TRUE)
#####     ## fill pmfeature, mmfeature, pm_qcmm tables    
#####     loadUnits.ngs(db, ndfdata[-controls,], isQc=FALSE)
#####     ## fill sequences table    
#####     seq.df <- data.frame(fid=ndfdata$fid, seq=ndfdata$PROBE_SEQUENCE, interrogation_position=ndfdata$POSITION)
#####     sql <- "INSERT INTO sequence (fid, seq, interrogation_position) VALUES (:fid, :seq, :interrogation_position)"
#####     dbBeginTransaction(db)
#####     dbGetPreparedQuery(db, sql, bind.data=seq.df)
#####     dbCommit(db)
##### }
##### 
##### loadNgd.ngs <- function(db, ngdfile){
#####   ## I had to fix ngdfile, there as a bad character somewhere and the whole table wasn't being loaded, specified comment column as character
#####   ngddata <- read.delim(ngdfile, sep="\t", as.is=TRUE, header=TRUE)
#####   ngddata <- as.data.frame(ngddata)
#####   colnames(ngddata) <- c("man_fsetid","comment")
#####   
#####   fset <- dbGetQuery(db, "SELECT fsetid, man_fsetid FROM featureSet")
#####   ngdIds <- match(fset[["man_fsetid"]],ngddata[["man_fsetid"]])
#####   ngddata <- cbind(fset,comment=ngddata[ngdIds,"comment"])
#####   
#####   db_cols <- c("comment")
#####   val_holders <- c(":comment")
#####   exprs <- paste(db_cols, " = ", val_holders, sep="", collapse=", ")
#####   sql <- paste("update featureSet set ", exprs,
#####                "where man_fsetid = :man_fsetid")
#####   
#####   dbBeginTransaction(db)
#####   dbGetPreparedQuery(db, sql, bind.data=ngddata)
#####   dbCommit(db)
##### }
##### 
##### ## TODO: fix for tiling arrays
##### ## DOING
##### loadPos.ngs <- function(db, posfile){
#####   posdataFull <- read.delim(posfile, as.is=TRUE, header=TRUE)
#####   posdata <- posdataFull[, c("SEQ_ID", "PROBE_ID", "POSITION")]
#####   posdata[["id"]] <- paste(posdata[["SEQ_ID"]], posdata[["PROBE_ID"]], sep=":::")
#####   posdata <- posdata[order(posdata[["id"]]),]
#####   df1 <- dbGetQuery(db, "SELECT fid, man_fsetid, probe_id FROM pmfeature_tmp, featureSet WHERE pmfeature_tmp.fsetid = featureSet.fsetid")
#####   df1[["id"]] <- paste(df1[["man_fsetid"]], df1[["probe_id"]], sep=":::")
#####   idx <- match(posdata[["id"]], df1[["id"]])
#####   df1 <- df1[idx,]
#####   tmp1 <- data.frame(fid=as.integer(df1$fid), position=as.integer(posdata$POSITION))
#####   sql <- "UPDATE pmfeature_tmp SET position = :position WHERE fid = :fid"
#####   dbBeginTransaction(db)
#####   dbGetPreparedQuery(db, sql, bind.data=tmp1)
#####   dbCommit(db)
#####     
#####   df1 <- dbGetQuery(db, "SELECT fid, man_fsetid, probe_id FROM mmfeature_tmp, featureSet WHERE mmfeature_tmp.fsetid = featureSet.fsetid")
#####   if (nrow(df1) > 0){
#####     df1[["id"]] <- paste(df1[["man_fsetid"]], df1[["probe_id"]], sep=":::")
#####     idx <- match(posdata[["id"]], df1[["id"]])
#####     df1 <- df1[idx,]
#####     tmp1 <- data.frame(fid=as.integer(df1$fid), position=as.integer(posdata$POSITION))
#####     sql <- "UPDATE mmfeature_tmp SET position = :position WHERE fid = :fid"
#####     dbBeginTransaction(db)
#####     dbGetPreparedQuery(db, sql, bind.data=tmp1)
#####     dbCommit(db)
#####   }
#####   
#####   ## update featureSet with chrom
#####   print("CHROMOSOME" %in% names(posdataFull))
#####   dup <- duplicated(posdataFull[["SEQ_ID"]])
#####   posdata <- posdataFull[!dup, c("SEQ_ID", "CHROMOSOME")]
#####   names(posdata) <- c("man_fsetid", "chrom")
#####   sql <- "UPDATE featureSet SET chrom = :chrom WHERE man_fsetid = :man_fsetid"
#####   dbBeginTransaction(db)
#####   dbGetPreparedQuery(db, sql, bind.data=posdata)
#####   dbCommit(db)
##### }
##### 
##### 
##### ## buildPdInfoDb.ngs <- function(ndfFile, posFile, csvSeqFile, dbFile, matFile,
##### ##                           batch_size=10000, verbose=FALSE) {
##### "buildPdInfoDb.ngsExprs" <- 
##### function(ndfdata, ngdFile, dbFile,
#####                           batch_size=10000, verbose=FALSE) {
##### ## debug
##### #ndfdata
##### #ngdFile = object@ngdFile
##### #dbFile = dbFilePath 
##### #batch_size=batch_size
##### #verbose=!quiet
##### ## end debug
#####     ST <- system.time
#####     printTime <- function(msg, t) {
#####         if (verbose) {
#####             m <- paste(msg, "took %.2f sec\n")
#####             cat(sprintf(m, t))
#####         }
#####     }
#####   
#####     # Initialize DB, create tables
#####     db <- initDb.ngs(dbFile)
#####     # Load primary feature information and sequence data    
#####     t <- ST(loadUnitsByBatch.ngs(db, ndfdata, batch_size=batch_size))
#####     printTime("loadUnitsByBatch", t[3])
#####     # load 
#####     if (file.exists(ngdFile)){
#####         t <- ST(loadNgd.ngs(db, ngdFile))
#####         printTime("loadNgd.ngs", t[3])
#####     }
#####     t <- ST({
#####         sortFeatureTables.ngs(db)
#####         createIndicesDb.ngs(db)
#####         createTableInfoTable.ngs(db)
#####         createFeatureTableInfo.ngs(db)
#####     })
#####     printTime("DB sort, index creation", t[3])
#####     closeDb(db)
##### }
##### 
##### ## TODO: Fix for tiling arrays
##### ## DOING
##### buildPdInfoDb.ngsTiling <- function(ndfdata, posFile, dbFile,
#####                                     batch_size=10000, verbose=FALSE) {
#####     ST <- system.time
#####     printTime <- function(msg, t) {
#####       if (verbose) {
#####         m <- paste(msg, "took %.2f sec\n")
#####         cat(sprintf(m, t))
#####       }
#####     }
#####    
#####     db <- initDb.ngs(dbFile)
#####     t <- ST(loadUnitsByBatch.ngs(db, ndfdata, batch_size=batch_size))
#####     printTime("loadUnitsByBatch", t[3])
#####     t <- ST(loadPos.ngs(db, posFile))
#####     printTime("loadAffySeqCsv", t[3])
#####     t <- ST({
#####         sortFeatureTables.ngs(db)
#####         createIndicesDb.ngs(db)
#####         createTableInfoTable.ngs(db)
#####         createFeatureTableInfo.ngs(db)
#####     })
#####     printTime("DB sort, index creation", t[3])
#####     closeDb(db)
##### }
##### 
##### 
##### ################################################
##### ## ngs.schema.R
##### ################################################
##### 
##### ## incorporate blocks
##### ## feature set table
##### createNgsFeatureSetSql <- ('
##### create table featureSet (
#####     fsetid integer primary key,
#####     man_fsetid text,
#####     alignment text,
#####     gene_symbol text,
#####     chrom text,
#####     ensembl text,
#####     strand integer,
#####     comment text)
##### ')
##### 
##### ## primary tables for pm and mm information
##### createNgsFeatureSql <- ('
##### create table %s (
#####     fid integer primary key,
#####     container text,
#####     unit_id integer,
#####     match_index integer,
#####     atom integer,
#####     probe_id text,
#####     x integer,
#####     y integer,
#####     position integer,
#####     fsetid integer not null references "featureSet" ("fsetid"))
##### ')
##### 
##### ##
##### createNgsPm_MmSql <- ('
##### create table %s (
#####     pm_fid integer primary key references "pmfeature" ("fid"),
#####     mm_fid integer references "mmfeature" ("fid"))
##### ')
##### 
##### ## create table sequence
##### createNgsSequenceSql <- ('
##### create table sequence (
#####     fid integer primary key,
#####     interrogation_position integer,
#####     seq text)
##### ')
##### 
##### 
##### 
##### ################################################
##### ## readNimblegen.R
##### ################################################
##### 
##### ## modified by Matt Settles June 2,2008
##### 
##### ## Nimblegen Header
##### ## (s)	"software": 
##### ## (s)	"version":
##### ## (s)   "imagefile":
##### ## (s)	"designfile":
##### ## (s)	"designname":
##### ##	(s)	"designid":
##### ## (d)	"date":
##### ## (i)	"border":
##### ## (i)	"ul_x":
##### ## (i)	"ul_y":
##### ## (i)   "ur_x":
##### ## (i)	"ur_y":
##### ## (i)	"lr_x":
##### ## (i)	"lr_y":
##### ## (i)	"ll_x":
##### ## (i)   "ll_y":
##### ## (i)	"score":
##### ## (i)	"qcscore":
##### ##	(b)	"locallyaligned":
##### ## (b)	"correctAstig":
##### ## (?)	"Knots":
##### ## (b)	"auto":
##### 
##### ## readPairHeader works with xys files or pair files
##### readPairHeader <- function (file){
#####   firstfield <- scan(file, what = "", sep = "\t", quote = "\"", #"
#####                      nlines = 100, flush = TRUE, quiet = TRUE, blank.lines.skip = FALSE,
#####                      multi.line = FALSE, allowEscape = FALSE)
#####   NHeaderRecords <- grep("# software=", firstfield)
#####   txt <- scan(file, what = "", sep = "\t", quote = "\"", nlines = NHeaderRecords, #"
#####               quiet = TRUE, allowEscape = FALSE)
#####   substring(text=txt[1],first=3)
#####   txt <- strsplit(txt,split="=")
#####     txt <- data.frame(Value=sapply(txt,function(x) x[2]),row.names=sapply(txt,function(x) x[1]),stringsAsFactors=FALSE)
#####   out <- list(NHeaderRecords = NHeaderRecords, BeginRawData = NHeaderRecords)
#####   out$Version <- txt["version",]
#####   out$Date <- txt["date",]
#####   out$DesignName <- txt["designname",]
#####   out$DesignId <- txt["designid",]
#####   out
##### }#readNimbleGenHeader
##### 
##### # XYS Contents in CAPITAL - contains Header 
##### ##	(x)	"X":
##### ##	(x)	"Y":
##### ##	(x)	"SIGNAL":
##### ##	(x)	"COUNT":
##### 
##### # PAIR Contents in CAPITAL - contains Header
##### ## (x)	"IMAGE_ID": 
##### ## (x)	"GENE_EXPR_OPTION": 
##### ## (x)	"SEQ_ID": 
##### ## (x)   "PROBE_ID":
##### ## (x) 	"POSITION": 
##### ##	(x)	"X":
##### ## (x)   "Y":
##### ## (x)   "MATCH_INDEX": 
##### ## (x)	"SEQ_URL": 
##### ## (x)	"PM":
##### ## (x)	"MM":
##### 
##### ## NDF Contents in CAPITAL - no Header
##### ## (x)   "PROBE_DESIGN_ID": not needed
##### ## (f)   "CONTAINER": Feature Table
##### ## (x)   "DESIGN_NOTE": not needed
##### ## (x)   "SELECTION_CRITERIA": not needed
##### ## (S)   "SEQ_ID": FeatureSet table (man_fsetid)
##### ## (s)   "PROBE_SEQUENCE": Sequence table (seq)
##### ## (x)   "MISMATCH": Will point to what table (pm/mm)Feature
##### ## (f)   "MATCH_INDEX": pair PM/MM - Feature Table
##### ## (f)   "FEATURE_ID": unit (4:9) - (unit_id)
##### ## (x)   "ROW_NUM": not required
##### ## (x)   "COL_NUM": not requied
##### ## (x)   "PROBE_CLASS": internal / not required (fiducial/control/blabla)
##### ## (f)   "PROBE_ID": Feature Table
##### ## (f)   "POSITION": Feature Table
##### ## (x)   "DESIGN_ID": not required
##### ## (f)   "X": Feature Table
##### ## (f)   "Y": Feature Table
##### ## (?)   "DMD"
##### 
##### ## NGD Contents in CAPITAL - no Header
##### ## (S)	"SEQ_ID": 
##### ## (S) 	"COMMENT": 
##### 
##### ## CALL Contents in CAPITAL - no Header
##### ##	(x)	"IMAGE_ID":
##### ## (x)	"SEQ_ID":
##### ## (x)	"PROBE_PAIRS":
##### ##	(x)	"FILTERED_PROBE_PAIRS":
##### ## (x)	"PM_AVG":
##### ##	(x)	"PM_CV":
##### ##	(x)	"MM_AVG":
##### ##	(x)	"MM_CV":
##### ##	(x)	"DIFF_AVG":
##### ##	(x)	"DIFF_CV":
##### ##	(x)	"GENE_EXPR_OPTION":
##### ##	(x)	"GENE_INFO":
##### 
##### ## POS Contents in CAPITAL - no Header
##### ##	(x)	"PROBE_ID":
##### ##	(x)	"SEQ_ID":
##### ##	(x)	"CHROMOSOME":
##### ##	(x)	"POSITION":
##### ##	(x)	"COUNT":
##### ##	(x)	"LENGTH":
##### ## (x)	"GC":
##### 
##### ## GFF Contents - contains Header ##gff-version	3 9 columns - no column header
##### ##	(x)	"seqid":
##### ##	(x)	"source":
##### ##	(x)	"type":
##### ##	(x)	"start":
##### ##	(x)	"end":
##### ##	(x)	"score":
##### ##	(x)	"strand":
##### ##	(x)	"phase":
##### ##	(x)	"attributes":


################################################
################################################
## FROM AllClasses.R
################################################
################################################

## setClass("AffyExpressionPDInfoPkgSeed",
##          contains="PDInfoPkgSeed",
##          representation=representation(
##                  cdfFile="ScalarCharacter",
##                  csvAnnoFile="ScalarCharacter",
##                  tabSeqFile="ScalarCharacter"),
##          prototype=prototype(
##                  manufacturer="Affymetrix",
##                  cdfFile=mkScalar(as.character(NA)),
##                  csvAnnoFile=mkScalar(as.character(NA)),
##                  tabSeqFile=mkScalar(as.character(NA))))
##

## modified by Matt Settles June 2,2008
## setClass("NgsExpressionPDInfoPkgSeed",
##          contains="NgsPDInfoPkgSeed",
##          representation=representation(
##            ndfFile="ScalarCharacter",
##            xysFile="ScalarCharacter", ## BC, 11/18/08
##            pairFile="ScalarCharacter",
##            ngdFile="ScalarCharacter"
##            ))
			
## modified by Matt Settles June 2,2008
## setClass("NgsTilingPDInfoPkgSeed",
##          contains="NgsPDInfoPkgSeed",
##          representation=representation(
##            ndfFile="ScalarCharacter",
##            xysFile="ScalarCharacter", ## BC, 11/18/08
##            pairFile="ScalarCharacter",
##            posFile="ScalarCharacter"
##            ))
## 
## validNgsTilingPDInfoPkgSeed <- function(object){
##   ndf <-  nchar(slot(object, "ndfFile")) > 0
##   xys <-  nchar(slot(object, "xysFile")) > 0
##   pair <- nchar(slot(object, "pairFile"))> 0
##   if (pair & xys) stop("Specify only one: XYS or PAIR file.")
##   if (!ndf) stop("NDF must be given.")
##   stopifnot(file.exists(slot(object, "ndfFile")))
##   if (xys) stopifnot(file.exists(slot(object, "xysFile")))
##   if (pair) stopifnot(file.exists(slot(object, "pairFile")))
##   TRUE
## }
## setValidity("NgsTilingPDInfoPkgSeed", validNgsTilingPDInfoPkgSeed)





################################################
## From chip-methods.R
################################################

## ## modified by Matt Settles June 2,2008
## #setMethod("chipName", "NgsPDInfoPkgSeed",
## #          function(object) {
## #            strsplit(tolower(object@ndfFile), ".ndf")[[1]]
## #          })
## 
## ## modified by Matt Settles June 2,2008
## setMethod("chipName", "NgsExpressionPDInfoPkgSeed",
##           function(object) {
##               ## compute chip name from the Pair file
##               header <- readPairHeader(object@pairFile)
##               header$DesignName
##           })
## 
## ## modified by Matt Settles June 2,2008


################################################
## From makePdInfoPackage.R
################################################

          
## modified by Matt Settles June 2,2008
## setMethod("makePdInfoPackage", "NgsExpressionPDInfoPkgSeed",
##     function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {
##   
## ## debug
## #object = pkg
## #destDir=targ
## #batch_size=10000
## #quiet=FALSE
## ## end debug
##     validInput <- function(x, destPath) {
##         msg <- NULL
##         ok <- sapply(c("pairFile", "xysFile"),
##             function(slt) file.exists(slot(x, slt)))
##         if (all(!ok))
##             msg <- "Must specify one of: pairFile or xysFile"
##         else if(all(ok))
##             msg <- "Must specify ONLY one of: pairFile or xysFile"
##         ok <- sapply(c("ndfFile"),
##             function(slt) file.exists(slot(x, slt)))
##         if (!all(ok))
##             msg <- paste(msg,"\n","missing file(s):",
##                 paste(sapply(names(ok[!ok]), function(slt) slt),
##                     "='",
##                     sapply(names(ok[!ok]), function(slt) slot(x, slt)),
##                     "'",
##                     collapse=", ", sep=""))
##         if (file.exists(destPath))
##             msg <- c(msg,"\n",
##                 paste("destination exists, remove or rename: '",
##                     destPath, "'", sep=""))
##         if (is.null(msg)) TRUE else msg
##     }
## 
##     chip <- chipName(object) ## only place that uses xys or pair file
##     pkgName <- cleanPlatformName(chip)
##     valid <- validInput(object, file.path(destDir, pkgName))
##     if (!identical(valid, TRUE))
##         stop(paste(valid, collapse="\n  "))
##     extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
##     dbFileName <- paste(pkgName, "sqlite", sep=".")
##     dbFilePath <- file.path(extdataDir, dbFileName)
##     ## get geometry from ndf File - MS
##     ndfdata <- read.delim(object@ndfFile, as.is=TRUE, header=TRUE)
##     geometry <- paste(max(ndfdata$Y), max(ndfdata$X), sep=";")
##     nx <- max(ndfdata$X)
##     
##     ## I don't see the point to all below
##     #ndf.idx <- as.integer(ndfdata$X+(ndfdata$Y-1)*nx)
##     ### do we really need any of this
##     #if(file.exists(object@pairFile)) { # pairFile used
##     #    pairdata <- read.delim(object@pairFile, as.is=TRUE, header=TRUE, comment.char="#")
##     #    pair.idx <- as.integer(pairdata$X+(pairdata$Y-1)*nx)
##     #    idx <- match(pair.idx, ndf.idx) ## match ndf to pair, removes Control probes 
##     #    ndfdata <- ndfdata[idx,]
##     #    rm(pairdata,pair.idx,idx)
##     #} else if (file.exists(object@xysfile)) { #xysFile used
##     #    xysdata <- read.delim(object@xysFile, as.is=TRUE, header=TRUE, comment.char="#")
##     #    xys.idx <- as.integer(xysdata$X+(xysdata$Y-1)*nx)
##     #    idx <- match(xys.idx, ndf.idx) ## match ndf to pair              
##     #    ndfdata <- ndfdata[idx,]
##     #    idx <- which(is.na(xysdata$SIGNAL)) 
##     #    if(length(idx) > 0) ndfdata[idx, "CONTAINER"] <- "OLIGO_CONTROL" # added if statement - MS                
##     #    rm(xysdata,xys.idx,idx) 
##     #}
##     ## in case one removes  H_CODE,  NGS_CONTROLS and V_CODE from the NDF file
##     ## in case two renames  H_CODE,  NGS_CONTROLS and V_CODE to OLIGO_CONTROL
##     # rm(ndf.idx, idx);gc();gc();
## 
##     ndfdata <- cbind(fid=as.integer(ndfdata$X+(ndfdata$Y-1)*nx) , ndfdata) 
##     ndfdata <- as.data.frame(ndfdata)
##     ndfdata <- ndfdata[order(ndfdata[["fid"]]),]
##     ndfdata[["fid"]] <- 1:nrow(ndfdata)  
##     syms <- list(MANUF=object@manufacturer,
##         VERSION=object@version,
##         GENOMEBUILD=object@genomebuild,
##         AUTHOR=object@author,
##         AUTHOREMAIL=object@email,
##         LIC=object@license,
##         DBFILE=dbFileName,
##         CHIPNAME=chip,
##         PKGNAME=pkgName,
##         PDINFONAME=pkgName,
##         PDINFOCLASS="NgsExpressionPDInfo",
##         GEOMETRY=geometry)
##     templateDir <- system.file("pd.PKG.template", package="pdInfoBuilder")
##     createPackage(pkgname=pkgName, destinationDir=destDir,
##     originDir=templateDir, symbolValues=syms, quiet=quiet)
##     dir.create(extdataDir, recursive=TRUE)
##     
##     buildPdInfoDb.ngsExprs(ndfdata,object@ngdFile, dbFilePath, # changed to buildPdInfoDb.ngsExprs - MS
##         batch_size=batch_size, verbose=!quiet)
## })

## TODO: fix for tiling arrays
## modified by Matt Settles June 2,2008
## setMethod("makePdInfoPackage", "NgsTilingPDInfoPkgSeed",
##           function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {
##             ## the validInput that was once here was moved to setValidity
##             chip <- chipName(object)
##             pkgName <- cleanPlatformName(chip)
##             extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
##             dbFileName <- paste(pkgName, "sqlite", sep=".")
##             dbFilePath <- file.path(extdataDir, dbFileName)
##             ## without geometry in header file this is the best alternative - MS
##             ## FIXME!!! need geometry 
##             ## Awww... it would be great if NGS had a header describing the array.
##             ## ndfdata was supposed to appear only on loadByBatch
##             ## appears here to avoid reading the same file twice
##             ndfdata <- read.delim(object@ndfFile, as.is=TRUE, header=TRUE)
##             geometry <- paste(max(ndfdata$Y), max(ndfdata$X), sep=";")
##             nx <- max(ndfdata$X)
##             ndf.idx <- as.integer(ndfdata$X+(ndfdata$Y-1)*nx)
##             if(file.exists(object@pairFile)) { # pairFile used
##               pairdata <- read.delim(object@pairFile, as.is=TRUE, header=TRUE, comment.char="#")
##               pair.idx <- as.integer(pairdata$X+(pairdata$Y-1)*nx)
##               idx <- match(pair.idx, ndf.idx) ## match ndf to pair
##               ndfdata <- ndfdata[idx,]
##               idx <- which(is.na(pairdata$PM))
##               rm(pairdata,pair.idx) 
##             }else if (file.exists(object@xysFile)) { #xysFile used
##               xysdata <- read.delim(object@xysFile, as.is=TRUE, header=TRUE, comment.char="#")
##               xys.idx <- as.integer(xysdata$X+(xysdata$Y-1)*nx)
##               idx <- match(xys.idx, ndf.idx) ## match ndf to xys
##               ndfdata <- ndfdata[idx,]
##               idx <- which(is.na(xysdata$SIGNAL))
##               rm(xysdata,xys.idx) 
##             }
##             if(length(idx) > 0) ndfdata[idx, "CONTAINER"] <- "OLIGO_CONTROL" # added if statement - MS                
##             ndfdata <- cbind(fid=as.integer(ndfdata$X+(ndfdata$Y-1)*nx), ndfdata)
##             ndfdata <- as.data.frame(ndfdata)
##             ndfdata <- ndfdata[order(ndfdata[["fid"]]),]
##             ndfdata[["fid"]] <- 1:nrow(ndfdata)
##             rm(ndf.idx, idx);gc();gc()
##             syms <- list(MANUF=object@manufacturer,
##                          VERSION=object@version,
##                          GENOMEBUILD=object@genomebuild,
##                          AUTHOR=object@author,
##                          AUTHOREMAIL=object@email,
##                          LIC=object@license,
##                          DBFILE=dbFileName,
##                          CHIPNAME=chip,
##                          PKGNAME=pkgName,
##                          PDINFONAME=pkgName,
##                          PDINFOCLASS="NgsTilingPDInfo",
##                          GEOMETRY=geometry)
##             templateDir <- system.file("pd.PKG.template",
##                                        package="pdInfoBuilder")
##             createPackage(pkgname=pkgName, destinationDir=destDir,
##                           originDir=templateDir, symbolValues=syms,
##                           quiet=quiet)
##             dir.create(extdataDir, recursive=TRUE)
##             buildPdInfoDb.ngsTiling(ndfdata, object@posFile, dbFilePath, ## changed to buildPdInfoDb.ngsTiling - MS
##                                     batch_size=batch_size, verbose=!quiet) 
##           })

## setMethod("makePdInfoPackage", "AffyExpressionPDInfoPkgSeed",
##     function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {
##         valid <- validInput(object,c("cdfFile"), c("csvAnnoFile", "tabSeqFile"))
##         if (!identical(valid, TRUE))
##             stop(paste(valid, collapse="\n  "))
## 
##         if(is.na(object@chipName)) object@chipName <- chipName(object)
##         pkgName <- cleanPlatformName(object@chipName)
##         
##         dbFileName <- paste(pkgName, "sqlite", sep=".")        
##         setupPackage(object,pkgName,destDir,dbFileName,unlink,quiet)
##         
##         extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
##         dbFilePath <- file.path(extdataDir, dbFileName)
##         seqMatFile <- file.path(extdataDir, "seqMat.rda")
##         
##         buildPdInfoDb.affyExpr(
##                 object@cdfFile, 
##                 object@csvAnnoFile,
##                 object@tabSeqFile, 
##                 dbFilePath, 
##                 seqMatFile,
##                 batch_size=batch_size, 
##                 verbose=!quiet)
## })

## Commented out to use the new builder. BC 03-11-09          
## setMethod("makePdInfoPackage", "AffyTilingPDInfoPkgSeed",
##     function(object, destDir=".", batch_size=1, quiet=FALSE, unlink=FALSE) {
## ## debug
## #destDir="."; batch_size=1; quiet=FALSE; unlink=TRUE
## ## end debug
##         valid <- validInput(object,c("bpmapFile", "celFile"))
##         if (!identical(valid, TRUE))
##             stop(paste(valid, collapse="\n  "))
##               
##         if(is.na(object@chipName)) object@chipName <- chipName(object)
##         pkgName <- cleanPlatformName(object@chipName)
##               
##         dbFileName <- paste(pkgName, "sqlite", sep=".")        
##         setupPackage(object,pkgName,destDir,dbFileName,unlink,quiet)
##         
##         extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
##         dbFilePath <- file.path(extdataDir, dbFileName)
##         seqMatFile <- file.path(extdataDir, "seqMat.rda")
##               
##         buildPdInfoDb.affyTiling(
##                 object@bpmapFile, 
##                 object@celFile,
##                 dbFilePath, 
##                 seqMatFile,
##                 batch_size=batch_size, 
##                 verbose=!quiet)
## })


##### setMethod("makePdInfoPackage", "AffySTPDInfoPkgSeed",
#####     function(object, destDir=".", batch_size=1, quiet=FALSE, unlink=FALSE) {
##### ## debug
##### #destDir="."; batch_size=1; quiet=FALSE; unlink=TRUE
##### ## end debug
##### 
#####         valid <- pdInfoBuilder:::validInput(object,c("pgfFile", "clfFile"),c("probeFile", "transFile"))
#####         if (!identical(valid, TRUE))
#####             stop(paste(valid, collapse="\n  "))
#####       
#####         if(is.na(object@chipName)) object@chipName <- chipName(object)
#####         pkgName <- cleanPlatformName(object@chipName)
#####       
#####         dbFileName <- paste(pkgName, "sqlite", sep=".")        
#####         pdInfoBuilder:::setupPackage(object,pkgName,destDir,dbFileName,unlink,quiet)
#####             
#####         extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
#####         dbFilePath <- file.path(extdataDir, dbFileName)
#####         seqMatFile <- file.path(extdataDir, "seqMat.rda")
#####         
##### #        buildPdInfoDb.affyST(
##### #                object@pgfFile, 
##### #                object@clfFile,
##### #                object@probeFile,
##### #                object@transFile, 
##### #                dbFilePath,
##### #                seqMatFile,
##### #                batch_size=batch_size,
##### #                verbose=!quiet)
#####             db <- pdInfoBuilder:::connectDb(dbFilePath)
#####             
#####             clf <- readClfEnv(object@clfFile)
#####             pgf <- readPgfEnv(object@pgfFile)
#####             
#####             t <- ST(loadUnits.affyST(db, pgf, clf))
#####             if (!quiet) printTime("loadUnitsByBatch", t[3])
#####             if(!is.na(slot(object,"transFile"))){
#####                 t <- ST(loadAffyCsv.affyST(db, object@transFile, batch_size=batch_size))
#####                 if (!quiet) printTime("loadAffyCsv", t[3])
#####             }
#####             if(!is.na(slot(object,"probeFile"))){
#####             t <- ST(loadAffySeqCsv.affyST(db, object@probeFile, pgf, batch_size=batch_size))
#####             if (!quiet) printTime("loadAffySeqCsv", t[3])
#####             }
#####             t <- ST({
#####                         sortFeatureTables.affyST(db)
#####                         createIndicesDb.affyST(db)
#####                         createTableInfoTable.affyST(db)
#####                         createFeatureTableInfo.affyST(db)
#####                     })
#####             if (!quiet) printTime("DB sort, index creation", t[3])
#####             
#####             ##     t <- ST({
#####             ##         seqMat <- createSeqMat(db)
#####             ##         save(seqMat, file=matFile, compress=TRUE)
#####             ##     })
#####             ##     printTime("sequence matrix", t[3])
#####             closeDb(db)
##### })

##### ### COPIED FROM ST - THIS WILL CHANGE
##### setMethod("makePdInfoPackage", "AffyExonPDInfoPkgSeed",
#####     function(object, destDir=".", batch_size=1, quiet=FALSE, unlink=FALSE) {
#####       valid <- pdInfoBuilder:::validInput(object, c("pgfFile", "clfFile"), c("probeFile", "transFile"))
#####       if (!identical(valid, TRUE))
#####         stop(paste(valid, collapse="\n  "))
#####       
#####       if(is.na(object@chipName)) object@chipName <- chipName(object)
#####       pkgName <- cleanPlatformName(object@chipName)
#####       
#####       dbFileName <- paste(pkgName, "sqlite", sep=".")        
#####       pdInfoBuilder:::setupPackage(object,pkgName,destDir,dbFileName,unlink,quiet)
#####       
#####       extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
#####       dbFilePath <- file.path(extdataDir, dbFileName)
#####       seqMatFile <- file.path(extdataDir, "seqMat.rda")
#####       db <- pdInfoBuilder:::connectDb(dbFilePath)
#####       
#####       clf <- readClfEnv(object@clfFile)
#####       pgf <- readPgfEnv(object@pgfFile)
#####       
#####       t <- ST(loadUnits.affyST(db, pgf, clf))
#####       if (!quiet) printTime("loadUnitsByBatch", t[3])
#####       if(!is.na(slot(object,"transFile"))){
#####         t <- ST(loadAffyCsv.affyST(db, object@transFile, batch_size=batch_size))
#####         if (!quiet) printTime("loadAffyCsv", t[3])
#####       }
#####       if(!is.na(slot(object,"probeFile"))){
#####         t <- ST(loadAffySeqCsv.affyST(db, object@probeFile, pgf, batch_size=batch_size))
#####         if (!quiet) printTime("loadAffySeqCsv", t[3])
#####       }
#####       t <- ST({
#####         sortFeatureTables.affyST(db)
#####         createIndicesDb.affyST(db)
#####         createTableInfoTable.affyST(db)
#####         createFeatureTableInfo.affyST(db)
#####       })
#####       if (!quiet) printTime("DB sort, index creation", t[3])
#####       closeDb(db)
##### })
##### 
##### ### COPIED FROM ST - THIS WILL CHANGE
##### setMethod("makePdInfoPackage", "AffyGenePDInfoPkgSeed",
#####     function(object, destDir=".", batch_size=1, quiet=FALSE, unlink=FALSE) {
#####       valid <- pdInfoBuilder:::validInput(object, c("pgfFile", "clfFile"), c("probeFile", "transFile"))
#####       if (!identical(valid, TRUE))
#####         stop(paste(valid, collapse="\n  "))
#####       
#####       if(is.na(object@chipName)) object@chipName <- chipName(object)
#####       pkgName <- cleanPlatformName(object@chipName)
#####       
#####       dbFileName <- paste(pkgName, "sqlite", sep=".")        
#####       pdInfoBuilder:::setupPackage(object,pkgName,destDir,dbFileName,unlink,quiet)
#####       
#####       extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
#####       dbFilePath <- file.path(extdataDir, dbFileName)
#####       seqMatFile <- file.path(extdataDir, "seqMat.rda")
#####       db <- pdInfoBuilder:::connectDb(dbFilePath)
#####       
#####       clf <- readClfEnv(object@clfFile)
#####       pgf <- readPgfEnv(object@pgfFile)
#####       
#####       t <- ST(loadUnits.affyST(db, pgf, clf))
#####       if (!quiet) printTime("loadUnitsByBatch", t[3])
#####       if(!is.na(slot(object,"transFile"))){
#####         t <- ST(loadAffyCsv.affyST(db, object@transFile, batch_size=batch_size))
#####         if (!quiet) printTime("loadAffyCsv", t[3])
#####       }
#####       if(!is.na(slot(object,"probeFile"))){
#####         t <- ST(loadAffySeqCsv.affyST(db, object@probeFile, pgf, batch_size=batch_size))
#####         if (!quiet) printTime("loadAffySeqCsv", t[3])
#####       }
#####       t <- ST({
#####         sortFeatureTables.affyST(db)
#####         createIndicesDb.affyST(db)
#####         createTableInfoTable.affyST(db)
#####         createFeatureTableInfo.affyST(db)
#####       })
#####       if (!quiet) printTime("DB sort, index creation", t[3])
#####       closeDb(db)
##### })
##### 


################################################
## From initialize-methods
################################################

## 
## setMethod("initialize", "AffyExpressionPDInfoPkgSeed",
##           function(.Object, cdfFile, csvAnnoFile, tabSeqFile, ...) {
##               .Object@cdfFile <- new("ScalarCharacter", cdfFile)
##               ## tabSeqFile and csvAnnoFile are optional
##               if (!missing(csvAnnoFile))
##                   .Object@csvAnnoFile <- new("ScalarCharacter", csvAnnoFile)
##               if (!missing(tabSeqFile))
##                   .Object@tabSeqFile <- new("ScalarCharacter", tabSeqFile)
##               .Object <- callNextMethod(.Object, ...)
##               .Object
##           })
## 

## modified by Matt Settles June 2,2008
## must specify one of pairFile or xysFile
## setMethod("initialize", "NgsExpressionPDInfoPkgSeed",
##           function(.Object, ndfFile, pairFile="", xysFile="", ngdFile="", ...) {
##               .Object@ndfFile <- new("ScalarCharacter", ndfFile)
##               .Object@pairFile <- new("ScalarCharacter", pairFile)
##               .Object@xysFile <- new("ScalarCharacter", xysFile)
##               .Object@ngdFile <- new("ScalarCharacter", ngdFile)
##               .Object <- callNextMethod(.Object, ...)
##               .Object
##           })

## modified by Matt Settles June 2,2008 ## do you combine multiple arrays sets?
## setMethod("initialize", "NgsTilingPDInfoPkgSeed",
##           function(.Object, ndfFile, posFile, xysFile, pairFile, ...) {
##               .Object@ndfFile <- new("ScalarCharacter", ndfFile)
##               .Object@posFile <- new("ScalarCharacter", posFile)
##               if (!missing(xysFile))
##                 .Object@xysFile <- new("ScalarCharacter", xysFile) ## BC 11/18/2008
##               if (!missing(pairFile))
##                 .Object@pairFile <- new("ScalarCharacter", pairFile)
##               .Object <- callNextMethod(.Object, ...)
##               .Object
##           })


################################################
################################################
################################################
################################################
