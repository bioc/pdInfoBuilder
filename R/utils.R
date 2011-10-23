## Utility functions for producing PDInfo packages
ST <- system.time

setPageSizeSql <- ('
pragma page_size = 32768;
')

setSyncOff <- ('
pragma synchronous = OFF;
')

increaseDbPerformance <- function(conn){
  dbGetQuery(conn, setPageSizeSql)
  dbGetQuery(conn, setSyncOff)
}


printTime <- function(msg, t) {
        m <- paste(msg, "took %.2f sec\n")
        cat(sprintf(m, t))
}

validInput <- function(object, required=c(), optional=c()) {
    msg <- NULL
    ok <- sapply(required, 
            function(slt) file.exists(slot(object, slt)))
    if (!all(ok))
        msg <- paste("missing required file(s):", 
                paste(sapply(names(ok[!ok]), function(slt) slt), "='", 
                        sapply(names(ok[!ok]), 
                                function(slt) slot(object, slt)), "'", collapse=", ", sep=""))
    ok <- sapply(optional, 
            function(slt) ifelse(is.na(slot(object, slt)), TRUE, file.exists(slot(object, slt))))
    if(!all(ok))
        msg <- paste(msg, "\n", "missing optional file(s):", 
                paste(sapply(names(ok[!ok]), function(slt) slt), "='", 
                        sapply(names(ok[!ok]), 
                                function(slt) slot(object, slt)), "'", collapse=", ", sep=""))
    
    if (is.null(msg)) TRUE else msg
} # end validInput

setupPackage <- function(object, pkgName, destDir, dbFileName, unlink, quiet){
    geometry  <- getGeometry(object)
    oligoc_dbi_version =installed.packages()['oligoClasses', 'Version']
    syms <- list(
            ## present in DESCRIPTION
            PKGNAME      =pkgName,
            MANUF        =object@manufacturer,
            CHIPNAME     =object@chipName,
            PKGVERSION   =object@version,
            AUTHOR       =object@author,
            AUTHOREMAIL  =object@email,
            MANUFURL     =object@url,
            LIC          =object@license,
            ORGANISM     =object@organism,
            GENOMEBUILD  =object@genomebuild,
            SPECIES      =object@species,
            BIOCVIEWS    =object@biocViews,
            OLIGOCVERSION=oligoc_dbi_version,
            ## present in namespace
            PDINFONAME   =pkgName,
            ## present in all.R
            ## not sure where these are used
            GEOMETRY     =paste(geometry$nrows, geometry$ncols, sep=";"), ## remove in future
            DBFILE       =dbFileName,
            PDINFOCLASS  = sub("PkgSeed", "", class(object)))
    
    templateDir <- system.file("pd.PKG.template", package="pdInfoBuilder")
    createPackage(pkgname=pkgName, destinationDir=destDir,
            originDir=templateDir, symbolValues=syms,
            unlink = unlink, quiet=quiet)
    file.rename(file.path(destDir, pkgName, "man/template.Rd"),
                file.path(destDir, pkgName, "man", paste(pkgName, "Rd", sep=".")))
}              

connectDb <- function(dbfile) {
    require("RSQLite")
    db <- dbConnect(SQLite(), dbname=dbfile, cache_size=6400, synchronous=0)
    sql <- ('
            pragma page_size = 8192;                
            ')
    sqliteQuickSQL(db, sql)
    db
}

closeDb <- function(db){
    dbDisconnect(db)
}

###################
## create chrom dictionary
## the function takes a vector with all chromosome
## values and returns a data frame with an integer id
## for each value.
## The integer ID is used in the pmfeature-like tables
##################
createChrDict <- function(x){
  possible <- as.character(na.omit(unique(x)))
  possible <- gsub("chr([1-9])$", "chr0\\1", possible)
  possible <- gsub("chr([1-9]_)", "chr0\\1", possible)
  dataSplit <- strsplit(possible, "_")
  len <- sapply(dataSplit, length)
  idx <- which(len == 1)
  basic <- unlist(dataSplit[idx])
  if (length(idx) < length(len)){
    suffixes <- sapply(dataSplit[-idx], function(x) paste(x[-1], collapse="_"))
    suffixes <- sort(unique(suffixes))
    out <- list()
    out[[1]] <- basic
    for (i in 1:length(suffixes))
      out[[i+1]] <- paste(basic, suffixes[i], sep="_")
    out <- unlist(out)
  }else{
    out <- sort(basic)
  }
  out <- gsub("chr0([1-9])", "chr\\1", out)
  data.frame(chrom=as.integer(1:length(out)),
             chrom_id=out,
             stringsAsFactors=FALSE)
}

### helpers
### table creation
dbCreateTableInfo <- function(db, verbose=FALSE) {
    tables <- dbListTables(db)
    counts <- integer(length(tables))
    sql <- "select count(*) from %s"
    for (i in seq(along=counts)) {
        if (verbose) message("Counting rows in ", tables[i])
        counts[i] <- dbGetQuery(db, sprintf(sql, tables[i]))[[1]][1]
    }

    df <- data.frame(tbl=tables, row_count=counts,
                     stringsAsFactors=FALSE)
    dbWriteTable(db, "table_info", df, row.names=FALSE)
}

dbCreateTable <- function(conn, tablename, col2type, col2key){
    col2type[names(col2key)] <- paste(col2type[names(col2key)], col2key, sep=" ")
    sql <- paste(names(col2type), col2type, sep=" ")
    sql <- paste(sql, collapse=", ")
    sql <- paste("CREATE TABLE ", tablename, " (", sql, ")", sep="")
    dbGetQuery(conn, sql)
}

dbInsertDataFrame <- function(conn, tablename, data, col2type, verbose=FALSE){
  cols <- names(col2type)
  if (!identical(sort(names(data)), sort(cols)))
    stop("cols in data frame 'data' don't match cols in table \"", tablename, "\"")
  values_template <- paste(paste(":", cols, sep=""), collapse=", ")
  sql_template <- paste("INSERT INTO ", tablename, " VALUES (", values_template, ")", sep="")
  if (verbose)
    simpleMessage("Inserting ", nrow(data), " rows into table ", tablename, "... ")
  dbBeginTransaction(conn)
  on.exit(dbCommit(conn))
  dbGetPreparedQuery(conn, sql_template, bind.data=data)
  if (verbose) msgOK()
}

dbCreateIndex <- function(conn, idxname, tblname, fieldname, unique=TRUE, verbose=TRUE){
  if (verbose) simpleMessage("Creating index ", idxname, " on ", tblname, "... ")
  sql <- paste("CREATE", ifelse(unique, "UNIQUE", ""),
               "INDEX", idxname, "ON", tblname,
               paste("(", fieldname, ")", sep=""))
  dbGetQuery(conn, sql)
  if (verbose) msgOK()
  NULL
}

dbCreateIndicesBg <- function(conn, verbose=TRUE){
  dbCreateIndex(conn, "idx_bgfsetid", "bgfeature", "fsetid", FALSE, verbose=verbose)
  dbCreateIndex(conn, "idx_bgfid", "bgfeature", "fid", TRUE, verbose=verbose)
}

dbCreateIndicesBgTiling <- function(conn, verbose=TRUE){
  dbCreateIndex(conn, "idx_bgfid", "bgfeature", "fid", TRUE, verbose=verbose)
}

dbCreateIndicesPm <- function(conn, verbose=TRUE){
  dbCreateIndex(conn, "idx_pmfsetid", "pmfeature", "fsetid", FALSE, verbose=verbose)
  dbCreateIndex(conn, "idx_pmfid", "pmfeature", "fid", TRUE, verbose=verbose)
}

dbCreateIndicesMm <- function(conn, verbose=TRUE){
  dbCreateIndex(conn, "idx_mmfid", "mmfeature", "fid", FALSE, verbose=verbose)
  dbCreateIndex(conn, "idx_mmpmfid", "mmfeature", "fidpm", FALSE, verbose=verbose)
}

dbCreateIndicesSnpPm <- function(conn, verbose=TRUE){
  dbCreateIndex(conn, "idx_pmfsetid", "pmfeature", "fsetid", FALSE, verbose=verbose)
  dbCreateIndex(conn, "idx_pmfid", "pmfeature", "fid", TRUE, verbose=verbose)
}

dbCreateIndicesCnvPm <- function(conn, verbose=TRUE){
  dbCreateIndex(conn, "idx_pmfsetidcnv", "pmfeatureCNV", "fsetid", FALSE, verbose=verbose)
  dbCreateIndex(conn, "idx_pmfidcnv", "pmfeatureCNV", "fid", TRUE, verbose=verbose)
}

dbCreateIndicesPmTiling <- function(conn, verbose=TRUE){
  dbCreateIndex(conn, "idx_pmfid", "pmfeature", "fid", TRUE, verbose=verbose)
}

dbCreateIndicesFs <- function(conn, verbose=TRUE){
  dbCreateIndex(conn, "idx_fsfsetid", "featureSet", "fsetid", TRUE, verbose=verbose)
}

dbCreateIndicesFsCnv <- function(conn, verbose=TRUE){
  dbCreateIndex(conn, "idx_fscnvfsetid", "featureSetCNV", "fsetid", TRUE, verbose=verbose)
}

## Messages

msgParsingFile <- function(fname)
  simpleMessage("Parsing file: ", basename(fname), "... ")

msgOK <- function() cat("OK\n")

msgBar <- function(){
  n <- options()[["width"]]
  cat(paste(c(rep("=", n), "\n"), collapse=""))
}

##
createSeqMat <- function(db) {
    allseq <- dbGetQuery(db, "select seq from sequence")[[1]]
    m <- .Call("PIB_25mers_to_mat", allseq, PACKAGE="pdInfoBuilder")
    remove(allseq)
    allfid <- dbGetQuery(db, "select fid from sequence")[[1]]
    rownames(m) <- allfid
    remove(allfid)
    m
}

createSeqMat.cnv <- function(db) {
    allseq <- dbGetQuery(db, "select seq from sequenceCNV")[[1]]
    m <- .Call("PIB_25mers_to_mat", allseq, PACKAGE="pdInfoBuilder")
    remove(allseq)
    allfid <- dbGetQuery(db, "select fid from sequence")[[1]]
    rownames(m) <- allfid
    remove(allfid)
    m
}


seqToMat <- function(seq) {
    .Call("PIB_25mers_to_mat", seq, PACKAGE="pdInfoBuilder")
}

annot2fdata <- function(csv){
    annot <- read.csv(csv, header=TRUE, comment.char="#",
                      na.strings="---", stringsAsFactors=FALSE)
    nms <- tolower(names(annot))
    nms <- gsub("_", "", nms)
    nms <- gsub("\\.", "", nms)
    names(annot) <- nms
    stopifnot('probesetid' %in% nms)
    rownames(annot) <- annot[['probesetid']]
    annot[['probeset_id']] <- NULL
    new("AnnotatedDataFrame", data=annot)
}


## insert in fragmentLengthTable

insertInFragmentLengthTable <- function(db, tblTarget, fragCol,
                                        manfsetid){
    tblref <- ifelse(length(grep("CNV", tblTarget)) == 0,
                     'featureSet', 'featureSetCNV')
    ref <- dbGetQuery(db, paste('SELECT man_fsetid, fsetid FROM', tblref))
    idx <- manfsetid %in% ref[['man_fsetid']]
    manfsetid <- manfsetid[idx]
    fragCol <- fragCol[idx]
    idx <- match(manfsetid, ref[['man_fsetid']])
    idx <- ref[idx, 'fsetid']
    rm(ref)
    flTable <- getFragmentLengthTable(fragCol, idx)
    rm(idx)
    dbInsertDataFrame(db, tablename=tblTarget, data=flTable,
                      col2type=c(fsetid='INTEGER', enzyme='TEXT',
                      length='INTEGER', start='INTEGER', stop='INTEGER'),
                      verbose=TRUE)
}
