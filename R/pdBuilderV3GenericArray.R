## This is the draft for a Generic Array Builder
## - Avoid genomic information (position / chromosome / gene symbols)
## - Use a V2Exon-like schema as template (ie., we always have meta-probesets)
##   * One table for probes
##   * Possibly many tables for probesets (at different levels)


##########################################
#### EXPRESSION AFFYMETRIX ARRAYS
####
#### probes: fid, fsetid, x, y (pm / mm)
#### probeset: fsetid, man_fsetid
####
#### ps1: one probe can belong to multiple probesets
#### ps2:

cdfUnits2table <- function(cdfFile, units){
    cdf <- readCdf(filename=cdfFile, units=units, readXY=TRUE, readBases=FALSE,
                   readIndexpos=FALSE, readAtoms=TRUE, readUnitType=TRUE,
                   readUnitDirection=TRUE, readUnitNumber=FALSE,
                   readGroupAtomNumbers=FALSE, readGroupDirection=FALSE,
                   readIndices=TRUE, readIsPm=TRUE)
    foreach (i=1:length(cdf), .combine=rbind, .inorder=FALSE) %do% {
        unit <- cdf[[i]]
        ## unit: x, y, indices, atom, ispm (in this sequence)
        nm <- names(cdf[i])
        tmp <- as.data.frame(unit$groups[[1]])
        tmp$ispm <- ifelse(tmp$ispm, 'pm', 'mm')
        tmp$ispm <- as.character(tmp$ispm)
        names(tmp) <- c('x', 'y', 'fid', 'atom', 'probe_type')
        tmp <- cbind(fs1_man_fsetid=nm, fs1_type=unit$unittype,
                     fs1_direction=unit$unitdirection,
                     tmp, stringsAsFactors=FALSE)
        tmp
    }
}

cdf2table <- function(cdfFile){
    units <- readCdfUnitNames(cdfFile)
    nworkers <- getDoParWorkers()
    grps <- sort(rep(1:nworkers, length.out=length(units)))
    unitNumbers <- split(1:length(units), grps)
    rm(units, nworkers, grps)
    foreach (unitLst=unitNumbers, .combine=rbind, .inorder=FALSE) %dopar% {
        cdfUnits2table(cdfFile, unitLst)
    }
}

sequenceParser <- function(seqFile){
    ## tab-delimited file
    ## assumes that the first 5 columns are (in order)
    ## name / x / y / position / sequence
    ## needs Biostrings
    seqtbl <- read.table(seqFile, header=TRUE, sep='\t', stringsAsFactors=FALSE)[, 1:5]
    names(seqtbl) <- c('name', 'x', 'y', 'position', 'sequence')
    unique(seqtbl[, -1])
}

getTotalFSetTables <- function(probetbl)
    max(as.integer(substr(grep('^fs([[:digit:]])', names(probetbl), value=TRUE), 3, 3)))

getFSetMpsTables <- function(probetbl, fsetID){
    ## fsetID is simply the "featureSet table number"
    theMax <- getTotalFSetTables(probetbl)
    if (fsetID > theMax) stop('The probe table has information for ', theMax, ' featureSet table(s)')
    regex <- paste0('^fs', fsetID, '_')
    flds <- grep(regex, names(probetbl), value=TRUE)
    tbl <- probetbl[, c(flds, 'fid', 'probe_type')]
    names(tbl) <- gsub(regex, '', names(tbl))
    flds <- gsub(regex, '', flds)
    tbl <- cbind(fsetid=as.integer(factor(tbl$man_fsetid)), tbl)
    out <- lapply(split(tbl[, c('fsetid', 'fid')], tbl$probe_type), unique)
    names(out) <- paste0('mps', fsetID, names(out))
    fset <- list(unique(tbl[, c('fsetid', flds)]))
    names(fset) <- paste0('featureSet', fsetID)
    rm(tbl)
    if (any(duplicated(fset[[1]]$man_fsetid)))
        stop('There are duplicates in man_fsetid for featureSet', fsetID)
    out <- append(fset, out)
    out
}

getAllFSetMpsTables <- function(probetbl){
    n <- getTotalFSetTables(probetbl)
    out <- foreach(i=1:n) %dopar% {
        getFSetMpsTables(probetbl, i)
    }
    unlist(out, recursive=FALSE)
}

getUniqueProbesTables <- function(probetbl){
    flds <- setdiff(names(probetbl), grep('^fs\\d\\_', names(probetbl), value=TRUE))
    if ( !('probe_type' %in% flds) ) stop('Probe table does not contain probe_type field')
    out <- lapply(split(probetbl[, setdiff(flds, 'probe_type')], probetbl$probe_type), unique)
    names(out) <- paste0(names(out), 'feature')
    out
}

## MUST HAVE:
##
## featureSetK
##  - fsetid: integer (PRIMARY KEY)
##  - man_fsetid: character
##
## xxfeature (xx = {pm, mm, ...})
##  - fid: integer (PRIMARY KEY)
##  - x: integer
##  - y: integer
##  - sequence: character
##
## mpsKxx
##  - fsetid: integer (REFERENCE TO featureSetK(fsetid))
##  - fid: integer (REFERENCE TO xxfeature(fid))

getSchema <- function(tbl){
    col2type <- toupper(unlist(lapply(tbl, storage.mode)))
    ## SQLite datatypes: INTEGER, REAL, TEXT, BLOB, NULL
    col2type[col2type=='DOUBLE'] <- 'REAL'
    col2type[col2type=='CHARACTER'] <- 'TEXT'
    col2type[col2type=='LOGICAL'] <- 'INTEGER'
    col2type
}

getKeys <- function(tblname){
    tblTypes <- c('featureSet', 'mps', 'feature')
    isFeatureSet <- length(grep('^featureSet', tblname)) == 1L
    isMPSpm <- length(grep('^mps[[:digit:]]pm$', tblname)) == 1L
    isMPSmm <- length(grep('^mps[[:digit:]]mm$', tblname)) == 1L
    isFeature <- length(grep('feature$', tblname)) == 1L
    if (isFeatureSet + isMPSpm + isMPSmm + isFeature != 1L) stop('Impossible to detect table type')
    if (isFeatureSet){
        col2key <- c(fsetid='PRIMARY KEY')
    } else if (isMPSpm){
        digit <- substr(tblname, 4, 4)
        col2key <- c(fid='REFERENCES pmfeature(fid)', fsetid=paste0('REFERENCES featureSet', digit, '(fsetid)'))
    } else if (isMPSmm){
        digit <- substr(tblname, 4, 4)
        col2key <- c(fid='REFERENCES mmfeature(fid)', fsetid=paste0('REFERENCES featureSet', digit, '(fsetid)'))
    } else if (isFeature){
        col2key <- c(fid='PRIMARY KEY')
    }
    col2key
}

getSchemaAndKeys <- function(lst2sql){
    tbls <- names(lst2sql)
    out <- vector('list', length(tbls))
    names(out) <- tbls
    for (tbl in tbls)
        out[[tbl]] <- list(col2type=getSchema(lst2sql[[tbl]]),
                           col2type=getKeys(tbl))
    out
}

setMethod("makePdInfoPackage", "GenericPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {

              ## Probe input table
              ## - fsK_man_fsetid (... fsK_fieldName) / fid / x / y / sequence / probe_type* [pm/mm/bg/...] / (... fields)
              
              ## Requirements for SQL (required) tables
              ## - featureSetN [the fields are extracted from the fsK_ prefixed columns]
              ##   * fsetid [automatically created] / man_fsetid / fieldName
              ## - pmfeature / mmfeature [[pm/mm/bg/...]feature tables are created based on the values in probe_type]
              ##   * fid / x / y / sequence / fields [everything not prefixed with fsK_ comes to this table]
              ## - mpsN (links probes to [meta-]probesets
              ##   * pmfeature(fid) / featureSetN(fsetidN)
              
              msgBar()
              message("Building annotation for Generic Array")
              msgBar()

              #######################################################################
              ## Part i) get array info (chipName, pkgName, dbname)
              #######################################################################
              pkgName <- object@pkgName
              extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
              dbFileName <- paste(pkgName, "sqlite", sep=".")
              dbFilePath <- file.path(extdataDir, dbFileName)

              #######################################################################
              ## Part ii) parse data. This should return a list of data.frames.
              ##          The names of the elements in the list are table names.
              #######################################################################
              geometry <- apply(object@table[, c('x', 'y')], 2, max)
              geometry <- paste(geometry, collapse=';')
              fsetTbls <- getAllFSetMpsTables(object@table)
              probeTbls <- getUniqueProbesTables(object@table)
              toSQL <- append(fsetTbls, probeTbls)
              rm(fsetTbls, probeTbls)
              allSchema <- getSchemaAndKeys(toSQL)

              #######################################################################
              ## Part iii) Create package from template
              #######################################################################
              syms <- list(MANUF=object@manufacturer,
                           VERSION=object@version,
                           GENOMEBUILD=object@genomebuild,
                           AUTHOR=object@author,
                           AUTHOREMAIL=object@email,
                           LIC=object@license,
                           DBFILE=dbFileName,
                           CHIPNAME=object@chipName,
                           PKGNAME=pkgName,
                           PDINFONAME=pkgName,
                           PDINFOCLASS='GenericPDInfo',
                           GEOMETRY=geometry)
              templateDir <- system.file("pd.PKG.template",
                                         package="pdInfoBuilder")
              createPackage(pkgname=pkgName, destinationDir=destDir,
                            originDir=templateDir, symbolValues=syms,
                            quiet=quiet)
              dir.create(extdataDir, recursive=TRUE)

              #######################################################################
              ## Part iv) Create SQLite database
              #######################################################################
              conn <- dbConnect(dbDriver("SQLite"), dbname=dbFilePath)
              for (tbl in names(allSchema)){
                  dbCreateTable(conn, tbl, allSchema[[tbl]][['col2type']], allSchema[[tbl]][['col2key']])
                  dbInsertDataFrame(conn, tbl, toSQL[[tbl]], allSchema[[tbl]][['col2type']], !quiet)
              }
              dbGetQuery(conn, "VACUUM")
              dbCreateTableInfo(conn, !quiet)
              dbDisconnect(conn)
          })
