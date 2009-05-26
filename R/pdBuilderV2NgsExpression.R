#######################################################################
## SECTION A - Db Schema
#######################################################################
ngsExprFeatureSetSchema <- list(col2type=c(
                                  fsetid="INTEGER",
                                  man_fsetid="TEXT"),
                                col2key=c(
                                  fsetid="PRIMARY KEY"
                                  ))
ngsExprPmFeatureSchema <- list(col2type=c(
                                 fid="INTEGER",
                                 fsetid="INTEGER",
                                 position="INTEGER",
                                 x="INTEGER",
                                 y="INTEGER"
                                 ),
                               col2key=c(
                                 fid="PRIMARY KEY"
                                 ))
## TODO: ngsExprMmFeatureSchema
ngsExprBgFeatureSchema <- list(col2type=c(
                                 fid="INTEGER",
                                 fsetid="INTEGER",
                                 x="INTEGER",
                                 y="INTEGER"
                                 ),
                               col2key=c(
                                 fid="PRIMARY KEY"
                                 ))


#######################################################################
## SECTION C - Parser specific for NGS Tiled Regions
##             This will take NDF/XYS pair and process (in memory)
##             the data, generating data.frames for each table
##             (featureSet, pmfeature, mmfeature*, bgfeature**)
##             to be created in the db.
#######################################################################

parseNgsPair <- function(ndfFile, xysFile, verbose=TRUE){
  stopifnot(!missing(ndfFile), !missing(xysFile))

  #######################################################################
  ## Step 1: Parse NDF
  #######################################################################
  if (verbose) msgParsingFile(ndfFile)
  ndfdata <- read.delim(ndfFile, stringsAsFactors=FALSE)
  if (verbose) msgOK()
  ndfdata[["fsetid"]] <- as.integer(as.factor(ndfdata[["SEQ_ID"]]))

  #######################################################################
  ## Step 3.1: Get XYS files and remove all controls (ie, NA in XYS)
  #######################################################################
  if (verbose) msgParsingFile(xysFile)
  xysdata <- read.delim(xysFile, comment="#")
  if (verbose) msgOK()
  xysdata[["fid"]] <- 1:nrow(xysdata)
  if (verbose) cat("Merging NDF and XYS files ...")
  ndfdata <- merge(ndfdata, xysdata, by.x=c("X", "Y"), by.y=c("X", "Y"))
  if (verbose) cat("OK\n")
  controls <- which(is.na(ndfdata[["SIGNAL"]]))
  ndfdata <- ndfdata[-controls,]
  rm(controls)
  
  #######################################################################
  ## Step 4: Prepare contents for featureSet table
  ## Fields featureSet: man_fsetid, chrom, start, end, type
  #######################################################################
  if (verbose) cat("Preparing contents for featureSet table ...")
  colsFS <- c("fsetid", "SEQ_ID")
  featureSet <- ndfdata[, colsFS]
  rm(colsFS)
  names(featureSet) <- c("fsetid", "man_fsetid")
  ok <- which(!duplicated(featureSet[["man_fsetid"]]))
  featureSet <- featureSet[ok,]
  rm(ok)
  if (verbose) cat("OK\n")

  
  #######################################################################
  ## Step 5: Prepare contents for pmfeature, mmfeature and bgfeature
  ## Fields pmfeature: fid, fsetid, position, x, y
  ##        mmfeature: fid, fsetid, position, x, y, fid(PM)
  ##        bgfeature: fid, fsetid, x, y
  ## FIX ME: Need to check the selection for BG probes
  ##         Implement mmfeature
  ## Comments: ideally, the tables would be ordered by chrom+position;
  ##           but 'fid' is the PRIMARY KEY, and the table will be
  ##           ordered by that
  #######################################################################
  features <- ndfdata[, c("X", "Y", "fsetid", "POSITION", "MISMATCH",
                          "MATCH_INDEX", "PROBE_SEQUENCE", "fid", "SEQ_ID")]
  names(features) <- c("x", "y", "fsetid", "position", "mismatch",
                       "match_index", "sequence", "fid", "man_fsetid")

  geometry <- paste(max(ndfdata[["Y"]]), max(ndfdata[["X"]]), sep=";")
  rm(xysdata, ndfdata)
  ## FIX ME: Double check ordering
  ## FIX ME: This should be passed by function

  if (verbose) cat("Preparing contents for bgfeature table ...")
  bgidx <- grep("RANDOM", features[["man_fsetid"]])
  bgFeatures <- bgSequence <- NULL
  if (length(bgidx) > 0){
    bgFeatures <- features[bgidx, c("fid", "fsetid", "x", "y")]
    bgSequence <- features[bgidx, c("fid", "sequence")]
    bgSequence <- bgSequence[order(bgSequence[["fid"]]),]
    bgSequence <- DataFrame(fid=bgSequence[["fid"]],
                            sequence=DNAStringSet(bgSequence[["sequence"]]))
    features <- features[-bgidx,]
  }
  rm(bgidx)
  if (verbose) cat("OK\n")

  if (verbose) cat("Preparing contents for pmfeature table ...")
  pmidx <- which(features[["mismatch"]] == 0)
  pmFeatures <- features[pmidx, c("fid", "fsetid", "position", "x", "y")]
  pmSequence <- features[pmidx, c("fid", "sequence")]
  pmSequence <- pmSequence[order(pmSequence[["fid"]]),]
  pmSequence <- DataFrame(fid=pmSequence[["fid"]],
                          sequence=DNAStringSet(pmSequence[["sequence"]]))
  mmFeatures <- features[-pmidx,]
  rm(pmidx)
  if (verbose) cat("OK\n")

  
  ## add mmSequence
  if (any(mmFeatures[["mismatch"]] >= 10000))
    stop("Control probe possibly identified as Experimental")
  if (nrow(mmFeatures) > 0)
    stop("Add methods for MMs")
  rm(mmFeatures)
  rm(features)

  return(list(featureSet=featureSet, pmFeatures=pmFeatures,
              bgFeatures=bgFeatures, geometry=geometry,
              pmSequence=pmSequence, bgSequence=bgSequence))
}


#######################################################################
## SECTION D - Package Maker
##             This shouldn't be extremely hard.
##             The idea is to: i) get array info (name, pkgname, dbname)
##             ii) parse data; iii) create pkg from template;
##             iv) dump the database
#######################################################################
setMethod("makePdInfoPackage", "NgsExpressionPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {

            msgBar()
            cat("Building annotation package for Nimblegen Expression Array\n")
            cat("NDF: ", basename(object@ndfFile), "\n")
            cat("XYS: ", basename(object@xysFile), "\n")
            msgBar()
            
            #######################################################################
            ## Part i) get array info (chipName, pkgName, dbname)
            #######################################################################
            chip <- chipName(object)
            pkgName <- cleanPlatformName(chip)
            extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
            dbFileName <- paste(pkgName, "sqlite", sep=".")
            dbFilePath <- file.path(extdataDir, dbFileName)

            #######################################################################
            ## Part ii) parse data. This should return a list of data.frames.
            ##          The names of the elements in the list are table names.
            #######################################################################
            parsedData <- parseNgsPair(object@ndfFile, object@xysFile, verbose=!quiet)

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
                         CHIPNAME=chip,
                         PKGNAME=pkgName,
                         PDINFONAME=pkgName,
                         PDINFOCLASS="NgsExpressionPDInfo",
                         GEOMETRY=parsedData[["geometry"]])
            templateDir <- system.file("pd.PKG.template",
                                       package="pdInfoBuilder")
            createPackage(pkgname=pkgName, destinationDir=destDir,
                          originDir=templateDir, symbolValues=syms,
                          quiet=quiet)
            dir.create(extdataDir, recursive=TRUE)

            #######################################################################
            ## Part iv) Create SQLite database
            ## FIX ME: Trusting tiledRegionFeatureSetSchema will be visible from
            ##         inside the method;
            ##         Fix ordering of the tables
            #######################################################################
            conn <- dbConnect(dbDriver("SQLite"), dbname=dbFilePath)
            increaseDbPerformance(conn)
            dbCreateTable(conn,
                          "featureSet",
                          ngsExprFeatureSetSchema[["col2type"]],
                          ngsExprFeatureSetSchema[["col2key"]])
            dbCreateTable(conn,
                          "pmfeature",
                          ngsExprPmFeatureSchema[["col2type"]],
                          ngsExprPmFeatureSchema[["col2key"]])
            containsMm <- "mmFeatures" %in% names(parsedData)
            if (containsMm)
              dbCreateTable(conn,
                            "mmfeature",
                            ngsExprMmFeatureSchema[["col2type"]],
                            ngsExprMmFeatureSchema[["col2key"]])
            containsBg <- !is.null(parsedData[["bgFeatures"]])
            if (containsBg)
              dbCreateTable(conn,
                            "bgfeature",
                            ngsExprBgFeatureSchema[["col2type"]],
                            ngsExprBgFeatureSchema[["col2key"]])

            dbInsertDataFrame(conn, "featureSet", parsedData[["featureSet"]],
                              ngsExprFeatureSetSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                              ngsExprPmFeatureSchema[["col2type"]], !quiet)
            if (containsMm)
              dbInsertDataFrame(conn, "mmfeature", parsedData[["mmFeatures"]],
                                ngsExprMmFeatureSchema[["col2type"]], !quiet)
            if (containsBg)
              dbInsertDataFrame(conn, "bgfeature", parsedData[["bgFeatures"]],
                                ngsExprBgFeatureSchema[["col2type"]], !quiet)
            dbGetQuery(conn, "VACUUM")

            dbCreateTableInfo(conn, !quiet)

            ## Create indices
            if (containsBg)
              dbCreateIndicesBg(conn, !quiet)
            dbCreateIndicesPm(conn, !quiet)
            dbCreateIndicesFs(conn, !quiet)
            
            
            dbDisconnect(conn)
            
            #######################################################################
            ## Part v) Save sequence DataFrames
            ## FIX ME: Fix ordering of the tables to match xxFeature tables
            #######################################################################
            datadir <- file.path(destDir, pkgName, "data")
            dir.create(datadir)

            if (!quiet) cat("Saving DataFrame object for PM.\n")
            pmSequence <- parsedData[["pmSequence"]]
            pmSeqFile <- file.path(datadir, "pmSequence.rda")
            save(pmSequence, file=pmSeqFile)

            if (containsBg){
              if (!quiet) cat("Saving DataFrame object for BG.\n")
              bgSequence <- parsedData[["bgSequence"]]
              bgSeqFile <- file.path(datadir, "bgSequence.rda")
              save(bgSequence, file=bgSeqFile)
            }
            if (!quiet) cat("Done.")
          })
