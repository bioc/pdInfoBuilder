## The following tables will be available:
## pmfeature: fid, fsetid, chrom, location, x, y
## featureSet: fsetid (region name), chrom, start, end
## mmfeature (opt): fid, fid (matching pm), x, y
## bgfeature: fid, x, y


#######################################################################
## SECTION A - Db Schema
#######################################################################
tiledRegionFeatureSetSchema <- list(col2type=c(
                                      fsetid="INTEGER",
                                      man_fsetid="TEXT",
                                      chrom="TEXT",
                                      start="INTEGER",
                                      end="INTEGER",
                                      type="TEXT"),
                                    col2key=c(
                                      fsetid="PRIMARY KEY"
                                      ))
tiledRegionPmFeatureSchema <- list(col2type=c(
                                     fid="INTEGER",
                                     fsetid="INTEGER",
                                     position="INTEGER",
                                     x="INTEGER",
                                     y="INTEGER"
                                     ),
                                   col2key=c(
                                     fid="PRIMARY KEY"
                                     ))
## TODO: tiledRegionMmFeatureSchema
tiledRegionBgFeatureSchema <- list(col2type=c(
                                     fid="INTEGER",
                                     fsetid="INTEGER",
                                     x="INTEGER",
                                     y="INTEGER"
                                     ),
                                   col2key=c(
                                     fid="PRIMARY KEY"
                                     ))


#######################################################################
## SECTION B - Utils - This should be moved from here
##             as everything in this section can be used on other cases
#######################################################################

#######################################################################
## SECTION C - Parser specific for NGS Tiled Regions
##             This will take NDF/POS/XYS trio and process (in memory)
##             the data, generating data.frames for each table
##             (featureSet, pmfeature, mmfeature*, bgfeature**)
##             to be created in the db.
#######################################################################
parseNgsTrio <- function(ndfFile, posFile, xysFile, verbose=TRUE){
  stopifnot(!missing(ndfFile), !missing(posFile), !missing(xysFile))

  #######################################################################
  ## Step 1: Parse NDF
  #######################################################################
  if (verbose) msgParsingFile(ndfFile)
  ndfdata <- read.delim(ndfFile, stringsAsFactors=FALSE)
  if (verbose) msgOK()
  ndfdata[["fsetid"]] <- as.integer(as.factor(ndfdata[["SEQ_ID"]]))

  #######################################################################
  ## Step 2: Parse POS
  #######################################################################
  if (verbose) msgParsingFile(posFile)
  posdata <- read.delim(posFile, stringsAsFactors=FALSE)
  posdata <- posdata[!duplicated(posdata),]
  rownames(posdata) <- NULL
  if (verbose) msgOK()
  
  #######################################################################
  ## Step 3: Match POS and NDF
  ##         and update positional info
  #######################################################################
  if (verbose) cat("Merging NDF and POS files... ")
  ndfdata <- merge(ndfdata, posdata, by.x=c("SEQ_ID", "PROBE_ID"), by.y=c("SEQ_ID", "PROBE_ID"), all.x=TRUE)
  if (all(c("POSITION.x", "POSITION.y") %in% names(ndfdata))){
    ndfdata[["POSITION.x"]] <- ndfdata[["POSITION.y"]]
    ndfdata[["POSITION.y"]] <- NULL
    ndfdata[["POSITION"]] <- ndfdata[["POSITION.x"]]
    ndfdata[["POSITION.x"]] <- NULL
  }
  rm(posdata)
  if (verbose) msgOK()
  
##   keyNdf <- paste(ndfdata[["SEQ_ID"]], ndfdata[["PROBE_ID"]], sep=":::")
##   keyPos <- paste(posdata[["SEQ_ID"]], posdata[["PROBE_ID"]], sep=":::")
##   idx <- match(keyPos, keyNdf)
##   rm(keyNdf, keyPos)
##   ndfdata[idx, "POSITION"] <- posdata[["POSITION"]]
##   ndfdata[["CHROMOSOME"]] <- NA
##   ndfdata[idx, "CHROMOSOME"] <- posdata[["CHROMOSOME"]]
##   rm(idx, posdata)

  #######################################################################
  ## Step 3.1: Get XYS files and remove all controls (ie, NA in XYS)
  #######################################################################
  if (verbose) msgParsingFile(xysFile)
  xysdata <- read.delim(xysFile, comment="#")
  if (verbose) msgOK()
  xysdata[["fid"]] <- 1:nrow(xysdata)
  ndfdata <- merge(ndfdata, xysdata, by.x=c("X", "Y"), by.y=c("X", "Y"))
  controls <- which(is.na(ndfdata[["SIGNAL"]]))
  ndfdata <- ndfdata[-controls,]
  rm(controls)
  
  #######################################################################
  ## Step 4: Prepare contents for featureSet table
  ## Fields featureSet: man_fsetid, chrom, start, end, type
  #######################################################################
  theMin <- aggregate(ndfdata[["POSITION"]], by=list(SEQ_ID=ndfdata[["SEQ_ID"]]), min)
  names(theMin) <- c("man_fsetid", "start")
  theMax <- aggregate(ndfdata[["POSITION"]], by=list(SEQ_ID=ndfdata[["SEQ_ID"]]), max)
  names(theMax) <- c("man_fsetid", "end")
  featureSet <- merge(theMin, theMax)
  rm(theMin, theMax)
  dups <- duplicated(ndfdata[["SEQ_ID"]])
  featureSet <- merge(featureSet,
                      ndfdata[!dups, c("SEQ_ID", "CHROMOSOME", "CONTAINER", "fsetid")],
                      by.x="man_fsetid", by.y="SEQ_ID")
  names(featureSet) <- c("man_fsetid", "start", "end", "chrom", "type", "fsetid")
  featureSet <- featureSet[, c("fsetid", "man_fsetid", "chrom", "start", "end", "type")]
  rm(dups)

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
  features <- ndfdata[, c("X", "Y", "fsetid", "POSITION", "MISMATCH", "MATCH_INDEX",
                          "CONTAINER", "CHROMOSOME", "PROBE_SEQUENCE", "fid", "PROBE_CLASS")]
  names(features) <- c("x", "y", "fsetid", "position", "mismatch",
                       "match_index", "container", "chromosome", "sequence", "fid", "class")

  geometry <- paste(max(ndfdata[["Y"]]), max(ndfdata[["X"]]), sep=";")
  rm(xysdata, ndfdata)
  ## FIX ME: Double check ordering
  features <- features[order(features[["chromosome"]], features[["position"]]),]

  ## FIX ME: This should be passed by function
  experimentalIDs <- "experimental"


  pmFeatures <- subset(features, mismatch == 0 & class %in% experimentalIDs)[, c("fid", "fsetid", "position", "x", "y")]
  pmSequence <- subset(features, mismatch == 0 & class %in% experimentalIDs)[, c("fid", "sequence")]
  pmSequence <- pmSequence[order(pmSequence[["fid"]]),]
  pmSequence <- DataFrame(fid=pmSequence[["fid"]],
                          sequence=DNAStringSet(pmSequence[["sequence"]]))
  mmFeatures <- subset(features, mismatch > 0 & class %in% experimentalIDs)
  ## add mmSequence
  if (any(mmFeatures[["mismatch"]] >= 10000))
    stop("Control probe possibly identified as Experimental")
  if (nrow(mmFeatures) > 0)
    stop("Add methods for MMs")
  rm(mmFeatures)
  bgFeatures <- subset(features, !(class %in% experimentalIDs))[, c("fid", "fsetid", "x", "y")]
  bgSequence <- subset(features, !(class %in% experimentalIDs))[, c("fid", "sequence")]
  bgSequence <- bgSequence[order(bgSequence[["fid"]]),]
  bgSequence <- DataFrame(fid=bgSequence[["fid"]],
                          sequence=DNAStringSet(bgSequence[["sequence"]]))
  rm(features, experimentalIDs)

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
setMethod("makePdInfoPackage", "NgsTilingPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {

            msgBar()
            cat("Building annotation package for Nimblegen Tiling Array\n")
            cat("NDF: ", basename(object@ndfFile), "\n")
            cat("POS: ", basename(object@posFile), "\n")
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
            parsedData <- parseNgsTrio(object@ndfFile,
                                       object@posFile,
                                       object@xysFile,
                                       verbose=!quiet)

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
                         PDINFOCLASS="NgsTilingPDInfo",
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
                          tiledRegionFeatureSetSchema[["col2type"]],
                          tiledRegionFeatureSetSchema[["col2key"]])
            dbCreateTable(conn,
                          "pmfeature",
                          tiledRegionPmFeatureSchema[["col2type"]],
                          tiledRegionPmFeatureSchema[["col2key"]])
            containsMm <- "mmFeatures" %in% names(parsedData)
            if (containsMm)
              dbCreateTable(conn,
                            "mmfeature",
                            tiledRegionMmFeatureSchema[["col2type"]],
                            tiledRegionMmFeatureSchema[["col2key"]])
            dbCreateTable(conn,
                          "bgfeature",
                          tiledRegionBgFeatureSchema[["col2type"]],
                          tiledRegionBgFeatureSchema[["col2key"]])

            dbInsertDataFrame(conn, "featureSet", parsedData[["featureSet"]],
                              tiledRegionFeatureSetSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                              tiledRegionPmFeatureSchema[["col2type"]], !quiet)
            if (containsMm)
              dbInsertDataFrame(conn, "mmfeature", parsedData[["mmFeatures"]],
                                tiledRegionMmFeatureSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "bgfeature", parsedData[["bgFeatures"]],
                              tiledRegionBgFeatureSchema[["col2type"]], !quiet)

            dbCreateTableInfo(conn, !quiet)

            ## Create indices
            dbCreateIndicesBg(conn, !quiet)
            dbCreateIndicesPm(conn, !quiet)
            dbCreateIndicesFs(conn, !quiet)
            
            dbGetQuery(conn, "VACUUM")
            dbDisconnect(conn)
            
            #######################################################################
            ## Part v) Save sequence DataFrames
            ## FIX ME: Fix ordering of the tables to match xxFeature tables
            #######################################################################
            datadir <- file.path(destDir, pkgName, "data")
            dir.create(datadir)
            pmSequence <- parsedData[["pmSequence"]]
            bgSequence <- parsedData[["bgSequence"]]
            pmSeqFile <- file.path(datadir, "pmSequence.rda")
            bgSeqFile <- file.path(datadir, "bgSequence.rda")
            if (!quiet) cat("Saving DataFrame object for PM.\n")
            save(pmSequence, file=pmSeqFile)
            if (!quiet) cat("Saving DataFrame object for BG.\n")
            save(bgSequence, file=bgSeqFile)
            if (!quiet) cat("Done.\n")
          })
