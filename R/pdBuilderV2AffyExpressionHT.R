#######################################################################
## SECTION A - Db Schema
#######################################################################
affyHTExpressionFeatureSetSchema <- list(col2type=c(
                                           fsetid="INTEGER",
                                           strand="INTEGER",
                                           man_fsetid="INTEGER"),
                                         col2key=c(
                                           fsetid="PRIMARY KEY"
                                           ))

affyHTExpressionPmFeatureSchema <- list(col2type=c(
                                          fid="INTEGER",
                                          fsetid="INTEGER",
                                          x="INTEGER",
                                          y="INTEGER"
                                          ),
                                        col2key=c(
                                          fid="PRIMARY KEY"
                                          ))

affyHTExpressionMmFeatureSchema <- list(col2type=c(
                                          fid="INTEGER",
                                          fsetid="INTEGER",
                                          x="INTEGER",
                                          y="INTEGER",
                                          fidpm="INTEGER"
                                          ),
                                        col2key=c(
                                          fid="PRIMARY KEY"
                                          ))

affyHTExpressionBgFeatureSchema <- list(col2type=c(
                                          fid="INTEGER",
                                          fsetid="INTEGER",
                                          x="INTEGER",
                                          y="INTEGER"
                                          ),
                                        col2key=c(
                                          fid="PRIMARY KEY"
                                          ))


parseCdfCelProbe <- function(cdfFile, celFile, probeFile, verbose=TRUE){
  if (verbose) msgParsingFile(cdfFile)
  cdf <- readCdf(cdfFile)
  if (verbose) msgOK()

  if (verbose) msgParsingFile(celFile)
  cel <- readCelHeader(celFile)
  if (verbose) msgOK()
  geometry <- c(cel[["rows"]], cel[["cols"]])
  rm(cel)

  if (verbose) msgParsingFile(probeFile)
  cols <- c("probe.x", "probe.y", "probe.sequence")
  probeSeq <- read.delim(probeFile, stringsAsFactors=FALSE)[, cols]
  rm(cols)
  names(probeSeq) <- c("x", "y", "sequence")
  if (verbose) msgOK()

  strands <- sapply(cdf, "[[", "unitdirection")
  strands <- ifelse(tolower(strands) == "sense",
                    as.integer(SENSE),
                    as.integer(ANTISENSE))
  if (verbose) cat("Getting information for featureSet table... ")
  featureSet <- data.frame(fsetid=1:length(strands),
                           man_fsetid=names(strands),
                           strand=strands,
                           stringsAsFactors=FALSE)
  rm(strands)
  if (verbose) msgOK()

  extractFromGroups <- function(x){
    ## x is a list and has "groups" as component
    ngroups <- length(x[["groups"]])
    natoms <- sapply(x[["groups"]], "[[", "natoms") * sapply(x[["groups"]], "[[", "ncellsperatom")
    mfsetid <- unlist(mapply('rep', names(x[["groups"]]), each=natoms))
    mfsetid <- as.character(mfsetid)
    probes <- lapply(x[["groups"]],
                     function(y){
                       data.frame(x=y[["x"]],
                                  y=y[["y"]],
                                  isPm=!(y[["tbase"]]==y[["pbase"]]),
                                  atom=y[["atom"]])
                     })
    probes <- do.call("rbind", probes)
    probes[["man_fsetid"]] <- mfsetid
    return(probes)
  }

  xy2i <- function(x, y, geom)
    as.integer(geom[1]*y+x+1)

  if (verbose) cat("Getting information for pm/mm feature tables ... ")
  allProbes <- lapply(cdf, extractFromGroups)
  allProbes <- do.call("rbind", allProbes)
  allProbes[["fid"]] <- xy2i(allProbes[["x"]], allProbes[["y"]], geometry)
  allProbes[["fsetid"]] <- featureSet[match(allProbes[["man_fsetid"]],
                                            featureSet[["man_fsetid"]]),
                                      "fsetid"]
  if (verbose) cat("OK\n")
  if (verbose) cat("Combining probe information with sequence information ... ")
  allProbes <- merge(allProbes, probeSeq,
                     by.x=c("x", "y"),
                     by.y=c("x", "y"),
                     all.x=TRUE)
  rm(probeSeq)
  if (verbose) msgOK()

  ## AFFX probes have sequence for only 1 probe of the probeset...
  if (verbose) cat("Getting sequence information for AFFX probes ...")
  missSeq <- which(is.na(allProbes[["sequence"]]) & allProbes[["isPm"]])
  missPS <- sort(unique(allProbes[missSeq, "fsetid"]))
  for (i in missPS){
    idx <- which(allProbes[["fsetid"]] == i & allProbes[["isPm"]])
    seq <- subset(allProbes, fsetid == i & !is.na(sequence) & isPm, sequence, drop=TRUE)
##     seq <- subset(allProbes, fsetid == i & !is.na(sequence) & isPm)[, "sequence"]
##     toGet <- (allProbes[["fsetid"]] == i) & (!is.na(allProbes[["sequence"]])) & (allProbes[["isPm"]])
##     seq <- allProbes[toGet, "sequence"]
    seq <- unique(seq)
    stopifnot(length(seq) == 1)
    allProbes[idx, "sequence"] <- seq
    rm(idx, seq)
  }
  rm(missSeq, missPS, i)
  
  idx <- grep("nonspecific", tolower(allProbes[["man_fsetid"]]))
  bgFeatures <- allProbes[idx,]
  allProbes <- allProbes[-idx,]
  if (any(!bgFeatures[["isPm"]]))
    warning("MM Background probes were ignored.")
  bgFeatures <- bgFeatures[bgFeatures[["isPm"]],]
  rm(idx)
  bgSequence <- XDataFrame(fid=bgFeatures[["fid"]],
                           sequence=bgFeatures[["sequence"]])
  bgFeatures <- bgFeatures[, c("fid", "fsetid", "x", "y")]
  if (verbose) msgOK()

  geometry <- paste(geometry, collapse=";")
  cols <- c("fid", "fsetid", "x", "y", "atom")
  cols2 <- c("fid", "sequence")
  pmidx <- which(allProbes[["isPm"]])
  pmFeatures <- allProbes[pmidx, cols]
  pmSequence <- allProbes[pmidx, cols2]

  missSeq <- which(is.na(pmSequence[["sequence"]]))
  if (any(is.na(pmSequence[["sequence"]])))
    cat("Problem with sequences. Check pmSequence for missing values.")

  pmSequence <- XDataFrame(fid=pmSequence[["fid"]],
                           sequence=DNAStringSet(pmSequence[["sequence"]]))
  mmFeatures <- allProbes[-pmidx, cols]
  rm(pmidx, allProbes)

  cols1 <- c("fsetid", "atom")
  cols2 <- c("fsetid", "fid", "atom")
  matchpm <- merge(mmFeatures[, cols2], pmFeatures[, cols2],
               by.x=cols1, by.y=cols1)[, c("fid.x", "fid.y")]
  rm(cols1, cols2)
  names(matchpm) <- c("fid", "fidpm")

  mmFeatures <- merge(mmFeatures, matchpm, by.x="fid", by.y="fid")
  mmFeatures[["atom"]] <- NULL
  pmFeatures[["atom"]] <- NULL
  return(list(featureSet=featureSet,
              pmSequence=pmSequence,
              pmFeatures=pmFeatures,
              mmFeatures=mmFeatures,
              geometry=geometry,
              bgFeatures=bgFeatures,
              bgSequence=bgSequence))
}

#######################################################################
## SECTION D - Package Maker
##             This shouldn't be extremely hard.
##             The idea is to: i) get array info (name, pkgname, dbname)
##             ii) parse data; iii) create pkg from template;
##             iv) dump the database
#######################################################################

setMethod("makePdInfoPackage", "AffyExpressionPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {

            msgBar()
            cat("Building annotation package for Affymetrix Expression array\n")
            cat("CDF...............: ", basename(object@cdfFile), "\n")
            cat("CEL...............: ", basename(object@celFile), "\n")
            cat("Sequence TAB-Delim: ", basename(object@tabSeqFile), "\n")
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
            parsedData <- parseCdfCelProbe(object@cdfFile,
                                           object@celFile,
                                           object@tabSeqFile,
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
                         PDINFOCLASS="AffyExpressionPDInfo",
                         GEOMETRY=parsedData[["geometry"]])
            templateDir <- system.file("pd.PKG.template",
                                       package="pdInfoBuilder")
            createPackage(pkgname=pkgName, destinationDir=destDir,
                          originDir=templateDir, symbolValues=syms,
                          quiet=quiet)
            dir.create(extdataDir, recursive=TRUE)

            #######################################################################
            ## Part iv) Create SQLite database
            ## FIX ME: Fix ordering of the tables
            #######################################################################
            conn <- dbConnect(dbDriver("SQLite"), dbname=dbFilePath)

            dbCreateTable(conn,
                          "featureSet",
                          affyHTExpressionFeatureSetSchema[["col2type"]],
                          affyHTExpressionFeatureSetSchema[["col2key"]])
            
            dbCreateTable(conn,
                          "pmfeature",
                          affyHTExpressionPmFeatureSchema[["col2type"]],
                          affyHTExpressionPmFeatureSchema[["col2key"]])

            dbCreateTable(conn,
                          "mmfeature",
                          affyHTExpressionMmFeatureSchema[["col2type"]],
                          affyHTExpressionMmFeatureSchema[["col2key"]])
            dbCreateTable(conn,
                          "bgfeature",
                          affyHTExpressionBgFeatureSchema[["col2type"]],
                          affyHTExpressionBgFeatureSchema[["col2key"]])

            dbInsertDataFrame(conn, "featureSet", parsedData[["featureSet"]],
                              affyHTExpressionFeatureSetSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                              affyHTExpressionPmFeatureSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "mmfeature", parsedData[["mmFeatures"]],
                              affyHTExpressionMmFeatureSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "bgfeature", parsedData[["bgFeatures"]],
                              affyHTExpressionBgFeatureSchema[["col2type"]], !quiet)
            dbGetQuery(conn, "VACUUM")

            dbCreateTableInfo(conn, !quiet)

            ## Create indices
            dbCreateIndicesBg(conn, !quiet)
            dbCreateIndicesPm(conn, !quiet)
            dbCreateIndicesFs(conn, !quiet)
            
            dbDisconnect(conn)
            
            #######################################################################
            ## Part v) Save sequence XDataFrames
            ## FIX ME: Fix ordering of the tables to match xxFeature tables
            #######################################################################
            datadir <- file.path(destDir, pkgName, "data")
            dir.create(datadir)
            pmSequence <- parsedData[["pmSequence"]]
            bgSequence <- parsedData[["bgSequence"]]
            pmSeqFile <- file.path(datadir, "pmSequence.rda")
            bgSeqFile <- file.path(datadir, "bgSequence.rda")
            if (!quiet) cat("Saving XDataFrame object for PM.\n")
            save(pmSequence, file=pmSeqFile)
            if (!quiet) cat("Saving XDataFrame object for BG.\n")
            save(bgSequence, file=bgSeqFile)
            if (!quiet) cat("Done.")
          })
