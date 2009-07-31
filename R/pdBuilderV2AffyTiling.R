#######################################################################
## SECTION A - Db Schema
#######################################################################

affyTilingPmFeatureSchema <- list(col2type=c(
                                    fid="INTEGER",
                                    chrom="INTEGER",
                                    position="INTEGER",
                                    x="INTEGER",
                                    y="INTEGER"),
                                  col2key=c(
                                    fid="PRIMARY KEY"
                                    ))

affyTilingMmFeatureSchema <- list(col2type=c(
                                    fid="INTEGER",
                                    fidpm="INTEGER",
                                    x="INTEGER",
                                    y="INTEGER"),
                                  col2key=c(
##                                    fid="PRIMARY KEY"
                                    ))

affyTilingBgFeatureSchema <- list(col2type=c(
                                    fid="INTEGER",
                                    x="INTEGER",
                                    y="INTEGER"),
                                  col2key=c(
                                    fid="PRIMARY KEY"
                                    ))

#######################################################################
## SECTION B - Utils - This should be moved from here
##             as everything in this section can be used on other cases
#######################################################################

#######################################################################
## SECTION C - Parser for BPMAP+CEL. The BPMAP file does not have
##             geometry information, that's why I use a CEL file.
##             The previous code used to use a CIF, but I believe
##             that the CEL file is more natural.
#######################################################################

parseBpmapCel <- function(bpmapFile, celFile, verbose=TRUE){
  if (verbose) msgParsingFile(bpmapFile)
  bpmap <- readBpmap(bpmapFile)
  if (verbose) msgOK()
  if (verbose) cat("Getting geometry from CEL file... ")
  celHeader <- readCelHeader(celFile)
  if (verbose) msgOK()

  ## Assuming that the "experimental" probes will be flagged as "tiling"
  ## and that background will be "arbitrary"

  getField <- function(x, field) x[["seqInfo"]][[field]]
##  experimental <- grep("^chr", sapply(bpmap, getField, "name"))
##  background <- grep("^AffxCtrlBkGr", sapply(bpmap, getField, "groupname"))
  grpNames <- sapply(bpmap, getField, "groupname")
  controls <- background <- grep("^AffxC", grpNames)
  experimental <- (1:length(grpNames))[-controls]
  
  problem <- intersect(experimental, background)
  if (length(problem) > 0)
    stop("Probes were identified as both Experimental and Background controls")
  rm(problem)
  stopifnot(length(experimental) > 0, length(background) > 0)

  experimental <- bpmap[experimental]
  background <- bpmap[background]
  rm(bpmap)

  geometry <- c(celHeader[["rows"]], celHeader[["cols"]])
  xy2i <- function(x, y, geom)
    as.integer(geom[1]*y+x+1)

  if (verbose) cat("Getting PMs... ")
  ## there are pmmm and pmonly
  idx <- which(sapply(experimental, getField, "mapping") == "onlypm")
  idx2 <- which(sapply(experimental, getField, "mapping") == "pmmm")
  cols <- c("fid", "chrom", "startpos", "pmx", "pmy", "probeseq")
  pmFeatures <- lapply(experimental,
                   function(x){
                     x[["chrom"]] <- x[["seqInfo"]][["name"]]
                     x[["fid"]] <- xy2i(x[["pmx"]], x[["pmy"]], geometry)
                     as.data.frame(x[cols], stringsAsFactors=FALSE)
                   })
  pmFeatures <- do.call("rbind", pmFeatures)
  names(pmFeatures) <- c("fid", "chrom", "position", "x", "y", "sequence")
  rownames(pmFeatures) <- NULL
  if (verbose) msgOK()

  if (length(idx2) > 0){
    if (verbose) cat("Getting MMs... ")
    cols <- c("fid", "fidpm", "mmx", "mmy")
    mmFeatures <- lapply(experimental[idx2],
                         function(x){
                           x[["fid"]] <- xy2i(x[["mmx"]], x[["mmy"]], geometry)
                           x[["fidpm"]] <- xy2i(x[["pmx"]], x[["pmy"]], geometry)
                           as.data.frame(x[cols], stringsAsFactors=FALSE)
                         })
    mmFeatures <- do.call("rbind", mmFeatures)
    names(mmFeatures) <- c("fid", "fidpm", "x", "y")
    if (verbose) msgOK()
  }
  rm(experimental, idx)

  if (verbose) cat("Getting background probes... ")
  cols <- c("fid", "pmx", "pmy", "probeseq")
  bgFeatures <- lapply(background,
                       function(x){
                         x[["fid"]] <- xy2i(x[["pmx"]], x[["pmy"]], geometry)
                         as.data.frame(x[cols], stringsAsFactors=FALSE)
                       })
  bgFeatures <- do.call("rbind", bgFeatures)
  names(bgFeatures) <- c("fid", "x", "y", "sequence")
  rownames(bgFeatures) <- NULL
  rm(background, cols)
  if (verbose) msgOK()
  
  chrom_dict <- createChrDict(pmFeatures[["chrom"]])
  
  pmFeatures[["chrom"]] <- match(pmFeatures[["chrom"]], chrom_dict[["chrom_id"]])

  ## Sequences
  if (verbose) cat("Getting sequences...")
  pmFeatures <- pmFeatures[order(pmFeatures[["fid"]]),]
  rownames(pmFeatures) <- NULL
  if (length(idx2) > 0){
    mmFeatures <- mmFeatures[order(mmFeatures[["fid"]]),]
    rownames(mmFeatures) <- NULL
  }else{
    mmFeatures <- NULL
  }
  bgFeatures <- bgFeatures[order(bgFeatures[["fid"]]),]
  rownames(bgFeatures) <- NULL
  pmSequence <- pmFeatures[, c("fid", "sequence")]
  pmSequence <-pmSequence[order(pmSequence[["fid"]]),]
  pmSequence <- DataFrame(fid=pmSequence[["fid"]],
                          sequence=DNAStringSet(pmSequence[["sequence"]]))
  bgSequence <- bgFeatures[, c("fid", "sequence")]
  bgSequence <- bgSequence[order(bgSequence[["fid"]]),]
  bgSequence <- DataFrame(fid=bgSequence[["fid"]],
                          sequence=DNAStringSet(bgSequence[["sequence"]]))

  pmFeatures[["sequence"]] <- NULL
  bgFeatures[["sequence"]] <- NULL
  if (verbose) msgOK()
  
  geometry <- paste(geometry, collapse=";")

  list(pmFeatures=pmFeatures, mmFeatures=mmFeatures,
       bgFeatures=bgFeatures, pmSequence=pmSequence,
       bgSequence=bgSequence, geometry=geometry,
       chrom_dict=chrom_dict)
}

#######################################################################
## SECTION D - Package Maker
##             This shouldn't be extremely hard.
##             The idea is to: i) get array info (name, pkgname, dbname)
##             ii) parse data; iii) create pkg from template;
##             iv) dump the database
#######################################################################

setMethod("makePdInfoPackage", "AffyTilingPDInfoPkgSeed",          

          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {

            msgBar()
            cat("Building annotation package for Affymetrix Tiling array\n")
            cat("BPMAP: ", basename(object@bpmapFile), "\n")
            cat("CEL..: ", basename(object@celFile), "\n")
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
            parsedData <- parseBpmapCel(object@bpmapFile,
                                        object@celFile,
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
                         PDINFOCLASS="AffyTilingPDInfo",
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
            increaseDbPerformance(conn)

            ## Adding new tables
            dbCreateTable(conn,
                          "chrom_dict",
                          chromDictTable[["col2type"]],
                          chromDictTable[["col2key"]])
            ## end adding
            
            dbCreateTable(conn,
                          "pmfeature",
                          affyTilingPmFeatureSchema[["col2type"]],
                          affyTilingPmFeatureSchema[["col2key"]])

            if (!is.null(parsedData[["mmFeatures"]]))
              dbCreateTable(conn,
                            "mmfeature",
                            affyTilingMmFeatureSchema[["col2type"]],
                            affyTilingMmFeatureSchema[["col2key"]])

            dbCreateTable(conn,
                          "bgfeature",
                          affyTilingBgFeatureSchema[["col2type"]],
                          affyTilingBgFeatureSchema[["col2key"]])

            ## Inserting data in new tables
            dbInsertDataFrame(conn, "chrom_dict", parsedData[["chrom_dict"]],
                              chromDictTable[["col2type"]], !quiet)
            ## end inserting

            dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                              affyTilingPmFeatureSchema[["col2type"]], !quiet)
            if (!is.null(parsedData[["mmFeatures"]]))
              dbInsertDataFrame(conn, "mmfeature", parsedData[["mmFeatures"]],
                                affyTilingMmFeatureSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "bgfeature", parsedData[["bgFeatures"]],
                              affyTilingBgFeatureSchema[["col2type"]], !quiet)
            dbGetQuery(conn, "VACUUM")

            dbCreateTableInfo(conn, !quiet)

            ## Create indices
            dbCreateIndicesBgTiling(conn, !quiet)
            dbCreateIndicesPmTiling(conn, !quiet)
            if (!is.null(parsedData[["mmFeatures"]]))
              dbCreateIndicesMm(conn, !quiet)
            
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
            if (!quiet) cat("Done.")
          })
