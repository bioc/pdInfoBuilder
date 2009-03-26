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
                                    fid="PRIMARY KEY"
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
  if (verbose) cat("Reading in ", bpmapFile, "... ")
  bpmap <- readBpmap(bpmapFile)
  if (verbose) cat("OK\nGetting geometry from CEL file... ")
  celHeader <- readCelHeader(celFile)
  if (verbose) cat("OK\n")

  ## Assuming that the "experimental" probes will be flagged as "tiling"
  ## and that background will be "arbitrary"

  getField <- function(x, field) x[["seqInfo"]][[field]]
  experimental <- grep("^chr", sapply(bpmap, getField, "name"))
  background <- grep("^AffxCtrlBkGr", sapply(bpmap, getField, "groupname"))
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

  if (verbose) cat("Getting PMs...")
  ## there are pmmm and pmonly
  idx <- which(sapply(experimental, getField, "mapping") == "onlypm")
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

  if (verbose) cat("OK\nGetting MMs...")
  cols <- c("fid", "fidpm", "mmx", "mmy")
  mmFeatures <- lapply(experimental[-idx],
                   function(x){
                     x[["fid"]] <- xy2i(x[["mmx"]], x[["mmy"]], geometry)
                     x[["fidpm"]] <- xy2i(x[["pmx"]], x[["pmy"]], geometry)
                     as.data.frame(x[cols], stringsAsFactors=FALSE)
                   })
  mmFeatures <- do.call("rbind", mmFeatures)
  names(mmFeatures) <- c("fid", "fidpm", "x", "y")
  rm(experimental, idx)

  if (verbose) cat("OK\nGetting background probes...")
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
  if (verbose) cat("OK\n")
  
  chrom_dict <- createChrDict(pmFeatures[["chrom"]])
  
  pmFeatures[["chrom"]] <- match(pmFeatures[["chrom"]], chrom_dict[["chrom_id"]])

  ## Sequences
  if (verbose) cat("Getting sequences...")
  pmFeatures <- pmFeatures[order(pmFeatures[["fid"]]),]
  rownames(pmFeatures) <- NULL
  mmFeatures <- mmFeatures[order(mmFeatures[["fid"]]),]
  rownames(mmFeatures) <- NULL
  bgFeatures <- bgFeatures[order(bgFeatures[["fid"]]),]
  rownames(bgFeatures) <- NULL
  pmSequence <- XDataFrame(fid=pmFeatures[["fid"]],
                           sequence=DNAStringSet(pmFeatures[["sequence"]]))
  bgSequence <- XDataFrame(fid=bgFeatures[["fid"]],
                           sequence=DNAStringSet(bgFeatures[["sequence"]]))

  pmFeatures[["sequence"]] <- NULL
  bgFeatures[["sequence"]] <- NULL
  if (verbose) cat("OK\n")
  
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


            message("============================================================")
            message("Building annotation package for Affymetrix Tiling array")
            message("BPMAP: ", basename(object@bpmapFile))
            message("CEL..: ", basename(object@celFile))
            message("============================================================")

            
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
            dbInsertDataFrame(conn, "mmfeature", parsedData[["mmFeatures"]],
                              affyTilingMmFeatureSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "bgfeature", parsedData[["bgFeatures"]],
                              affyTilingBgFeatureSchema[["col2type"]], !quiet)
            dbGetQuery(conn, "VACUUM")

            dbCreateTableInfo(conn, !quiet)
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
            if (!quiet) message("Saving XDataFrame object for PM.")
            save(pmSequence, file=pmSeqFile)
            if (!quiet) message("Saving XDataFrame object for BG.")
            save(bgSequence, file=bgSeqFile)
            if (!quiet) message("Done.")
          })