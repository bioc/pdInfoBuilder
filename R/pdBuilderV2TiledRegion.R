## This is a major change to pdInfoBuilder
## We're going to remove any biological metadata
## Such info will come from other packages

## Please don't change.

## For the moment, this file contains the
## initialization struct for the db to be
## used with NimbleGen arrays ("Tiled Regions")

## The following tables will be available:
## pmfeature: fid, fsetid, chr, location, x, y
## featureSet: fsetid (region name), chr, start, end
## mmfeature (opt): fid, fid (matching pm), x, y
## bgfeature: fid, x, y


#######################################################################
## SECTION A - Db Schema
#######################################################################
tiledRegionFeatureSetSchema <- list(col2type=c(
                                      fsetid="INTEGER",
                                      man_fsetid="TEXT",
                                      chr="TEXT",
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
    cat("Inserting ", nrow(data), " rows into table \"", tablename, "\"... ", sep="")
  dbBeginTransaction(conn)
  on.exit(dbCommit(conn))
  dbGetPreparedQuery(conn, sql_template, bind.data=data)
  if (verbose)
    cat("OK\n")
}

#######################################################################
## SECTION C - Parser specific for NGS Tiled Regions
##             This will take NDF/POS/XYS trio and process (in memory)
##             the data, generating data.frames for each table
##             (featureSet, pmfeature, mmfeature*, bgfeature**)
##             to be created in the db.
#######################################################################
parseNgsTrio <- function(ndfFile, posFile, xysFile){
  stopifnot(!missing(ndfFile), !missing(posFile), !missing(xysFile))

  #######################################################################
  ## Step 1: Parse NDF
  #######################################################################
  ndfdata <- read.delim(ndfFile, stringsAsFactors=FALSE)
  ndfdata[["fsetid"]] <- as.integer(as.factor(ndfdata[["SEQ_ID"]]))

  #######################################################################
  ## Step 2: Parse POS
  #######################################################################
  posdata <- read.delim(posFile, stringsAsFactors=FALSE)
  
  #######################################################################
  ## Step 3: Match POS and NDF
  ##         and update positional info
  #######################################################################
  keyNdf <- paste(ndfdata[["SEQ_ID"]], ndfdata[["PROBE_ID"]], sep=":::")
  keyPos <- paste(posdata[["SEQ_ID"]], posdata[["PROBE_ID"]], sep=":::")
  idx <- match(keyPos, keyNdf)
  rm(keyNdf, keyPos)
  ndfdata[idx, "POSITION"] <- posdata[["POSITION"]]
  ndfdata[["CHROMOSOME"]] <- NA
  ndfdata[idx, "CHROMOSOME"] <- posdata[["CHROMOSOME"]]
  rm(idx, posdata)

  #######################################################################
  ## Step 4: Prepare contents for featureSet table
  ## Fields featureSet: man_fsetid, chr, start, end, type
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
  names(featureSet) <- c("man_fsetid", "start", "end", "chr", "type", "fsetid")
  featureSet <- featureSet[, c("fsetid", "man_fsetid", "chr", "start", "end", "type")]
  rm(dups)

  #######################################################################
  ## Step 5: Prepare contents for pmfeature, mmfeature and bgfeature
  ## Fields pmfeature: fid, fsetid, position, x, y
  ##        mmfeature: fid, fsetid, position, x, y, fid(PM)
  ##        bgfeature: fid, fsetid, x, y
  ## FIX ME: Need to check the selection for BG probes
  ##         Implement mmfeature
  ## Comments: ideally, the tables would be ordered by chr+position;
  ##           but 'fid' is the PRIMARY KEY, and the table will be
  ##           ordered by that
  #######################################################################
  xysdata <- read.delim(xysFile, comment="#")
  xysdata[["fid"]] <- 1:nrow(xysdata)
  features <- merge(ndfdata[, c("X", "Y", "fsetid", "POSITION", "MISMATCH",
                                "MATCH_INDEX", "CONTAINER", "CHROMOSOME")],
                    xysdata[, c("X", "Y", "fid")],
                    by.x=c("X", "Y"), by.y=c("X", "Y"))
  names(features) <- c("x", "y", "fsetid", "position", "mismatch",
                       "match_index", "container", "chromosome", "fid")
  geometry <- paste(max(ndfdata$Y), max(ndfdata$X), sep=";")
  rm(xysdata, ndfdata)
  features <- features[order(features[["chromosome"]], features[["position"]]),]
  pmFeatures <- subset(features, mismatch == 0)[, c("fid", "fsetid", "position", "x", "y")]
  mmFeatures <- subset(features, mismatch > 0 & mismatch < 10000)
  if (nrow(mmFeatures) > 0)
    stop("Add methods for MMs")
  rm(mmFeatures)
  bgFeatures <- subset(features, mismatch >= 10000)[, c("fid", "fsetid", "x", "y")]
  rm(features)
  return(list(featureSet=featureSet, pmFeatures=pmFeatures, bgFeatures=bgFeatures, geometry=geometry))
}

#######################################################################
## SECTION D - Package Maker
##             This shouldn't be extremely hard.
##             The idea is to: i) get array info (name, pkgname, dbname)
##             ii) parse data; iii) create pkg from template;
##             iv) dump the database
#######################################################################
setMethod("makePdInfoPackage", "NgsTiledRegionPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {
            
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
                                       object@xysFile)

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
            dbCreateTable(conn, "featureSet", tiledRegionFeatureSetSchema[["col2type"]], tiledRegionFeatureSetSchema[["col2key"]])
            dbCreateTable(conn, "pmfeature", tiledRegionPmFeatureSchema[["col2type"]], tiledRegionPmFeatureSchema[["col2key"]])
            containsMm <- "mmFeatures" %in% names(parsedData)
            if (containsMm)
              dbCreateTable(conn, "mmfeature", tiledRegionMmFeatureSchema[["col2type"]], tiledRegionMmFeatureSchema[["col2key"]])
            dbCreateTable(conn, "bgfeature", tiledRegionBgFeatureSchema[["col2type"]], tiledRegionBgFeatureSchema[["col2key"]])

            dbInsertDataFrame(conn, "featureSet", parsedData[["featureSet"]],
                              tiledRegionFeatureSetSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                              tiledRegionPmFeatureSchema[["col2type"]], !quiet)
            if (containsMm)
              dbInsertDataFrame(conn, "mmfeature", parsedData[["mmFeatures"]],
                                tiledRegionMmFeatureSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "bgfeature", parsedData[["bgFeatures"]],
                              tiledRegionBgFeatureSchema[["col2type"]], !quiet)
            dbGetQuery(conn, "VACUUM")
            dbDisconnect(conn)
          })
