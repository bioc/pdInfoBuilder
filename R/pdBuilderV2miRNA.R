miRNAFeatureSetSchema <- list(col2type=c(
                              fsetid="INTEGER",
                              man_fsetid="TEXT",
                              type="INTEGER"),
                              col2key=c(
                              fsetid="PRIMARY KEY",
                              type ="REFERENCES type_dict(type)"
                              ))


getTypeSchema <- function(probesets) {
     tab <- table(probesets[["probeset_type"]])
     mainloc <- grep("main", names(tab))
     data.frame(type = seq_len(length(tab)),
                type_id = names(tab)[c(mainloc, seq_len(length(tab))[-mainloc])],
                stringsAsFactors = FALSE)
 }

getLevelSchema <- function()
    data.frame(level=as.integer(1:5),
               level_id=c("core", "extended", "full", "free", "ambiguous"),
               stringsAsFactors=FALSE)

setMethod("makePdInfoPackage", "AffyMiRNAPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {

            msgBar()
            message("Building annotation package for Affymetrix miRNA Array")
            message("PGF.........: ", basename(object@pgfFile))
            message("CLF.........: ", basename(object@clfFile))
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
            parsedData <- combinePgfClfProbesetsMps(object@pgfFile, object@clfFile, WT=FALSE)

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
                         PDINFOCLASS='AffyExpressionPDInfo',
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
                          "type_dict",
                          typeDictTable[["col2type"]],
                          typeDictTable[["col2key"]])
            ## end adding

            dbCreateTable(conn,
                          "featureSet",
                          miRNAFeatureSetSchema[["col2type"]],
                          miRNAFeatureSetSchema[["col2key"]])

            ## Identical structures as Gene ST
            dbCreateTable(conn,
                          "pmfeature",
                          genePmFeatureSchema[["col2type"]],
                          genePmFeatureSchema[["col2key"]])
            
            containsMm <- nrow(parsedData[["mmFeatures"]]) > 0
            if (containsMm)
                dbCreateTable(conn,
                              "mmfeature",
                              genePmFeatureSchema[["col2type"]],
                              genePmFeatureSchema[["col2key"]])
            ## end creating tables
            

            ## Inserting data in new tables
            dbInsertDataFrame(conn, "type_dict", parsedData[["type_dict"]],
                              typeDictTable[["col2type"]], !quiet)

            dbInsertDataFrame(conn, "featureSet", parsedData[["featureSet"]],
                              miRNAFeatureSetSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                              genePmFeatureSchema[["col2type"]], !quiet)
            if (containsMm)
                dbInsertDataFrame(conn, "mmfeature", parsedData[["mmFeatures"]],
                                  genePmFeatureSchema[["col2type"]], !quiet)
            ## end inserting

            dbCreateTableInfo(conn, !quiet)

            ## Create indices
            ## dbCreateIndicesPm(conn, !quiet)
            dbCreateIndex(conn, "idx_pmfsetid", "pmfeature", "fsetid", FALSE, verbose=!quiet)
            dbCreateIndex(conn, "idx_pmfid", "pmfeature", "fid", FALSE, verbose=!quiet)
            if (containsMm) {
                dbCreateIndex(conn, "idx_mmfsetid", "pmfeature", "fsetid", FALSE, verbose=!quiet)
                dbCreateIndex(conn, "idx_mmfid", "pmfeature", "fid", FALSE, verbose=!quiet)
            }
##            dbCreateIndicesMm(conn, !quiet)
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
            pmSeqFile <- file.path(datadir, "pmSequence.rda")
            if (!quiet) message("Saving DataFrame object for PM.")
            save(pmSequence, file=pmSeqFile, compress='xz')
            mmSequence <- parsedData[["mmSequence"]]
            mmSeqFile <- file.path(datadir, "mmSequence.rda")
            if (!quiet) message("Saving DataFrame object for MM.")
            save(mmSequence, file=mmSeqFile, compress='xz')

            #######################################################################
            ## Part vi) Save NetAffx Annotation to extdata
            #######################################################################
            ## miRNA arrays do not have NetAffx Annotation files

            if (!quiet) message("Done.")
          })
