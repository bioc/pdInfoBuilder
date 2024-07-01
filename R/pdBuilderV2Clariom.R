setMethod("makePdInfoPackage", "AffyClariomSPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {

    msgBar()
    message("Building annotation package for Affymetrix Clariom S Array")
    message("PGF.........: ", basename(object@pgfFile))
    message("CLF.........: ", basename(object@clfFile))
    message("MPS.........: ", basename(object@coreMps))
    message("Transcript..: ", basename(object@transFile))
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

    ## Affy is mixing up the man_fsetid and fsetids, so fix
    mps <- mpsParser(object@coreMps)
    matcher <- match(mps$fsetid, parsedData$featureSet$fsetid)
    parsedData$featureSet$man_fsetid[matcher] <- mps$transcript_cluster_id
    

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

    dbExecute(conn, "VACUUM")
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
    if (!quiet) message("Saving NetAffx Annotation... ", appendLF=FALSE)
    netaffxTranscript <- annot2fdata(object@transFile)
    save(netaffxTranscript, file=file.path(extdataDir,
                                           'netaffxTranscript.rda'), compress='xz')
    if (!quiet) msgOK()

    

    if (!quiet) message("Done.")
})
