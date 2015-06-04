## tables
htaMpsSchema <- list(col2type=c(
                         meta_fsetid="TEXT",
                         transcript_cluster_id="TEXT",
                         fsetid="INTEGER"
                         ),
                     col2key=c(
                         fsetid="REFERENCES featureSet(fsetid)"
                         ))


htaFeatureSetSchema <- list(col2type=c(
                                fsetid="INTEGER",
                                man_fsetid="TEXT",
                                strand="INTEGER",
                                start="INTEGER",
                                stop="INTEGER",
                                transcript_cluster_id="INTEGER",
                                exon_id="INTEGER",
                                crosshyb_type="INTEGER",
                                level="INTEGER",
                                junction_start_edge="INTEGER",
                                junction_stop_edge="INTEGER",
                                junction_sequence="TEXT",
                                has_cds="INTEGER",
                                chrom="INTEGER",
                                type="INTEGER"),
                            col2key=c(
                                fsetid="PRIMARY KEY",
                                chrom="REFERENCES chrom_dict(chrom_id)",
                                level="REFERENCES level_dict(level_id)",
                                type ="REFERENCES type_dict(type_id)"
                                ))

###

parseHtaProbesetCSV <- function(probeFile, verbose=TRUE){
  ## Variables
  SENSE <- as.integer(0)
  ANTISENSE <- as.integer(1)

###################################################
## TABLES TO ADD
###################################################
  if (verbose) simpleMessage("Creating dictionaries... ")
  ## chromosome dictionary moved to after reading
  ##  the CSV file

  ## level_schema
  level_schema <- getLevelSchema()

  ## type_schema
  ## type_schema <- getTypeSchema()
  ## do this on the fly, below
  if (verbose) msgOK()


  ## the "probesets" df is to be the featureSet table
  if (verbose) msgParsingFile(probeFile)
  probesets <- read.csv(probeFile, comment.char="#",
                        stringsAsFactors=FALSE, na.strings="---")
  if (verbose) msgOK()

  ## added: "junction_start_edge", "junction_stop_edge",
  ## "junction_sequence", "has_cds"
  cols <- c("probeset_id", "seqname", "strand", "start", "stop",
            "transcript_cluster_id", "exon_id",
            "crosshyb_type", "level", "probeset_type",
            "junction_start_edge", "junction_stop_edge",
            "junction_sequence", "has_cds")
  probesets <- probesets[, cols]
  cols[1] <- "man_fsetid"
  names(probesets) <- cols
  rm(cols)
##  probesets$fsetid <- as.integer(factor(probesets$man_fsetid))

  ## type_schema
  type_schema <- getTypeSchema(probesets)

  chromosome_schema <- createChrDict(probesets[["seqname"]])

  probesets[["chrom"]] <- match(probesets[["seqname"]], chromosome_schema[["chrom_id"]])
  probesets[["seqname"]] <- NULL
  probesets[["strand"]] <- ifelse(probesets[["strand"]] == "-", ANTISENSE, SENSE)
  probesets[["level"]] <- match(tolower(probesets[["level"]]), level_schema[["level_id"]])
  probesets[["type"]] <- match(probesets[["probeset_type"]], type_schema[["type_id"]])
  probesets[["probeset_type"]] <- NULL

  list(probesets=probesets, level=level_schema, chromosome=chromosome_schema, type=type_schema)
}

combinePgfClfProbesetsMpsHTA <- function(pgfFile, clfFile, probeFile,
                                         coreMps, verbose=TRUE){
    tmp <- parsePgfClf(pgfFile=pgfFile, clfFile=clfFile, verbose=verbose)
    probes.table <- tmp[["probes.table"]]
    geom <- tmp[["geometry"]]
    rm(tmp)
    
    probesetInfo <- parseHtaProbesetCSV(probeFile, verbose=verbose)
    
    ## levels table
    ## id
    ## desc
    level_dict <- probesetInfo[["level"]]
    
    ## chromosome table
    ## id
    ## chrom_id
    chrom_dict <- probesetInfo[["chromosome"]]
    
    ## types table
    ## id
    ## type_id
    type_dict <- probesetInfo[["type"]]
    
    ## featureSet table - Fields
    ## probeset_id
    ## strand
    ## start
    ## stop
    ## transcript_cluster_id
    ## exon_id
    ## crosshyb_type
    ## level
    ## junction_start_edge
    ## junction_stop_edge
    ## junction_sequence
    ## chrom
    ## type

    ## Starting with the MTA 2.0 array, there are (811) PSR and JUC probesets
    ## that are part of the design, but because they are targeting regions that
    ## are too repetitive to synthesize a probe, or too short for a probe (25-mer) to fit
    ## they don't actually have these probesets on the array. But they kept them in the
    ## design, so they show up in the csv files. This causes problems when we create the featureSet
    ## data.frame, because we are doing a merge with all = TRUE, so we get NA values for the
    ## fsetids that don't exist, but for which we do have man_fsetids. These NA values
    ## are problematic when we insert into the database, as the fsetids have to be unique.
    ## So we remove them here.

    if(any(!probesetInfo$probesets$man_fsetid %in% probes.table$man_fsetid)){
        ind <- probesetInfo$probesets$man_fsetid %in% probes.table$man_fsetid
        probesetInfo[["probesets"]] <- probesetInfo[["probesets"]][ind,]
    }

    featureSet <- merge(probesetInfo[["probesets"]],
                        unique(probes.table[, c('man_fsetid', 'fsetid')]),
                        by='man_fsetid', all=TRUE)

    missFeatureSet <- setdiff(unique(probes.table[["man_fsetid"]]),
                              unique(featureSet[["man_fsetid"]]))
    if (length(missFeatureSet) > 0){
        missFS = data.frame(man_fsetid=missFeatureSet)
        cols <- names(featureSet)
        cols <- cols[cols != "man_fsetid"]
        for (i in cols)
            missFS[[i]] <- NA
        missFS <- missFS[, names(featureSet)]
        featureSet <- rbind(featureSet, missFS)
        rm(missFS, cols, i)
    }
    rm(missFeatureSet)
    
    ## pmfeature table - Fields
    ##  fid
    ##  fsetid
    ##  chr (NA)
    ##  location (NA)
    ##  x
    ##  y
    ## IMPORTANT:
    ##    ignoring strand
    ##    keeping atom to match with MM's
    pmFeatures <- subset(probes.table,
                         substr(probes.table[["ptype"]], 1, 2) == "pm",
                         select=c("fid", "fsetid", "atom", "x", "y", "sequence"))
    
    pmSequence <- pmFeatures[, c("fid", "sequence")]
    pmFeatures[["sequence"]] <- NULL
    pmSequence <- pmSequence[order(pmSequence[["fid"]]),]
    pmSequence <- DataFrame(fid=pmSequence[["fid"]],
                            sequence=DNAStringSet(pmSequence[["sequence"]]))
    
    ## mmfeature table - Fields
    ##  fid
    ##  fid of matching pm
    ##  x
    ##  y
    ## IMPORTANT:
    ##    ignoring strand
    ##    keeping atom to match with MM's
    ##    ADD sequence for MM
    mmFeatures <- subset(probes.table, substr(probes.table$ptype, 1, 2) =="mm",
                         select=c("fid", "fsetid", "atom", "x", "y", "sequence"))
    if (nrow(mmFeatures) > 0){
        mmSequence <- mmFeatures[, c("fid", "sequence")]
        mmFeatures[["sequence"]] <- NULL
        mmSequence <- mmSequence[order(mmSequence[["fid"]]),]
        mmSequence <- DataFrame(fid=mmSequence[["fid"]],
                                sequence=DNAStringSet(mmSequence[["sequence"]]))
    }else{
        mmFeatures <- data.frame()
        mmSequence <- data.frame()
    }
    
    ## IMPORTANT: for the moment, bgfeature will contain everything (that is PM) but 'main'
    ## bgfeature table - Fields
    ##  fid
    ##  x
    ##  y
    ##  fs_type: featureSet type: genomic/antigenomic
    ##  f_type: pm/mm at/st
    ## old code:
    ## subset using cols
    ## cols <- c("fid", "fsetid", "pstype", "ptype", "x", "y", "sequence")
    rm(probes.table)
    
    core <- mpsParser(coreMps, verbose=verbose)
        
    ## Here we should have the following tables available:
    ##  featureSet: fsetid, type
    ##  pmfeature: fid, fsetid, atom, x, y
    ##  bgfeature: fid, fsetid, fs_type, f_type, x, y  - NOT ANYMORE
    ##  pmSequence: fid, sequence
    ##  bgSequence: fid, sequence  - NOT ANYMORE
    ##  core, extended, full: meta_fsetid, trancript_cluster_id, fsetid
    ##  mmfeatures/mmSequence
    
    out <- list(featureSet=featureSet, pmFeatures=pmFeatures,
                mmFeatures=mmFeatures, geometry=geom,
                pmSequence=pmSequence, mmSequence=mmSequence,
                chrom_dict=chrom_dict, level_dict=level_dict,
                type_dict=type_dict, core=core)
    return(out)
}

#######################################################################
## SECTION D - Package Maker
##             This shouldn't be extremely hard.
##             The idea is to: i) get array info (name, pkgname, dbname)
##             ii) parse data; iii) create pkg from template;
##             iv) dump the database
#######################################################################

setMethod("makePdInfoPackage", "AffyHTAPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {
            msgBar()
            message("Building annotation package for Affymetrix HTA Array")
            message("PGF.........: ", basename(object@pgfFile))
            message("CLF.........: ", basename(object@clfFile))
            message("Probeset....: ", basename(object@probeFile))
            message("Transcript..: ", basename(object@transFile))
            message("Core MPS....: ", basename(object@coreMps))
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
            parsedData <- combinePgfClfProbesetsMpsHTA(object@pgfFile,
                                                       object@clfFile,
                                                       object@probeFile,
                                                       object@coreMps,
                                                       verbose=!quiet)

            #######################################################################
            ## Part iii) Create package from template
            #######################################################################
            pdInfoClass <- "AffyHTAPDInfo"
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
                         PDINFOCLASS=pdInfoClass,
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
            dbCreateTable(conn,
                          "level_dict",
                          levelDictTable[["col2type"]],
                          levelDictTable[["col2key"]])
            dbCreateTable(conn,
                          "type_dict",
                          typeDictTable[["col2type"]],
                          typeDictTable[["col2key"]])
            dbCreateTable(conn,
                          "core_mps",
                          htaMpsSchema[["col2type"]],
                          htaMpsSchema[["col2key"]])
            ## end adding

            dbCreateTable(conn,
                          "featureSet",
                          htaFeatureSetSchema[["col2type"]],
                          htaFeatureSetSchema[["col2key"]])

            dbCreateTable(conn, "pmfeature",
                          genePmFeatureSchema[["col2type"]],
                          genePmFeatureSchema[["col2key"]])

            containsMm <- nrow(parsedData[["mmFeatures"]]) > 0
            if (containsMm)
              dbCreateTable(conn,
                            "mmfeature",
                            exonTranscriptionMmFeatureSchema[["col2type"]],
                            exonTranscriptionMmFeatureSchema[["col2key"]])

            ## Inserting data in new tables
            dbInsertDataFrame(conn, "chrom_dict", parsedData[["chrom_dict"]],
                              chromDictTable[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "level_dict", parsedData[["level_dict"]],
                              levelDictTable[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "type_dict", parsedData[["type_dict"]],
                              typeDictTable[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "core_mps", parsedData[["core"]],
                              mpsSchema[["col2type"]], !quiet)
            ## end inserting

            dbInsertDataFrame(conn, "featureSet", parsedData[["featureSet"]],
                              htaFeatureSetSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                              genePmFeatureSchema[["col2type"]], !quiet)
            if (containsMm)
              dbInsertDataFrame(conn, "mmfeature", parsedData[["mmFeatures"]],
                                exonTranscriptionMmFeatureSchema[["col2type"]], !quiet)

            dbCreateTableInfo(conn, !quiet)

            ## Create indices
            dbCreateIndex(conn, "idx_pmfsetid", "pmfeature", "fsetid", FALSE, verbose=!quiet)
            dbCreateIndex(conn, "idx_pmfid", "pmfeature", "fid", FALSE, verbose=!quiet)
            dbCreateIndicesFs(conn, !quiet)
            dbCreateIndex(conn, "idx_core_meta_fsetid", "core_mps", "meta_fsetid", FALSE, verbose=!quiet)
            dbCreateIndex(conn, "idx_core_fsetid", "core_mps", "fsetid", FALSE, verbose=!quiet)

            if (containsMm){
              dbCreateIndex(conn, "idx_mmfsetid", "mmfeature", "fsetid", FALSE, verbose=!quiet)
              dbCreateIndex(conn, "idx_mmfid", "mmfeature", "fid", FALSE, verbose=!quiet)
            }

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
            if (containsMm){
              mmSequence <- parsedData[["mmSequence"]]
              mmSeqFile <- file.path(datadir, "mmSequence.rda")
              if (!quiet) message("Saving DataFrame object for MM.")
              save(mmSequence, file=mmSeqFile, compress='xz')
            }


            #######################################################################
            ## Part vi) Save NetAffx Annotation to extdata
            #######################################################################
            if (!quiet) message("Saving NetAffx Annotation... ", appendLF=FALSE)
            netaffxProbeset <- annot2fdata(object@probeFile)
            save(netaffxProbeset, file=file.path(extdataDir,
                                  'netaffxProbeset.rda'), compress='xz')
            netaffxTranscript <- annot2fdata(object@transFile)
            save(netaffxTranscript, file=file.path(extdataDir,
                                    'netaffxTranscript.rda'), compress='xz')
            if (!quiet) msgOK()

            if (!quiet) message("Done.")
          })
