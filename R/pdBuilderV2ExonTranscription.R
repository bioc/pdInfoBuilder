## Define schema

#######################################################################
## SECTION A - Db Schema
#######################################################################
chromDictTable <- list(col2type=c(
                         chrom="INTEGER",
                         chrom_id="TEXT"),
                       col2key=c(
                         chrom="PRIMARY KEY"
                         ))

levelDictTable <- list(col2type=c(
                         level="INTEGER",
                         level_id="TEXT"),
                       col2key=c(
                         level="PRIMARY KEY"
                         ))

typeDictTable <- list(col2type=c(
                        type="INTEGER",
                        type_id="TEXT"),
                      col2key=c(
                        type="PRIMARY KEY"
                        ))

fset2gene <- list(col2type=c(
                    fsetid="INTEGER",
                    gid="INTEGER"),
                  col2key=c(
                    fsetid="REFERENCES exonTranscriptionFeatureSetSchema(fsetid)",
                    gid="REFERENCES gene(gid)"
                    ))

gene <- list(col2type=c(
               accession="TEXT",
               symbol="TEXT",
               gid="INTEGER"),
             col2key=c(
               gid="PRIMARY KEY"
               ))

exonTranscriptionFeatureSetSchema <- list(col2type=c(
                                            fsetid="INTEGER",
                                            strand="INTEGER",
                                            start="INTEGER",
                                            stop="INTEGER",
                                            transcript_cluster_id="INTEGER",
                                            exon_id="INTEGER",
                                            crosshyb_type="INTEGER",
                                            level="INTEGER",
                                            chrom="INTEGER",
                                            type="INTEGER"),
                                          col2key=c(
                                            fsetid="PRIMARY KEY",
                                            chrom="REFERENCES chrom_dict(chrom_id)",
                                            level="REFERENCES level_dict(level_id)",
                                            type ="REFERENCES type_dict(type_id)"
                                            ))

exonTranscriptionPmFeatureSchema <- list(col2type=c(
                                           fid="INTEGER",
                                           fsetid="INTEGER",
                                           atom="INTEGER",
                                           x="INTEGER",
                                           y="INTEGER"
                                           ),
                                         col2key=c(
                                           fid="PRIMARY KEY"
                                           ))
## TODO: exonTranscriptionMmFeatureSchema
exonTranscriptionBgFeatureSchema <- list(col2type=c(
                                           fid="INTEGER",
                                           fsetid="INTEGER",
                                           fs_type="TEXT",
                                           f_type="TEXT",
                                           x="INTEGER",
                                           y="INTEGER"
                                           ),
                                         col2key=c(
                                           fid="PRIMARY KEY"
                                           ))


#######################################################################
## SECTION B - Utils - Already done.
#######################################################################
createChrDict <- function(x){
  possible <- as.character(na.omit(unique(x)))
  possible <- gsub("chr([1-9])$", "chr0\\1", possible)
  possible <- gsub("chr([1-9]_)", "chr0\\1", possible)
  dataSplit <- strsplit(possible, "_")
  len <- sapply(dataSplit, length)
  idx <- which(len == 1)
  basic <- unlist(dataSplit[idx])
  if (length(idx) < length(len)){
    suffixes <- sapply(dataSplit[-idx], function(x) paste(x[-1], collapse="_"))
    suffixes <- sort(unique(suffixes))
    out <- list()
    out[[1]] <- basic
    for (i in 1:length(suffixes))
      out[[i+1]] <- paste(basic, suffixes[i], sep="_")
    out <- unlist(out)
  }else{
    out <- sort(basic)
  }
  out <- gsub("chr0([1-9])", "chr\\1", out)
  data.frame(chrom=as.integer(1:length(out)),
             chrom_id=out,
             stringsAsFactors=FALSE)
}


#######################################################################
## SECTION C - Parser specific for Affy Exon array
##             This will take PGF/CLF and process (in memory)
##             the data, generating data.frames for each table
##             (featureSet, pmfeature, mmfeature*, bgfeature**)
##             to be created in the db.
#######################################################################

.idxToIdx <- function(pgfFrom, pgfTo, ids) {
  starts <- pgfFrom[ids]
  ends <- pgfFrom[ids+1] - 1
  ends[is.na(ends)] <- length(pgfTo)
  mapply(":", starts, ends, SIMPLIFY=FALSE)
}

.combineIdx <- function(psToAtomIdx, atomToProbeIdx) {
  probesPerAtom <- with(atomToProbeIdx,
                        sapply(split(atomIdx, atomIdx), length))
  cbind(probesetIdx=rep(psToAtomIdx$probesetIdx, probesPerAtom),
        atomToProbeIdx)
}

probesetIdxToAtomIdx <- function(pgf, probesetIdx) {
  atoms <- .idxToIdx(pgf[["probesetStartAtom"]],
                     pgf[["atomId"]], probesetIdx)
  data.frame(probesetIdx=rep(probesetIdx, sapply(atoms, length)),
             atomIdx=unlist(atoms))
}

atomIdxToProbeIdx <- function(pgf, atomIdx) {
  probes <- .idxToIdx(pgf[["atomStartProbe"]],
                      pgf[["probeId"]], atomIdx)
  data.frame(atomIdx=rep(atomIdx, sapply(probes, length)),
             probeIdx=unlist(probes))
}

probesetIdxToTripletIdx <- function(pgf, probesetIdx) {
  df1 <- probesetIdxToAtomIdx(pgf, probesetIdx)
  df2 <- atomIdxToProbeIdx(pgf, df1$atomIdx)
  .combineIdx(df1, df2)
}

parseProbesetCSV <- function(probeFile, verbose=TRUE){
  ## Variables
  SENSE <- as.integer(0)
  ANTISENSE <- as.integer(1)
  
###################################################
  ## TABLES TO ADD
###################################################
  if (verbose) cat("Creating dictionaries... ")
  ## chromosome dictionary moved to after reading
  ##  the CSV file
  
  ## level_schema
  level_schema <- data.frame(level=as.integer(1:5),
                             level_id=c("core", "extended", "full", "free", "ambiguous"),
                             stringsAsFactors=FALSE)
  
  ## type_schema
  type_schema <- data.frame(type=as.integer(1:8),
                            type_id=c("main", "control->affx",
                                      "control->chip",
                                      "control->bgp->antigenomic",
                                      "control->bgp->genomic",
                                      "normgene->exon",
                                      "normgene->intron",
                                      "rescue->FLmRNA->unmapped"),
                            stringsAsFactors=FALSE)
  if (verbose) cat("OK\n")
  
###################################################
  
  ## the "probesets" df is to be the featureSet table
  if (verbose) msgParsingFile(probeFile)
  probesets <- read.csv(probeFile, comment.char="#",
                        stringsAsFactors=FALSE, na.strings="---")
  if (verbose) msgOK()
  cols <- c("probeset_id", "seqname", "strand", "start", "stop",
            "transcript_cluster_id", "exon_id", "gene_assignment",
            "crosshyb_type", "level", "probeset_type")
  probesets <- probesets[, cols]
  rm(cols)

  chromosome_schema <- createChrDict(probesets[["seqname"]])
  
  probesets[["chrom"]] <- match(probesets[["seqname"]], chromosome_schema[["chrom_id"]])
  probesets[["seqname"]] <- NULL
  probesets[["strand"]] <- ifelse(probesets[["strand"]] == "-", ANTISENSE, SENSE)
  probesets[["level"]] <- match(tolower(probesets[["level"]]), level_schema[["level_id"]])
  probesets[["type"]] <- match(probesets[["probeset_type"]], type_schema[["type_id"]])
  probesets[["probeset_type"]] <- NULL
  
  ## probesets won't have Gene information due to multiplicity
  ## must remove "gene_assignment" later
  
  ## create one table (ps2genes) which will map probeset -> genes
  ## ps  |  gid
  ##  1  |  1
  ##  1  |  2
  ##  2  |  9 ... etc

  if (verbose) cat("Creating probeset -> gene table... ")
  ps2genes <- strsplit(probesets[["gene_assignment"]], " /// ")
  names(ps2genes) <- probesets[["probeset_id"]]
  probesets[["gene_assignment"]] <- NULL
  
  idx <- which(sapply(ps2genes, function(x) is.na(x[1])))
  ps2genes <- ps2genes[-idx]
  psids <- unlist(mapply(rep, as.integer(names(ps2genes)), sapply(ps2genes, length)))
  rm(idx)
  
  f <- function(x)
    do.call("rbind", strsplit(x, " // "))
  
  tmp <- do.call("rbind", lapply(ps2genes, f))
  ps2genes <- data.frame(probeset_id=psids, accession=tmp[,1], symbol=tmp[,2], stringsAsFactors=FALSE)
  rm(psids, tmp)
  ps2genes[["key"]] <- paste(ps2genes[["accession"]], ps2genes[["symbol"]], sep=":")
  rownames(ps2genes) <- NULL
  if (verbose) cat("OK\n")

  if (verbose) cat("Creating genes table... ")
  cols <- c("accession", "symbol", "key")
  onlyGenes <- ps2genes[!duplicated(ps2genes[["key"]]), cols]
  onlyGenes <- onlyGenes[order(onlyGenes[["accession"]], onlyGenes[["symbol"]]),]
  onlyGenes[["gid"]] <- 1:nrow(onlyGenes)
  rownames(onlyGenes) <- NULL
  rm(cols)
  
  ps2genes[["gid"]] <- match(ps2genes[["key"]], onlyGenes[["key"]])
  onlyGenes[["key"]] <- ps2genes[["accession"]] <- ps2genes[["symbol"]] <- ps2genes[["key"]] <- NULL

  tmp <- names(probesets)
  tmp[tmp == "probeset_id"] <- "fsetid"
  names(probesets) <- tmp
  tmp <- names(ps2genes)
  tmp[tmp == "probeset_id"] <- "fsetid"
  names(ps2genes) <- tmp
  rm(tmp)
  if (verbose) cat("OK\n")
  
  list(probesets=probesets, ps2genes=ps2genes, genes=onlyGenes,
       level=level_schema, chromosome=chromosome_schema, type=type_schema)
}

parsePgfClf <- function(pgfFile, clfFile, probeFile, geneArray=FALSE, verbose=TRUE){
  if (verbose) msgParsingFile(pgfFile)
  pgf <- readPgf(pgfFile)
  if (verbose) msgOK()
  if (verbose) msgParsingFile(clfFile)
  clf <- readClf(clfFile)
  if (verbose) msgOK()
  geometry <- paste(clf[["dims"]], collapse=";")
  triplet <- probesetIdxToTripletIdx(pgf, 1:length(pgf[["probesetId"]]))
  fid <- pgf[["probeId"]][triplet[["probeIdx"]]]
  i <- match(fid, pgf[["probeId"]])
  ii <- match(fid, clf[["id"]])
  if (verbose) cat("Creating initial table for probes...")
  probes.table <- data.frame(fid=fid,
                             fsetid=pgf[["probesetId"]][triplet[["probesetIdx"]]],
                             pstype=pgf[["probesetType"]][triplet[["probesetIdx"]]],
                             atom=pgf[["atomId"]][triplet[["atomIdx"]]],
                             x=clf[["x"]][ii],
                             y=clf[["y"]][ii],
                             ptype=pgf[["probeType"]][i],
                             sequence=pgf[["probeSequence"]][i],
                             stringsAsFactors=FALSE)
  rm(i, ii, triplet, fid, pgf, clf)
  if (verbose) msgOK()
  probesetInfo <- parseProbesetCSV(probeFile, verbose=verbose)

  ## probesets not in the CSV file
  ## remove them
  toKeep <- which(probes.table[["fsetid"]] %in% probesetInfo[["probesets"]][["fsetid"]])
  probes.table <- probes.table[toKeep,]
  rm(toKeep)

  ## probesets -> genes table
  ## probeset_ids
  ## gid
  ps2genes <- probesetInfo[["ps2genes"]]

  ## genes table
  ## accession
  ## symbol
  ## gid
  genes <- probesetInfo[["genes"]]

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
  ##  fsetid
  ##  type
  ## IMPORTANT: I'm filtering types to "main" and "bpg" probes...
  featureSet <- probesetInfo[["probesets"]]

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
  idx <- probes.table[["pstype"]] == "main" & substr(probes.table[["ptype"]], 1, 2) == "pm"
  idx <- which(idx)
  cols <- c("fid", "fsetid", "atom", "x", "y")
  pmFeatures <- probes.table[idx, cols]
  rm(cols)

  ## if it is a Gene ST array, must split the pmFeatures table
  if (geneArray){
    cols2 <- c("fid", "fsetid", "atom")
    f2fset <- pmFeatures[, cols2]
    keys <- apply(f2fset, 1, paste, collapse=":")
    dups <- !duplicated(keys)
    f2fset <- f2fset[dups,]
    pmFeatures[["fsetid"]] <- NULL
    pmFeatures[["atom"]] <- NULL
    dups <- !duplicated(pmFeatures[["fid"]])
    pmFeatures <- pmFeatures[dups,]
    rm(cols2, dups)
  }

  pmSequence <- probes.table[idx, c("fid", "sequence")]
  pmSequence <-pmSequence[order(pmSequence[["fid"]]),]
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
  mmfeatures <- subset(probes.table, pstype=="main" & substr(ptype, 1, 2) =="mm")
  cols <- c("fid", "fsetid", "atom", "x", "y")
  mmfeatures <- mmfeatures[, cols]
  rm(cols)
  if (nrow(mmfeatures) > 0)
    stop("Add tables for MM")
  rm(mmfeatures)

  ## bgfeature table - Fields
  ##  fid
  ##  x
  ##  y
  ##  fs_type: featureSet type: genomic/antigenomic
  ##  f_type: pm/mm at/st
  idx <- grep("control->bgp", probes.table[["pstype"]])
  cols <- c("fid", "fsetid", "x", "y", "pstype", "ptype")
  bgFeatures <- probes.table[idx, cols]
  names(bgFeatures) <- c("fid", "fsetid", "x", "y", "fs_type", "f_type")
  bgFeatures <- bgFeatures[, c("fid", "fsetid", "fs_type", "f_type", "x", "y")]
  bgSequence <- probes.table[idx, c("fid", "sequence")]
  bgSequence <-bgSequence[order(bgSequence[["fid"]]),]
  bgSequence <- DataFrame(fid=bgSequence[["fid"]],
                          sequence=DNAStringSet(bgSequence[["sequence"]]))
  rm(idx, cols, probes.table)
  
  ## Here we should have the following tables available:
  ##  featureSet: fsetid, type
  ##  pmfeature: fid, fsetid, atom, x, y
  ##  bgfeature: fid, fsetid, fs_type, f_type, x, y
  ##  pmSequence: fid, sequence
  ##  bgSequence: fid, sequence

  out <- list(featureSet=featureSet, pmFeatures=pmFeatures,
              bgFeatures=bgFeatures, geometry=geometry,
              pmSequence=pmSequence, bgSequence=bgSequence,
              chrom_dict=chrom_dict, level_dict=level_dict,
              type_dict=type_dict, fset2gene=ps2genes,
              gene=genes)

  if (geneArray)
    out[["f2fset"]] <- f2fset
  
  return(out)
}

#######################################################################
## SECTION D - Package Maker
##             This shouldn't be extremely hard.
##             The idea is to: i) get array info (name, pkgname, dbname)
##             ii) parse data; iii) create pkg from template;
##             iv) dump the database
#######################################################################

## setMethod("makePdInfoPackage", "ExonTranscriptionPDInfoPkgSeed",
setMethod("makePdInfoPackage", "AffySTPDInfoPkgSeed",          
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {
            geneArray <- object@geneArray
            stopifnot(geneArray %in% c(TRUE, FALSE))
            if (geneArray){
              msg <- "Building annotation package for Affymetrix Gene ST Array"
            }else{
              msg <- "Building annotation package for Affymetrix Exon ST Array"
            }

            msgBar()
            cat(msg, "\n")
            cat("PGF.....: ", basename(object@pgfFile), "\n")
            cat("CLF.....: ", basename(object@clfFile), "\n")
            cat("Probeset: ", basename(object@probeFile), "\n")
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
            parsedData <- parsePgfClf(object@pgfFile,
                                      object@clfFile,
                                      object@probeFile,
                                      verbose=!quiet,
                                      geneArray=geneArray)
            
            #######################################################################
            ## Part iii) Create package from template
            #######################################################################
            pdInfoClass <- ifelse(geneArray, "AffyGenePDInfo", "AffyExonPDInfo")
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
                          "fset2gene",
                          fset2gene[["col2type"]],
                          fset2gene[["col2key"]])
            dbCreateTable(conn,
                          "gene",
                          gene[["col2type"]],
                          gene[["col2key"]])
            ## end adding
            
            dbCreateTable(conn,
                          "featureSet",
                          exonTranscriptionFeatureSetSchema[["col2type"]],
                          exonTranscriptionFeatureSetSchema[["col2key"]])

            if (geneArray){
              dbCreateTable(conn, "pmfeature",
                            geneStPmFeatureSchema[["col2type"]],
                            geneStPmFeatureSchema[["col2key"]])
              dbCreateTable(conn, "f2fset",
                            f2fsetSchema[["col2type"]],
                            f2fsetSchema[["col2key"]])
            }else{
              dbCreateTable(conn,
                            "pmfeature",
                            exonTranscriptionPmFeatureSchema[["col2type"]],
                            exonTranscriptionPmFeatureSchema[["col2key"]])
            }
            containsMm <- "mmFeatures" %in% names(parsedData)
            if (containsMm)
              dbCreateTable(conn,
                            "mmfeature",
                            exonTranscriptionMmFeatureSchema[["col2type"]],
                            exonTranscriptionMmFeatureSchema[["col2key"]])
            dbCreateTable(conn,
                          "bgfeature",
                          exonTranscriptionBgFeatureSchema[["col2type"]],
                          exonTranscriptionBgFeatureSchema[["col2key"]])

            ## Inserting data in new tables
            dbInsertDataFrame(conn, "chrom_dict", parsedData[["chrom_dict"]],
                              chromDictTable[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "level_dict", parsedData[["level_dict"]],
                              levelDictTable[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "type_dict", parsedData[["type_dict"]],
                              typeDictTable[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "fset2gene", parsedData[["fset2gene"]],
                              fset2gene[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "gene", parsedData[["gene"]],
                              gene[["col2type"]], !quiet)
            ## end inserting

            dbInsertDataFrame(conn, "featureSet", parsedData[["featureSet"]],
                              exonTranscriptionFeatureSetSchema[["col2type"]], !quiet)
            if (geneArray){
              dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                                geneStPmFeatureSchema[["col2type"]], !quiet)
              dbInsertDataFrame(conn, "f2fset", parsedData[["f2fset"]],
                                f2fsetSchema[["col2type"]], !quiet)
            }else{
              dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                                exonTranscriptionPmFeatureSchema[["col2type"]], !quiet)
            }
            if (containsMm)
              dbInsertDataFrame(conn, "mmfeature", parsedData[["mmFeatures"]],
                                exonTranscriptionMmFeatureSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "bgfeature", parsedData[["bgFeatures"]],
                              exonTranscriptionBgFeatureSchema[["col2type"]], !quiet)
            dbGetQuery(conn, "VACUUM")

            dbCreateTableInfo(conn, !quiet)

            ## Create indices
            dbCreateIndicesBg(conn, !quiet)
            if (geneArray){
              dbCreateIndicesPmTiling(conn, !quiet)
              dbCreateIndex(conn, "idx_f2fsfid", "f2fset", "fid", FALSE, verbose=!quiet)
              dbCreateIndex(conn, "idx_f2fsfsetid", "f2fset", "fsetid", FALSE, verbose=!quiet)
            }else{
              dbCreateIndicesPm(conn, !quiet)
            }
            dbCreateIndicesFs(conn, !quiet)
            dbCreateIndex(conn, "idx_fs2gfsetid", "fset2gene", "fsetid", FALSE, verbose=!quiet)
            dbCreateIndex(conn, "idx_fs2ggid", "fset2gene", "gid", FALSE, verbose=!quiet)
            dbCreateIndex(conn, "idx_genegid", "gene", "gid", FALSE, verbose=!quiet)
            
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
