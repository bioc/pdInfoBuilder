#######################################################################
## SECTION A - Db Schema
#######################################################################

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

parseCdfSeqAnnotSnp <- function(cdfFile, probeseqFileSNP, annotFileSNP,
                                axiom=FALSE, verbose=TRUE){
  if (verbose) msgParsingFile(cdfFile)
  cdf <- readCdf(cdfFile, readIndices=TRUE, stratifyBy="pm", readBases=FALSE)
  geometry <- paste(unlist(readCdfHeader(cdfFile)[c("nrows", "ncols")]),
                    collapse=";", sep="")
  if (verbose) msgOK()
  
  theTypes <- sapply(cdf, "[[", "unittype")
  
  ## SNP
  if (verbose) simpleMessage("Getting SNP probes... ")

  if (!axiom){
      idx <- which(theTypes == "genotyping")
  }else{
      idx <- setdiff(1:length(cdf), grep("^AFFX", names(cdf)))
  }
  if (verbose) msgOK()
  
  if (verbose) simpleMessage("Organizing PM probes for SNPs... ")
  lens <- sapply(cdf[idx], function(x)
                 sum(sapply(x[["groups"]], function(y)
                            length(y[["indices"]]))))
  pmfeatureSNP <- do.call("rbind", lapply(cdf[idx], readCdfUnitToMat, verify.pmmm=FALSE))
  rm(idx)
  pmfeatureSNP <- as.data.frame(pmfeatureSNP)
  pmfeatureSNP[["man_fsetid"]] <- rep(names(lens), lens)
  rm(lens)
  pmfeatureSNP[["fsetid"]] <- match(pmfeatureSNP[["man_fsetid"]], names(cdf))
  
  ## pmfeature SNP
  ## columns: x y fid fsetid allele strand
  cols <- c("man_fsetid", "fsetid", "indices", "x", "y", "strand", "allele", "atom")
  pmfeatureSNP <- pmfeatureSNP[, cols]
  rm(cols)
  names(pmfeatureSNP) <- c("man_fsetid", "fsetid", "fid", "x", "y", "strand", "allele", "atom")

  if (axiom){
      g1 <- c('strand', 'allele', 'atom')
      g2 <- setdiff(names(pmfeatureSNP), g1)
      pmfeatureSNP <- aggregate(pmfeatureSNP[, g1],
                                by=pmfeatureSNP[, g2], mean)
      pmfeatureSNP <- pmfeatureSNP[order(pmfeatureSNP[['fsetid']],
                                         pmfeatureSNP[['fid']],
                                         pmfeatureSNP[['allele']]),]
      rownames(pmfeatureSNP) <- NULL
      easy <- which(pmfeatureSNP[['allele']] == .5)
      pmfeatureSNP[easy, 'allele'] <- 2
      rm(easy)
      pmfeatureSNP[['allele']] <- as.integer(pmfeatureSNP[['allele']])
  }
  
  if (verbose) msgOK()
  
  ## featureSet SNP
  ## columns: fsetid man_fsetid chr location rsid
  if (verbose) simpleMessage("Getting SNP information... ")
  featureSetSNP <- data.frame(man_fsetid=pmfeatureSNP[["man_fsetid"]],
                              fsetid=pmfeatureSNP[["fsetid"]],
                              stringsAsFactors=FALSE)
  dups <- duplicated(pmfeatureSNP[["fsetid"]])
  featureSetSNP <- featureSetSNP[!dups,]
  rownames(featureSetSNP) <- NULL
  rm(dups)
  pmfeatureSNP[["man_fsetid"]] <- NULL
  if (verbose) msgOK()

  ## Sequence files
  if (verbose) simpleMessage("Getting sequences for SNPs... ")
  probeseqSNP <- parseProbeSequenceFile(probeseqFileSNP, axiom=axiom)
  if (verbose) msgOK()
  
  if (verbose) simpleMessage("Merging sequence information for SNPs... ")
  if (!axiom){
      cols <- c("x", "y")
      pmSequenceSNP <- merge(pmfeatureSNP[, c(cols, "fid")],
                             probeseqSNP[, c(cols, "sequence")],
                             by.x=cols, by.y=cols)[, c("fid",
                                        "sequence")]
      rm(cols)
  }else{
      pmSequenceSNP <- merge(pmfeatureSNP[, c('fid', 'fsetid', 'allele')],
                             merge(probeseqSNP, featureSetSNP))[,c('fid', 'sequence', 'count')]
  }
  pmSequenceSNP <- pmSequenceSNP[order(pmSequenceSNP[["fid"]]),]
  rownames(pmSequenceSNP) <- NULL
  if (verbose) msgOK()
  
  rm(probeseqSNP)
  if (verbose) simpleMessage("Creating Biostrings objects... ")
  if (!axiom){
      pmSequenceSNP <- DataFrame(fid=pmSequenceSNP[["fid"]],
                                 sequence=DNAStringSet(pmSequenceSNP[["sequence"]]))
  }else{
      pmSequenceSNP <- DataFrame(fid=pmSequenceSNP[["fid"]],
                                 sequence=DNAStringSet(pmSequenceSNP[["sequence"]]),
                                 count=pmSequenceSNP[['count']])
  }
  if (verbose) msgOK()
  
  ## Annotation files
  if (verbose) msgParsingFile(annotFileSNP)
  annotSNP <- parseAnnotFile(annotFileSNP, snp=TRUE, axiom=axiom)
  if (verbose) msgOK()
  if (verbose) simpleMessage("Merging information... ")
  featureSetSNP <- merge(featureSetSNP, annotSNP, all.x=TRUE)
  rm(annotSNP)
  if (verbose) msgOK()
  
  out <- list(featureSet=featureSetSNP,
              pmFeatures=pmfeatureSNP,
              pmSequenceSNP=pmSequenceSNP,
              geometry=geometry)
  return(out)
}

#######################################################################
## SECTION D - Package Maker
##             This shouldn't be extremely hard.
##             The idea is to: i) get array info (name, pkgname, dbname)
##             ii) parse data; iii) create pkg from template;
##             iv) dump the database
#######################################################################
setMethod("makePdInfoPackage", "AffySNPPDInfoPkgSeed2",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {

            msgBar()
            message("Building annotation package for Affymetrix SNP Array\n")
            message("CDF...........: ", basename(object@cdfFile))
            message("SNP Annotation: ", basename(object@csvAnnoFile))
            message("SNP Sequence..: ", basename(object@csvSeqFile))
            message("Axiom Chip....: ", object@axiom)
            msgBar()
            
            #######################################################################
            ## Part i) get array info (chipName, pkgName, dbname)
            #######################################################################
            chip <- chipName(object)
            pkgName <- cleanPlatformName(chip)
            humanchips <- c("pd.mapping50k.xba240", "pd.mapping50k.hind240", "pd.mapping250k.nsp", "pd.mapping250k.sty")
            if (pkgName %in% humanchips)
              warning("The package ", pkgName,
                      " *IS* available on BioConductor.",
                      " Use that instead.")
            
            extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
            dbFileName <- paste(pkgName, "sqlite", sep=".")
            dbFilePath <- file.path(extdataDir, dbFileName)

            #######################################################################
            ## Part ii) parse data. This should return a list of data.frames.
            ##          The names of the elements in the list are table names.
            #######################################################################
            parsedData <- parseCdfSeqAnnotSnp(object@cdfFile,
                                              object@csvSeqFile,
                                              object@csvAnnoFile,
                                              axiom=object@axiom,
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
                         PDINFOCLASS="AffySNPPDInfo",
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
                          affySnpFeatureSetSchema[["col2type"]],
                          affySnpFeatureSetSchema[["col2key"]])
##            if (!object@axiom){
                dbCreateTable(conn,
                              "pmfeature",
                              affySnpPmFeatureSchema[["col2type"]],
                              affySnpPmFeatureSchema[["col2key"]])
##            }else{
##                dbCreateTable(conn,
##                              "pmfeature",
##                              affyAxiomSnpPmFeatureSchema[["col2type"]],
##                              affyAxiomSnpPmFeatureSchema[["col2key"]])
##            }                

            dbInsertDataFrame(conn, "featureSet", parsedData[["featureSet"]],
                              affySnpFeatureSetSchema[["col2type"]],
                              !quiet)
##            if (!object@axiom){
                dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                                  affySnpPmFeatureSchema[["col2type"]],
                                  !quiet)
##            }else{
##                dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
##                                  affyAxiomSnpPmFeatureSchema[["col2type"]],
##                                  !quiet)
##            }

            dbCreateTableInfo(conn, !quiet)

            ## Create indices
            dbCreateIndicesSnpPm(conn, !quiet)
            dbCreateIndicesFs(conn, !quiet)
            dbGetQuery(conn, "VACUUM")
            dbDisconnect(conn)
            
            #######################################################################
            ## Part v) Save sequence DataFrames
            ## FIX ME: Fix ordering of the tables to match xxFeature tables
            #######################################################################
            datadir <- file.path(destDir, pkgName, "data")
            dir.create(datadir)
            pmSequence <- parsedData[["pmSequenceSNP"]]
            pmSeqFile <- file.path(datadir, "pmSequence.rda")
            if (!quiet) cat("Saving DataFrame object for SNPs / PM.\n")
            save(pmSequence, file=pmSeqFile)
            if (!quiet) cat("Done.\n")
          })
