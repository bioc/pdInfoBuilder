#######################################################################
## SECTION A - Db Schema
#######################################################################
affySnpFeatureSetSchema <- list(col2type=c(
                                  man_fsetid="TEXT",
                                  fsetid="INTEGER",
                                  rsid="TEXT",
                                  chrom="TEXT",
                                  position="INTEGER",
                                  final="INTEGER"),
                                col2key=c(
                                  fsetid="PRIMARY KEY"
                                  ))

affySnpPmFeatureSchema <- list(col2type=c(
                                 fsetid="INTEGER",
                                 fid="INTEGER",
                                 x="INTEGER",
                                 y="INTEGER",
                                 strand="INTEGER",
                                 allele="INTEGER",
                                 atom="INTEGER"),
                               col2key=c(
                                 fid="PRIMARY KEY"
                                 ))

## affyCnvPmFeatureSchema <- list(col2type=c(
##                                  man_fsetid="TEXT",
##                                  fid="INTEGER",
##                                  strand="INTEGER",
##                                  x="INTEGER",
##                                  y="INTEGER",
##                                  fsetid="INTEGER",
##                                  chrom="TEXT",
##                                  position="INTEGER",
##                                  final="INTEGER"),
##                                col2key=c(
##                                  fid="PRIMARY KEY"
##                                  ))


## make chrom to be integer
affyCnvFeatureSetSchema <- list(col2type=c(
                                  man_fsetid="TEXT",
                                  fsetid="INTEGER",
                                  chrom="TEXT",
                                  position="INTEGER",
                                  strand="INTEGER",
                                  final="INTEGER"),
                                col2key=c(
                                  fsetid="PRIMARY KEY"
                                  ))

affyCnvPmFeatureSchema <- list(col2type=c(
                                 fsetid="INTEGER",
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

parseAnnotFile <- function(annotFile, snp=TRUE){
  annot <- read.csv(annotFile, comment.char="#", stringsAsFactors=FALSE, na.strings="---")
  names(annot) <- gsub("\\.", "", tolower(names(annot)))
  if (snp){
    cols <- c("probesetid", "dbsnprsid", "chromosome", "physicalposition", "infinallist")
  }else{
    cols <- c("probesetid", "chromosome", "chromosomestart", "infinallist")
  }
  ok <- all(cols %in% names(annot))
  if (!ok)
    stop("Annotation file must have '", paste(cols, collapse="', '"), "' columns")
  rm(ok)
  annot <- annot[, cols]
  if (snp){
    names(annot) <- c("man_fsetid", "rsid", "chrom", "position", "final")
  }else{
    names(annot) <- c("man_fsetid", "chrom", "position", "final")
  }
  annot[["final"]] <- as.integer(tolower(annot[["final"]]) == "yes")
  return(annot)
}

parseProbeSequenceFile <- function(probeseqFile, snp=TRUE, axiom=FALSE){
  probeseqTab <- read.delim(probeseqFile, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  if (!axiom){
    cols <- c("PROBESET_ID", "PROBE_X_POS", "PROBE_Y_POS", "PROBE_SEQUENCE")
    colsOut <- c("man_fsetid", "x", "y", "sequence")
  }else{
    cols <- c("PROBESET_ID", "PROBE_SEQUENCE", "ALLELE_A", "ALLELE_B")
    colsOut <- c("man_fsetid", "sequence", "allele_a", "allele_b")
  }
  checkFields(cols, names(probeseqTab))
  probeseqTab <- probeseqTab[, cols]
  names(probeseqTab) <- colsOut

  ## Convert to proper IUPAC code
  seq <- probeseqTab[["sequence"]]
  seq <- gsub("\\[G\\/T\\]", "K", seq, ignore.case=TRUE)
  seq <- gsub("\\[T\\/G\\]", "K", seq, ignore.case=TRUE)
  seq <- gsub("\\[A\\/G\\]", "R", seq, ignore.case=TRUE)
  seq <- gsub("\\[G\\/A\\]", "R", seq, ignore.case=TRUE)
  seq <- gsub("\\[C\\/T\\]", "Y", seq, ignore.case=TRUE)
  seq <- gsub("\\[C\\/T\\]", "Y", seq, ignore.case=TRUE)
  seq <- gsub("\\[G\\/C\\]", "S", seq, ignore.case=TRUE)
  seq <- gsub("\\[C\\/G\\]", "S", seq, ignore.case=TRUE)
  seq <- gsub("\\[A\\/T\\]", "W", seq, ignore.case=TRUE)
  seq <- gsub("\\[T\\/A\\]", "W", seq, ignore.case=TRUE)
  seq <- gsub("\\[A\\/C\\]", "M", seq, ignore.case=TRUE)
  seq <- gsub("\\[C\\/A\\]", "M", seq, ignore.case=TRUE)
  probeseqTab[["sequence"]] <- seq
  
  rm(cols, colsOut, seq)
  return(probeseqTab)
}

#######################################################################
## SECTION C - Parser specific for NGS Tiled Regions
##             This will take NDF/POS/XYS trio and process (in memory)
##             the data, generating data.frames for each table
##             (featureSet, pmfeature, mmfeature*, bgfeature**)
##             to be created in the db.
#######################################################################
simpleMessage <- function(...)
  message(..., appendLF=FALSE)

parseCdfSeqAnnot <- function(cdfFile, probeseqFileSNP, probeseqFileCNV, annotFileSNP, annotFileCNV, verbose=TRUE){
  if (verbose) msgParsingFile(cdfFile)
  cdf <- readCdf(cdfFile, readIndices=TRUE)
  geometry <- paste(unlist(readCdfHeader(cdfFile)[c("nrows", "ncols")]),
                    collapse=";", sep="")
  if (verbose) msgOK()
  
  theTypes <- sapply(cdf, "[[", "unittype")
  
  ## SNP
  if (verbose) simpleMessage("Getting SNP probes... ")
  idx <- which(theTypes == "genotyping")
  if (verbose) msgOK()
  
  ## assumes no MM

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
  names(pmfeatureSNP) <- c("man_fsetid", "fsetid", "fid", "x", "y", "strand", "allele", "atom")
  rm(cols)
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
  
  ## CNV
  if (verbose) simpleMessage("Organizing PM probes for CNVs... ")
  idx <- which(theTypes == "copynumber")
  rm(theTypes)
  
  ## pmfeature CNV
  ## columns: x y fid fsetid  location
  lens <- sapply(cdf[idx], function(x)
                 sum(sapply(x[["groups"]], function(y)
                            length(y[["indices"]]))))
  pmfeatureCNV <- do.call("rbind", lapply(cdf[idx], readCdfUnitToMat.cnv))
  rm(idx)
  pmfeatureCNV <- as.data.frame(pmfeatureCNV)
  pmfeatureCNV[["man_fsetid"]] <- rep(names(lens), lens)
  rm(lens)
  pmfeatureCNV[["fsetid"]] <- match(pmfeatureCNV[["man_fsetid"]], names(cdf))
  names(pmfeatureCNV) <- c("fid", "strand", "x", "y", "man_fsetid", "fsetid")
  rm(cdf)
  if (verbose) msgOK()

  ## Sequence files
  if (verbose) simpleMessage("Getting sequences for SNPs... ")
  probeseqSNP <- parseProbeSequenceFile(probeseqFileSNP)
  if (verbose) msgOK()
  if (verbose) simpleMessage("Getting sequences for CNVs... ")
  probeseqCNV <- parseProbeSequenceFile(probeseqFileCNV, snp=FALSE)
  if (verbose) msgOK()
  
  cols <- c("x", "y")
  if (verbose) simpleMessage("Merging sequence information for SNPs... ")
  pmSequenceSNP <- merge(pmfeatureSNP[, c(cols, "fid")],
                         probeseqSNP[, c(cols, "sequence")],
                         by.x=cols, by.y=cols)[, c("fid", "sequence")]
  pmSequenceSNP <- pmSequenceSNP[order(pmSequenceSNP[["fid"]]),]
  if (verbose) msgOK()
  if (verbose) simpleMessage("Merging sequence information for CNVs... ")
  pmSequenceCNV <- merge(pmfeatureCNV[, c(cols, "fid")],
                         probeseqCNV[, c(cols, "sequence")],
                         by.x=cols, by.y=cols)[, c("fid", "sequence")]
  pmSequenceCNV <- pmSequenceCNV[order(pmSequenceCNV[["fid"]]),]
  if (verbose) msgOK()
  
  rm(cols, probeseqSNP, probeseqCNV)
  if (verbose) simpleMessage("Creating Biostrings objects... ")
  pmSequenceSNP <- DataFrame(fid=pmSequenceSNP[["fid"]],
                             sequence=DNAStringSet(pmSequenceSNP[["sequence"]]))
  pmSequenceCNV <- DataFrame(fid=pmSequenceCNV[["fid"]],
                             sequence=DNAStringSet(pmSequenceCNV[["sequence"]]))
  if (verbose) msgOK()
  
  ## Annotation files
  if (verbose) msgParsingFile(annotFileSNP)
  annotSNP <- parseAnnotFile(annotFileSNP)
  if (verbose) msgOK()
  if (verbose) simpleMessage("Merging information... ")
  featureSetSNP <- merge(featureSetSNP, annotSNP, all.x=TRUE)
  rm(annotSNP)
  if (verbose) msgOK()
  
  if (verbose) msgParsingFile(annotFileCNV)
  annotCNV <- parseAnnotFile(annotFileCNV, snp=FALSE)
  if (verbose) msgOK()
  if (verbose) simpleMessage("Merging information... ")
  pmfeatureCNV <- merge(pmfeatureCNV, annotCNV, all.x=TRUE)
  rm(annotCNV)
  if (verbose) msgOK()

  if (verbose) simpleMessage("Creating featureSetCNV table... ")
  featureSetCNV <- pmfeatureCNV[, c("man_fsetid", "fsetid", "chrom", "position", "strand", "final")]
  featureSetCNV <- unique(featureSetCNV)
  rownames(featureSetCNV) <- NULL
  pmfeatureCNV <- pmfeatureCNV[, c("fsetid", "fid", "x", "y")]
  if (verbose) msgOK()
  
  out <- list(featureSet=featureSetSNP,
              pmFeatures=pmfeatureSNP,
              featureSetCNV=featureSetCNV,
              pmFeaturesCNV=pmfeatureCNV,
              pmSequenceSNP=pmSequenceSNP,
              pmSequenceCNV=pmSequenceCNV,
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
setMethod("makePdInfoPackage", "AffySNPCNVPDInfoPkgSeed2",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {

            msgBar()
            cat("Building annotation package for Affymetrix SNP/CNV Array\n")
            cat("CDF...........: ", basename(object@cdfFile), "\n")
            cat("SNP Annotation: ", basename(object@csvAnnoFile), "\n")
            cat("CNV Annotation: ", basename(object@csvAnnoFileCnv), "\n")
            cat("SNP Sequence..: ", basename(object@csvSeqFile), "\n")
            cat("CNV Sequence..: ", basename(object@csvSeqFileCnv), "\n")
            msgBar()
            
            #######################################################################
            ## Part i) get array info (chipName, pkgName, dbname)
            #######################################################################
            chip <- chipName(object)
            pkgName <- cleanPlatformName(chip)
            humanchips <- c("pd.genomewidesnp.5", "pd.genomewidesnp.6")
            if (pkgName %in% humanchips)
              warning("The package ", pkgName,
                      " *IS* available on BioConductor.",
                      " This one does *NOT* provide the data required by CRLMM.",
                      " If you have ", pkgName,
                      " downloaded/installed directly from BioConductor,",
                      " this one might overwrite the BioConductor one",
                      " and CRLMM will fail to work.")
            extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
            dbFileName <- paste(pkgName, "sqlite", sep=".")
            dbFilePath <- file.path(extdataDir, dbFileName)

            #######################################################################
            ## Part ii) parse data. This should return a list of data.frames.
            ##          The names of the elements in the list are table names.
            #######################################################################
            parsedData <- parseCdfSeqAnnot(object@cdfFile,
                                           object@csvSeqFile,
                                           object@csvSeqFileCnv,
                                           object@csvAnnoFile,
                                           object@csvAnnoFileCnv,
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
                         PDINFOCLASS="AffySNPCNVPDInfo",
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
            dbCreateTable(conn,
                          "featureSetCNV",
                          affyCnvFeatureSetSchema[["col2type"]],
                          affyCnvFeatureSetSchema[["col2key"]])
            dbCreateTable(conn,
                          "pmfeature",
                          affySnpPmFeatureSchema[["col2type"]],
                          affySnpPmFeatureSchema[["col2key"]])
            dbCreateTable(conn,
                          "pmfeatureCNV",
                          affyCnvPmFeatureSchema[["col2type"]],
                          affyCnvPmFeatureSchema[["col2key"]])
            
            dbInsertDataFrame(conn, "featureSet", parsedData[["featureSet"]],
                              affySnpFeatureSetSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "featureSetCNV", parsedData[["featureSetCNV"]],
                              affyCnvFeatureSetSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                              affySnpPmFeatureSchema[["col2type"]], !quiet)
            dbInsertDataFrame(conn, "pmfeatureCNV", parsedData[["pmFeaturesCNV"]],
                              affyCnvPmFeatureSchema[["col2type"]], !quiet)
            
            dbCreateTableInfo(conn, !quiet)

            ## Create indices
            dbCreateIndicesSnpPm(conn, !quiet)
            dbCreateIndicesCnvPm(conn, !quiet)
            dbCreateIndicesFs(conn, !quiet)
            dbCreateIndicesFsCnv(conn, !quiet)
            
            dbGetQuery(conn, "VACUUM")
            dbDisconnect(conn)
            
            #######################################################################
            ## Part v) Save sequence DataFrames
            ## FIX ME: Fix ordering of the tables to match xxFeature tables
            #######################################################################
            datadir <- file.path(destDir, pkgName, "data")
            dir.create(datadir)
            pmSequence <- parsedData[["pmSequenceSNP"]]
            pmSequenceCNV <- parsedData[["pmSequenceCNV"]]
            pmSeqFile <- file.path(datadir, "pmSequence.rda")
            pmSeqFileCNV <- file.path(datadir, "pmSequenceCNV.rda")
            if (!quiet) cat("Saving DataFrame object for SNPs / PM.\n")
            save(pmSequence, file=pmSeqFile)
            if (!quiet) cat("Saving DataFrame object for CNV / PM.\n")
            save(pmSequenceCNV, file=pmSeqFileCNV)
            if (!quiet) cat("Done.\n")
          })
