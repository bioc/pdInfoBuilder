setClass("PkgSeed",
         representation=representation(
           name="character",
           version="character",
           license="character",
           author="character",
           email="character",
           url="character",
           biocViews="character"),
         prototype=prototype(
           license="Artistic License, Version 2.0"),
         validity=function(object) .isValidPkgSeed(object))

setClass("PDInfoPkgSeed",
         contains="PkgSeed",
         representation=representation(
           chipName="character",
           manufacturer="character",
           genomebuild="character",
           pdInfoObjectName="character"))

setClass("AffySNPPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         representation=representation(
           cdfFile="character",
           csvAnnoFile="character",
           csvSeqFile="character"
           ),
         prototype=list(manufacturer="Affymetrix"))

## FIXME: we should offer to determine packageName based on header of
## CDF file and oligo::cleanPlatformName.

setGeneric("chipName", function(object) standardGeneric("chipName"))

setMethod("chipName", "AffySNPPDInfoPkgSeed",
          function(object) {
              ## compute chip name from the CDF file
              header <- readCdfHeader(object@cdfFile)
              header$chiptype
          })


setGeneric("makePdInfoPackage", signature=c("object"),
           function(object, destDir, batch_size=10000, quiet=FALSE) {
               standardGeneric("makePdInfoPackage")
           })


setMethod("makePdInfoPackage", "AffySNPPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE) {
              chip <- chipName(object)
              pkgName <- cleanPlatformName(chip)
              extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
              dbFileName <- paste(pkgName, "sqlite", sep=".")
              dbFilePath <- file.path(extdataDir, dbFileName)
              seqMatFile <- file.path(extdataDir, "seqMat.rda")
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
                           PDINFOCLASS="AffySNPPDInfo")

              templateDir <- system.file("pd.PKG.template",
                                         package="pdInfoBuilder")
              createPackage(pkgname=pkgName, destinationDir=destDir,
                            originDir=templateDir, symbolValues=syms,
                            quiet=quiet)
              dir.create(extdataDir, recursive=TRUE)
              buildPdInfoDb(object@cdfFile, object@csvAnnoFile,
                            object@csvSeqFile, dbFilePath, seqMatFile,
                            batch_size=batch_size, verbose=!quiet)
          })


.isValidPkgSeed <- function(object) {
    email <- object@email
    if (length(email) != 1 || grep("@", email) != 1)
      return("invalid email address")
    TRUE
}
