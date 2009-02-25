
setClass("PkgSeed",
         representation=representation(
                name="character",
                version="character",
                license="character",
                author="character",
                email="character",
                biocViews="character"),
        prototype=prototype(
                name="The Name",
                version="0.0.1",
                license="Artistic License, Version 2.0",
                author="My Name",
                email="my@email.com",
                biocViews="AnnotationData"),
        validity=function(object) .isValidPkgSeed(object))

setClass("PDInfoPkgSeed",
        contains="PkgSeed",
        representation=representation(
                chipName="character",
                manufacturer="character",
                url="character",
                genomebuild="character",
                organism = "character",
                species = "character"),
        prototype=prototype(
                chipName="The Chip Name",
                manufacturer="The Manufacturer's Name",
                url="http://www.manufacturer.com",
                genomebuild="The Genome Build",
                organism="Organism",
                species="Species")
)


setClass("AffySNPPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         representation=representation(
           cdfFile="ScalarCharacter",
           csvAnnoFile="ScalarCharacter",
           csvSeqFile="ScalarCharacter",
           splineParamFile="ScalarCharacter",
           crlmmInfoFile="ScalarCharacter",
           referenceDistFile="ScalarCharacter"
           ),
         prototype=list(manufacturer="Affymetrix"))

setClass("AffySNPCNVPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         representation=representation(
           cdfFile="ScalarCharacter",
           csvAnnoFile="ScalarCharacter",
           csvSeqFile="ScalarCharacter",
           csvAnnoFileCnv="ScalarCharacter",
           csvSeqFileCnv="ScalarCharacter",
           splineParamFile="ScalarCharacter",
           crlmmInfoFile="ScalarCharacter",
           referenceDistFile="ScalarCharacter"
           ),
         prototype=list(manufacturer="Affymetrix"))


.isValidPkgSeed <- function(object) {
    email <- object@email
    if (length(email) != 1 || grep("@", email) != 1)
      return("invalid email address")
    TRUE
}

setClass("NgsPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         prototype=list(manufacturer="NimbleGen"))

## modified by Matt Settles June 2,2008
setClass("NgsExpressionPDInfoPkgSeed",
         contains="NgsPDInfoPkgSeed",
         representation=representation(
           ndfFile="ScalarCharacter",
           xysFile="ScalarCharacter", ## BC, 11/18/08
           pairFile="ScalarCharacter",
           ngdFile="ScalarCharacter"
           ))
			
## modified by Matt Settles June 2,2008
setClass("NgsTilingPDInfoPkgSeed",
         contains="NgsPDInfoPkgSeed",
         representation=representation(
           ndfFile="ScalarCharacter",
           xysFile="ScalarCharacter", ## BC, 11/18/08
           pairFile="ScalarCharacter",
           posFile="ScalarCharacter"
           ))


validNgsTilingPDInfoPkgSeed <- function(object){
  ndf <-  nchar(slot(object, "ndfFile")) > 0
  xys <-  nchar(slot(object, "xysFile")) > 0
  pair <- nchar(slot(object, "pairFile"))> 0
  if (pair & xys) stop("Specify only one: XYS or PAIR file.")
  if (!ndf) stop("NDF must be given.")
  stopifnot(file.exists(slot(object, "ndfFile")))
  if (xys) stopifnot(file.exists(slot(object, "xysFile")))
  if (pair) stopifnot(file.exists(slot(object, "pairFile")))
  TRUE
}
setValidity("NgsTilingPDInfoPkgSeed", validNgsTilingPDInfoPkgSeed)



setClass("AffyExpressionPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         representation=representation(
                 cdfFile="ScalarCharacter",
                 csvAnnoFile="ScalarCharacter",
                 tabSeqFile="ScalarCharacter"),
         prototype=prototype(
                 manufacturer="Affymetrix",
                 cdfFile=mkScalar(as.character(NA)),
                 csvAnnoFile=mkScalar(as.character(NA)),
                 tabSeqFile=mkScalar(as.character(NA))))

 ## changed cif file to cel file, cif doesn't seem to really provide anything, expanted prototype
 setClass("AffyTilingPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         representation=representation(
                 bpmapFile="ScalarCharacter",
                 celFile="ScalarCharacter"),
         prototype=prototype(
                 manufacturer="Affymetrix",
                 bpmapFile=mkScalar(as.character(NA)),
                 celFile=mkScalar(as.character(NA))))

setClass("AffySTPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         representation=representation(
                 pgfFile="ScalarCharacter",
                 clfFile="ScalarCharacter",
                 probeFile="ScalarCharacter",
                 transFile="ScalarCharacter"),
         prototype=prototype(
                 manufacturer="Affymetrix",
                 pgfFile=mkScalar(as.character(NA)),
                 clfFile=mkScalar(as.character(NA)),
                 probeFile=mkScalar(as.character(NA)),
                 transFile=mkScalar(as.character(NA))))
 
 setClass("AffyExonPDInfoPkgSeed", contains="AffySTPDInfoPkgSeed")
 setClass("AffyGenePDInfoPkgSeed", contains="AffySTPDInfoPkgSeed")
 
##### PD Info v2 by BC
##### please don't modify
setClass("NgsTiledRegionPDInfoPkgSeed",
         contains="NgsPDInfoPkgSeed",
         representation=representation(
           ndfFile="ScalarCharacter",
           xysFile="ScalarCharacter",
           posFile="ScalarCharacter"
           ),
         validity=function(object) file.exists(object@ndfFile) & file.exists(object@xysFile) & file.exists(object@posFile))
