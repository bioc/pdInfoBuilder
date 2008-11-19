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
           pdInfoObjectName="character",
           geometry="integer"))

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
         prototype=list(manufacturer="Affymetrix"))

setClass("AffyTilingPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         representation=representation(
           bpmapFile="ScalarCharacter",
           cifFile="ScalarCharacter"),
         prototype=list(manufacturer="Affymetrix"))

setClass("AffyGenePDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         representation=representation(
           pgfFile="ScalarCharacter",
           clfFile="ScalarCharacter",
           probeFile="ScalarCharacter",
           transFile="ScalarCharacter"),
         prototype=list(manufacturer="Affymetrix"))
