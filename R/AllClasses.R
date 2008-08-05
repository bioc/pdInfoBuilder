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
## at this point only using pair/xys file for chipName in header
## can also get chipName from ndf filename 
## modified by Matt Settles June 2,2008
setClass("NgsPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         representation=representation(
           ndfFile="ScalarCharacter"         
           ),
         prototype=list(manufacturer="NimbleGen"))

## modified by Matt Settles June 2,2008
setClass("NgsExpressionPDInfoPkgSeed",
         contains="NgsPDInfoPkgSeed",
         representation=representation(
           xysFile="ScalarCharacter", ## backwards capatablilty
           pairFile="ScalarCharacter",
           ngdFile="ScalarCharacter"
           ))
			
## modified by Matt Settles June 2,2008
setClass("NgsTilingPDInfoPkgSeed",
         contains="NgsPDInfoPkgSeed",
         representation=representation(
           pairFile="ScalarCharacter",
           posFile="ScalarCharacter"
           ))


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
