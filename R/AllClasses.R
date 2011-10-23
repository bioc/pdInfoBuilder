#################################################################
## Base class
#################################################################
setClass("PDInfoPkgSeed",
         representation=representation(
           version="character",
           license="character",
           author="character",
           email="character",
           biocViews="character",
           chipName="character",
           manufacturer="character",
           url="character",
           genomebuild="character",
           organism="character",
           species="character"),
         prototype=prototype(
           version="0.0.1",
           license="Artistic-2.0",
           author="My Name",
           email="my@email.com",
           biocViews="AnnotationData",
           chipName="The Chip Name",
           manufacturer="The Manufacturer's Name",
           url="http://www.manufacturer.com",
           genomebuild="The Genome Build",
           organism="Organism",
           species="Species"),
         validity=
         function(object){
           email <- object@email
           if (length(email) != 1 || grep("@", email) != 1)
             return("invalid email address")
           TRUE
         })

#################################################################
## Manufacturer-specific classes: Affymetrix and NimbleGen are
## supported for the moment
#################################################################
setClass("AffymetrixPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         prototype=list(
           manufacturer="Affymetrix",
           url="http://www.affymetrix.com"
           ))

setClass("NimbleGenPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         prototype=list(
           manufacturer="NimbleGen",
           url="http://www.nimblegen.com"
           ))

#################################################################
## Affymetrix seeds
#################################################################
setClass("AffySNPPDInfoPkgSeed2",
         contains="AffymetrixPDInfoPkgSeed",
         representation=representation(
           cdfFile="ScalarCharacter",
           csvAnnoFile="ScalarCharacter",
           csvSeqFile="ScalarCharacter",
           axiom="logical"
         ))

setClass("AffySNPPDInfoPkgSeed",
         contains="AffySNPPDInfoPkgSeed2",
         representation=representation(
           splineParamFile="ScalarCharacter",
           crlmmInfoFile="ScalarCharacter",
           referenceDistFile="ScalarCharacter"))

setClass("AffySNPCNVPDInfoPkgSeed2",
         contains="AffySNPPDInfoPkgSeed2",
         representation=representation(
           csvAnnoFileCnv="ScalarCharacter",
           csvSeqFileCnv="ScalarCharacter"))

setClass("AffySNPCNVPDInfoPkgSeed",
         contains="AffySNPCNVPDInfoPkgSeed2",
         representation=representation(
           splineParamFile="ScalarCharacter",
           crlmmInfoFile="ScalarCharacter",
           referenceDistFile="ScalarCharacter"))

setClass("AffyTilingPDInfoPkgSeed",
         contains="AffymetrixPDInfoPkgSeed",
         representation=representation(
           bpmapFile="ScalarCharacter",
           celFile="ScalarCharacter"))

setClass("AffySTPDInfoPkgSeed",
         contains="AffymetrixPDInfoPkgSeed",
         representation=representation(
           pgfFile="ScalarCharacter",
           clfFile="ScalarCharacter",
           probeFile="ScalarCharacter",
           transFile="ScalarCharacter",
           coreMps="ScalarCharacter",
           fullMps="ScalarCharacter",
           extendedMps="ScalarCharacter",
           geneArray="logical"))

setClass("AffyExonPDInfoPkgSeed",
         contains="AffySTPDInfoPkgSeed",
         prototype=list(geneArray=FALSE))

setClass("AffyGenePDInfoPkgSeed",
         contains="AffySTPDInfoPkgSeed",
         prototype=list(geneArray=TRUE))

setClass("AffyExpressionPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         representation=representation(
           cdfFile="ScalarCharacter",
           celFile="ScalarCharacter",
           tabSeqFile="ScalarCharacter"))

#################################################################
## NimbleGen seeds
#################################################################
setClass("NgsTilingPDInfoPkgSeed",
         contains="NimbleGenPDInfoPkgSeed",
         representation=representation(
           ndfFile="ScalarCharacter",
           xysFile="ScalarCharacter",
           posFile="ScalarCharacter"
           ))

setClass("NgsExpressionPDInfoPkgSeed",
         contains="NimbleGenPDInfoPkgSeed",
         representation=representation(
           ndfFile="ScalarCharacter",
           xysFile="ScalarCharacter"
           ))
