#################################################################
## Base class
#################################################################
setClass("fileName", contains="ScalarCharacter",
         validity=function(object){
           res <- file.exists(object)
           if (!res)
             message("File '", object, "' does not exist.")
           res
         })

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
           license="Artistic License, Version 2.0",
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
           cdfFile="fileName",
           csvAnnoFile="fileName",
           csvSeqFile="fileName"))

setClass("AffySNPPDInfoPkgSeed",
         contains="AffySNPPDInfoPkgSeed2",
         representation=representation(
           splineParamFile="fileName",
           crlmmInfoFile="fileName",
           referenceDistFile="fileName"))

setClass("AffySNPCNVPDInfoPkgSeed2",
         contains="AffySNPPDInfoPkgSeed2",
         representation=representation(
           csvAnnoFileCnv="fileName",
           csvSeqFileCnv="fileName"))

setClass("AffySNPCNVPDInfoPkgSeed",
         contains="AffySNPPDInfoPkgSeed2",
         representation=representation(
           splineParamFile="fileName",
           crlmmInfoFile="fileName",
           referenceDistFile="fileName"))

## setClass("AffySNPPDInfoPkgSeed",
##          contains="AffymetrixPDInfoPkgSeed",
##          representation=representation(
##            cdfFile="fileName",
##            csvAnnoFile="fileName",
##            csvSeqFile="fileName",
##            splineParamFile="fileName",
##            crlmmInfoFile="fileName",
##            referenceDistFile="fileName"))
## 
## setClass("AffySNPCNVPDInfoPkgSeed",
##          contains="AffySNPPDInfoPkgSeed",
##          representation=representation(
##            csvAnnoFileCnv="fileName",
##            csvSeqFileCnv="fileName"))

setClass("AffyTilingPDInfoPkgSeed",
         contains="AffymetrixPDInfoPkgSeed",
         representation=representation(
           bpmapFile="fileName",
           celFile="fileName"))

setClass("AffySTPDInfoPkgSeed",
         contains="AffymetrixPDInfoPkgSeed",
         representation=representation(
           pgfFile="fileName",
           clfFile="fileName",
           probeFile="fileName",
           transFile="fileName",
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
           cdfFile="fileName",
           celFile="fileName",
           tabSeqFile="fileName"))

#################################################################
## NimbleGen seeds
#################################################################
setClass("NgsTilingPDInfoPkgSeed",
         contains="NimbleGenPDInfoPkgSeed",
         representation=representation(
           ndfFile="fileName",
           xysFile="fileName",
           posFile="fileName"
           ))

setClass("NgsExpressionPDInfoPkgSeed",
         contains="NimbleGenPDInfoPkgSeed",
         representation=representation(
           ndfFile="fileName",
           xysFile="fileName"
           ))
