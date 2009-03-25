#################################################################
## Base class
#################################################################

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
         validity=
         function(object){
           email <- object@email
           if (length(email) != 1 || grep("@", email) != 1)
             return("invalid email address")
           TRUE
         })

setClass("PDInfoPkgSeed",
         contains="PkgSeed",
         representation=representation(
           chipName="character",
           manufacturer="character",
           url="character",
           genomebuild="character",
           organism="character",
           species="character"),
         prototype=prototype(
           chipName="The Chip Name",
           manufacturer="The Manufacturer's Name",
           url="http://www.manufacturer.com",
           genomebuild="The Genome Build",
           organism="Organism",
           species="Species")
         )

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

setClass("AffySNPPDInfoPkgSeed",
         contains="AffymetrixPDInfoPkgSeed",
         representation=representation(
           cdfFile="ScalarCharacter",
           csvAnnoFile="ScalarCharacter",
           csvSeqFile="ScalarCharacter",
           splineParamFile="ScalarCharacter",
           crlmmInfoFile="ScalarCharacter",
           referenceDistFile="ScalarCharacter"))

setClass("AffySNPCNVPDInfoPkgSeed",
         contains="AffySNPPDInfoPkgSeed",
         representation=representation(
           csvAnnoFileCnv="ScalarCharacter",
           csvSeqFileCnv="ScalarCharacter"))

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
           ),
         validity=function(object) file.exists(object@ndfFile) & file.exists(object@xysFile) & file.exists(object@posFile))

setClass("NgsExpressionPDInfoPkgSeed",
         contains="NimbleGenPDInfoPkgSeed",
         representation=representation(
           ndfFile="ScalarCharacter",
           xysFile="ScalarCharacter"
           ))
