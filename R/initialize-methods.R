############## Affymetrix Arrays ###############
setMethod("initialize", "AffyExpressionPDInfoPkgSeed",
          function(.Object, cdfFile="TheCdfFile", celFile="TheCelFile", tabSeqFile="TheTabSeqFile", ...) {
              .Object@cdfFile <- new("ScalarCharacter", cdfFile)
              .Object@celFile <- new("ScalarCharacter", celFile)
              .Object@tabSeqFile <- new("ScalarCharacter", tabSeqFile)
              callNextMethod(.Object, ...)
          })

setMethod("initialize", "AffySNPPDInfoPkgSeed2",
          function(.Object, cdfFile="TheCdfFile", csvAnnoFile="TheAnnotationFile", csvSeqFile="TheSequenceFile", ...){
            .Object@cdfFile <- new("ScalarCharacter", cdfFile)
            .Object@csvAnnoFile <- new("ScalarCharacter", csvAnnoFile)
            .Object@csvSeqFile <- new("ScalarCharacter", csvSeqFile)
            callNextMethod(.Object, ...)
          })

setMethod("initialize", "AffySNPPDInfoPkgSeed",
          function(.Object, splineParamFile="TheSplineParamFile", crlmmInfoFile="TheCrlmmInfoFile", referenceDistFile="TheReferenceDistFile", ...){
            .Object@splineParamFile <- new("ScalarCharacter", splineParamFile)
            .Object@crlmmInfoFile <- new("ScalarCharacter", crlmmInfoFile)
            .Object@referenceDistFile <- new("ScalarCharacter", referenceDistFile)
            callNextMethod(.Object, ...)
          })

setMethod("initialize", "AffySNPCNVPDInfoPkgSeed2",
          function(.Object, csvAnnoFileCnv="TheAnnoFileCnv", csvSeqFileCnv="TheCsvSeqFileCnv", ...){
            .Object@csvAnnoFileCnv <- new("ScalarCharacter", csvAnnoFileCnv)
            .Object@csvSeqFileCnv <- new("ScalarCharacter", csvSeqFileCnv)
            callNextMethod(.Object, ...)
          })

setMethod("initialize", "AffySNPCNVPDInfoPkgSeed",
          function(.Object, splineParamFile="TheSplineParamFile", crlmmInfoFile="TheCrlmmInfoFile", referenceDistFile="TheReferenceDistFile", ...){
            .Object@splineParamFile <- new("ScalarCharacter", splineParamFile)
            .Object@crlmmInfoFile <- new("ScalarCharacter", crlmmInfoFile)
            .Object@referenceDistFile <- new("ScalarCharacter", referenceDistFile)
            callNextMethod(.Object, ...)
          })

#### setMethod("initialize", "AffySNPPDInfoPkgSeed",
####           function(.Object, cdfFile, csvAnnoFile, csvSeqFile,
####                    splineParamFile="", crlmmInfoFile="",
####                    referenceDistFile="", ...) {
####               .Object@cdfFile <- new("ScalarCharacter", cdfFile)
####               .Object@csvAnnoFile <- new("ScalarCharacter", csvAnnoFile)
####               .Object@csvSeqFile <- new("ScalarCharacter", csvSeqFile)
####               .Object@splineParamFile <- new("ScalarCharacter", splineParamFile)
####               .Object@crlmmInfoFile <- new("ScalarCharacter", crlmmInfoFile)
####               .Object@referenceDistFile <- new("ScalarCharacter", referenceDistFile)
####               .Object <- callNextMethod(.Object, ...)
####               .Object
####           })
#### 
#### setMethod("initialize", "AffySNPCNVPDInfoPkgSeed",
####           function(.Object, cdfFile, csvAnnoFile, csvSeqFile, csvAnnoFileCnv, csvSeqFileCnv,
####                    splineParamFile="", crlmmInfoFile="",
####                    referenceDistFile="", ...) {
####               .Object@cdfFile <- new("ScalarCharacter", cdfFile)
####               .Object@csvAnnoFile <- new("ScalarCharacter", csvAnnoFile)
####               .Object@csvSeqFile <- new("ScalarCharacter", csvSeqFile)
####               .Object@csvAnnoFileCnv <- new("ScalarCharacter", csvAnnoFileCnv)
####               .Object@csvSeqFileCnv <- new("ScalarCharacter", csvSeqFileCnv)
####               .Object@splineParamFile <- new("ScalarCharacter", splineParamFile)
####               .Object@crlmmInfoFile <- new("ScalarCharacter", crlmmInfoFile)
####               .Object@referenceDistFile <- new("ScalarCharacter", referenceDistFile)
####               .Object <- callNextMethod(.Object, ...)
####               .Object
####           })
  
setMethod("initialize", "AffyTilingPDInfoPkgSeed",
          function(.Object, bpmapFile="TheBpmapFile", celFile="TheCelFile", ...) {
              .Object@bpmapFile <- new("ScalarCharacter", bpmapFile)
              .Object@celFile <- new("ScalarCharacter", celFile)
              callNextMethod(.Object, ...)
          })
  
setMethod("initialize", "AffySTPDInfoPkgSeed",
          function(.Object, pgfFile="ThePgfFile", clfFile="TheClfFile",
          probeFile="TheProbeFile", transFile="TheTranscriptFile", coreMps="coreMps", fullMps="fullMps", extendedMps="extendedMps", ...) {
            .Object@pgfFile <- new("ScalarCharacter", pgfFile)
            .Object@clfFile <- new("ScalarCharacter", clfFile)
            .Object@probeFile <- new("ScalarCharacter", probeFile)
            .Object@transFile <- new("ScalarCharacter", transFile)
            .Object@coreMps <- new("ScalarCharacter", coreMps)
            .Object@fullMps <- new("ScalarCharacter", fullMps)
            .Object@extendedMps <- new("ScalarCharacter", extendedMps)
            callNextMethod(.Object, ...)
          })
  
########### Nimblegen Arrays ###############
setMethod("initialize", "NgsExpressionPDInfoPkgSeed",
          function(.Object, ndfFile="TheNdfFile", xysFile="TheXysFile", ...) {
            .Object@ndfFile <- new("ScalarCharacter", ndfFile)
            .Object@xysFile <- new("ScalarCharacter", xysFile)
            callNextMethod(.Object, ...)
        })

setMethod("initialize", "NgsTilingPDInfoPkgSeed",
          function(.Object, ndfFile="TheNdfFile", xysFile="TheXysFile", posFile="ThePosFile", ...){
            .Object@ndfFile <- new("ScalarCharacter", ndfFile)
            .Object@xysFile <- new("ScalarCharacter", xysFile)
            .Object@posFile <- new("ScalarCharacter", posFile)
            callNextMethod(.Object, ...)
          })
