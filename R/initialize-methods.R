############## Affymetrix Arrays ###############
setMethod("initialize", "AffyExpressionPDInfoPkgSeed",
          function(.Object, cdfFile, celFile, tabSeqFile, ...) {
              .Object@cdfFile <- new("ScalarCharacter", cdfFile)
              .Object@celFile <- new("ScalarCharacter", celFile)
              .Object@tabSeqFile <- new("ScalarCharacter", tabSeqFile)
              .Object <- callNextMethod(.Object, ...)
              .Object
          })

setMethod("initialize", "AffySNPPDInfoPkgSeed",
          function(.Object, cdfFile, csvAnnoFile, csvSeqFile,
                   splineParamFile="", crlmmInfoFile="",
                   referenceDistFile="", ...) {
              .Object@cdfFile <- new("ScalarCharacter", cdfFile)
              .Object@csvAnnoFile <- new("ScalarCharacter", csvAnnoFile)
              .Object@csvSeqFile <- new("ScalarCharacter", csvSeqFile)
              .Object@splineParamFile <- new("ScalarCharacter", splineParamFile)
              .Object@crlmmInfoFile <- new("ScalarCharacter", crlmmInfoFile)
              .Object@referenceDistFile <- new("ScalarCharacter", referenceDistFile)
              .Object <- callNextMethod(.Object, ...)
              .Object
          })

setMethod("initialize", "AffySNPCNVPDInfoPkgSeed",
          function(.Object, cdfFile, csvAnnoFile, csvSeqFile, csvAnnoFileCnv, csvSeqFileCnv,
                   splineParamFile="", crlmmInfoFile="",
                   referenceDistFile="", ...) {
              .Object@cdfFile <- new("ScalarCharacter", cdfFile)
              .Object@csvAnnoFile <- new("ScalarCharacter", csvAnnoFile)
              .Object@csvSeqFile <- new("ScalarCharacter", csvSeqFile)
              .Object@csvAnnoFileCnv <- new("ScalarCharacter", csvAnnoFileCnv)
              .Object@csvSeqFileCnv <- new("ScalarCharacter", csvSeqFileCnv)
              .Object@splineParamFile <- new("ScalarCharacter", splineParamFile)
              .Object@crlmmInfoFile <- new("ScalarCharacter", crlmmInfoFile)
              .Object@referenceDistFile <- new("ScalarCharacter", referenceDistFile)
              .Object <- callNextMethod(.Object, ...)
              .Object
          })
  
setMethod("initialize", "AffyTilingPDInfoPkgSeed",
          function(.Object, bpmapFile, celFile, ...) {
              .Object@bpmapFile <- new("ScalarCharacter", bpmapFile)
              .Object@celFile <- new("ScalarCharacter", celFile)
              .Object <- callNextMethod(.Object, ...)
              .Object
          })
  
setMethod("initialize", "AffySTPDInfoPkgSeed",
          function(.Object, pgfFile, clfFile, probeFile, transFile, ...) {
              .Object@pgfFile <- new("ScalarCharacter", pgfFile)
              .Object@clfFile <- new("ScalarCharacter", clfFile)
              ## probeFile and transFile are optional
              if (!missing(probeFile))
                  .Object@probeFile <- new("ScalarCharacter", probeFile)
              if (!missing(transFile))
                  .Object@transFile <- new("ScalarCharacter", transFile)
              .Object <- callNextMethod(.Object, ...)
              .Object
          })
  

########### Nimblegen Arrays ###############
setMethod("initialize", "NgsExpressionPDInfoPkgSeed",
          function(.Object, ndfFile, xysFile, pairFile, ngdFile, ...) {
            .Object@ndfFile <- new("ScalarCharacter", ndfFile)
            if (!missing(xysFile))
              .Object@xysFile <- new("ScalarCharacter", xysFile)
            callNextMethod(.Object, ...)
            .Object <- callNextMethod(.Object, ...)
            .Object
        })

setMethod("initialize", "NgsTilingPDInfoPkgSeed",
          function(.Object, ndfFile, xysFile, posFile, ...){
            .Object@ndfFile <- new("ScalarCharacter", ndfFile)
            .Object@xysFile <- new("ScalarCharacter", xysFile)
            .Object@posFile <- new("ScalarCharacter", posFile)
            .Object <- callNextMethod(.Object, ...)
            .Object
          })
