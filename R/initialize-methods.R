############## Affymetrix Arrays ###############
setMethod("initialize", "AffyExpressionPDInfoPkgSeed",
          function(.Object, cdfFile, celFile, tabSeqFile, ...) {
              .Object@cdfFile <- new("fileName", cdfFile)
              .Object@celFile <- new("fileName", celFile)
              .Object@tabSeqFile <- new("fileName", tabSeqFile)
              .Object <- callNextMethod(.Object, ...)
              .Object
          })

setMethod("initialize", "AffySNPPDInfoPkgSeed2",
          function(.Object, cdfFile, csvAnnoFile, csvSeqFile, ...){
            message("I'm in 'initialize' for AffySNPPDInfoPkgSeed2... ")
            .Object@cdfFile <- new("fileName", cdfFile)
            .Object@csvAnnoFile <- new("fileName", csvAnnoFile)
            .Object@csvSeqFile <- new("fileName", csvSeqFile)
            callNextMethod()
          })

setMethod("initialize", "AffySNPPDInfoPkgSeed",
          function(.Object, splineParamFile, crlmmInfoFile, referenceDistFile, ...){
            .Object@splineParamFile <- new("fileName", splineParamFile)
            .Object@crlmmInfoFile <- new("fileName", crlmmInfoFile)
            .Object@referenceDistFile <- new("fileName", referenceDistFile)
            callNextMethod()
          })

setMethod("initialize", "AffySNPCNVPDInfoPkgSeed2",
          function(.Object, csvAnnoFileCnv, csvSeqFileCnv, ...){
            .Object@csvAnnoFileCnv <- new("fileName", csvAnnoFileCnv)
            .Object@csvSeqFileCnv <- new("fileName", csvSeqFileCnv)
            callNextMethod()
          })

setMethod("initialize", "AffySNPCNVPDInfoPkgSeed",
          function(.Object, splineParamFile, crlmmInfoFile, referenceDistFile, ...){
            .Object@splineParamFile <- new("fileName", splineParamFile)
            .Object@crlmmInfoFile <- new("fileName", crlmmInfoFile)
            .Object@referenceDistFile <- new("fileName", referenceDistFile)
            callNextMethod()
          })

#### setMethod("initialize", "AffySNPPDInfoPkgSeed",
####           function(.Object, cdfFile, csvAnnoFile, csvSeqFile,
####                    splineParamFile="", crlmmInfoFile="",
####                    referenceDistFile="", ...) {
####               .Object@cdfFile <- new("fileName", cdfFile)
####               .Object@csvAnnoFile <- new("fileName", csvAnnoFile)
####               .Object@csvSeqFile <- new("fileName", csvSeqFile)
####               .Object@splineParamFile <- new("fileName", splineParamFile)
####               .Object@crlmmInfoFile <- new("fileName", crlmmInfoFile)
####               .Object@referenceDistFile <- new("fileName", referenceDistFile)
####               .Object <- callNextMethod(.Object, ...)
####               .Object
####           })
#### 
#### setMethod("initialize", "AffySNPCNVPDInfoPkgSeed",
####           function(.Object, cdfFile, csvAnnoFile, csvSeqFile, csvAnnoFileCnv, csvSeqFileCnv,
####                    splineParamFile="", crlmmInfoFile="",
####                    referenceDistFile="", ...) {
####               .Object@cdfFile <- new("fileName", cdfFile)
####               .Object@csvAnnoFile <- new("fileName", csvAnnoFile)
####               .Object@csvSeqFile <- new("fileName", csvSeqFile)
####               .Object@csvAnnoFileCnv <- new("fileName", csvAnnoFileCnv)
####               .Object@csvSeqFileCnv <- new("fileName", csvSeqFileCnv)
####               .Object@splineParamFile <- new("fileName", splineParamFile)
####               .Object@crlmmInfoFile <- new("fileName", crlmmInfoFile)
####               .Object@referenceDistFile <- new("fileName", referenceDistFile)
####               .Object <- callNextMethod(.Object, ...)
####               .Object
####           })
  
setMethod("initialize", "AffyTilingPDInfoPkgSeed",
          function(.Object, bpmapFile, celFile, ...) {
              .Object@bpmapFile <- new("fileName", bpmapFile)
              .Object@celFile <- new("fileName", celFile)
              .Object <- callNextMethod(.Object, ...)
              .Object
          })
  
setMethod("initialize", "AffySTPDInfoPkgSeed",
          function(.Object, pgfFile, clfFile, probeFile, transFile, ...) {
              .Object@pgfFile <- new("fileName", pgfFile)
              .Object@clfFile <- new("fileName", clfFile)
              ## probeFile and transFile are optional
              if (!missing(probeFile))
                  .Object@probeFile <- new("fileName", probeFile)
              if (!missing(transFile))
                  .Object@transFile <- new("fileName", transFile)
              .Object <- callNextMethod(.Object, ...)
              .Object
          })
  

########### Nimblegen Arrays ###############
setMethod("initialize", "NgsExpressionPDInfoPkgSeed",
          function(.Object, ndfFile, xysFile, pairFile, ngdFile, ...) {
            .Object@ndfFile <- new("fileName", ndfFile)
            if (!missing(xysFile))
              .Object@xysFile <- new("fileName", xysFile)
            callNextMethod(.Object, ...)
            .Object <- callNextMethod(.Object, ...)
            .Object
        })

setMethod("initialize", "NgsTilingPDInfoPkgSeed",
          function(.Object, ndfFile, xysFile, posFile, ...){
            .Object@ndfFile <- new("fileName", ndfFile)
            .Object@xysFile <- new("fileName", xysFile)
            .Object@posFile <- new("fileName", posFile)
            .Object <- callNextMethod(.Object, ...)
            .Object
          })
