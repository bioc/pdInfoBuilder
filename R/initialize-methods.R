############## Affymetrix Arrays ###############
setMethod("initialize", "AffyExpressionPDInfoPkgSeed",
          function(.Object, cdfFile, celFile, tabSeqFile, ...) {
              if (!missing(cdfFile))
                .Object@cdfFile <- new("fileName", cdfFile)
              if (!missing(celFile))
                .Object@celFile <- new("fileName", celFile)
              if (!missing(tabSeqFile))
                .Object@tabSeqFile <- new("fileName", tabSeqFile)
              callNextMethod(.Object, ...)
          })

setMethod("initialize", "AffySNPPDInfoPkgSeed2",
          function(.Object, cdfFile, csvAnnoFile, csvSeqFile, ...){
            if (!missing(cdfFile))
              .Object@cdfFile <- new("fileName", cdfFile)
            if (!missing(csvAnnoFile))
              .Object@csvAnnoFile <- new("fileName", csvAnnoFile)
            if (!missing(csvSeqFile))
              .Object@csvSeqFile <- new("fileName", csvSeqFile)
            callNextMethod(.Object, ...)
          })

setMethod("initialize", "AffySNPPDInfoPkgSeed",
          function(.Object, splineParamFile, crlmmInfoFile, referenceDistFile, ...){
            if (!missing(splineParamFile))
              .Object@splineParamFile <- new("fileName", splineParamFile)
            if (!missing(crlmmInfoFile))
              .Object@crlmmInfoFile <- new("fileName", crlmmInfoFile)
            if (!missing(referenceDistFile))
              .Object@referenceDistFile <- new("fileName", referenceDistFile)
            callNextMethod(.Object, ...)
          })

setMethod("initialize", "AffySNPCNVPDInfoPkgSeed2",
          function(.Object, csvAnnoFileCnv, csvSeqFileCnv, ...){
            if (!missing(csvAnnoFileCnv))
              .Object@csvAnnoFileCnv <- new("fileName", csvAnnoFileCnv)
            if (!missing(csvSeqFileCnv))
              .Object@csvSeqFileCnv <- new("fileName", csvSeqFileCnv)
            callNextMethod(.Object, ...)
          })

setMethod("initialize", "AffySNPCNVPDInfoPkgSeed",
          function(.Object, splineParamFile, crlmmInfoFile, referenceDistFile, ...){
            if (!missing(splineParamFile))
              .Object@splineParamFile <- new("fileName", splineParamFile)
            if (!missing(crlmmInfoFile))
              .Object@crlmmInfoFile <- new("fileName", crlmmInfoFile)
            if (!missing(referenceDistFile))
              .Object@referenceDistFile <- new("fileName", referenceDistFile)
            callNextMethod(.Object, ...)
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
              if (!missing(bpmapFile))
                .Object@bpmapFile <- new("fileName", bpmapFile)
              if (!missing(celFile))
                .Object@celFile <- new("fileName", celFile)
              callNextMethod(.Object, ...)
          })
  
setMethod("initialize", "AffySTPDInfoPkgSeed",
          function(.Object, pgfFile, clfFile, probeFile, transFile, ...) {
              if (!missing(pgfFile))
                .Object@pgfFile <- new("fileName", pgfFile)
              if (!missing(clfFile))
                .Object@clfFile <- new("fileName", clfFile)
              ## probeFile and transFile are optional
              if (!missing(probeFile))
                .Object@probeFile <- new("fileName", probeFile)
              if (!missing(transFile))
                .Object@transFile <- new("fileName", transFile)
              callNextMethod(.Object, ...)
          })
  

########### Nimblegen Arrays ###############
setMethod("initialize", "NgsExpressionPDInfoPkgSeed",
          function(.Object, ndfFile, xysFile, pairFile, ngdFile, ...) {
            if (!missing(ndfFile))
              .Object@ndfFile <- new("fileName", ndfFile)
            if (!missing(xysFile))
              .Object@xysFile <- new("fileName", xysFile)
            if (!missing(pairFile))
              .Object@pairFile <- new("fileName", pairFile)
            if (!missing(ngdFile))
              .Object@ngdFile <- new("fileName", ngdFile)
            callNextMethod(.Object, ...)
        })

setMethod("initialize", "NgsTilingPDInfoPkgSeed",
          function(.Object, ndfFile, xysFile, posFile, ...){
            if (!missing(ndfFile))
              .Object@ndfFile <- new("fileName", ndfFile)
            if (!missing(xysFile))
              .Object@xysFile <- new("fileName", xysFile)
            if (!missing(posFile))
              .Object@posFile <- new("fileName", posFile)
            callNextMethod(.Object, ...)
          })
