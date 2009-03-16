############## Affymetrix Arrays ###############
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

setMethod("initialize", "AffyExpressionPDInfoPkgSeed",
          function(.Object, cdfFile, csvAnnoFile, tabSeqFile, ...) {
              .Object@cdfFile <- new("ScalarCharacter", cdfFile)
              ## tabSeqFile and csvAnnoFile are optional
              if (!missing(csvAnnoFile))
                  .Object@csvAnnoFile <- new("ScalarCharacter", csvAnnoFile)
              if (!missing(tabSeqFile))
                  .Object@tabSeqFile <- new("ScalarCharacter", tabSeqFile)
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
##             if (!missing(pairFile))
##                 .Object@pairFile <- new("ScalarCharacter", pairFile)
##             .Object@ngdFile <- new("ScalarCharacter", ngdFile)
            callNextMethod(.Object, ...)
            .Object <- callNextMethod(.Object, ...)
            .Object
        })


setMethod("initialize", "NgsTilingPDInfoPkgSeed",
        function(.Object, ndfFile, xysFile, pairFile, posFile, ...) {
            .Object@ndfFile <- new("ScalarCharacter", ndfFile)
            if (!missing(xysFile))
                .Object@xysFile <- new("ScalarCharacter", xysFile)
            if (!missing(pairFile))
                .Object@pairFile <- new("ScalarCharacter", pairFile)
            .Object@posFile <- new("ScalarCharacter", posFile)
            callNextMethod(.Object, ...)
        })


## modified by Matt Settles June 2,2008
## must specify one of pairFile or xysFile
## setMethod("initialize", "NgsExpressionPDInfoPkgSeed",
##           function(.Object, ndfFile, pairFile="", xysFile="", ngdFile="", ...) {
##               .Object@ndfFile <- new("ScalarCharacter", ndfFile)
##               .Object@pairFile <- new("ScalarCharacter", pairFile)
##               .Object@xysFile <- new("ScalarCharacter", xysFile)
##               .Object@ngdFile <- new("ScalarCharacter", ngdFile)
##               .Object <- callNextMethod(.Object, ...)
##               .Object
##           })

## modified by Matt Settles June 2,2008 ## do you combine multiple arrays sets?
## setMethod("initialize", "NgsTilingPDInfoPkgSeed",
##           function(.Object, ndfFile, posFile, xysFile, pairFile, ...) {
##               .Object@ndfFile <- new("ScalarCharacter", ndfFile)
##               .Object@posFile <- new("ScalarCharacter", posFile)
##               if (!missing(xysFile))
##                 .Object@xysFile <- new("ScalarCharacter", xysFile) ## BC 11/18/2008
##               if (!missing(pairFile))
##                 .Object@pairFile <- new("ScalarCharacter", pairFile)
##               .Object <- callNextMethod(.Object, ...)
##               .Object
##           })

