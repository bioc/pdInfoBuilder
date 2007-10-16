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

setMethod("initialize", "NgsPDInfoPkgSeed",
          function(.Object, ndfFile, posFile, xysFile, ...) {
              .Object@ndfFile <- new("ScalarCharacter", ndfFile)
              .Object@posFile <- new("ScalarCharacter", posFile)
              .Object@xysFile <- new("ScalarCharacter", xysFile)
              .Object <- callNextMethod(.Object, ...)
              .Object
          })

setMethod("initialize", "AffyExpressionPDInfoPkgSeed",
          function(.Object, cdfFile, csvAnnoFile, tabSeqFile, ...) {
              .Object@cdfFile <- new("ScalarCharacter", cdfFile)
              .Object@csvAnnoFile <- new("ScalarCharacter", csvAnnoFile)
              .Object@tabSeqFile <- new("ScalarCharacter", tabSeqFile)
              .Object <- callNextMethod(.Object, ...)
              .Object
          })

