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