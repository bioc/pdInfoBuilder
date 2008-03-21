setMethod("chipName", "AffySNPPDInfoPkgSeed",
          function(object) {
              ## compute chip name from the CDF file
              header <- readCdfHeader(object@cdfFile)
              header$chiptype
          })

setMethod("chipName", "AffySNPCNVPDInfoPkgSeed",
          function(object) {
              ## compute chip name from the CDF file
              header <- readCdfHeader(object@cdfFile)
              header$chiptype
          })

setMethod("chipName", "NgsPDInfoPkgSeed",
          function(object) {
            strsplit(tolower(object@ndfFile), ".ndf")[[1]]
          })

setMethod("chipName", "AffyExpressionPDInfoPkgSeed",
          function(object) {
              ## compute chip name from the CDF file
              header <- readCdfHeader(object@cdfFile)
              header$chiptype
          })

setMethod("chipName", "AffyTilingPDInfoPkgSeed",
          function(object) {
              ## compute chip name from the CDF file
            strsplit(tolower(object@bpmapFile), ".bpmap")[[1]]
          })

setMethod("chipName", "AffyGeneSTPDInfoPkgSeed",
          function(object) {
              ## compute chip name from the PGF file
            readPgf(object@pgfFile, readBody = FALSE)$header$chip_type
          })
