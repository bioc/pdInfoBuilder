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
