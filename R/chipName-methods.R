########## Affymetrix Arrays ###########
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

setMethod("chipName", "AffySNPPDInfoPkgSeed2",
          function(object) {
              ## compute chip name from the CDF file
              header <- readCdfHeader(object@cdfFile)
              header$chiptype
          })

setMethod("chipName", "AffySNPCNVPDInfoPkgSeed2",
          function(object) {
              ## compute chip name from the CDF file
              header <- readCdfHeader(object@cdfFile)
              header$chiptype
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
              readCelHeader(object@celFile)$chiptype
          })
  
setMethod("chipName", "AffySTPDInfoPkgSeed",
          function(object) {
              ## compute chip name from the PGF file
              readPgfHeader(object@pgfFile)$header$chip_type[[1]]
          })

############ Nimblegen Arrays #############

setMethod("chipName", "NimbleGenPDInfoPkgSeed",
        function(object) {
          if ("pairFile" %in% slotNames(object))
            if(nchar(slot(object, "pairFile"))> 0)
                return(readPairHeader(object@pairFile)[["DesignName"]])
          return(strsplit(tolower(basename(object@ndfFile)), ".ndf")[[1]])
        })
