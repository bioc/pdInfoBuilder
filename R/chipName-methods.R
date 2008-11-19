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

## ## modified by Matt Settles June 2,2008
## #setMethod("chipName", "NgsPDInfoPkgSeed",
## #          function(object) {
## #            strsplit(tolower(object@ndfFile), ".ndf")[[1]]
## #          })
## 
## ## modified by Matt Settles June 2,2008
## setMethod("chipName", "NgsExpressionPDInfoPkgSeed",
##           function(object) {
##               ## compute chip name from the Pair file
##               header <- readPairHeader(object@pairFile)
##               header$DesignName
##           })
## 
## ## modified by Matt Settles June 2,2008

setMethod("chipName", "NgsPDInfoPkgSeed",
          function(object) {
            if(nchar(slot(object, "pairFile"))> 0){
              ## compute chip name from the Pair file
              header <- readPairHeader(object@pairFile)
              return(header$DesignName)
            }else{
              strsplit(tolower(object@ndfFile), ".ndf")[[1]]
            }
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

setMethod("chipName", "AffyGenePDInfoPkgSeed",
          function(object) {
              ## compute chip name from the PGF file
            readPgfEnv(object@pgfFile, readBody = FALSE)$header$lib_set_name
          })
