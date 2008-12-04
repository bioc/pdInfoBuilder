setMethod("chipName", "AffySNPPDInfoPkgSeed",
          function(object) {
              return(paste(readCdfHeader(object@cdfFile)[c("nrows", "ncols")], collapse=";"))
          })

setMethod("getGeometry", "AffySNPCNVPDInfoPkgSeed",
          function(object) {
              return(paste(readCdfHeader(object@cdfFile)[c("nrows", "ncols")], collapse=";"))
          })

setMethod("getGeometry", "NgsPDInfoPkgSeed",
          function(object) {
              ndfdata <- read.delim(object@ndfFile, as.is=TRUE, header=TRUE)
              return(paste(max(ndfdata$Y), max(ndfdata$X), sep=";"))
          })

setMethod("getGeometry", "AffyExpressionPDInfoPkgSeed",
          function(object) {
              return(paste(readCdfHeader(object@cdfFile)[c("nrows", "ncols")], collapse=";"))
          })

setMethod("getGeometry", "AffyTilingPDInfoPkgSeed",
          function(object) {
              ## get geometry, but why do we need here?
              celh <- readCelHeader(celFile)
              nx <- as.integer(celh$cols)
              ny <- as.integer(celh$rows)
              return(paste(ny, nx, sep=";"))
          })
  
setMethod("getGeometry", "AffySTPDInfoPkgSeed",
          function(object) {
              ## compute chip name from the PGF file
              readPgfHeader(object@pgfFile)$header$chip_type[[1]]
          })
