######### Affymetrix Arrays ##########
setMethod("getGeometry", "AffySNPPDInfoPkgSeed",
          function(object) {
              cdfh <- readCdfHeader(object@cdfFile)
              return(list(nrows=cdfh$nrows,ncols=cdfh$ncols))
          })

setMethod("getGeometry", "AffySNPCNVPDInfoPkgSeed",
          function(object) {
              cdfh <- readCdfHeader(object@cdfFile)
              return(list(nrows=cdfh$nrows,ncols=cdfh$ncols))
          })

setMethod("getGeometry", "AffyExpressionPDInfoPkgSeed",
          function(object) {
              cdfh <- readCdfHeader(object@cdfFile)
              return(list(nrows=cdfh$nrows,ncols=cdfh$ncols))
          })

setMethod("getGeometry", "AffyTilingPDInfoPkgSeed",
          function(object) {
              celh <- readCelHeader(object@celFile)
              return(list(nrows=as.integer(celh$rows), ncols=as.integer(celh$cols)))
          })
  
setMethod("getGeometry", "AffySTPDInfoPkgSeed",
          function(object) {
              clfh <- readClfHeader(object@clfFile)
              return(list(nrows=clfh$header$rows,ncols=clfh$header$cols))
          })

########### Nimblegen Arrays ############
setMethod("getGeometry", "NimbleGenPDInfoPkgSeed",
        function(object) {
            ndfdata <- read.delim(object@ndfFile, as.is=TRUE, header=TRUE)
            return(list(nrows=max(ndfdata$Y),ncols= max(ndfdata$X)))
        })

