setGeneric("chipName", function(object) standardGeneric("chipName"))

setGeneric("getGeometry", function(object) standardGeneric("getGeometry"))

setGeneric("makePdInfoPackage",
           function(object, destDir, batch_size=10000, quiet=FALSE, unlink=FALSE){
             standardGeneric("makePdInfoPackage")
           })

## setGeneric("makePdInfoPackage", signature=c("object"),
##            function(object, destDir, batch_size=10000, quiet=FALSE) {
##                standardGeneric("makePdInfoPackage")
##            })
