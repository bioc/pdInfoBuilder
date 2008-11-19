setGeneric("chipName", function(object) standardGeneric("chipName"))

setGeneric("makePdInfoPackage",
           function(object, destDir, batch_size=10000, quiet=FALSE){
             standardGeneric("makePdInfoPackage")
           })

## setGeneric("makePdInfoPackage", signature=c("object"),
##            function(object, destDir, batch_size=10000, quiet=FALSE) {
##                standardGeneric("makePdInfoPackage")
##            })
