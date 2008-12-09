## Utility functions for producing PDInfo packages
ST <- system.time

printTime <- function(msg, t) {
        m <- paste(msg, "took %.2f sec\n")
        cat(sprintf(m, t))
}

validInput <- function(object, required=c(), optional=c()) {
    msg <- NULL
    ok <- sapply(required,
            function(slt) file.exists(slot(object, slt)))
    if (!all(ok))
        msg <- paste("missing required file(s):",
                paste(sapply(names(ok[!ok]), function(slt) slt), "='",
                        sapply(names(ok[!ok]),
                                function(slt) slot(object, slt)),"'",collapse=", ", sep=""))
    ok <- sapply(optional,
            function(slt) ifelse(is.na(slot(object,slt)),TRUE, file.exists(slot(object, slt))))
    if(!all(ok))
        msg <- paste(msg,"\n","missing optional file(s):",
                paste(sapply(names(ok[!ok]), function(slt) slt), "='",
                        sapply(names(ok[!ok]),
                                function(slt) slot(object, slt)),"'",collapse=", ", sep=""))
    
    if (is.null(msg)) TRUE else msg
} # end validInput

setupPackage <- function(object,pkgName,destDir,dbFileName,unlink,quiet){
    geometry  <- getGeometry(object)
    oligoc_dbi_version =installed.packages()['oligoClasses','Version']
    syms <- list(
            ## present in DESCRIPTION
            PKGNAME      =pkgName,
            MANUF        =object@manufacturer,
            CHIPNAME     =object@chipName,
            PKGVERSION   =object@version,
            AUTHOR       =object@author,
            AUTHOREMAIL  =object@email,
            MANUFURL     =object@url,
            LIC          =object@license,
            ORGANISM     =object@organism,
            GENOMEBUILD  =object@genomebuild,
            SPECIES      =object@species,                
            BIOCVIEWS    =object@biocViews,
            OLIGOCVERSION=oligoc_dbi_version,
            ## present in namespace
            PDINFONAME   =pkgName,
            ## present in all.R
            ## not sure where these are used
            GEOMETRY     =paste(geometry$nrows,geometry$ncols, sep=";"), ## remove in future
            DBFILE       =dbFileName,
            PDINFOCLASS  = sub("PkgSeed","", class(object)))
    
    templateDir <- system.file("pdInfo.PKG.template", package="pdInfoBuilder")
    createPackage(pkgname=pkgName, destinationDir=destDir,
            originDir=templateDir, symbolValues=syms,
            unlink = unlink,quiet=quiet)
}              

"connectDb" <- function(dbfile) {
    if (!file.exists(dbfile))
        stop("DB file '", dbfile, "' not found")
    require("RSQLite")
    dbConnect(SQLite(), dbname=dbfile, cache_size=6400, page_size = 8192, synchronous=0)
}

"closeDb" <- function(db){
    dbDisconnect(db)
}
