datacache <- new.env(hash=TRUE, parent=emptyenv())

datacache$DEBUG <- FALSE

## setup the path at package level so that DB can be accessed
## during package install/lazyload db creation.
##
## We reset the DB_PATH in .onLoad since we need to
## get the right one based on libpath
#datacache$DB_PATH <- system.file("extdata", "@DBFILE@",
#                               package="@PKGNAME@")
#if (nchar(datacache$DB_PATH) == 0)
#  stop("Unable to locate DB file")

"initDbConnection" <- function() {
    datacache$dbCon <- dbConnect(dbDriver("SQLite"), dbname=datacache$DB_PATH)
    datacache$dbCon
}

"getDb"  <- function() {
    if (!is.null(datacache$dbCon) && isIdCurrent(datacache$dbCon))
      return(datacache$dbCon)
    initDbConnection()
}

"closeDb" <- function() {
    ## FIXME: check for valid connection?
    sapply(dbListResults(datacache$dbCon), dbClearResult)
    dbDisconnect(datacache$dbCon)
    remove(dbCon, envir=datacache)
}

.onLoad <- function(libname, pkgname) {
    require("methods", quietly=TRUE)
    datacache$DB_PATH <- system.file("extdata", "@DBFILE@",
                                   package=pkgname,
                                   lib.loc=libname)
    if (nchar(datacache$DB_PATH) == 0)
      stop("Unable to locate DB file")
    ## Establish a connection to the SQLite DB
    initDbConnection()
}

.onUnload <- function(libpath) {
    closeDb()
}

@PDINFONAME@ <- new("@PDINFOCLASS@",
                    genomebuild="@GENOMEBUILD@",
                    getdb=getDb,
                    geometry=as.integer(strsplit("@GEOMETRY@", ";")[[1]])) ## modify in future
