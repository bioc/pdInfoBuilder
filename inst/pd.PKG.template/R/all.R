globals <- new.env(hash=TRUE, parent=emptyenv())

globals$DEBUG <- FALSE

globals$DB_PATH <- system.file("extdata", "@DBFILE@",
                               package="@PKGNAME@")
if (nchar(globals$DB_PATH) == 0)
  stop("Unable to locate DB file")

initDbConnection <- function() {
    globals$dbCon <- dbConnect(dbDriver("SQLite"), dbname=globals$DB_PATH)
    globals$dbCon
}


getDb  <- function() {
    if (!is.null(globals$dbCon) && isIdCurrent(globals$dbCon))
      return(globals$dbCon)
    initDbConnection()
}


closeDb <- function() {
    ## FIXME: check for valid connection?
    sapply(dbListResults(globals$dbCon), dbClearResult)
    dbDisconnect(globals$dbCon)
    remove(dbCon, envir=globals)
}



.onLoad <- function(libname, pkgname) {
    require("methods", quietly=TRUE)
    ## Establish a connection to the SQLite DB
    initDbConnection()
}

.onUnload <- function(libpath) {
    closeDb()
}

@PDINFONAME@ <- new("@PDINFOCLASS@",
                    genomebuild="@GENOMEBUILD@",
                    getdb=getDb,
                    geometry=as.integer(strsplit("@GEOMETRY@", ";")[[1]]))

