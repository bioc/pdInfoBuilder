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
    if (!is.null(globals$dbCon))
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

.onAttach <- function(libname, pkgname) {
    ## FIXME: once we've made a more careful pass over the data, this message
    ## should go away!
    cat("\n\n")
    cat("Welcome to the pd.@CHIPNAME@ prototype pdInfo package\n")
    cat("WARNING: DO NOT USE THIS PACKAGE FOR ANY ANALYSIS.\n")
    cat("THIS PACKAGE IS FOR INTERFACE PROTOTYPE USE ONLY!\n")
    cat("THE DATA HAS NOT BEEN VALIDATED AND LIKELY HAS ERRORS.\n")
    cat("Have fun!\n\n")
}

.onUnload <- function(libpath) {
    closeDb()
}

@PDINFONAME@ <- new("@PDINFOCLASS@",
                    genomebuild="@GENOMEBUILD@",
                    getdb=getDb,
                    geometry=as.integer(strsplit("@GEOMETRY@", ";")[[1]]))

