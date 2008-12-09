datacache <- new.env(hash=TRUE, parent=emptyenv())

datacache$DEBUG <- FALSE

"initDbConnection" <- function(dbfile) {
    if (!file.exists(dbfile))
        stop("DB file '", dbfile, "' not found")
    require("RSQLite")
    dbConnect(SQLite(), dbname=dbfile, cache_size=6400, synchronous=0)
}

"closeDb" <- function(dbconn) {
    dbDisconnect(dbconn)
}

"getDb"  <- function() {
    if (!is.null(get("dbCon",datacache)) && isIdCurrent(get("dbCon",datacache)))
        return(get("dbCon",datacache))
    initDbConnection(get("DB_PATH",datacache))
}



@PDINFONAME@ <- new("@PDINFOCLASS@",
        genomebuild="@GENOMEBUILD@",
        getdb=getDb,
        geometry=as.integer(strsplit("@GEOMETRY@", ";")[[1]])) ## modify in future


.onLoad <- function(libname, pkgname) {
    require("methods", quietly=TRUE)
    DB_PATH <- system.file("extdata", "@DBFILE@",
                                   package=pkgname,
                                   lib.loc=libname)
    assign("DB_PATH", DB_PATH, envir=datacache)
    ## Establish a connection to the SQLite DB
    dbCon <- initDbConnection(DB_PATH)
    assign("dbCon",dbCon,envir=datacache)
}

.onUnload <- function(libpath) {
    closeDb(get(dbCon,envir=datacache))
}
