initDb <- function(dbname) {
    db <- dbConnect(dbDriver("SQLite"), dbname)

    ## Set page size
    dbGetQuery(db, setPageSizeSql)
    
    ## Create tables
    dbGetQuery(db, createFeatureSetSql)
    
    dbGetQuery(db, sprintf(createFeatureSql, "pmfeature"))
    
    dbGetQuery(db, sprintf(createFeatureSql, "mmfeature"))

    ## Create index
    ## NOTE: might be more efficient to create this after loading,
    ## but current perf is ok.
    sql <- 'create index man_fsetid_idx on featureSet ("man_fsetid")'
    dbGetQuery(db, sql)
    db
}

closeDb <- function(db) dbDisconnect(db)
