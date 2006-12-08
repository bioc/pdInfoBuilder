initDb <- function(dbname) {
    db <- dbConnect(dbDriver("SQLite"), dbname)

    ## Set page size
    dbGetQuery(db, setPageSizeSql)
    
    ## Create tables
    ## BC: Soon we need to add a table for the control probes
    dbGetQuery(db, createFeatureSetSql)
    
    dbGetQuery(db, sprintf(createFeatureSql, "pmfeature_tmp"))
    
    dbGetQuery(db, sprintf(createFeatureSql, "mmfeature_tmp"))

    dbGetQuery(db, createPm_MmSql)

    dbGetQuery(db, createSequenceSql)

    ## Create index
    ## NOTE: might be more efficient to create this after loading,
    ## but current perf is ok.
    sql <- 'create index man_fsetid_idx on featureSet ("man_fsetid")'
    dbGetQuery(db, sql)
    db
}

createIndicesDb <- function(db) {
    makeIndex <- function(name, t, cols) {
        sql <- paste("create index", name, "on", t,
                     paste("(", paste(cols, collapse=","), ")"))
        dbGetQuery(db, sql)
    }
    
    ## Create DB indices and fix ordering
    dbGetQuery(db, sprintf(createFeatureSql, "pmfeature"))
    dbGetQuery(db, sprintf(createFeatureSql, "mmfeature"))

    ## Reorder XXfeature tables
    fillSql <- paste("insert into %s select * from %s order by",
                     "strand, allele, fsetid")
    dbBeginTransaction(db)
    dbGetQuery(db, sprintf(fillSql, "pmfeature", "pmfeature_tmp"))
    dbGetQuery(db, sprintf(fillSql, "mmfeature", "mmfeature_tmp"))
    dbCommit(db)
    ## drop temp tables
    dbGetQuery(db, "drop table pmfeature_tmp")
    dbGetQuery(db, "drop table mmfeature_tmp")

    makeIndex("pmf_idx_fsetid", "pmfeature", "fsetid")
    makeIndex("mmf_idx_fsetid", "mmfeature", "fsetid")

    makeIndex("fset_idx_chrom", "featureSet", "chrom")
    
}

closeDb <- function(db) dbDisconnect(db)
