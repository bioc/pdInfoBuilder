initDb <- function(dbname) {
    db <- dbConnect(dbDriver("SQLite"), dbname)

    ## Set page size
    dbGetQuery(db, setPageSizeSql)

    ## Create tables
    ## BC: Soon we need to add a table for the control probes
    dbGetQuery(db, createSnpFeatureSetSql)

    dbGetQuery(db, sprintf(createSnpFeatureSql, "pmfeature_tmp"))

    dbGetQuery(db, sprintf(createSnpFeatureSql, "mmfeature_tmp"))

    dbGetQuery(db, sprintf(createSnpFeatureSql, "qcpmfeature"))
    dbGetQuery(db, sprintf(createSnpFeatureSql, "qcmmfeature"))

    dbGetQuery(db, sprintf(createSnpPm_MmSql, "pm_mm"))
    dbGetQuery(db, sprintf(createSnpPm_MmSql, "qcpm_qcmm"))

    dbGetQuery(db, createSnpSequenceSql)

    ## Create index
    ## NOTE: might be more efficient to create this after loading,
    ## but current perf is ok.
    sql <- 'create index man_fsetid_idx on featureSet ("man_fsetid")'
    dbGetQuery(db, sql)
    db
}

sortFeatureTables <- function(db) {
    dbGetQuery(db, sprintf(createSnpFeatureSql, "pmfeature"))
    dbGetQuery(db, sprintf(createSnpFeatureSql, "mmfeature"))

    ## Reorder XXfeature tables
    fillSql <- paste("insert into %s select * from %s order by",
                     "fsetid, allele, strand, pos")
    dbBeginTransaction(db)
    dbGetQuery(db, sprintf(fillSql, "pmfeature", "pmfeature_tmp"))
    dbGetQuery(db, sprintf(fillSql, "mmfeature", "mmfeature_tmp"))
    dbCommit(db)
    ## drop temp tables
    dbBeginTransaction(db)
    dbGetQuery(db, "drop table pmfeature_tmp")
    dbGetQuery(db, "drop table mmfeature_tmp")
    dbCommit(db)
}


createIndicesDb <- function(db) {
    makeIndex <- function(name, t, cols) {
        sql <- paste("create index", name, "on", t,
                     paste("(", paste(cols, collapse=","), ")"))
        dbGetQuery(db, sql)
    }

    ## Create DB indices and fix ordering
    makeIndex("pmf_idx_fsetid", "pmfeature", "fsetid")
    makeIndex("mmf_idx_fsetid", "mmfeature", "fsetid")

    makeIndex("fset_idx_chrom", "featureSet", "chrom")
    makeIndex("fset_idx_fsetid", "featureSet", "fsetid")

    ## finally, run analyze (SQLite specific?)
    dbGetQuery(db, "analyze")
}


createTableInfoTable <- function(db, verbose=FALSE) {
    tables <- dbListTables(db)
    counts <- integer(length(tables))
    sql <- "select count(*) from %s"
    for (i in seq(along=counts)) {
        if (verbose)
          cat("counting rows in ", tables[i], "\n")
        counts[i] <- dbGetQuery(db, sprintf(sql, tables[i]))[[1]][1]
    }

    df <- data.frame(tbl=tables, row_count=counts,
                     stringsAsFactors=FALSE)
    dbWriteTable(db, "table_info", df, row.names=FALSE)
}


createFeatureTableInfo <- function(db, tname) {
    return(FALSE)
    ## FIXME: add code to determine offsets of sorted
    ## strand and allele
}


closeDb <- function(db) dbDisconnect(db)
