initDb.affyST <- function(dbname) {
    db <- dbConnect(dbDriver("SQLite"), dbname)

    ## Set page size
    dbGetQuery(db, setPageSizeSql)

    ## Create tables

    dbGetQuery(db, createAffySTFeatureSetSql)
    dbGetQuery(db, sprintf(createAffySTFeatureSql, "pmfeature_tmp"))
    dbGetQuery(db, sprintf(createAffySTFeatureSql, "mmfeature_tmp"))
    dbGetQuery(db, sprintf(createAffySTFeatureSql, "qcpmfeature"))
    dbGetQuery(db, sprintf(createAffySTFeatureSql, "qcmmfeature"))
    dbGetQuery(db, sprintf(createAffySTPm_MmSql, "pm_mm"))
    dbGetQuery(db, sprintf(createAffySTPm_MmSql, "qcpm_qcmm"))
    dbGetQuery(db, createAffySTSequenceSql)

    ## Create index
    ## NOTE: might be more efficient to create this after loading,
    ## but current perf is ok.

    sql <- 'create index man_fsetid_idx on featureSet ("man_fsetid")'
    dbGetQuery(db, sql)
    db
}

sortFeatureTables.affyST <- function(db) {
    dbGetQuery(db, sprintf(createAffySTFeatureSql, "pmfeature"))
    dbGetQuery(db, sprintf(createAffySTFeatureSql, "mmfeature"))

    ## Reorder XXfeature tables

    ## on other arrays, fid is primary key (therefore, tbls ordered by
    ## fid, no matter what)
    
    ## fillSql <- paste("insert into %s select * from %s order by",
    ##                  "fsetid, atom")
    fillSql <- paste("insert into %s select * from %s order by fid")
    
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


createIndicesDb.affyST <- function(db) {
    makeIndex <- function(name, t, cols) {
        sql <- paste("create index", name, "on", t,
                     paste("(", paste(cols, collapse=","), ")"))
        dbGetQuery(db, sql)
    }

    ## Create DB indices and fix ordering
    makeIndex("pmf_idx_fsetid", "pmfeature", "fsetid")
    makeIndex("mmf_idx_fsetid", "mmfeature", "fsetid")

##    makeIndex("fset_idx_chrom", "featureSet", "chrom")
    makeIndex("fset_idx_fsetid", "featureSet", "fsetid")

    ## finally, run analyze (SQLite specific?)
    dbGetQuery(db, "analyze")
}


createTableInfoTable.affyST <- function(db, verbose=FALSE) {
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


createFeatureTableInfo.affyST <- function(db, tname) {
    return(FALSE)
    ## FIXME: add code to determine offsets of sorted
    ## strand and allele
}


closeDb <- function(db) dbDisconnect(db)
