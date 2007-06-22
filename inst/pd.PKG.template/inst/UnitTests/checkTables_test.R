pdinfo <- @PDINFONAME@

test_table_set <- function() {
    ##
    ## this test assumes we know what tables have been defined.  if we
    ## change the schema, this test needs to be updated
    ##
    alltabs.sort <- sort(c("featureSet", "mmfeature", "pm_mm", "pmfeature",
                           "qcmmfeature", "qcpm_qcmm", "qcpmfeature",
                           "sequence", "sqlite_stat1", "table_info"))
    dd  <- pdinfo@getdb()
    ll  <- sort(dbListTables(dd))
    checkEquals(alltabs.sort, ll)
}

test_table_info <- function() {
    ##
    ## this test verifies that the table_info table knows the correct
    ## number of rows for featureSet, and that we can get the featureSet
    ## table into R to count them for 50K chip this is reasonable --
    ## actually a little slow but maybe worth it
    ##
    dd  <- pdinfo@getdb()
    tinfo  <- dbGetQuery(dd, "select * from table_info")
    nrfs  <- as.numeric(tinfo[ tinfo$tbl == "featureSet", "row_count" ])
    fstab  <- dbGetQuery(dd, "select * from featureSet")
    checkEquals(nrfs, nrow(fstab))
}

