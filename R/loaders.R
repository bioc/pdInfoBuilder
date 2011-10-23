loadUnitNames <- function(db, unames) {
    dbBeginTransaction(db)
    ## To use an auto-incrementing field via RSQLite, you need
    ## to be careful to pass integer NA's
    df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
    values <- "(:id, :name)"
    sql <- "insert into featureSet (fsetid, man_fsetid) values"
    dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
    dbCommit(db)
}

loadUnits <- function(db, batch, isQc=FALSE) {
    pmfeature <- "pmfeature_tmp"
    mmfeature <- "mmfeature_tmp"
    pmmm <- "pm_mm"
    if (isQc) {
        pmfeature <- "qcpmfeature"
        mmfeature <- "qcmmfeature"
        pmmm <- "qcpm_qcmm"
    }
    batchMat <- do.call(rbind, lapply(batch, readCdfUnitToMat))

    loadUnitNames(db, names(batch))

    ## Find internal featureSet IDs for these features
    batchLens <- sapply(batch, function(x)
                        sum(sapply(x$groups, function(y)
                                   length(y$indices))))
    sql <- paste("select fsetid from featureSet where man_fsetid in (",
                 paste('"', names(batch), '"', sep="", collapse=","),
                 ") order by fsetid")
    batchIds <- dbGetQuery(db, sql)[[1]]
    batchIds <- rep(batchIds, batchLens)
    batchMat <- cbind(batchMat, fsetid=batchIds)

    ## Insert pm and mm into respective tables
    isPm <- as.logical(batchMat[, "ispm"])
    values <- "(:indices, :strand, :allele, :fsetid, :indexpos, :x, :y)"
    sql <- paste("insert into", pmfeature, "values", values)
    dbBeginTransaction(db)
    rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[isPm, ]))
    dbClearResult(rset)

    sql <- paste("insert into", mmfeature, "values", values)
    rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[!isPm, ]))
    dbClearResult(rset)
    dbCommit(db)

    ## Insert pm <--> mm link
    values <- "(:pmi, :mmi)"
    sql <- paste("insert into", pmmm, "values", values)
    dbBeginTransaction(db)
    rset <- dbSendPreparedQuery(db, sql,
                                data.frame(pmi=batchMat[isPm, "indices"],
                                           mmi=batchMat[!isPm, "indices"]))
    dbClearResult(rset)
    dbCommit(db)
}


loadUnitsByBatch <- function(db, cdfFile, batch_size=10000,
                             max_units=NULL, verbose=FALSE) {
    unames <- readCdfUnitNames(cdfFile)
    if (is.null(max_units))
      max_units <- length(unames)
    offset <- 1
    whQc <- grep("^AFFX", unames)
    if (length(whQc)) {                 # load all QC at once
        qcunits <- readCdf(cdfFile, units=whQc, readGroupDirection=TRUE,
                           readIndices=TRUE, readIsPm=TRUE)
        loadUnits(db, qcunits, isQc=TRUE)
        offset <- max(whQc) + 1
    }
    extra <- (length(unames)-length(whQc)) %% batch_size
    num_batches <- (length(unames)-length(whQc)) %/% batch_size
    if (extra != 0)
      num_batches <- num_batches + 1
    done <- FALSE
    while (!done) {
        end <- min(offset + batch_size, max_units)
        if (end == max_units)
          done <- TRUE
        wanted <- seq.int(offset, end)
        offset <- offset + batch_size + 1
        vvunits <- readCdf(cdfFile, units=wanted, readGroupDirection=TRUE,
                           readIndices=TRUE, readIsPm=TRUE)
        loadUnits(db, vvunits)
    }
}

loadAffySeqCsv <- function(db, csvFile, cdfFile, batch_size=5000) {
    cdfHeader <- readCdfHeader(cdfFile)
    ncol <- cdfHeader$ncol
    xy2i <- function(x, y) x + 1 + y * ncol

    complementBase <- function(x, special=FALSE){
      bases <- c("A", "C", "G", "T")
      if (!special)
        comp <- c("T", "G", "C", "A")
      else
        comp <- c("G", "T", "A", "C")
      comp[match(x, bases)]
    }

    con <- file(csvFile, open="r")
    on.exit(close(con))
    header <- c("fset.name", "x", "y", "offset", "seq", "tstrand", "type",
                "tallele")
    done <- FALSE
    while (!done) {
        pmdf <- read.table(con, sep="\t", stringsAsFactors=FALSE,
                           nrows=batch_size, na.strings="---",
                           header=FALSE)
        if (nrow(pmdf) < batch_size) {
            done <- TRUE
            if (nrow(pmdf) == 0)
              break
        }
        names(pmdf) <- header
        pmdf[["fid"]] <- xy2i(pmdf[["x"]], pmdf[["y"]])
        N <- nrow(pmdf)

        ## process MMs
        mmSql <- paste("select mm_fid, pm_fid from pm_mm where pm_fid in (",
                       paste(pmdf[["fid"]], collapse=","), ")")
        pairedIds <- dbGetQuery(db, mmSql)
        foundIdIdx <- match(pairedIds[["pm_fid"]], pmdf[["fid"]], 0)
        mmdf <- pmdf[foundIdIdx, ]
        mmdf[["fid"]] <-  pairedIds[["mm_fid"]]

        ## Assuming 25mers
        midbase <- substr(mmdf$seq, 13, 13)
        types <- apply(table(mmdf$fset.name, mmdf$tallele)>0, 1,
                       function(v) paste(c("A", "C", "G", "T")[v], collapse=""))
        types <- rep(types, as.integer(table(mmdf$fset.name)))
        isSpecial <- (types == "AT" | types == "CG") & mmdf$offset == 0
        rm(types)
        midbase[isSpecial] <- complementBase(midbase[isSpecial], T)
        midbase[!isSpecial] <- complementBase(midbase[!isSpecial])
        rm(isSpecial)
        mmdf$seq <- paste(substr(mmdf$seq, 1, 12), midbase,
                          substr(mmdf$seq, 14, 25), sep="")
        rm(midbase)
        ## end MM seq

        values <- "(:fid, :offset, :tstrand, :tallele, :seq)"
        sql <- paste("insert into sequence values", values)
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, sql, bind.data=pmdf)
        dbGetPreparedQuery(db, sql, bind.data=mmdf)
        dbCommit(db)
    }

}


buildPdInfoDb <- function(cdfFile, csvFile, csvSeqFile, dbFile, matFile,
                          batch_size=10000, verbose=FALSE) {

    ST <- system.time
    printTime <- function(msg, t) {
        if (verbose) {
            m <- paste(msg, "took %.2f sec\n")
            cat(sprintf(m, t))
        }
    }

    db <- initDb(dbFile)

    t <- ST(loadUnitsByBatch(db, cdfFile, batch_size=batch_size))
    printTime("loadUnitsByBatch", t[3])
    t <- ST(loadAffyCsv(db, csvFile, batch_size=batch_size))
    printTime("loadAffyCsv", t[3])
    t <- ST(loadAffySeqCsv(db, csvSeqFile, cdfFile, batch_size=batch_size))
    printTime("loadAffySeqCsv", t[3])
    t <- ST({
        sortFeatureTables(db)
        createIndicesDb(db)
        createTableInfoTable(db)
        createFeatureTableInfo(db)
    })
    printTime("DB sort, index creation", t[3])

    t <- ST({
        seqMat <- createSeqMat(db)
        save(seqMat, file=matFile, compress='xz')
    })
    printTime("sequence matrix", t[3])
    closeDb(db)
}

loadAffyCsv <- function(db, csvFile, batch_size=5000) {
    con <- file(csvFile, open="r")
    on.exit(close(con))

    ## This is the header of annotation files NA24
    ## The needed fields have (*)
    ## 01 - (*) - Probe Set ID
    ## 02 - (*) - Affy SNP ID
    ## 03 - (*) - dbSNP RS ID
    ## 04 - (*) - Chromosome
    ## 05 - (*) - Physical Position
    ## 06 - (*) - Strand
    ## 07 - ( ) - ChrX pseudo-autosomal region 1
    ## 08 - (*) - Cytoband
    ## 09 - ( ) - Flank
    ## 10 - (*) - Allele A
    ## 11 - (*) - Allele B
    ## 12 - (*) - Associated Gene
    ## 13 - ( ) - Genetic Map
    ## 14 - ( ) - Microsatellite
    ## 15 - (*) - Fragment Enzyme Length Start Stop
    ## 16 - ( ) - Allele Frequencies
    ## 17 - ( ) - Heterozygous Allele Frequencies
    ## 18 - ( ) - Number of individuals/Number of chromosomes
    ## 19 - ( ) - In Hapmap
    ## 20 - (*) - Strand Versus dbSNP
    ## 21 - (*) - Copy Number Variation
    ## 22 - ( ) - Probe Count
    ## 23 - ( ) - ChrX pseudo-autosomal region 2
    ## 24 - ( ) - In Final List
    ## 25 - ( ) - Minor Allele
    ## 26 - ( ) - Minor Allele Frequency
    ## wantedCols <- c(1, 2, 3, 4, 5, 6, 8, 10, 11, 12, 15, 20, 21)

    ## This is the header of annotation files NA29
    ##############################################################
    ## CHANGES:
    ##     A) Affymetrix got rid of 'Affy SNP ID'
    ##     B) Affymetrix changed the format of
    ##        "Fragment Enzyme Length Start Stop" to
    ##        "Fragment Enzyme Type Length Start Stop"
    ##     C) Affymetrix added GC and OMIM
    ##############################################################
    ## ACTIONS: 13 / OCT / 2009
    ##     A) Remove 'Affy SNP ID' from SQL schema
    ##     B) Change function to get the right fields
    ##     C) None needed
    ##############################################################
    
    ## The needed fields have (*)
    ## 01 - (*) - Probe Set ID
    ## 02 - (*) - dbSNP RS ID
    ## 03 - (*) - Chromosome
    ## 04 - (*) - Physical Position
    ## 05 - (*) - Strand
    ## 06 - ( ) - ChrX pseudo-autosomal region 1
    ## 07 - (*) - Cytoband
    ## 08 - ( ) - Flank
    ## 09 - (*) - Allele A
    ## 10 - (*) - Allele B
    ## 11 - (*) - Associated Gene
    ## 12 - ( ) - Genetic Map
    ## 13 - ( ) - Microsatellite
    ## 14 - (*) - Fragment Enzyme Type Length Start Stop
    ## 15 - ( ) - Allele Frequencies
    ## 16 - ( ) - Heterozygous Allele Frequencies
    ## 17 - ( ) - Number of individuals/Number of chromosomes
    ## 18 - ( ) - In Hapmap
    ## 19 - (*) - Strand Versus dbSNP
    ## 20 - (*) - Copy Number Variation
    ## 21 - ( ) - Probe Count
    ## 22 - ( ) - ChrX pseudo-autosomal region 2
    ## 23 - ( ) - In Final List
    ## 24 - ( ) - Minor Allele
    ## 25 - ( ) - Minor Allele Frequency
    ## 26 - ( ) - % GC
    ## 27 - ( ) - OMIM
    wantedCols <- c(1:5, 7, 9:11, 14, 19:20)

    df <- read.table(con, sep=",", stringsAsFactors=FALSE,
                     na.strings="---", header=TRUE)[, wantedCols]
    header <- gsub(".", "_", names(df), fixed=TRUE)
    names(df) <- header
    df[["Strand_Versus_dbSNP"]] <- as.integer(df[["Strand_Versus_dbSNP"]] == "same")
    df[["Copy_Number_Variation"]] <- as.character(df[["Copy_Number_Variation"]])
    df[["Associated_Gene"]] <- as.character(df[["Associated_Gene"]])

    insertInFragmentLengthTable(db, 'fragmentLength',
                                df[['Fragment_Enzyme_Type_Length_Start_Stop']],
                                df[['Probe_Set_ID']])
    df[['Fragment_Enzyme_Type_Length_Start_Stop']] <- NULL

    db_cols <- c("dbsnp_rs_id", "chrom", "physical_pos", "strand",
                 "cytoband", "allele_a", "allele_b", "gene_assoc",
                 "dbsnp", "cnv")

    val_holders <- c(":dbSNP_RS_ID", ":Chromosome",
                     ":Physical_Position", ":Strand", ":Cytoband",
                     ":Allele_A", ":Allele_B", ":Associated_Gene",
                     ":Strand_Versus_dbSNP", ":Copy_Number_Variation")

    exprs <- paste(db_cols, " = ", val_holders, sep="", collapse=", ")
    sql <- paste("update featureSet set ", exprs,
                 "where man_fsetid = :Probe_Set_ID")

    dbBeginTransaction(db)
    dbGetPreparedQuery(db, sql, bind.data=df)
    dbCommit(db)
}

