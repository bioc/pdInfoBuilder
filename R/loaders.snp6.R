snp6.loadUnitNames <- function(db, unames) {
    dbBeginTransaction(db)
    ## To use an auto-incrementing field via RSQLite, you need
    ## to be careful to pass integer NA's
    df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
    values <- "(:id, :name)"
    sql <- "insert into featureSet (fsetid, man_fsetid) values"
    dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
    dbCommit(db)
}

snp6.loadUnitNames.cnv <- function(db, unames) {
    dbBeginTransaction(db)
    ## To use an auto-incrementing field via RSQLite, you need
    ## to be careful to pass integer NA's
    df <- data.frame(id=rep(as.integer(NA), length(unames)), name=unames)
    values <- "(:id, :name)"
    sql <- "insert into featureSetCNV (fsetid, man_fsetid) values"
    dbGetPreparedQuery(db, paste(sql, values), bind.data=df)
    dbCommit(db)
}

snp6.loadUnits.snp <- function(db, batch, isQc=FALSE) {
  pmfeature <- "pmfeature_tmp"

  ## Don't check PM/MM matching.
  ## SNP 5/6 is PM-only
  batchMat <- do.call(rbind, lapply(batch, readCdfUnitToMat, verify.pmmm=FALSE))

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
  
  ## Insert pm 
  isPm <- as.logical(batchMat[, "ispm"])
  values <- "(:indices, :strand, :allele, :fsetid, :indexpos, :x, :y)"
  sql <- paste("insert into", pmfeature, "values", values)
  dbBeginTransaction(db)
  rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[isPm, ]))
  dbClearResult(rset)
  dbCommit(db)
}

readCdfUnitToMat.cnv <- function(u){
  mat <- do.call("rbind",
                 lapply(u$groups,
                        function(vv){
                          theDir <- groupDirectionToInt(vv$groupdirection)
                          cbind(indices=vv$indices, strand=theDir, x=vv$x, y=vv$y, ispm=vv$ispm)
                        }))
  return(mat)
}


readCdfUnitToMat.cnv.old <- function(u){
  mat <- t(sapply(u$groups,
                  function(vv){
                    theDir <- groupDirectionToInt(vv$groupdirection)
                    c(indices=vv$indice, strand=theDir, x=vv$x, y=vv$y, ispm=vv$ispm)
                  }))
}

snp6.loadUnits.cnv2 <- function(db, batch, isQc=FALSE) {
  pmfeature <- "pmfeatureCNV_tmp"

  ## Don't check PM/MM matching.
  ## SNP 5/6 is PM-only
  batchMat <- do.call(rbind, lapply(batch, readCdfUnitToMat.cnv))

  snp6.loadUnitNames.cnv(db, names(batch))

  ## Find internal featureSet IDs for these features
  batchLens <- sapply(batch, function(x)
                      sum(sapply(x$groups, function(y)
                                 length(y$indices))))
  sql <- paste("select fsetid from featureSetCNV where man_fsetid in (",
               paste('"', names(batch), '"', sep="", collapse=","),
               ") order by fsetid")

  batchIds <- dbGetQuery(db, sql)[[1]]
  batchIds <- rep(batchIds, batchLens)
  batchMat <- cbind(batchMat, fsetid=batchIds)
  
  ## Insert pm 
  isPm <- as.logical(batchMat[, "ispm"])
  values <- "(:indices, :strand, :fsetid, :x, :y)"
  sql <- paste("insert into", pmfeature, "values", values)
  dbBeginTransaction(db)
  rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[isPm, ]))
  dbClearResult(rset)
  dbCommit(db)
}

snp6.loadUnits.cnv <- function(db, batch, isQc=FALSE) {
  pmfeature <- "pmfeatureCNV_tmp"

  ## Don't check PM/MM matching.
  ## SNP 5/6 is PM-only
  batchMat <- do.call(rbind, lapply(batch, readCdfUnitToMat.cnv))

  snp6.loadUnitNames.cnv(db, names(batch))

  ## Find internal featureSet IDs for these features
  batchLens <- sapply(batch, function(x)
                      sum(sapply(x$groups, function(y)
                                 length(y$indices))))
  theNames <- names(batch)
  sets <- split(theNames, rep(1:length(theNames), each=1000, length.out=length(theNames)))
  batchIds <- NULL
  for (i in 1:length(sets)){
    sql <- paste("select fsetid from featureSetCNV where man_fsetid in (",
                 paste('"', sets[[i]], '"', sep="", collapse=","),
                 ") order by fsetid")
    batchIds <- c(batchIds, dbGetQuery(db, sql)[[1]])
  }
##  batchIds <- dbGetQuery(db, sql)[[1]]
  batchIds <- rep(batchIds, batchLens)
  batchMat <- cbind(batchMat, fsetid=batchIds)
  
  ## Insert pm 
  isPm <- as.logical(batchMat[, "ispm"])
  values <- "(:indices, :strand, :fsetid, :x, :y)"
  sql <- paste("insert into", pmfeature, "values", values)
  dbBeginTransaction(db)
  rset <- dbSendPreparedQuery(db, sql, as.data.frame(batchMat[isPm, ]))
  dbClearResult(rset)
  dbCommit(db)
}

snp6.loadUnitsByBatch <- function(db, cdfFile, batch_size=10000,
                                  max_units=NULL, verbose=FALSE) {
  unames <- readCdfUnitNames(cdfFile)
  if (is.null(max_units))
    max_units <- length(unames)

  snp.probes <- grep("^SNP", unames)
  cnv.probes <- grep("^CN", unames)
  whQc <- (1:length(unames))[-c(snp.probes, cnv.probes)]

  ## SNP probes
  done <- FALSE
  while(!done){
    if (length(snp.probes) >= batch_size){
      wanted <- snp.probes[1:batch_size]
    }else{
      wanted <- snp.probes
      done <- TRUE
    }
    vvunits <- readCdf(cdfFile, units=wanted, readGroupDirection=TRUE,
                       readIndices=TRUE, readIsPm=TRUE)
    snp6.loadUnits.snp(db, vvunits)
    if (!done)
      snp.probes <- snp.probes[-(1:batch_size)]
  }
  
  ## CNV probes
  done <- FALSE
  while(!done){
    if (length(snp.probes) >= batch_size){
      wanted <- cnv.probes[1:batch_size]
    }else{
      wanted <- cnv.probes
      done <- TRUE
    }
    vvunits <- readCdf(cdfFile, units=wanted, readGroupDirection=TRUE,
                       readIndices=TRUE, readIsPm=TRUE)
    snp6.loadUnits.cnv(db, vvunits)
    if (!done)
      cnv.probes <- cnv.probes[-(1:batch_size)]
  }
}

snp6.loadAffySeqCsv <- function(db, csvFile, cdfFile, batch_size=5000) {
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
    pmdf <- read.table(con, sep="\t", stringsAsFactors=FALSE, nrows=1, header=FALSE)
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

        values <- "(:fid, :offset, :tstrand, :tallele, :seq)"
        sql <- paste("insert into sequence values", values)
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, sql, bind.data=pmdf)
        dbCommit(db)
    }

}


snp6.buildPdInfoDb <- function(cdfFile, csvFile, csvSeqFile, csvFileCnv,
                               csvSeqFileCnv, dbFile, matFile, matFileCnv,
                               batch_size=10000, verbose=FALSE) {
  ST <- system.time
  printTime <- function(msg, t) {
    if (verbose) {
      m <- paste(msg, "took %.2f sec\n")
      cat(sprintf(m, t))
    }
  }

  db <- snp6.initDb(dbFile)
  t <- ST(snp6.loadUnitsByBatch(db, cdfFile, batch_size=batch_size))
  printTime("loadUnitsByBatch", t[3])
  t <- ST(snp6.loadAffyCsv(db, csvFile, batch_size=batch_size))
  printTime("loadAffyCsv", t[3])
  t <- ST(snp6.loadAffySeqCsv(db, csvSeqFile, cdfFile, batch_size=batch_size))
  printTime("loadAffySeqCsv", t[3])
  t <- ST(snp6.loadAffyCsv.cnv(db, csvFileCnv, batch_size=batch_size))
  printTime("loadAffyCsv-CNV", t[3])
  t <- ST(snp6.loadAffySeqCsv.cnv(db, csvSeqFileCnv, cdfFile, batch_size=batch_size))
  printTime("loadAffySeqCsv-CNV", t[3])
  t <- ST({
    snp6.sortFeatureTables(db)
    snp6.createIndicesDb(db)
    snp6.createTableInfoTable(db)
    snp6.createFeatureTableInfo(db)
  })
  printTime("DB sort, index creation", t[3])
  
  t <- ST({
    seqMat <- createSeqMat(db)
    save(seqMat, file=matFile, compress=TRUE)
    seqMatCnv <- createSeqMat(db)
    save(seqMatCnv, file=matFileCnv, compress=TRUE)
  })
  printTime("sequence matrix", t[3])
  closeDb(db)
}

getFragLength.na29.old <- function(v){
  tmp <- strsplit(v, " /// ")
  t(sapply(tmp,
           function(y){
             at.snp <- strsplit(y, " // ")
             n <- length(at.snp)
             enzymes <- toupper(sapply(at.snp, "[[", 1))
             if(all(c(n==2, enzymes == "---")))
               stop("2 enzymes missing IDs")
             out <- rep(NA, 2)
             i <- grep("NSP", enzymes)
             if (length(i) > 0) out[1] <- as.integer(at.snp[[i]][3])
             i <- grep("STY", enzymes)
             if (length(i) > 0) out[2] <- as.integer(at.snp[[i]][3])
             i <- grep("---", enzymes)
             if (length(i) > 0) out[1] <- as.integer(at.snp[[i]][3])
             out
           }))
}

getFragLength.na29 <- function(v){
  tmp <- strsplit(v, " /// ")
  lens <- sapply(tmp, length)
  if (max(lens) > 2) warning("For, at least, one SNP there is more than 1 record by enzyme.")
  t(sapply(tmp,
           function(y){
             at.snp <- strsplit(y, " // ")
             n <- length(at.snp)
             enzymes <- toupper(sapply(at.snp, "[[", 1))
             if(all(c(n==2, enzymes == "---")))
               stop("2 enzymes missing IDs")
             out <- rep(NA, 2)
             i <- grep("NSP", enzymes)
             if (length(i) > 0) out[1] <- mean(as.integer(sapply(at.snp[i], "[", 3)))
             i <- grep("STY", enzymes)
             if (length(i) > 0) out[2] <- mean(as.integer(sapply(at.snp[i], "[", 3)))

###              i <- grep("---", enzymes)
###              if (length(i) > 0) out[1] <- as.integer(mean(as.integer(sapply(at.snp[i], "[", 3))))

             out
           }))
}


snp6.loadAffyCsv <- function(db, csvFile, batch_size=5000) {
  con <- file(csvFile, open="r")
  on.exit(close(con))

  getFragLength <- function(v){
    tmp <- sapply(strsplit(v, " // "), function(obj) obj[[1]])
    tmp[tmp == "---"] <- NA
    as.integer(tmp)
  }
    
##    wantedCols <- c(1,2,3,4,7,8,10,12,13,14,17) # added 10/14

  ## NA24
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

  ## NA29
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
  
  df <- read.table(con, sep=",", stringsAsFactors=FALSE, nrows=10,
                   na.strings="---", header=TRUE)[, wantedCols]
  header <- gsub(".", "_", names(df), fixed=TRUE)
  names(df) <- header
  df[["Strand_Versus_dbSNP"]] <- as.integer(df[["Strand_Versus_dbSNP"]] == "same")
  df[["Copy_Number_Variation"]] <- as.character(df[["Copy_Number_Variation"]])
  df[["Associated_Gene"]] <- as.character(df[["Associated_Gene"]])
  
  FRAG_COL <- "Fragment_Enzyme_Type_Length_Start_Stop"
  df[[FRAG_COL]] <- as.character(df[[FRAG_COL]])
  tmp.length <- getFragLength.na29(df[[FRAG_COL]])
  
  df[["Strand"]] <- as.integer(ifelse(df[["Strand"]] == "+", SENSE, ANTISENSE))
  df[["frag1"]] <- tmp.length[,1]
  df[["frag2"]] <- tmp.length[,2]
  df[[FRAG_COL]] <- NULL

  db_cols <- c("dbsnp_rs_id", "chrom", "physical_pos", "strand",
               "cytoband", "allele_a", "allele_b", "gene_assoc",
               "dbsnp", "cnv", "fragment_length", "fragment_length2")

  val_holders <- c(":dbSNP_RS_ID", ":Chromosome",
                   ":Physical_Position", ":Strand", ":Cytoband",
                   ":Allele_A", ":Allele_B", ":Associated_Gene",
                   ":Strand_Versus_dbSNP", ":Copy_Number_Variation",
                   ":frag1", ":frag2")

  exprs <- paste(db_cols, " = ", val_holders, sep="", collapse=", ")
  sql <- paste("update featureSet set ", exprs,
               "where man_fsetid = :Probe_Set_ID")

  dbBeginTransaction(db)
  dbGetPreparedQuery(db, sql, bind.data=df)
  dbCommit(db)

  ## Now do the rest in batches
  done <- FALSE
  while (!done) {
    df <- read.table(con, sep=",", stringsAsFactors=FALSE,
                     nrows=batch_size, na.strings="---",
                     header=FALSE)[, wantedCols]
    if (nrow(df) < batch_size) {
      done <- TRUE
      if (nrow(df) == 0)
        break
    }
    names(df) <- header
    df[["Strand_Versus_dbSNP"]] <- as.integer(df[["Strand_Versus_dbSNP"]] == "same")
    df[["Copy_Number_Variation"]] <- as.character(df[["Copy_Number_Variation"]])
    df[["Associated_Gene"]] <- as.character(df[["Associated_Gene"]])
    
    df[[FRAG_COL]] <- as.character(df[[FRAG_COL]])
    tmp.length <- getFragLength.na29(df[[FRAG_COL]])
    df[["Strand"]] <- as.integer(ifelse(df[["Strand"]] == "+", SENSE, ANTISENSE))
    df[["frag1"]] <- tmp.length[,1]
    df[["frag2"]] <- tmp.length[,2]
    
    dbBeginTransaction(db)
    dbGetPreparedQuery(db, sql, bind.data=df)
    dbCommit(db)

  }
}

snp6.loadAffyCsv.cnv <- function(db, csvFile, batch_size=5000) {
  getPAR <- function(theDF){
    ## PAR: Pseudo-Autosomal Region
    ##   0: No / 1: PAR1 / 2: PAR2
    PAR <- rep(as.integer(0), nrow(theDF))
    PAR[theDF[["ChrX_pseudo_autosomal_region_1"]] == 1] <- 1
    PAR[theDF[["ChrX_pseudo_autosomal_region_2"]] == 2] <- 2
    theDF[["XPAR"]] <- as.integer(PAR)
    theDF[["ChrX_pseudo_autosomal_region_1"]] <- NULL
    theDF[["ChrX_pseudo_autosomal_region_2"]] <- NULL
    theDF
  }

  ## CNV probes
  con <- file(csvFile, open="r")
  on.exit(close(con))
  
  wantedCols <- c(1, 2, 3, 4, 5, 6, 7, 8, 10, 13, 11)
  df <- read.table(con, sep=",", stringsAsFactors=FALSE, nrows=10,
                   na.strings="---", header=TRUE)[, wantedCols]
  header <- gsub(".", "_", names(df), fixed=TRUE)
  names(df) <- header

  df[["Associated_Gene"]] <- as.character(df[["Associated_Gene"]])
  
  FRAG_COL <- "Fragment_Enzyme_Type_Length_Start_Stop"
  
  df <- getPAR(df)
  df[["Strand"]] <- as.integer(ifelse(df[["Strand"]] == "+", SENSE, ANTISENSE))
  df[[FRAG_COL]] <- as.character(df[[FRAG_COL]])
  tmp.length <- getFragLength.na29(df[[FRAG_COL]])
  df[["frag1"]] <- tmp.length[,1]
  df[["frag2"]] <- tmp.length[,2]

  df[["Copy_Number_Variation"]] <- as.character(df[["Copy_Number_Variation"]])

  
  db_cols <- c("chrom", "chrom_start", "chrom_stop", "strand",
               "cytoband", "gene_assoc", "xpar", "fragment_length",
               "fragment_length2", "cnv")
  
  val_holders <- c(":Chromosome", ":Chromosome_Start",
                   ":Chromosome_Stop", ":Strand", ":Cytoband",
                   ":Associated_Gene", ":XPAR", ":frag1", ":frag2",
                   ":Copy_Number_Variation")
  
  exprs <- paste(db_cols, " = ", val_holders, sep="", collapse=", ")
  sql <- paste("update featureSetCNV set ", exprs,
               "where man_fsetid = :Probe_Set_ID")

  dbBeginTransaction(db)
  dbGetPreparedQuery(db, sql, bind.data=df)
  dbCommit(db)

  ## Now do the rest in batches
  done <- FALSE
  while (!done) {
    df <- read.table(con, sep=",", stringsAsFactors=FALSE,
                     nrows=batch_size, na.strings="---",
                     header=FALSE)[, wantedCols]
    if (nrow(df) < batch_size) {
      done <- TRUE
      if (nrow(df) == 0)
        break
    }
    names(df) <- header
    df[["Associated_Gene"]] <- as.character(df[["Associated_Gene"]])

    df <- getPAR(df)
    df[["Strand"]] <- as.integer(ifelse(df[["Strand"]] == "+", SENSE, ANTISENSE))
    df[[FRAG_COL]] <- as.character(df[[FRAG_COL]])
    tmp.length <- getFragLength.na29(df[[FRAG_COL]])
    df[["frag1"]] <- tmp.length[,1]
    df[["frag2"]] <- tmp.length[,2]

    df[["Copy_Number_Variation"]] <- as.character(df[["Copy_Number_Variation"]])

    dbBeginTransaction(db)
    dbGetPreparedQuery(db, sql, bind.data=df)
    dbCommit(db)
  }


}

snp6.loadAffySeqCsv.cnv <- function(db, csvFile, cdfFile, batch_size=5000) {
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
    header <- c("fset.name", "x", "y", "offset", "seq", "tstrand",
                "type", "chromosome", "strand", "probe_start_pos",
                "in_xpar1", "in_xpar2", "nsp_frag_start", "nsp_frag_end",
                "sty_frag_start", "sty_frag_end")
    wanted <- 1:7
    header <- header[wanted]
    done <- FALSE
    pmdf <- read.table(con, sep="\t", stringsAsFactors=FALSE, nrows=1, header=FALSE)[, wanted]
    while (!done) {
        pmdf <- read.table(con, sep="\t", stringsAsFactors=FALSE,
                           nrows=batch_size, na.strings="---",
                           header=FALSE)[, wanted]
        if (nrow(pmdf) < batch_size) {
            done <- TRUE
            if (nrow(pmdf) == 0)
              break
        }
        names(pmdf) <- header
        pmdf[["fid"]] <- xy2i(pmdf[["x"]], pmdf[["y"]])
        N <- nrow(pmdf)
        values <- "(:fid, :offset, :tstrand, :seq)"
        sql <- paste("insert into sequenceCNV values", values)
        dbBeginTransaction(db)
        dbGetPreparedQuery(db, sql, bind.data=pmdf)
        dbCommit(db)
    }

}

