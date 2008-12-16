#############################################################################
#############################################################################
###
### Functions for importing Human Exon Array CSV annotation files into a
### database.
###
### Affymetrix provides 2 files:
###   1. The Probe Set CSV file (HuEx-1_0-st-v2.na21.hg18.probeset.csv)
###   2. The Transcript CSV file (HuEx-1_0-st-v2.na21.hg18.transcript.csv)
###
### This file is divided in 7 sections:
###   A. Dictionary manipulation.
###   B. General purpose low-level SQL helper functions.
###   C. CSV field description and DB schema.
###   D. Objects and functions shared by E. and F.
###   E. Importation of the Transcript CSV file (167M, aka the "small" file).
###   F. Importation of the Probe Set CSV file (475M, aka the "big" file).
###   G. Importation of the 2 CSV files (Transcript + Probe Set).
### 
### WARNING: This is a WORK IN PROGRESS!!! (if pdInfoBuilder had a NAMESPACE,
### nothing should be exported from this file for now)
###
#############################################################################



### =========================================================================
### A. Dictionary manipulation (should probably go somewhere else).
### -------------------------------------------------------------------------

### Create a new dictionary.
new.dict <- function(keys)
{
    dict <- new.env(hash=TRUE, parent=emptyenv())
    for (key in keys)
        assign(key, NULL, envir=dict)
    dict
}

dict.keys <- function(dict)
{
    ls(dict, all.names=TRUE)
}

dict.rm <- function(dict, keys)
{
    remove(list=keys, envir=dict)
}

dict.lengths <- function(dict)
{
    ## Slow
    #sapply(dict.keys(dict), function(key) length(dict[[key]]))
    ## 2x faster
    lengths <- unlist(eapply(dict, length, all.names=TRUE))
    lengths[dict.keys(dict)]
}

## Return a named character vector.
dict.flatten <- function(dict, length1vals.only=FALSE)
{
    ans <- character(0)
    for (key in dict.keys(dict)) {
        val <- dict[[key]]
        if (is.null(val) || (length1vals.only && length(val) != 1))
            next
        if (!is.character(val))
            val <- as.character(val)
        names(val) <- rep(key, length.out=length(val))
        ans <- c(ans, val)
    }
    ans
}

dict.toString <- function(dict)
{
    ans <- sapply(dict.keys(dict),
                  function(key)
                      paste(key, ":", paste(dict[[key]], collapse=","), sep=""))
    paste("[", paste(ans, collapse="|"), "]", sep="")
}



### =========================================================================
### B. General purpose low-level SQL helper functions (should probably go
###    somewhere else).
### -------------------------------------------------------------------------

### 'vals' must be a character vector with NAs for the SQL NULL.
### Return a character vector of the same length as 'vals' (and with the
### same names if 'vals' has names).
toSQLValues <- function(vals, col2type)
{
    nvals <- length(vals)
    if (nvals == 0)
        stop("'vals' is empty")
    cols <- names(vals)
    if (is.null(cols)) {
        if (nvals != length(col2type)) {
            cat("          vals  = ", vals, "\n", sep="|")
            cat("names(col2type) = ", names(col2type), "\n", sep="|")
            cat("      col2type  = ", col2type, "\n", sep="|")
            stop("when unamed, 'vals' must be of the same length as 'col2type'")
        }
        types <- col2type
    } else {
        if (any(is.na(cols)) || any(cols == ""))
            stop("'vals' has invalid names")
        if (any(duplicated(cols)))
            stop("'vals' has duplicated names")
        types <- col2type[cols]
        if (any(is.na(types)))
            stop("no type found for some of the values in 'vals'")
    }
    ## From here, 'vals' and 'types' have the same length
    ## and the mapping between them is positional.
    for (i in seq_len(nvals)) {
        val <- vals[i]
        if (is.na(val)) {
            vals[i] <- "NULL"
            next
        }
        if (types[i] == "INTEGER") {
            x <- as.numeric(val)
            xx <- as.integer(x)
            if (is.na(xx) || xx != x)
                stop("vals[", i, "]=\"", val, "\" not an integer")
            vals[i] <- xx
            next
        }
        vals[i] <- paste("'", gsub("'", "''", val, fixed=TRUE), "'", sep="")
    }
    vals
}

### ID generation.
### WARNING: Those functions are unsafe in a "multi db" context because when
### several db connections are open at the same time (e.g. 'conn1' and 'conn2')
### then 'tablename' is not a unique key anymore.

.DB.CURRENT.IDS <- new.env(hash=TRUE, parent=emptyenv())

.db.set.id <- function(tablename, id=0)
{
    assign(tablename, id, envir=.DB.CURRENT.IDS)
}

.db.next.id <- function(tablename)
{
    if (!exists(tablename, envir=.DB.CURRENT.IDS))
        stop("no current id for table \"", tablename, "\"")
    id <- get(tablename, envir=.DB.CURRENT.IDS)
    id <- id + 1
    .db.set.id(tablename, id)
    id
}

### 'conn' must be a DBIConnection object, a filename or a file object
.dbGetQuery <- function(conn, sql)
{
    if (is(conn, "DBIConnection"))
        dbGetQuery(conn, sql)
    else
        cat(sql, ";\n", file=conn, sep="", append=TRUE)
}

dbCreateTable <- function(conn, tablename, col2type, col2key)
{
    col2type[names(col2key)] <- paste(col2type[names(col2key)], col2key, sep=" ")
    sql <- paste(names(col2type), col2type, sep=" ")
    sql <- paste(sql, collapse=", ")
    sql <- paste("CREATE TABLE ", tablename, " (", sql, ")", sep="")
    .dbGetQuery(conn, sql)
    .db.set.id(tablename)
}

### 'row' must be a character vector with or without names.
dbInsertRow <- function(conn, tablename, row, col2type, is.fullrow=TRUE)
{
    sqlvals <- try(toSQLValues(row, col2type), silent=TRUE)
    if (is(sqlvals, "try-error"))
        stop("table=\"", tablename, "\": ", sqlvals)
    sql <- paste(sqlvals, collapse=", ")
    sql <- paste("VALUES (", sql, ")", sep="")
    cols <- names(row)
    if (!is.null(cols) && !identical(cols, names(col2type))) {
        if (is.fullrow)
            stop("table=\"", tablename, "\": 'names(row)' and 'names(col2type)' are not identical")
        cols <- paste(cols, collapse=", ")
        sql <- paste("(", cols, ") ", sql, sep="")
    }
    sql <- paste("INSERT INTO ", tablename, " ", sql, sep="")
    .dbGetQuery(conn, sql)
}

### 'row' must be a named character vector.
### Return NULL if 'row' is not found in 'tablename'.
### Note that it wouldn't make sense to pass a row such that 'row[unique_col]'
### is NA.
dbGetThisRow <- function(conn, tablename, unique_col, row, col2type)
{
    if (is.null(names(row)))
        stop("'row' must be a named character vector")
    unique_sqlval <- try(toSQLValues(row[unique_col], col2type), silent=TRUE)
    if (is(unique_sqlval, "try-error"))
        stop("table=\"", tablename, "\": ", unique_sqlval)
    sql <- paste("SELECT * FROM ", tablename, " WHERE ",
                 unique_col, "=", unique_sqlval, " LIMIT 1", sep="")
    data <- dbGetQuery(conn, sql)
    ## 'data' data frame can only have 0 or 1 row
    if (nrow(data) == 0)
        return(NULL)
    row0 <- unlist(data)
    for (col in names(row)) {
        if (is.na(row[col]))
            next
        if (row[col] != row0[col]) {
            cat("row = ", row, "\n", sep="|")
            cat("row0 = ", row0, "\n", sep="|")
            stop("row in table ", tablename, " has different values")
        }
    }
    row0
}

dbInsertDataFrame <- function(conn, tablename, data, col2type, verbose=FALSE)
{
    cols <- names(col2type)
    if (!identical(sort(names(data)), sort(cols)))
        stop("cols in data frame 'data' don't match cols in table \"", tablename, "\"")
    values_template <- paste(paste(":", cols, sep=""), collapse=", ")
    sql_template <- paste("INSERT INTO ", tablename, " VALUES (", values_template, ")", sep="")
    if (verbose)
        cat("Inserting ", nrow(data), " rows into table \"", tablename, "\"... ", sep="")
    dbBeginTransaction(conn)
    on.exit(dbCommit(conn))
    dbGetPreparedQuery(conn, sql_template, bind.data=data)
    if (verbose)
        cat("OK\n")
}



### =========================================================================
### C. CSV field description and DB schema.
###
###    3 sub-sections:
###      C.a. Transcript CSV file: sub-fields names for each multipart field.
###      C.b. Probe Set CSV file: sub-fields names for each multipart field.
###      C.c. The DB schema.
###
###    The original field description is provided in
###      HuEx-1_0-st-v2.na21.hg18.AFFX_README.NetAffx-CSV-Files.txt
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### C.a. Transcript CSV file: sub-fields names for each multipart field.
###

### In HuEx-1_0-st-v2.na21.hg18.transcript.csv, the "gene_assignment" multipart
### field is NA ("---") 89% of the time (for 278540 transcript cluster IDs).
### The rest of the time (33828 lines), it has between 1 and 254 parts with an
### average of 3.26 parts (254 parts for transcript cluster ID "3581637").
### Across all lines, the total number of unique values per sub-field is:
###   - 62220 for "accession"
###   - 27033 for "entrez_gene_id"
### The "accession" sub-field is always unique within the parts of the same
### field value except for 34 lines. Among those 34 lines, 30 contain only
### 1 duplicated accession. The 4 lines with more than 1 duplicated are:
###
###   transcript_cluster_ID | duplicated accessions in the "gene_assignment" field
###   ----------------------|-----------------------------------------------------
###   3078948               | XM_936040 (2), XM_934915 (2), XM_935763 (2)
###   3287928               | AY704155 (2), AY704147 (2), AY704156 (2), AY704150 (2)
###   3584443               | NR_001294 (3), NR_001292 (2), NR_001290 (2)
###   3655450               | AK122733 (2), AK094769 (2), AK097386 (2)
###
TRsubfields.gene_assignment <- c(
    "accession",
    "gene_symbol",
    "gene_title",
    "cytoband",
    "entrez_gene_id"
)

### In HuEx-1_0-st-v2.na21.hg18.transcript.csv, the "mrna_assignment" multipart
### field is NA ("---") 76% of the time (for 238561 transcript cluster IDs).
### The rest of the time (73807 lines), it has between 1 and 5096 parts with an
### average of 4.8 parts (5096 parts for transcript cluster ID "3581637", same
### as for the max number of parts in "gene_assignment").
### Across all lines, the total number of unique values in the "accession"
### sub-field is 237713.
### The "accession" sub-field is always unique within the parts of the same
### field value except for 12 lines. In each of those 12 lines, only 1
### accession is repeated ("NC_001807" in 10 of them + "AK097625" and
### "GENSCAN00000038476" in the other 2 lines). Among those 12 lines there are
### only 5 lines where the number of occurences of this repeated accession
### is >= 3: for transcript cluster IDs "2315301" (4), "2315315" (3),
### "3531553" (3), "4037638" (3) and "4037708" (6). In those 5 lines, the
### repeated accession is always "NC_001807".
###
TRsubfields.mrna_assignment <- c(
    "accession",
    "source_name",
    "description",
    "assignment_seqname",
    "assignment_score",
    "assignment_coverage",
    "direct_probes",
    "possible_probes",
    "xhyb"
)

TRsubfields.swissprot <- c(
    "accession",
    "swissprot_accession"
)

TRsubfields.unigene <- c(
    "accession",
    "unigene_id",
    "unigene_expr"
)

TRsubfields.GO_biological_process <- c(
    "accession",
    "GO_id",
    "GO_term",
    "GO_evidence"   # can be "---" (-> NA in R, -> NULL in SQL)
)

TRsubfields.pathway <- c(
    "accession",
    "source",
    "pathway_name"
)

TRsubfields.protein_domains <- c(
    "accession",
    "source",
    "accession_or_domain_name",
    "domain_description"
)

TRsubfields.protein_families <- c(
    "accession",
    "source",
    "family_accession",
    "family_description"
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### C.b. Probe Set CSV file: sub-fields names for each multipart field.
###

PBSsubfields.gene_assignment <- c(
    "accession",
    "gene_symbol"
)

PBSsubfields.mrna_assignment <- c(
    "accession",
    "assignment_seqname",
    "assignment_score",
    "direct_probes",
    "possible_probes",
    "xhyb"
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### C.c. The DB schema.
###
### The DB schema used for importing the annotations for the Affymetrix Human
### Exon Array probe sets and transcript clusters is called AFFYHUEX_DB.
### This sub-section provides a definition of this schema at the R level.
###

### The "transcript_cluster" table.
transcript_cluster_desc <- list(
    col2type=c(
        transcript_cluster_ID="INTEGER",    # PRIMARY KEY
        seqname="TEXT",
        strand="CHAR(1)",
        start="INTEGER",
        stop="INTEGER",
        total_probes="INTEGER"
    ),
    col2key=c(
        transcript_cluster_ID="PRIMARY KEY"
    )
)

### The "gene" table.
###
### Note that we can't put a UNIQUE constraint on "gene_symbol".
### For example in HuEx-1_0-st-v2.na21.hg18.transcript.csv, line 10430
### (transcript_cluster_ID=2564225), the "gene_assignment" field (multipart)
### contains the following parts:
###
###      accession | gene_symbol | gene_title                    | cytoband | entrez_gene_id
###   -------------|-------------|-------------------------------|----------|---------------
###      NM_002339 | LSP1        | lymphocyte-specific protein 1 |  11p15.5 |           4046
###      XM_946374 | LSP1        | lymphocyte-specific protein 1 |   2p11.1 |         654342
###   NM_001013254 | LSP1        | lymphocyte-specific protein 1 |  11p15.5 |           4046
###   NM_001013255 | LSP1        | lymphocyte-specific protein 1 |  11p15.5 |           4046
###   NM_001013253 | LSP1        | lymphocyte-specific protein 1 |  11p15.5 |           4046
###       AK093859 | LSP1        | lymphocyte-specific protein 1 |  11p15.5 |           4046
###       AK092071 | LSP1        | lymphocyte-specific protein 1 |  11p15.5 |           4046
###
gene_desc <- list(
    col2type=c(
        entrez_gene_id="INTEGER",   # PRIMARY KEY
        gene_symbol="TEXT",
        gene_title="TEXT",
        cytoband="TEXT"
    ),
    col2key=c(
        entrez_gene_id="PRIMARY KEY"
    )
)

### The "mrna" table.
mrna_desc <- list(
    col2type=c(
        `_mrna_id`="INTEGER",       # internal id (PRIMARY KEY)
        accession="VARCHAR(100)"
    ),
    col2key=c(
        `_mrna_id`="PRIMARY KEY",
        accession="UNIQUE"
    )
)

### The "mrna2gene" table.
mrna2gene_desc <- list(
    col2type=c(
        `_mrna_id`="INTEGER",       # REFERENCES mrna(_mrna_id)
        entrez_gene_id="INTEGER"    # REFERENCES gene(entrez_gene_id)
    ),
    col2key=c(
        `_mrna_id`="REFERENCES mrna(_mrna_id)",
        entrez_gene_id="REFERENCES gene(entrez_gene_id)"
    )
)

### The "TR2mrna" table.
### TODO: Add a UNIQUE constraint on (transcript_cluster_ID, _mrna_id).
TR2mrna_desc <- list(
    col2type=c(
        `_TR2mrna_id`="INTEGER",
        transcript_cluster_ID="INTEGER",    # REFERENCES transcript_cluster(transcript_cluster_ID)
        `_mrna_id`="INTEGER"        # REFERENCES mrna(_mrna_id)
    ),
    col2key=c(
        `_TR2mrna_id`="PRIMARY KEY",
        transcript_cluster_ID="NOT NULL REFERENCES transcript_cluster(transcript_cluster_ID)",
        `_mrna_id`="NOT NULL REFERENCES mrna(_mrna_id)"
    )
)

### The "TR2mrna_details" table.
TR2mrna_details_desc <- list(
    col2type=c(
        source_name="TEXT",
        description="TEXT",
        assignment_seqname="TEXT",
        assignment_score="INTEGER",
        assignment_coverage="INTEGER",
        direct_probes="INTEGER",
        possible_probes="INTEGER",
        xhyb="INTEGER",
        `_TR2mrna_id`="INTEGER"         # REFERENCES TR2mrna(_TR2mrna_id)
    ),
    col2key=c(
        `_TR2mrna_id`="NOT NULL REFERENCES TR2mrna(_TR2mrna_id)"
    )
)

### The "swissprot" table.
### TODO: Add a UNIQUE constraint on (swissprot_accession, _mrna_id).
swissprot_desc <- list(
    col2type=c(
        swissprot_accession="TEXT", 
        `_mrna_id`="INTEGER"        # REFERENCES mrna(_mrna_id)
    ),
    col2key=c(
        `_mrna_id`="NOT NULL REFERENCES mrna(_mrna_id)"
    )
)

### The "unigene" table.
### TODO: Add a UNIQUE constraint on (unigene_id, _mrna_id).
unigene_desc <- list(
    col2type=c(
        unigene_id="TEXT",
        unigene_expr="TEXT",
        `_mrna_id`="INTEGER"        # REFERENCES mrna(_mrna_id)
    ),
    col2key=c(
        `_mrna_id`="NOT NULL REFERENCES mrna(_mrna_id)"
    )
)

### The "GO_biological_process" table.
GO_biological_process_desc <- list(
    col2type=c(
        GO_id="TEXT",
        GO_term="TEXT",
        GO_evidence_code="TEXT",    # can be NULL
        entrez_gene_id="INTEGER"    # REFERENCES gene(entrez_gene_id)
    ),
    col2key=c(
        `entrez_gene_id`="REFERENCES gene(entrez_gene_id)"
    )
)

### The "pathway" table.
### TODO: Add a UNIQUE constraint on (source, _mrna_id).
pathway_desc <- list(
    col2type=c(
        source="TEXT",
        pathway_name="TEXT",
        `_mrna_id`="INTEGER"        # REFERENCES mrna(_mrna_id)
    ),
    col2key=c(
        `_mrna_id`="NOT NULL REFERENCES mrna(_mrna_id)"
    )
)

### The "protein_domains" table.
### TODO: Add a UNIQUE constraint on (source, accession_or_domain_name).
protein_domains_desc <- list(
    col2type=c(
        `_protein_domains_id`="INTEGER",  # PRIMARY KEY
        source="TEXT",
        accession_or_domain_name="TEXT",  # can be "" (empty string)
        domain_description="TEXT"
    ),
    col2key=c(
        `_protein_domains_id`="PRIMARY KEY"
    )
)

### The "TR2mrna_2_protein_domains" table.
### TODO: Add a UNIQUE constraint on (_TR2mrna_id, _protein_domains_id).
TR2mrna_2_protein_domains_desc <- list(
    col2type=c(
        `_TR2mrna_id`="INTEGER",
        `_protein_domains_id`="INTEGER"
    ),
    col2key=c(
        `_TR2mrna_id`="NOT NULL REFERENCES TR2mrna(_TR2mrna_id)",
        `_protein_domains_id`="NOT NULL REFERENCES protein_domains(_protein_domains_id)"
    )
)

### The "protein_families" table.
### TODO: Add a UNIQUE constraint on (source, _mrna_id).
protein_families_desc <- list(
    col2type=c(
        source="TEXT",
        family_accession="TEXT",
        family_description="TEXT",
        `_mrna_id`="INTEGER"        # REFERENCES mrna(_mrna_id)
    ),
    col2key=c(
        `_mrna_id`="NOT NULL REFERENCES mrna(_mrna_id)"
    )
)

### The "probeset" table.
probeset_desc <- list(
    col2type=c(
        probeset_ID="INTEGER",              # PRIMARY KEY
        seqname="TEXT",
        strand="CHAR(1)",
        start="INTEGER",
        stop="INTEGER",
        probe_count="INTEGER",
        transcript_cluster_ID="INTEGER",    # REFERENCES transcript_cluster(transcript_cluster_ID)
        exon_id="INTEGER",
        psr_id="INTEGER",
        crosshyb_type="INTEGER",
        number_independent_probes="INTEGER",
        number_cross_hyb_probes="INTEGER",
        number_nonoverlapping_probes="INTEGER",
        level="TEXT",
        bounded="INTEGER",
        noBoundedEvidence="INTEGER",
        has_cds="INTEGER",
        fl="INTEGER",
        mrna="INTEGER",
        est="INTEGER",
        vegaGene="INTEGER",
        vegaPseudoGene="INTEGER",
        ensGene="INTEGER",
        sgpGene="INTEGER",
        exoniphy="INTEGER",
        twinscan="INTEGER",
        geneid="INTEGER",
        genscan="INTEGER",
        genscanSubopt="INTEGER",
        mouse_fl="INTEGER",
        mouse_mrna="INTEGER",
        rat_fl="INTEGER",
        rat_mrna="INTEGER",
        microRNAregistry="INTEGER",
        rnaGene="INTEGER",
        mitomap="INTEGER",
        probeset_type="TEXT"
    ),
    col2key=c(
        probeset_ID="PRIMARY KEY",
        transcript_cluster_ID="REFERENCES transcript_cluster(transcript_cluster_ID)"
    )
)

### The "PBS2mrna" table.
### TODO: Add a UNIQUE constraint on (probeset_ID, _mrna_id).
PBS2mrna_desc <- list(
    col2type=c(
        #`_PBS2mrna_id`="INTEGER",
        assignment_seqname="TEXT",
        assignment_score="INTEGER",
        direct_probes="INTEGER",
        possible_probes="INTEGER",
        xhyb="INTEGER",
        probeset_ID="INTEGER",          # REFERENCES probeset(probeset_ID)
        `_mrna_id`="INTEGER"            # REFERENCES mrna(_mrna_id)
    ),
    col2key=c(
        #`_PBS2mrna_id`="PRIMARY KEY",
        probeset_ID="NOT NULL REFERENCES probeset(probeset_ID)",
        `_mrna_id`="NOT NULL REFERENCES mrna(_mrna_id)"
    )
)

### Global schema (17 tables).
AFFYHUEX_DB_schema <- list(
    transcript_cluster=transcript_cluster_desc,
    gene=gene_desc,
    mrna=mrna_desc,
    mrna2gene=mrna2gene_desc,
    TR2mrna=TR2mrna_desc,
    TR2mrna_details=TR2mrna_details_desc,
    swissprot=swissprot_desc,
    unigene=unigene_desc,
    GO_biological_process=GO_biological_process_desc,
    GO_cellular_component=GO_biological_process_desc,
    GO_molecular_function=GO_biological_process_desc,
    pathway=pathway_desc,
    protein_domains=protein_domains_desc,
    TR2mrna_2_protein_domains=TR2mrna_2_protein_domains_desc,
    protein_families=protein_families_desc,
    probeset=probeset_desc,
    PBS2mrna=PBS2mrna_desc
)



### =========================================================================
### D. Objects and functions shared by E. and F. 
### -------------------------------------------------------------------------


dbCreateTables.AFFYHUEX_DB <- function(conn, tablenames=NULL)
{
    if (is.null(tablenames))
        tablenames <- names(AFFYHUEX_DB_schema)
    for (tablename in tablenames) {
        col2type <- AFFYHUEX_DB_schema[[tablename]]$col2type
        col2key <- AFFYHUEX_DB_schema[[tablename]]$col2key
        dbCreateTable(conn, tablename, col2type, col2key)
    }
}

### TODO: Implement this!
dbInitTableIds.AFFYHUEX_DB <- function(conn)
{
    stop("not ready")
    for (tablename in names(AFFYHUEX_DB_schema)) {
        #id0 :== "SELECT max(id) FROM tablename" + 1
        .db.set.id(tablename, id0)
    }
}

### Utilities for reporting errors/warnings in the CSV data

.CSVIMPORT.VARS <- new.env(hash=TRUE, parent=emptyenv())

.CSVimport.var <- function(objname, objval, defval)
{
    if (!is.null(objval))
        return(assign(objname, objval, envir=.CSVIMPORT.VARS))
    if (exists(objname, envir=.CSVIMPORT.VARS))
        return(get(objname, envir=.CSVIMPORT.VARS))
    defval
}

.CSVimport.verbose <- function(verbose=NULL)
{
    .CSVimport.var("verbose", verbose, defval=FALSE)
}

.CSVimport.fullmode <- function(fullmode=NULL)
{
    .CSVimport.var("fullmode", fullmode, defval=FALSE)
}

.CSVimport.infile <- function(file=NULL)
{
    if (.CSVimport.verbose() && !is.null(file))
        cat("- CSV_file=\"", file, "\"\n", sep="")
    .CSVimport.var("infile", file)
}

.CSVimport.dataline_nb <- function(nb=NULL)
{
    if (.CSVimport.verbose() && !is.null(nb))
        cat("-- dataline_nb=", nb, "\n", sep="")
    .CSVimport.var("dataline_nb", nb)
}

.CSVimport.line_ID <- function(line_ID=NULL)
{
    if (.CSVimport.verbose() && !is.null(line_ID))
        cat("-- line_ID=", line_ID, "\n", sep="")
    .CSVimport.var("line_ID", line_ID)
}

.CSVimport.field <- function(field=NULL)
{
    if (.CSVimport.verbose() && !is.null(field))
        cat("--- field=\"", field, "\"\n", sep="")
    .CSVimport.var("field", field)
}

.CSVimport.logfile <- function(file=NULL)
{
    .CSVimport.var("logfile", file)
}

csv_current_pos <- function()
{
    paste("CSV_file=\"", .CSVimport.infile(), "\"",
          ", dataline_nb=", .CSVimport.dataline_nb(), 
          ", [line_ID=", .CSVimport.line_ID(), "]",
          ", field=\"", .CSVimport.field(), "\"",
          sep="")
}
          
data_error <- function(msg)
{
    msg <- paste("*** <DATAERROR>\n",
                 "*** In ", csv_current_pos(), ":\n*** ", msg, "\n",
                 "*** </DATAERROR>\n", sep="")
    cat("\n", msg, "\n", sep="")
    stop(msg)
}

data_warning <- function(msg)
{
    msg <- paste("*** <DATAWARNING>\n",
                 "*** In ", csv_current_pos(), ":\n*** ", msg, "\n",
                 "*** </DATAWARNING>\n", sep="")
    cat("\n", msg, "\n", sep="")
    warning(msg)
}

getAffyHuExArrayAnnCSVHeader <- function(data)
{
    header <- names(data)
    header[header == "probeset_id"] <- "probeset_ID"
    header[header == "transcript_cluster_id"] <- "transcript_cluster_ID"
    header
}

### Convert a multipart field value to a matrix.
mvalToMat <- function(mval, subfields, min.nsubfields=length(subfields))
{
    ncol <- length(subfields)
    if (is.na(mval)) {
        mat <- matrix(data=character(0), ncol=ncol)
    } else {
        vals <- strsplit(mval, " /// ", fixed=TRUE)[[1]]
        rows <- strsplit(vals, " // ", fixed=TRUE)
        for (i in seq_len(length(rows))) {
            nsubfields <- length(rows[[i]])
            if (nsubfields < min.nsubfields || nsubfields > ncol) {
                msg <- "wrong number of subfields"
                data_error(msg)
            }
            if (nsubfields < ncol)
                length(rows[[i]]) <- ncol
        }
        data <- unlist(rows)
        data[data == "---"] <- NA
        mat <- matrix(data=data, ncol=ncol, byrow=TRUE)
    }
    colnames(mat) <- subfields
    mat
}

### The vectorized form of "mvalToMat".
mvalsToMats <- function(mvals, subfields)
{
    mvals[mvals == "---"] <- NA
    lapply(mvals, function(x) mvalToMat(x, subfields))
}

### Authority: http://www.geneontology.org/GO.evidence.shtml
GO_evidence_codes <- c(
    IC="Inferred by Curator",
    IDA="Inferred from Direct Assay",
    IEA="Inferred from Electronic Annotation",
    IEP="Inferred from Expression Pattern",
    IGC="Inferred from Genomic Context",
    IGI="Inferred from Genetic Interaction",
    IMP="Inferred from Mutant Phenotype",
    IPI="Inferred from Physical Interaction",
    ISS="Inferred from Sequence or Structural Similarity",
    NAS="Non-traceable Author Statement",
    ND="No biological Data available",
    RCA="inferred from Reviewed Computational Analysis",
    TAS="Traceable Author Statement",
    NR="Not Recorded"
)

GO_evidence_codes_lower <- tolower(GO_evidence_codes)

replaceGOEvidenceByCode <- function(mat)
{
    evidences <- mat[ , "GO_evidence"]
    notNA <- !is.na(evidences)
    notNA_evidences <- evidences[notNA]
    pos <- match(notNA_evidences, GO_evidence_codes_lower)
    if (any(is.na(pos))) {
        msg <- paste("unknown evidence '", notNA_evidences[is.na(pos)][1], "'", sep="")
        data_error(msg)
    }
    evidences[notNA] <- names(GO_evidence_codes_lower)[pos]
    mat[ , "GO_evidence"] <- evidences
    colnames(mat)[colnames(mat) == "GO_evidence"] <- "GO_evidence_code"
    mat
}

### Return a list of matrices.
### 'acc2id' must be a named character vector.
splitMatrix <- function(mat, acc2id)
{
    if (ncol(mat) < 2)
        stop("won't split a matrix with less than 2 cols")
    if (colnames(mat)[1] != "accession")
        stop("first 'mat' col name must be \"accession\"")
    if (!is.character(acc2id) || is.null(names(acc2id)))
        stop("'acc2id' is not a named character vector")
    accessions <- mat[ , 1]
    id2submat <- new.dict(unique(acc2id))
    for (accession in unique(accessions)) {
        submat <- mat[accessions == accession, -1 , drop=FALSE]
        id <- acc2id[accession]
        submat0 <- id2submat[[id]]
        if (is.null(submat0)) {
            id2submat[[id]] <- submat
            next
        }
        if (!identical(submat, submat0)) {
            show(mat)
            show(acc2id)
            data_warning("no way to properly split this matrix (see above)")
        }
    }
    id2submat
}

### Comparison of the data contained in a character matrix ('mat') and a data
### frame ('dat').
### 'mat' _must_ have col names. If 'dat' is not a 0-row data frame then it
### _must_  have col names too and they _must_ be exactly the same as 'mat'
### col names. It's an error if one of those "_must_" is not satisfied.
haveTheSameData <- function(mat, dat)
{
    if (is.null(colnames(mat)))
        stop("'mat' has no col names")
    if (nrow(dat) == 0)
        return(nrow(mat) == 0)
    if (!identical(colnames(mat), colnames(dat)))
        stop("'mat' and 'dat' have different col names")
    if (nrow(mat) != nrow(dat))
        return(FALSE)
    ## Convert 'dat' to a character matrix (we can't just use as.matrix
    ## because of the infamous "format" feature).
    mat2 <- do.call(cbind, args=lapply(dat, function(x) if (is.character(x)) x else as.character(x)))
    ## Find the order of the rows in the 2 matrices.
    ii1 <- do.call(order, args=lapply(colnames(mat), function(col) mat[ , col]))
    ii2 <- do.call(order, args=lapply(colnames(mat2), function(col) mat2[ , col]))
    ## Compare row by row.
    for (i in seq_len(nrow(mat)))
        if (!identical(mat[ii1[i], ], mat2[ii2[i], ]))
            return(FALSE)
    return(TRUE)
}

dbInsert_multipart_data <- function(conn, tablename, mat, insres)
{
    if (colnames(mat)[1] != "accession")
        stop("first col name in 'mat' must be \"accession\"")
    acc2ids <- insres$acc2ids

    ## Remove orphan parts
    accessions <- mat[ , 1]
    ignored_parts <- !(accessions %in% dict.keys(acc2ids))
    if (any(ignored_parts)) {
        msg <- paste(accessions[ignored_parts], collapse=",")
        msg <- paste("Ignoring orphan parts: ", msg, sep="")
        msg <- paste(msg, "\nacc2ids = ", dict.toString(acc2ids), sep="")
        data_warning(msg)
        mat <- mat[ignored_parts, , drop=FALSE]
    }
    if (length(acc2ids) == 0)
        return()

    ## Remove ambiguous parts
    accessions <- mat[ , 1] # might have changed
    nb_ids <- dict.lengths(acc2ids)
    multiid_accessions <- names(nb_ids)[nb_ids > 1]
    ignored_parts <- accessions %in% multiid_accessions
    if (any(ignored_parts)) {
        msg <- paste(accessions[ignored_parts], collapse=",")
        msg <- paste("Ignoring ambiguous parts: ", msg, sep="")
        msg <- paste(msg, "\nacc2ids = ", dict.toString(acc2ids), sep="")
        data_warning(msg)
        mat <- mat[ignored_parts, , drop=FALSE]
    }

    acc2id <- dict.flatten(acc2ids, length1vals.only=TRUE)
    id2submat <- splitMatrix(mat, acc2id)
    new_ids <- insres$new_ids
    col2type <- AFFYHUEX_DB_schema[[tablename]]$col2type
    cols0 <- paste(colnames(mat)[-1], collapse=",")
    link0 <- names(col2type)[length(col2type)]
    sql0 <- paste("SELECT ", cols0, " FROM ", tablename, " WHERE ", link0, "=", sep="")
    for (id in dict.keys(id2submat)) {
        submat <- id2submat[[id]]
        if (id %in% new_ids) {
            if (is.null(submat))
                next
            for (i in seq_len(nrow(submat))) {
                row1 <- c(submat[i, ], id)
                names(row1) <- names(col2type)
                res <- try(dbInsertRow(conn, tablename, row1, col2type), silent=TRUE)
                if (is(res, "try-error"))
                    stop("In ", csv_current_pos(), ":\n", res, "\n")
            }
        } else {
            sql <- paste(sql0, id, sep="")
            data0 <- dbGetQuery(conn, sql)
            if (is.null(submat)) {
                if (nrow(data0) == 0)
                    next
            } else {
                if (haveTheSameData(submat, data0))
                    next
            }
            msg <- paste("parts in CSV field don't match records already in ",
                         "table \"", tablename, "\" for ", link0, "=", id, sep="")
            data_error(msg)
        }
    }
}



### =========================================================================
### E. Importation of the Transcript CSV file (167M, aka the "small" file).
### -------------------------------------------------------------------------


dbInsertRows.gene <- function(conn, genes)
{
    accessions <- unique(genes[ , "accession"])
    acc2ids <- new.dict(accessions)
    new_ids <- character(0)
    col2type <- gene_desc$col2type
    for (i in seq_len(nrow(genes))) {
        gene <- genes[i, ]
        row1 <- gene[names(col2type)]
        id <- row1[["entrez_gene_id"]] # [[ ]] to get rid of the name
        row0 <- try(dbGetThisRow(conn, "gene", "entrez_gene_id", row1, col2type), silent=TRUE)
        if (is(row0, "try-error"))
            stop("In ", csv_current_pos(), ":\n", row0, "\n")
        if (is.null(row0)) {
            res <- try(dbInsertRow(conn, "gene", row1, col2type), silent=TRUE)
            if (is(res, "try-error"))
                stop("In ", csv_current_pos(), ":\n", res, "\n")
            new_ids <- c(new_ids, id)
        }
        accession <- gene["accession"]
        acc2ids[[accession]] <- c(acc2ids[[accession]], id)
    }
    list(acc2ids=acc2ids, new_ids=new_ids)
}

dbInsertRows.mrna <- function(conn, accessions)
{
    accessions <- unique(accessions)
    acc2ids <- new.dict(accessions)
    new_ids <- character(0)
    col2type <- mrna_desc$col2type
    for (accession in accessions) {
        row1 <- c(NA, accession)
        names(row1) <- names(col2type)
        row0 <- try(dbGetThisRow(conn, "mrna", "accession", row1, col2type), silent=TRUE)
        if (is(row0, "try-error"))
            stop("In ", csv_current_pos(), ":\n", row0, "\n")
        if (is.null(row0)) {
            id <- as.character(.db.next.id("mrna"))
            row1["_mrna_id"] <- id
            res <- try(dbInsertRow(conn, "mrna", row1, col2type), silent=TRUE)
            if (is(res, "try-error"))
                stop("In ", csv_current_pos(), ":\n", res, "\n")
            new_ids <- c(new_ids, id)
        } else {
            id <- row0[["_mrna_id"]] # [[ ]] to get rid of the name
        }
        acc2ids[[accession]] <- id
    }
    list(acc2ids=acc2ids, new_ids=new_ids)
}

dbInsertRows.mrna2gene <- function(conn, mrna_insres, gene_insres)
{
    mrna_acc2ids <- mrna_insres$acc2ids
    gene_acc2ids <- gene_insres$acc2ids
    col2type <- mrna2gene_desc$col2type
    for (accession in dict.keys(mrna_acc2ids)) {
        mrna.id <- mrna_acc2ids[[accession]] # always of length 1
        if (!(mrna.id %in% mrna_insres$new_ids))
            next
        for (entrez_gene_id in gene_acc2ids[[accession]]) {
            row1 <- c(mrna.id, entrez_gene_id)
            names(row1) <- names(col2type)
            res <- try(dbInsertRow(conn, "mrna2gene", row1, col2type), silent=TRUE)
            if (is(res, "try-error"))
                stop("In ", csv_current_pos(), ":\n", res, "\n")
        }
    }
}

dbInsertRows.TR2mrna <- function(conn, transcript_cluster_ID, mrna_acc2ids)
{
    acc2ids <- new.dict(dict.keys(mrna_acc2ids))
    col2type <- TR2mrna_desc$col2type
    for (accession in dict.keys(mrna_acc2ids)) {
        id <- .db.next.id("TR2mrna")
        row1 <- c(id, transcript_cluster_ID, mrna_acc2ids[[accession]])
        names(row1) <- names(col2type)
        res <- try(dbInsertRow(conn, "TR2mrna", row1, col2type), silent=TRUE)
        if (is(res, "try-error"))
            stop("In ", csv_current_pos(), ":\n", res, "\n")
        acc2ids[[accession]] <- id
    }
    acc2ids
}

dbInsertRows.TR2mrna_details <- function(conn, mrna_assignment, TR2mrna_acc2ids)
{
    col2type <- TR2mrna_details_desc$col2type
    for (i in seq_len(nrow(mrna_assignment))) {
        accession <- mrna_assignment[i, "accession"]
        TR2mrna.id <- TR2mrna_acc2ids[[accession]]
        row1 <- mrna_assignment[i, names(col2type)[1:8]]
        row1 <- c(row1, TR2mrna.id)
        names(row1) <- names(col2type)
        res <- try(dbInsertRow(conn, "TR2mrna_details", row1, col2type), silent=TRUE)
        if (is(res, "try-error"))
            stop("In ", csv_current_pos(), ":\n", res, "\n")
    }
}

dbImportLine.AFFYHUEX_DB.Transcript <- function(conn, dataline)
{
    transcript_cluster_ID <- dataline[["transcript_cluster_ID"]] # [[ ]] to get rid of the name
    .CSVimport.line_ID(transcript_cluster_ID)
    #if (!is(conn, "DBIConnection"))
    #    .dbGetQuery(conn, paste("-- transcript_cluster_ID ", transcript_cluster_ID, sep=""))
    dataline[dataline == "---"] <- NA

    ## Extract the simple fields

    row1 <- dataline[names(transcript_cluster_desc$col2type)]
    dbInsertRow(conn, "transcript_cluster", row1, transcript_cluster_desc$col2type)

    ## Extract and insert the "gene_assignment" and "mrna_assignment" data

    field <- "gene_assignment"
    .CSVimport.field(field)
    gene_assignment <- mvalToMat(dataline[field], TRsubfields.gene_assignment)
    gene_insres <- dbInsertRows.gene(conn, gene_assignment)

    field <- "mrna_assignment"
    .CSVimport.field(field)
    mrna_assignment <- mvalToMat(dataline[field], TRsubfields.mrna_assignment)
    accessions <- mrna_assignment[ , "accession"]
    mrna_insres <- dbInsertRows.mrna(conn, accessions)
    dbInsertRows.mrna2gene(conn, mrna_insres, gene_insres)
    TR2mrna_acc2ids <- dbInsertRows.TR2mrna(conn, transcript_cluster_ID, mrna_insres$acc2ids)
    dbInsertRows.TR2mrna_details(conn, mrna_assignment, TR2mrna_acc2ids)

    if (!.CSVimport.fullmode())
        return()

    ## Extract and insert the "swissprot" data

    field <- "swissprot"
    .CSVimport.field(field)
    mat <- mvalToMat(dataline[field], TRsubfields.swissprot)
    dbInsert_multipart_data(conn, field, mat, mrna_insres)

    ## Extract and insert the "unigene" data
    ## TODO: Like the "GO" stuff, the "unigene" stuff needs to be linked to the
    ## genes not to the mrnas. This requires a change in the db schema.

    field <- "unigene"
    .CSVimport.field(field)
    mat <- mvalToMat(dataline[field], TRsubfields.unigene, min.nsubfields=2)
    dbInsert_multipart_data(conn, field, mat, mrna_insres)

    ## Extract and insert the "GO" data
    ##
    ## Note that those GO annotations can be ambiguous.
    ## For example in HuEx-1_0-st-v2.na21.hg18.transcript.csv, line 16411
    ## (transcript_cluster_ID="2718075"), the "gene_assignment" field (multipart)
    ## contains the following parts:
    ##
    ##      accession | gene_symbol | gene_title                    | cytoband | entrez_gene_id
    ##   -------------|-------------|-------------------------------|----------|---------------
    ##   NM_001040448 | DEFB131     | defensin, beta 131            |   4p16.1 |         644414
    ##      XM_938410 | DEFB131     | defensin, beta 131            |   4p16.1 |         644414
    ##      XM_938410 | LOC649335   | similar to Beta-defensin 131..|          |         649335
    ##   ...
    ## When this happens, then the GO annotations are ambiguous (this is
    ## because Affymetrix has choosen to link GO terms to mrnas, not to genes).
    ## For example, on the same line (transcript_cluster_ID="2718075"), the
    ## "GO_biological_process" field (multipart) contains:
    ##
    ##      accession | GO_id | GO_term                                | GO_evidence_code
    ##   -------------|-------|----------------------------------------|-----------------
    ##   NM_001040448 |  9613 | response to pest, pathogen or parasite | IEA
    ##   NM_001040448 | 42742 | defense response to bacteria           | IEA
    ##      XM_938410 |  9613 | response to pest, pathogen or parasite | IEA
    ##      XM_938410 | 42742 | defense response to bacteria           | IEA
    ##   NM_001037804 |  9613 | response to pest, pathogen or parasite | IEA
    ##   NM_001037804 | 42742 | defense response to bacteria           | IEA
    ##
    ## 3rd and 4th parts are linked to which gene?

    field <- "GO_biological_process"
    .CSVimport.field(field)
    mat <- mvalToMat(dataline[field], TRsubfields.GO_biological_process)
    mat <- replaceGOEvidenceByCode(mat)
    dbInsert_multipart_data(conn, field, mat, gene_insres)

    field <- "GO_cellular_component"
    .CSVimport.field(field)
    mat <- mvalToMat(dataline[field], TRsubfields.GO_biological_process)
    mat <- replaceGOEvidenceByCode(mat)
    dbInsert_multipart_data(conn, field, mat, gene_insres)

    field <- "GO_molecular_function"
    .CSVimport.field(field)
    mat <- mvalToMat(dataline[field], TRsubfields.GO_biological_process)
    mat <- replaceGOEvidenceByCode(mat)
    dbInsert_multipart_data(conn, field, mat, gene_insres)

    ## Extract and insert the "pathway" data

    field <- "pathway"
    .CSVimport.field(field)
    mat <- mvalToMat(dataline[field], TRsubfields.pathway)
    dbInsert_multipart_data(conn, field, mat, mrna_insres)

    return() # that's all for now

    ## The code below is not ready...
    protein_domains <- mvalToMat(dataline["protein_domains"],
                                         TRsubfields.protein_domains, min.nsubfields=3)
    dbInsert_multipart_data(conn, "protein_domains", protein_domains,
                                   acc2id, new_accessions, transcript_cluster_ID)

    protein_families <- mvalToMat(dataline["protein_families"], TRsubfields.protein_families)
    dbInsert_multipart_data(conn, "protein_families", protein_families,
                                   acc2id, new_accessions, transcript_cluster_ID)
}

dbImportData.AFFYHUEX_DB.Transcript <- function(conn, csv_file, seqname, nrows=-1)
{
    if (!is.null(csv_file)) {
        .CSVimport.infile(csv_file)
        csv_con <- file(csv_file, open="r")
        on.exit(close(csv_con))
        dataline_nb <- 0
        while (nrows == -1 || dataline_nb < nrows) {
            if (dataline_nb == 0) {
                data <- read.table(csv_con, header=TRUE, sep=",", quote="\"",
                                   nrows=1, stringsAsFactors=FALSE)
                names(data) <- header <- getAffyHuExArrayAnnCSVHeader(data)
            } else {
                data <- read.table(csv_con, header=FALSE, sep=",", quote="\"",
                                   col.names=header, nrows=1, stringsAsFactors=FALSE)
            }
            if (nrow(data) == 0)
                break
            dataline_nb <- dataline_nb + 1
            .CSVimport.dataline_nb(dataline_nb)
            dataline <- unlist(data[1, ])
            dbImportLine.AFFYHUEX_DB.Transcript(conn, dataline)
        }
    } else {
        infile <- paste("tr_", seqname, ".rda", sep="")
        .CSVimport.infile(infile)
        tmp_envir <- new.env(parent=emptyenv())
        load(infile, envir=tmp_envir)
        data <- get("tr_table", envir=tmp_envir)
        for (dataline_nb in seq_len(nrow(data))) {
            .CSVimport.dataline_nb(dataline_nb)
            dataline <- unlist(data[dataline_nb, ])
            dbImportLine.AFFYHUEX_DB.Transcript(conn, dataline)
        }
    }
}



### =========================================================================
### F. Importation of the Probe Set CSV file (475M, aka the "big" file).
### -------------------------------------------------------------------------


### Affymetrix README file (HuEx-1_0-st-v2.na21.hg18.AFFX_README.NetAffx-CSV-Files.txt)
### contains this statement:
###   Only limited details are provided for each assigned mRNA. For more
###   information about these mRNAs, see the corresponding entry in the
###   transcript cluster CSV file (use the probe set's transcript cluster
###   ID to join).
### This strongly suggests that we should be able to assume the following:
###   mRNAs linked to a given probeset_ID are also linked to the
###   corresponding transcript_cluster_ID
### Unfortunately, we can't!
### For example, in HuEx-1_0-st-v2.na21.hg18.probeset.csv,
### probeset_ID="3935513" is linked to mRNAs "GENSCAN00000012720" and
### "GENSCAN00000040267". But in HuEx-1_0-st-v2.na21.hg18.transcript.csv,
### transcript_cluster_ID="3935512" (the transcript cluster corresponding
### to probeset_ID="3935513") is linked to mRNA "GENSCAN00000012720" only.
### Even worse: HuEx-1_0-st-v2.na21.hg18.transcript.csv contains no
### information about mRNA "GENSCAN00000040267"!
dbInsertRows.PBS2mrna <- function(conn, mrna_assignment, probeset_ID)
{
    col2type <- PBS2mrna_desc$col2type
    for (i in seq_len(nrow(mrna_assignment))) {
        accession <- mrna_assignment[i, "accession"]
        row1 <- c(NA, accession)
        names(row1) <- names(mrna_desc$col2type)
        row0 <- try(dbGetThisRow(conn, "mrna", "accession", row1, mrna_desc$col2type), silent=TRUE)
        if (is(row0, "try-error"))
            stop("In ", csv_current_pos(), ":\n", row0, "\n")
        if (is.null(row0)) {
            msg <- paste("table \"mrna\" has no record ",
                         "with accession ", accession, sep="")
            data_warning(msg)
            mrna.id <- .db.next.id("mrna")
            row1["_mrna_id"] <- mrna.id
            res <- try(dbInsertRow(conn, "mrna", row1, mrna_desc$col2type), silent=TRUE)
            if (is(res, "try-error"))
                stop("In ", csv_current_pos(), ":\n", res, "\n")
        } else {
            mrna.id <- row0["_mrna_id"]
        }
        row1 <- mrna_assignment[i, names(col2type)[1:5]]
        row1 <- c(row1, probeset_ID, mrna.id)
        names(row1) <- names(col2type)
        res <- try(dbInsertRow(conn, "PBS2mrna", row1, col2type), silent=TRUE)
        if (is(res, "try-error"))
            stop("In ", csv_current_pos(), ":\n", res, "\n")
    }
}

dbImportLine.AFFYHUEX_DB.ProbeSet <- function(conn, dataline)
{
    probeset_ID <- dataline[["probeset_ID"]] # [[ ]] to get rid of the name
    .CSVimport.line_ID(probeset_ID)
    #if (!is(conn, "DBIConnection"))
    #    .dbGetQuery(conn, paste("-- probeset_ID ", probeset_ID, sep=""))
    dataline[dataline == "---"] <- NA

    ## Extract the simple fields

    probeset_row <- dataline[names(probeset_desc$col2type)]
    dbInsertRow(conn, "probeset", probeset_row, probeset_desc$col2type)

    ## Extract and insert the "mrna_assignment" data

    field <- "mrna_assignment"
    .CSVimport.field(field)
    mrna_assignment <- mvalToMat(dataline[field], PBSsubfields.mrna_assignment)
    accessions <- mrna_assignment[ , "accession"]
    if (any(duplicated(accessions))) {
        msg <- paste(accessions, collapse=",")
        msg <- paste("duplicated accessions in ", msg, sep="")
        data_warning(msg)
    }
    dbInsertRows.PBS2mrna(conn, mrna_assignment, probeset_ID)
}

dbImportData.AFFYHUEX_DB.ProbeSet <- function(conn, csv_file, seqname, nrows=-1)
{
    if (!is.null(csv_file)) {
        .CSVimport.infile(csv_file)
        csv_con <- file(csv_file, open="r")
        on.exit(close(csv_con))
        dataline_nb <- 0
        while (nrows == -1 || dataline_nb < nrows) {
            if (dataline_nb == 0) {
                data <- read.table(csv_con, header=TRUE, sep=",", quote="\"",
                                   nrows=1, stringsAsFactors=FALSE)
                names(data) <- header <- getAffyHuExArrayAnnCSVHeader(data)
            } else {
                data <- read.table(csv_con, header=FALSE, sep=",", quote="\"",
                                   col.names=header, nrows=1, stringsAsFactors=FALSE)
            }
            if (nrow(data) == 0)
                break
            dataline_nb <- dataline_nb + 1
            .CSVimport.dataline_nb(dataline_nb)
            dataline <- unlist(data[1, ])
            dbImportLine.AFFYHUEX_DB.ProbeSet(conn, dataline)
        }
    } else {
        infile <- paste("pbs_", seqname, ".rda", sep="")
        .CSVimport.infile(infile)
        tmp_envir <- new.env(parent=emptyenv())
        load(infile, envir=tmp_envir)
        data <- get("pbs_table", envir=tmp_envir)
        for (dataline_nb in seq_len(nrow(data))) {
            .CSVimport.dataline_nb(dataline_nb)
            dataline <- unlist(data[dataline_nb, ])
            dbImportLine.AFFYHUEX_DB.ProbeSet(conn, dataline)
        }
    }
}



### =========================================================================
### G. Importation of the 2 CSV files (Transcript + Probe Set).
### -------------------------------------------------------------------------

buildTrDicts <- function(tr_file, safe=TRUE, nrows=-1)
{
    ## File "HuEx-1_0-st-v2.na21.hg18.transcript.csv" has 312368 lines and 17
    ## fields. Loading takes about 44 seconds on gladstone and 1 min on gopher6.
    cat("Loading the Transcript table from \"", tr_file, "\"... ", sep="")
    tr_table <- read.table(tr_file, header=TRUE, sep=",", quote="\"",
                           nrows=nrows, stringsAsFactors=FALSE)
    cat("OK (", nrow(tr_table), " lines loaded)\n", sep="")
    names(tr_table) <- getAffyHuExArrayAnnCSVHeader(tr_table)

    ## Extract the transcript_cluster table
    tr_cols <- names(transcript_cluster_desc$col2type)
    transcript_cluster_table <- tr_table[ , tr_cols]
    cat("Saving transcript_cluster table to transcript_cluster.rda... ")
    save(transcript_cluster_table, file="transcript_cluster.rda")
    cat("OK\n")

    mrna_assignment_list <- mvalsToMats(tr_table$mrna_assignment, TRsubfields.mrna_assignment)
    gene_assignment_list <- mvalsToMats(tr_table$gene_assignment, TRsubfields.gene_assignment)
    gene_cols <- names(gene_desc$col2type)

    ## Extract all accessions from the "mrna_assignment" col
    cat("Extract all accessions from the \"mrna_assignment\" col... ", sep="")
    allaccs <- unique(unlist(lapply(mrna_assignment_list, function(x) x[, 1])))
    cat("OK (", length(allaccs), " accessions found)\n", sep="")
    acc2id_dict <- new.env(hash=TRUE, parent=emptyenv(), size=length(allaccs))
    for (i in seq_len(length(allaccs))) acc2id_dict[[allaccs[i]]] <- i
    cat("Saving accession index to acc2id_dict.rda... ")
    save(acc2id_dict, file="acc2id_dict.rda")
    cat("OK\n")

    tr_ID_col <- as.character(tr_table$transcript_cluster_ID)

    TR2mrna_dict <- new.env(hash=TRUE, parent=emptyenv(), size=1000000L)
    TR2mrna_details_dict <- new.env(hash=TRUE, parent=emptyenv(), size=1000000L)
    gene_dict <- new.env(hash=TRUE, parent=emptyenv(), size=30000L)
    acc2genes_dict <- new.env(hash=TRUE, parent=emptyenv(), size=300000L)
    for (n in seq_len(nrow(tr_table))) {
        tr_ID <- tr_ID_col[n]
        cat(n, "/", length(tr_ID_col), ": tr_ID=", tr_ID, "\n", sep="")

        mrna_assignment <- mrna_assignment_list[[n]]
        df <- data.frame(mrna_assignment[ , -1, drop=FALSE],
                         check.names=FALSE, stringsAsFactors=FALSE)
        acc2details <- split(df, mrna_assignment[ , 1])
        accessions <- names(acc2details)
        ## Feed TR2mrna_dict and TR2mrna_details_dict
        for (acc in accessions) {
            ## With 1L instead of 1 then the arg to as.character() is integer
            ## so for 200000 it returns "200000", not "2e+05"!
            TR2mrna_id <- as.character(length(TR2mrna_dict) + 1L)
            TR2mrna_dict[[TR2mrna_id]] <- c(tr_ID, acc2id_dict[[acc]])
            TR2mrna_details_dict[[TR2mrna_id]] <- acc2details[[acc]]
        }

        gene_assignment <- gene_assignment_list[[n]]
        ## Feed gene_dict
        for (i in seq_len(nrow(gene_assignment))) {
            gene <- gene_assignment[i, gene_cols[-1]]
            gene_id <- gene_assignment[i, gene_cols[1]]
            gene0 <- gene_dict[[gene_id]]
            if (is.null(gene0)) {
                gene_dict[[gene_id]] <- gene
            } else {
                if (safe && !identical(gene, gene0))
                    stop("gene ", gene_id, " has unexpected sub-fields")
            }
        }
        ## Feed acc2genes_dict
        acc2genes <- split(gene_assignment[ , "entrez_gene_id"], gene_assignment[ , "accession"])
        for (acc in accessions) {
            genes <- acc2genes[[acc]]
            if (is.null(genes))
                genes <- character(0)
            genes0 <- acc2genes_dict[[acc]]
            if (is.null(genes0)) {
                acc2genes_dict[[acc]] <- genes
            } else {
                if (safe && !identical(genes, genes0))
                    stop("accession ", acc, " mapped to unexpected genes")
            }
        }
    }
    cat("Saving accession index to TR2mrna_dict.rda... ")
    save(TR2mrna_dict, file="TR2mrna_dict.rda")
    cat("OK\n")
    cat("Saving accession index to TR2mrna_details_dict.rda... ")
    save(TR2mrna_details_dict, file="TR2mrna_details_dict.rda")
    cat("OK\n")
    cat("Saving accession index to gene_dict.rda... ")
    save(gene_dict, file="gene_dict.rda")
    cat("OK\n")
    cat("Saving accession index to acc2genes_dict.rda... ")
    save(acc2genes_dict, file="acc2genes_dict.rda")
    cat("OK\n")
}

### 'dict' must be a dictionary where values are data frames or named atomic
### vectors. 'class' can be "data.frame" or one of the atomic types
### ("character", etc...).
dictToDataFrame <- function(dict, class, cols, col0)
{
    good_vals <- eapply(dict,
                        function(val)
                            class(val) == class && identical(names(val), cols),
                        all.names=TRUE
                 )
    if (!all(good_vals))
        stop("dictionary 'dict' contains values with bad class or names")
    data <- list()
    for (j in seq_len(length(cols))) {
        cat("Extracting col \"", cols[j], "\"... ", sep="")
        col_as_list <- eapply(dict, function(val) val[[j]], all.names=TRUE)
        ## Not sure I can trust eapply() to _always_ produce a list with the
        ## names in the same order when applied consecutively on the same
        ## environment! Hence the paranoid check...
        if (!identical(names(col_as_list), names(good_vals)))
            stop("eapply() didn't return a list with the expected names")
        if (class == "data.frame")
            data[[j]] <- unlist(col_as_list)
        else
            data[[j]] <- as(col_as_list, class)
        cat("OK\n")
    }
    names(data) <- cols
    cat("Extracting col \"", col0, "\"... ", sep="")
    if (class == "data.frame") {
        key2nrows <- eapply(dict, nrow, all.names=TRUE)
        if (!identical(names(key2nrows), names(good_vals)))
            stop("eapply() didn't return a list with the expected names")
        data[[col0]] <- rep(names(key2nrows), unlist(key2nrows))
    } else {
        data[[col0]] <- names(good_vals)
    }
    cat("OK\n")
    data.frame(data, check.names=FALSE, stringsAsFactors=FALSE)
}

importTrDicts <- function(db_file)
{
    if (is(db_file, "DBIConnection")) {
        is_new_db <- TRUE
        conn <- db_file
    } else {
        is_new_db <- !file.exists(db_file)
        conn <- dbConnect(dbDriver("SQLite"), dbname=db_file)
        on.exit(dbDisconnect(conn))
    }
    if (is_new_db) {
        tablenames <- c(
            "transcript_cluster",
            "gene",
            "mrna",
            "mrna2gene",
            "TR2mrna",
            "TR2mrna_details"
        )
        dbCreateTables.AFFYHUEX_DB(conn, tablenames)
    }

    cat("Loading the dictionaries... ")
    load("transcript_cluster.rda")
    load("acc2id_dict.rda")
    load("TR2mrna_dict.rda")
    load("TR2mrna_details_dict.rda")
    load("gene_dict.rda")
    load("acc2genes_dict.rda")
    cat("OK\n")

    cat("START IMPORTING THE DICTIONARIES...\n")

    tablename <- "transcript_cluster"
    col2type <- AFFYHUEX_DB_schema[[tablename]]$col2type
    data <- data.frame(transcript_cluster_table, check.names=FALSE, stringsAsFactors=FALSE)
    dbInsertDataFrame(conn, tablename, data, col2type, TRUE)

    tablename <- "mrna"
    col2type <- AFFYHUEX_DB_schema[[tablename]]$col2type
    if (!all(eapply(acc2id_dict, length, all.names=TRUE) == 1))
        stop("some accessions in acc2id_dict are mapped to 0 or > 1 ids")
    acc2id_list <- as.list(acc2id_dict, all.names=TRUE)
    acc2id <- sort(unlist(acc2id_list))
    data <- list()
    data[[1]] <- acc2id
    names(data[[1]]) <- NULL
    data[[2]] <- names(acc2id)
    names(data) <- names(col2type)
    data <- data.frame(data, check.names=FALSE, stringsAsFactors=FALSE)
    dbInsertDataFrame(conn, tablename, data, col2type, TRUE)

    tablename <- "TR2mrna"
    col2type <- AFFYHUEX_DB_schema[[tablename]]$col2type
    if (!all(eapply(TR2mrna_dict, length, all.names=TRUE) == 2))
        stop("some values in TR2mrna_dict have a length != 2")
    keys0 <- ls(TR2mrna_dict, all.names=TRUE)
    bad_keys <- setdiff(keys0, as.character(seq_len(length(keys0))))
    if (length(bad_keys) != 0)
        stop("TR2mrna_dict contains bad keys: \"", paste(bad_keys, collapse="\", \""), "\"")
    keys1 <- as.character(sort(as.integer(keys0)))
    data <- list()
    data[[1]] <- keys1
    data[[2]] <- sapply(keys1, function(key) TR2mrna_dict[[key]][[1]])
    data[[3]] <- sapply(keys1, function(key) TR2mrna_dict[[key]][[2]])
    names(data) <- names(col2type)
    data <- data.frame(data, check.names=FALSE, stringsAsFactors=FALSE)
    dbInsertDataFrame(conn, tablename, data, col2type, TRUE)

    ## Objects in TR2mrna_details_dict are data frames
    tablename <- "TR2mrna_details"
    col2type <- AFFYHUEX_DB_schema[[tablename]]$col2type
    keys2 <- ls(TR2mrna_details_dict, all.names=TRUE)
    if (!identical(keys2, keys0))
        stop("TR2mrna_details_dict and TR2mrna_dict don't have the same keys")
    cols <- names(col2type)[1:8]
    col0 <- names(col2type)[9]
    data <- dictToDataFrame(TR2mrna_details_dict, "data.frame", cols, col0)
    dbInsertDataFrame(conn, tablename, data, col2type, TRUE)

    ## Objects in gene_dict are named character vectors
    tablename <- "gene"
    col2type <- AFFYHUEX_DB_schema[[tablename]]$col2type
    cols <- names(col2type)[2:4]
    col0 <- names(col2type)[1]
    data <- dictToDataFrame(gene_dict, "character", cols, col0)
    dbInsertDataFrame(conn, tablename, data, col2type, TRUE)

    tablename <- "mrna2gene"
    col2type <- AFFYHUEX_DB_schema[[tablename]]$col2type
    ngenes <- eapply(acc2genes_dict, length, all.names=TRUE)
    entrez_gene_id <- as.list(acc2genes_dict, all.names=TRUE)
    ## It seems that eapply() and as.list() collect objects in the same order
    ## from the same environment but who knows for sure...
    if (!identical(names(ngenes), names(entrez_gene_id)))
        stop("eapply() and as.list() didn't return lists with the same names")
    data <- list()
    data[[1]] <- rep(names(ngenes), unlist(ngenes))
    names(entrez_gene_id) <- NULL
    data[[2]] <- unlist(entrez_gene_id)
    names(data) <- names(col2type)
    data <- data.frame(data, check.names=FALSE, stringsAsFactors=FALSE)
    dbInsertDataFrame(conn, tablename, data, col2type, TRUE)

    cat("IMPORT COMPLETED.\n")
}

### Typical use:
###   > tr_file <- "srcdata/HuEx-1_0-st-v2.na21.hg18.transcript.csv"
###   > pbs_file <- "srcdata/HuEx-1_0-st-v2.na21.hg18.probeset.csv"
###   > dbImport.AffyHuExArrayAnnCSV("test.sqlite", tr_file, pbs_file, verbose=TRUE)
### To import "tr_chr22.rda" and "pbs_chr22.rda":
###   > dbImport.AffyHuExArrayAnnCSV("test.sqlite", seqname="chr22", verbose=TRUE)
dbImport.AffyHuExArrayAnnCSV <- function(db_file, tr_file=NULL, pbs_file=NULL,
                                         fullmode=FALSE, seqname=NULL,
                                         tr_nrows=-1, pbs_nrows=-1, verbose=FALSE)
{
    .CSVimport.verbose(verbose)
    .CSVimport.fullmode(fullmode)
    if (is(db_file, "DBIConnection")) {
        is_new_db <- TRUE
        conn <- db_file
    } else {
        is_new_db <- !file.exists(db_file)
        conn <- dbConnect(dbDriver("SQLite"), dbname=db_file)
        on.exit(dbDisconnect(conn))
    }
    if (is_new_db)
        dbCreateTables.AFFYHUEX_DB(conn)
    else
        dbInitTableIds.AFFYHUEX_DB(conn)
    cat("START IMPORTING THE DATA...\n")
    dbImportData.AFFYHUEX_DB.Transcript(conn, tr_file, seqname, tr_nrows)
    if (!is.null(pbs_file) || !is.null(seqname))
        dbImportData.AFFYHUEX_DB.ProbeSet(conn, pbs_file, seqname, pbs_nrows)
    cat("IMPORTATION COMPLETED.\n")
}

