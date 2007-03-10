#############################################################################
#############################################################################
###
### Functions for importing the CSV files containing annotations for the
### Affymetrix Human Exon Array probe sets and transcript clusters.
###
### Affymetrix provides 2 files:
###   1. The Probe Set CSV file (HuEx-1_0-st-v2.na21.hg18.probeset.csv)
###   2. The Transcript CSV file (HuEx-1_0-st-v2.na21.hg18.transcript.csv)
###
### This file is divided in 6 sections:
###   A. General purpose low-level SQL helper functions.
###   B. CSV field description and DB schema.
###   C. Objects and functions shared by D. and E.
###   D. Importation of the Transcript CSV file (167M, aka the "small" file).
###   E. Importation of the Probe Set CSV file (475M, aka the "big" file).
###   F. Importation of the 2 CSV files (Transcript + Probe Set).
### 
### WARNING: This is a WORK IN PROGRESS!!! (if pdInfoBuilder had a NAMESPACE,
### nothing should be exported from this file for now)
###
#############################################################################



### =========================================================================
### A. General purpose low-level SQL helper functions (should probably go
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



### =========================================================================
### B. CSV field description and DB schema.
###
###    3 sub-sections:
###      B.a. CSV field description for the Transcript CSV file.
###      B.b. CSV field description for the Probe Set CSV file.
###      B.c. The DB schema.
###
###    The original field description is provided in
###      HuEx-1_0-st-v2.na21.hg18.AFFX_README.NetAffx-CSV-Files.txt
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### B.a. CSV field description for the Transcript CSV file.
###
### The character vectors below list the CSV sub-fields for each multipart
### field.
###

TRsubfields.gene_assignment <- c(
    "accession",
    "gene_symbol",
    "gene_title",
    "cytoband",
    "entrez_gene_id"
)

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
### B.b. CSV field description for the Probe Set CSV file.
###
### The character vectors below list the CSV sub-fields for each multipart
### field.
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
### B.c. The DB schema.
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
        accession="TEXT",
        entrez_gene_id="INTEGER"    # REFERENCES gene(entrez_gene_id)
    ),
    col2key=c(
        `_mrna_id`="PRIMARY KEY",
        accession="UNIQUE",
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

### Global schema (16 tables).
AFFYHUEX_DB_schema <- list(
    transcript_cluster=transcript_cluster_desc,
    gene=gene_desc,
    mrna=mrna_desc,
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
### C. Objects and functions shared by D. and E. 
### -------------------------------------------------------------------------


dbCreateTables.AFFYHUEX_DB <- function(conn)
{
    for (tablename in names(AFFYHUEX_DB_schema)) {
        col2type <- AFFYHUEX_DB_schema[[tablename]]$col2type
        col2key <- AFFYHUEX_DB_schema[[tablename]]$col2key
        dbCreateTable(conn, tablename, col2type, col2key)
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

.CSVimport.chr1.only <- function(chr1.only=NULL)
{
    .CSVimport.var("chr1.only", chr1.only, defval=FALSE)
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

multipartToMatrix <- function(multipart_val, subfields, min.nsubfields=length(subfields))
{
    ncol <- length(subfields)
    if (is.na(multipart_val)) {
        mat <- matrix(data=character(0), ncol=ncol)
    } else {
        vals <- strsplit(multipart_val, " /// ", fixed=TRUE)[[1]]
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
splitMatrix <- function(mat, acc2id)
{
    if (ncol(mat) < 2)
        stop("won't split a matrix with less than 2 cols")
    if (colnames(mat)[1] != "accession")
        stop("first 'mat' col name must be \"accession\"")
    accessions <- mat[ , 1]
    if (!all(accessions %in% names(acc2id)))
        data_error("'mat' contains invalid accessions")
    uids <- unique(acc2id)
    id2submat <- list()
    length(id2submat) <- length(uids)
    names(id2submat) <- uids
    for (accession in accessions) {
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
### 'mat' and 'dat' _must_ have exactly the same col names (it's an error if
### they don't).
haveTheSameData <- function(mat, dat)
{
    if (is.null(colnames(mat)))
        stop("'mat' has no col names")
    if (!identical(colnames(mat), colnames(dat)))
        stop("'mat' and 'dat' have different col names")
    if (nrow(mat) != nrow(dat))
        return(FALSE)
    ## Convert 'dat' to a character matrix (we can't just use as.matrix
    ## because of the infamous "format" feature).
    mat2 <- do.call("cbind", args=lapply(dat, function(x) if (is.character(x)) x else as.character(x)))
    ## Find the order of the rows in the 2 matrices.
    ii1 <- do.call("order", args=lapply(colnames(mat), function(col) mat[ , col]))
    ii2 <- do.call("order", args=lapply(colnames(mat2), function(col) mat2[ , col]))
    ## Compare row by row.
    for (i in seq_len(nrow(mat)))
        if (!identical(mat[ii1[i], ], mat2[ii2[i], ]))
            return(FALSE)
    return(TRUE)
}

dbInsert_multipart_data <- function(conn, tablename, mat, insres)
{
    if (colnames(mat)[1] != "accession")
        stop("first 'mat' col name must be \"accession\"")
    acc2id <- insres$acc2id
    id2submat <- splitMatrix(mat, acc2id)
    new_ids <- insres$new_ids
    col2type <- AFFYHUEX_DB_schema[[tablename]]$col2type
    cols0 <- paste(colnames(mat)[-1], collapse=",")
    link0 <- names(col2type)[length(col2type)]
    sql0 <- paste("SELECT ", cols0, " FROM ", tablename, " WHERE ", link0, "=", sep="")
    for (id in names(id2submat)) {
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
### D. Importation of the Transcript CSV file (167M, aka the "small" file).
### -------------------------------------------------------------------------


dbInsertRows.gene <- function(conn, genes)
{
    acc2id <- character(nrow(genes))
    names(acc2id) <- genes[ , "accession"]
    new_ids <- character(0)
    col2type <- gene_desc$col2type
    for (i in seq_len(nrow(genes))) {
        row1 <- genes[i, names(col2type)]
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
        acc2id[i] <- id
    }
    list(acc2id=acc2id, new_ids=new_ids)
}

dbInsertRows.mrna <- function(conn, accessions, gene_acc2id)
{
    acc2id <- character(length(accessions))
    names(acc2id) <- accessions
    new_ids <- character(0)
    col2type <- mrna_desc$col2type
    for (i in seq_len(length(accessions))) {
        accession <- accessions[i]
        entrez_gene_id <- gene_acc2id[accession]
        row1 <- c(NA, accession, entrez_gene_id)
        names(row1) <- names(col2type)
        row0 <- try(dbGetThisRow(conn, "mrna", "accession", row1, col2type), silent=TRUE)
        if (is(row0, "try-error"))
            stop("In ", csv_current_pos(), ":\n", row0, "\n")
        if (is.null(row0)) {
            id <- .db.next.id("mrna")
            row1["_mrna_id"] <- id
            res <- try(dbInsertRow(conn, "mrna", row1, col2type), silent=TRUE)
            if (is(res, "try-error"))
                stop("In ", csv_current_pos(), ":\n", res, "\n")
            new_ids <- c(new_ids, id)
        } else {
            id <- row0["_mrna_id"]
        }
        acc2id[i] <- id
    }
    list(acc2id=acc2id, new_ids=new_ids)
}

dbInsertRows.TR2mrna <- function(conn, transcript_cluster_ID, mrna_acc2id)
{
    acc2id <- mrna_acc2id
    acc2id[] <- ""
    col2type <- TR2mrna_desc$col2type
    for (i in seq_len(length(mrna_acc2id))) {
        id <- .db.next.id("TR2mrna")
        row1 <- c(id, transcript_cluster_ID, mrna_acc2id[i])
        names(row1) <- names(col2type)
        res <- try(dbInsertRow(conn, "TR2mrna", row1, col2type), silent=TRUE)
        if (is(res, "try-error"))
            stop("In ", csv_current_pos(), ":\n", res, "\n")
        acc2id[i] <- id
    }
    acc2id
}

dbInsertRows.TR2mrna_details <- function(conn, mrna_assignment, TR2mrna_acc2id)
{
    col2type <- TR2mrna_details_desc$col2type
    for (i in seq_len(nrow(mrna_assignment))) {
        accession <- mrna_assignment[i, "accession"]
        TR2mrna.id <- TR2mrna_acc2id[accession]
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

    ## Extract and insert the "gene_assignment" data

    field <- "gene_assignment"
    .CSVimport.field(field)
    gene_assignment <- multipartToMatrix(dataline[field], TRsubfields.gene_assignment)
    gene_insres <- dbInsertRows.gene(conn, gene_assignment)

    ## There should never be more than 1 part with the same accession in the
    ## "gene_assignment" field. Unfortunately this happens sometimes (very rarely though).
    ## For example in HuEx-1_0-st-v2.na21.hg18.transcript.csv, line 16411
    ## (transcript_cluster_ID=2718075), the "gene_assignment" field (multipart)
    ## contains the following parts:
    ##
    ##      accession | gene_symbol | gene_title                    | cytoband | entrez_gene_id
    ##   -------------|-------------|-------------------------------|----------|---------------
    ##   NM_001040448 | DEFB131     | defensin, beta 131            |   4p16.1 |         644414
    ##      XM_938410 | DEFB131     | defensin, beta 131            |   4p16.1 |         644414
    ##      XM_938410 | LOC649335   | similar to Beta-defensin 131..|          |         649335
    ##   ...
    ## This poses 2 problems: (1) the current DB schema can't handle this (a
    ## given mrna can only be linked to 1 or 0 gene), (2) when this happens,
    ## then the GO annotations are ambiguous (this is because Affymetrix has
    ## choosen to link GO terms to mrnas, not to genes).
    ## For example, on the same line (transcript_cluster_ID=2718075), the
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
    ##
    ## So we can't do this anymore:
    ##   if (any(duplicated(names(gene_insres$acc2id))))
    ##      stop("in CSV line for transcript_cluster_ID=", transcript_cluster_ID, ": ",
    ##           "\"gene_assignment\" has more than 1 part with the same accession")
    ## For now we ignore the duplicated and issue a warning
    dup_gene_acc2id <- duplicated(names(gene_insres$acc2id))
    if (any(dup_gene_acc2id)) {
        msg <- paste(names(gene_insres$acc2id), collapse=",")
        msg <- paste(msg, " [", paste(gene_insres$acc2id, collapse=","), "]", sep="")
        msg <- paste("\"gene_assignment\" has more than 1 part with the same accession: ",
                     msg, sep="")
        data_warning(msg)
        gene_insres$acc2id <- gene_insres$acc2id[!dup_gene_acc2id]
    }

    ## Extract and insert the "mrna_assignment" data

    field <- "mrna_assignment"
    .CSVimport.field(field)
    mrna_assignment <- multipartToMatrix(dataline[field], TRsubfields.mrna_assignment)
    accessions <- unique(mrna_assignment[ , "accession"])
    if (!all(names(gene_insres$acc2id) %in% accessions)) {
        msg <- "\"gene_assignment\" has unlinked parts"
        data_error(msg)
    }
    mrna_insres <- dbInsertRows.mrna(conn, accessions, gene_insres$acc2id)
    TR2mrna_acc2id <- dbInsertRows.TR2mrna(conn, transcript_cluster_ID, mrna_insres$acc2id)
    dbInsertRows.TR2mrna_details(conn, mrna_assignment, TR2mrna_acc2id)

    ## Extract and insert the "swissprot" data

    field <- "swissprot"
    .CSVimport.field(field)
    mat <- multipartToMatrix(dataline[field], TRsubfields.swissprot)
    dbInsert_multipart_data(conn, field, mat, mrna_insres)

    ## Extract and insert the "unigene" data

    field <- "unigene"
    .CSVimport.field(field)
    mat <- multipartToMatrix(dataline[field], TRsubfields.unigene, min.nsubfields=2)
    dbInsert_multipart_data(conn, field, mat, mrna_insres)

    ## Extract and insert the "GO" data

    field <- "GO_biological_process"
    .CSVimport.field(field)
    mat <- multipartToMatrix(dataline[field], TRsubfields.GO_biological_process)
    mat <- replaceGOEvidenceByCode(mat)
    dbInsert_multipart_data(conn, field, mat, gene_insres)

    field <- "GO_cellular_component"
    .CSVimport.field(field)
    mat <- multipartToMatrix(dataline[field], TRsubfields.GO_biological_process)
    mat <- replaceGOEvidenceByCode(mat)
    dbInsert_multipart_data(conn, field, mat, gene_insres)

    field <- "GO_molecular_function"
    .CSVimport.field(field)
    mat <- multipartToMatrix(dataline[field], TRsubfields.GO_biological_process)
    mat <- replaceGOEvidenceByCode(mat)
    dbInsert_multipart_data(conn, field, mat, gene_insres)

    ## Extract and insert the "pathway" data

    field <- "pathway"
    .CSVimport.field(field)
    mat <- multipartToMatrix(dataline[field], TRsubfields.pathway)
    dbInsert_multipart_data(conn, field, mat, mrna_insres)

    return() # that's all for now

    ## The code below is not ready...
    protein_domains <- multipartToMatrix(dataline["protein_domains"],
                                         TRsubfields.protein_domains, min.nsubfields=3)
    dbInsert_multipart_data(conn, "protein_domains", protein_domains,
                                   acc2id, new_accessions, transcript_cluster_ID)

    protein_families <- multipartToMatrix(dataline["protein_families"], TRsubfields.protein_families)
    dbInsert_multipart_data(conn, "protein_families", protein_families,
                                   acc2id, new_accessions, transcript_cluster_ID)
}

### File "HuEx-1_0-st-v2.na21.hg18.transcript.csv" has 312368 lines
### and 17 fields. To load the entire file at once:
###   > data <- read.table(csv_file, header=TRUE, sep=",", quote="\"", stringsAsFactors=FALSE)
### It takes about 1 min on gopher6.
###
dbImportData.AFFYHUEX_DB.Transcript <- function(conn, csv_file, nrows=-1)
{
    .CSVimport.infile(csv_file)
    csv_con <- file(csv_file, open="r")
    on.exit(close(csv_con))
    dataline_nb <- 0
    while (nrows == -1 || dataline_nb < nrows) {
        if (dataline_nb == 0) {
            data <- read.table(csv_con, header=TRUE, sep=",", quote="\"",
                               nrows=1, stringsAsFactors=FALSE)
            header <- names(data)
            header[header == "transcript_cluster_id"] <- "transcript_cluster_ID"
            names(data) <- header
        } else {
            data <- read.table(csv_con, header=FALSE, sep=",", quote="\"",
                               col.names=header, nrows=1, stringsAsFactors=FALSE)
        }
        if (nrow(data) == 0)
            break
        dataline_nb <- dataline_nb + 1
        .CSVimport.dataline_nb(dataline_nb)
        dataline <- unlist(data[1, ])
        if (.CSVimport.chr1.only() && dataline["seqname"] != "chr1")
            break
        dbImportLine.AFFYHUEX_DB.Transcript(conn, dataline)
    }
}



### =========================================================================
### E. Importation of the Probe Set CSV file (475M, aka the "big" file).
### -------------------------------------------------------------------------


dbInsertRows.PBS2mrna <- function(conn, mrna_assignment, probeset_ID)
{
    col2type <- PBS2mrna_desc$col2type
    for (i in seq_len(nrow(mrna_assignment))) {
        accession <- mrna_assignment[i, "accession"]
        row1 <- c(NA, accession, NA)
        names(row1) <- names(mrna_desc$col2type)
        row0 <- try(dbGetThisRow(conn, "mrna", "accession", row1, col2type), silent=TRUE)
        if (is(row0, "try-error"))
            stop("In ", csv_current_pos(), ":\n", row0, "\n")
        if (is.null(row0)) {
            msg <- paste("table \"", tablename, "\" has no record ",
                         "with accession ", accession, sep="")
            data_error(msg)
        }
        mrna.id <- row0["_mrna_id"]
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
    mrna_assignment <- multipartToMatrix(dataline[field], PBSsubfields.mrna_assignment)
    accessions <- mrna_assignment[ , "accession"]
    if (any(duplicated(accessions))) {
        msg <- paste(accessions, collapse=",")
        msg <- paste("duplicated accessions in ", msg, sep="")
        data_warning(msg)
    }
    dbInsertRows.PBS2mrna(conn, mrna_assignment, probeset_ID)
}

### File "HuEx-1_0-st-v2.na21.hg18.probeset.csv" has 1425647 lines
### and 39 fields. Trying to load the entire file at once with:
###   > data <- read.table(csv_file, header=TRUE, sep=",", quote="\"", stringsAsFactors=FALSE)
### takes 20 minutes on gladstone! (32G of RAM)
###
dbImportData.AFFYHUEX_DB.ProbeSet <- function(conn, csv_file, nrows=-1)
{
    .CSVimport.infile(csv_file)
    csv_con <- file(csv_file, open="r")
    on.exit(close(csv_con))
    dataline_nb <- 0
    while (nrows == -1 || dataline_nb < nrows) {
        if (dataline_nb == 0) {
            data <- read.table(csv_con, header=TRUE, sep=",", quote="\"",
                               nrows=1, stringsAsFactors=FALSE)
            header <- names(data)
            header[header == "probeset_id"] <- "probeset_ID"
            header[header == "transcript_cluster_id"] <- "transcript_cluster_ID"
            names(data) <- header
        } else {
            data <- read.table(csv_con, header=FALSE, sep=",", quote="\"",
                               col.names=header, nrows=1, stringsAsFactors=FALSE)
        }
        if (nrow(data) == 0)
            break
        dataline_nb <- dataline_nb + 1
        .CSVimport.dataline_nb(dataline_nb)
        dataline <- unlist(data[1, ])
        if (.CSVimport.chr1.only() && dataline["seqname"] != "chr1")
            break
        dbImportLine.AFFYHUEX_DB.ProbeSet(conn, dataline)
    }
}



### =========================================================================
### F. Importation of the 2 CSV files (Transcript + Probe Set).
### -------------------------------------------------------------------------


### Typical use:
###   > transcript_file <- "srcdata/HuEx-1_0-st-v2.na21.hg18.transcript.csv"
###   > probeset_file <- "srcdata/HuEx-1_0-st-v2.na21.hg18.probeset.csv"
###   > dbImport.AffyHuExArrayAnnCSV(transcript_file, probeset_file,
###                                  "test.sqlite", chr1.only=TRUE, verbose=TRUE)
### To skip importation of the "Probe Set" file:
###   > dbImport.AffyHuExArrayAnnCSV(transcript_file, "", "test.sqlite", , 20, verbose=TRUE)
###
dbImport.AffyHuExArrayAnnCSV <- function(transcript_file, probeset_file, db_file,
                                         chr1.only=FALSE,
                                         transcript_nrows=-1, probeset_nrows=-1,
                                         verbose=FALSE)
{
    .CSVimport.verbose(verbose)
    .CSVimport.chr1.only(chr1.only)
    conn <- dbConnect(dbDriver("SQLite"), dbname=db_file)
    on.exit(dbDisconnect(conn))
    dbCreateTables.AFFYHUEX_DB(conn)
    dbImportData.AFFYHUEX_DB.Transcript(conn, transcript_file, transcript_nrows)
    if (is.null(probeset_file) || is.na(probeset_file) || probeset_file == "")
        return()
    dbImportData.AFFYHUEX_DB.ProbeSet(conn, probeset_file, probeset_nrows)
}

