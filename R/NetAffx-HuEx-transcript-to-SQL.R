### =========================================================================
### Utilities for converting a NetAffx HuEx transcript CSV file to SQL
### 
### ... work in progress ... (nothing from this file is exported)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### General purpose low-level SQL helper functions (should probably go
### somewhere else).
###

### 'conn' must be a DBIConnection object, a filename or a file object
.dbGetQuery <- function(conn, sql)
{
    if (is(conn, "DBIConnection"))
        dbGetQuery(conn, sql)
    else
        cat(sql, ";\n", file=conn, sep="", append=TRUE)
}

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
            x <- as.integer(val)
            if (is.na(x) || x != val)
                stop("vals[", i, "]=\"", val, "\" not an integer")
            vals[i] <- x
            next
        }
        vals[i] <- paste("'", gsub("'", "''", val, fixed=TRUE), "'", sep="")
    }
    vals
}

dbCreateTable <- function(conn, tablename, col2type, col2key)
{
    col2type[names(col2key)] <- paste(col2type[names(col2key)], col2key, sep=" ")
    sql <- paste(names(col2type), col2type, sep=" ")
    sql <- paste(sql, collapse=", ")
    sql <- paste("CREATE TABLE ", tablename, " (", sql, ")", sep="")
    .dbGetQuery(conn, sql)
}

### 'row' must be a character vector with or without names.
dbInsertRow <- function(conn, tablename, row, col2type, is.fullrow=TRUE)
{
    sqlvals <- toSQLValues(row, col2type)
    sql <- paste(sqlvals, collapse=", ")
    sql <- paste("VALUES (", sql, ")", sep="")
    cols <- names(row)
    if (!is.null(cols) && !identical(cols, names(col2type))) {
        if (is.fullrow)
            stop("'names(row)' and 'names(col2type)' are not identical")
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
    unique_sqlval <- toSQLValues(row[unique_col], col2type)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ID generation.
###

.db.current.ids <- new.env(hash=TRUE)

set.id <- function(tablename, id=0)
{
    assign(tablename, id, envir=.db.current.ids, inherits=FALSE)
}

next.id <- function(tablename)
{
    if (!exists(tablename, envir=.db.current.ids, inherits=FALSE))
        stop("no current id for table \"", tablename, "\"")
    id <- get(tablename, envir=.db.current.ids, inherits=FALSE)
    id <- id + 1
    set.id(tablename, id)
    id
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions specific to NetAffx CSV files.
###

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
            if (nsubfields < min.nsubfields || nsubfields > ncol)
                stop("bad number of subfields")
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

### Return a list of matrices.
splitMatrix <- function(mat, acc2id)
{
    if (colnames(mat)[1] != "accession")
        stop("first 'mat' col name must be \"accession\"")
    accessions <- mat[ , 1]
    if (!all(accessions %in% names(acc2id)))
        stop("'mat' contains invalid accessions")
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
            stop("can't split 'mat'")
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### CSV sub-fields of the multipart fields.
###

gene_assignment_subfields <- c(
    "accession",
    "gene_symbol",
    "gene_title",
    "cytoband",
    "entrez_gene_id"
)

mrna_assignment_subfields <- c(
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

swissprot_subfields <- c(
    "accession",
    "swissprot_accession"
)

unigene_subfields <- c(
    "accession",
    "unigene_id",
    "unigene_expr"
)

GO_biological_process_subfields <- c(
    "accession",
    "GO_id",
    "GO_term",
    "GO_evidence"
)

pathway_subfields <- c(
    "accession",
    "source",
    "pathway_name"
)

protein_domains_subfields <- c(
    "accession",
    "source",
    "accession_or_domain_name",
    "domain_description"
)

protein_families_subfields <- c(
    "accession",
    "source",
    "family_accession",
    "family_description"
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### An R representation of the NETAFFX_HUEX_TRANSCRIPT_DB schema.
###

### The "probeset" table.
###
### TODO: Rename this table "featureSet". (This naming style maybe popular in
### the R culture (especially for functions/methods) but not really in the SQL
### culture (especially for table names)).
### Also rename cols: "probeset_id" -> "fsetid" and "seqname" -> "chrom".
###
probeset_desc <- list(
    col2type=c(
        probeset_id="INTEGER",      # PRIMARY KEY
        seqname="TEXT",
        strand="CHAR(1)",
        start="INTEGER",
        stop="INTEGER",
        total_probes="INTEGER"
    ),
    col2key=c(
        probeset_id="PRIMARY KEY"
    )
)

### The "gene" table.
###
### Note that we can't put a UNIQUE constraint on "gene_symbol".
### For example in HuEx-1_0-st-v2.na21.hg18.transcript.csv, line 10430
### (probeset_id=2564225), the "gene_assignment" field (multipart) contains
### the following parts:
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

### The "mrna_assignment" table.
### TODO: Add a UNIQUE constraint on (probeset_id, _mrna_id).
mrna_assignment_desc <- list(
    col2type=c(
        `_mrna_assignment_id`="INTEGER",
        probeset_id="INTEGER",      # REFERENCES probeset(probeset_id)
        `_mrna_id`="INTEGER"        # REFERENCES mrna(_mrna_id)
    ),
    col2key=c(
        `_mrna_assignment_id`="PRIMARY KEY",
        probeset_id="NOT NULL REFERENCES probeset(probeset_id)",
        `_mrna_id`="NOT NULL REFERENCES mrna(_mrna_id)"
    )
)

### The "mrna_assignment_details" table.
mrna_assignment_details_desc <- list(
    col2type=c(
        source_name="TEXT",
        description="TEXT",
        assignment_seqname="TEXT",
        assignment_score="INTEGER",
        assignment_coverage="INTEGER",
        direct_probes="INTEGER",
        possible_probes="INTEGER",
        xhyb="INTEGER",
        `_mrna_assignment_id`="INTEGER"        # REFERENCES mrna_assignment(_mrna_assignment_id)
    ),
    col2key=c(
        `_mrna_assignment_id`="NOT NULL REFERENCES mrna_assignment(_mrna_assignment_id)"
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
### TODO: Add a UNIQUE constraint on (GO_id, GO_evidence).
GO_biological_process_desc <- list(
    col2type=c(
        GO_id="TEXT",
        GO_term="TEXT",
        GO_evidence="TEXT",
        entrez_gene_id="INTEGER"        # REFERENCES gene(entrez_gene_id)
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

### The "mrna_assignment_2_protein_domains" table.
### TODO: Add a UNIQUE constraint on (_mrna_assignment_id, _protein_domains_id).
mrna_assignment_2_protein_domains_desc <- list(
    col2type=c(
        `_mrna_assignment_id`="INTEGER",
        `_protein_domains_id`="INTEGER"
    ),
    col2key=c(
        `_mrna_assignment_id`="NOT NULL REFERENCES mrna_assignment(_mrna_assignment_id)",
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

### Global schema (14 tables).
NETAFFX_HUEX_TRANSCRIPT_DB_schema <- list(
    probeset=probeset_desc,
    gene=gene_desc,
    mrna=mrna_desc,
    mrna_assignment=mrna_assignment_desc,
    mrna_assignment_details=mrna_assignment_details_desc,
    swissprot=swissprot_desc,
    unigene=unigene_desc,
    GO_biological_process=GO_biological_process_desc,
    GO_cellular_component=GO_biological_process_desc,
    GO_molecular_function=GO_biological_process_desc,
    pathway=pathway_desc,
    protein_domains=protein_domains_desc,
    mrna_assignment_2_protein_domains=mrna_assignment_2_protein_domains_desc,
    protein_families=protein_families_desc
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Create the "transcript" tables and insert the "transcript" data into
### them.
###

create_NetAffx_HuEx_transcript_tables <- function(conn)
{
    for (tablename in names(NETAFFX_HUEX_TRANSCRIPT_DB_schema)) {
        col2type <- NETAFFX_HUEX_TRANSCRIPT_DB_schema[[tablename]]$col2type
        col2key <- NETAFFX_HUEX_TRANSCRIPT_DB_schema[[tablename]]$col2key
        dbCreateTable(conn, tablename, col2type, col2key)
        set.id(tablename)
    }
}

insert_gene_data <- function(conn, genes)
{
    acc2id <- character(nrow(genes))
    names(acc2id) <- genes[ , "accession"]
    new_ids <- character(0)
    col2type <- gene_desc$col2type
    for (i in seq_len(nrow(genes))) {
        row1 <- genes[i, names(col2type)]
        id <- row1[["entrez_gene_id"]] # [[ ]] to get rid of the name
        row0 <- dbGetThisRow(conn, "gene", "entrez_gene_id", row1, col2type)
        if (is.null(row0)) {
            dbInsertRow(conn, "gene", row1, col2type)
            new_ids <- c(new_ids, id)
        }
        acc2id[i] <- id
    }
    list(acc2id=acc2id, new_ids=new_ids)
}

insert_mrna_data <- function(conn, accessions, gene_acc2id)
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
        row0 <- dbGetThisRow(conn, "mrna", "accession", row1, col2type)
        if (is.null(row0)) {
            id <- next.id("mrna")
            row1["_mrna_id"] <- id
            dbInsertRow(conn, "mrna", row1, col2type)
            new_ids <- c(new_ids, id)
        } else {
            id <- row0["_mrna_id"]
        }
        acc2id[i] <- id
    }
    list(acc2id=acc2id, new_ids=new_ids)
}

insert_mrna_assignment_data <- function(conn, probeset_id, mrna_acc2id)
{
    acc2id <- mrna_acc2id
    acc2id[] <- ""
    col2type <- mrna_assignment_desc$col2type
    for (i in seq_len(length(mrna_acc2id))) {
        id <- next.id("mrna_assignment")
        row1 <- c(id, probeset_id, mrna_acc2id[i])
        names(row1) <- names(col2type)
        dbInsertRow(conn, "mrna_assignment", row1, col2type)
        acc2id[i] <- id
    }
    acc2id
}

insert_mrna_assignment_details_data <- function(conn, mrna_assignment, mrna_assignment_acc2id)
{
    col2type <- mrna_assignment_details_desc$col2type
    for (i in seq_len(nrow(mrna_assignment))) {
        accession <- mrna_assignment[i, "accession"]
        mrna_assignment.id <- mrna_assignment_acc2id[accession]
        row1 <- mrna_assignment[i, names(col2type)[1:8]]
        row1 <- c(row1, mrna_assignment.id)
        names(row1) <- names(col2type)
        dbInsertRow(conn, "mrna_assignment_details", row1, col2type)
    }
}

insert_NetAffx_multipart_field <- function(conn, tablename, mat, insres, probeset_id, verbose=FALSE)
{
    if (colnames(mat)[1] != "accession")
        stop("first 'mat' col name must be \"accession\"")
    acc2id <- insres$acc2id
    id2submat <- splitMatrix(mat, acc2id)
    new_ids <- insres$new_ids
    col2type <- NETAFFX_HUEX_TRANSCRIPT_DB_schema[[tablename]]$col2type
    cols0 <- paste(colnames(submat), collapse=",")
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
                dbInsertRow(conn, tablename, row1, col2type)
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
            stop("in CSV line for probeset_id=", probeset_id, ": ",
                 accession, " parts in \"", tablename, "\" are not the expected ones")
        }
    }
}

insert_NetAffx_HuEx_transcript_data <- function(conn, data, verbose=FALSE)
{
    .mrna.id <- 0
    for (i in seq_len(nrow(data))) {
        #if (is(conn, "SQLiteConnection"))
        #    dbBeginTransaction(conn)
        row <- unlist(data[i, ])
        probeset_id <- row[["probeset_id"]] # [[ ]] to get rid of the name
        if (verbose)
            cat(i, ": probeset_id=", probeset_id, "\n", sep="")
        #if (!is(conn, "DBIConnection"))
        #    .dbGetQuery(conn, paste("-- probeset_id ", probeset_id, sep=""))
        row[row == "---"] <- NA

        ## Extract the simple fields
        probeset_row <- row[names(probeset_desc$col2type)]
        dbInsertRow(conn, "probeset", probeset_row, probeset_desc$col2type)

        ## Extract and insert the "gene_assignment" data
        gene_assignment <- multipartToMatrix(row["gene_assignment"], gene_assignment_subfields)
        gene_insres <- insert_gene_data(conn, gene_assignment)
        if (any(duplicated(names(gene_insres$acc2id))))
            stop("in CSV line for probeset_id=", probeset_id, ": ",
                 "\"gene_assignment\" has more than 1 part with the same accession")

        ## Extract and insert the "mrna_assignment" data
        mrna_assignment <- multipartToMatrix(row["mrna_assignment"], mrna_assignment_subfields)
        accessions <- unique(mrna_assignment[ , "accession"])
        if (!all(names(gene_insres$acc2id) %in% accessions))
            stop("in CSV line for probeset_id=", probeset_id, ": ",
                 "\"gene_assignment\" has unlinked parts")
        mrna_insres <- insert_mrna_data(conn, accessions, gene_insres$acc2id)
        mrna_assignment_acc2id <- insert_mrna_assignment_data(conn, probeset_id, mrna_insres$acc2id)
        insert_mrna_assignment_details_data(conn, mrna_assignment, mrna_assignment_acc2id)

        ## Extract and insert the "swissprot" data
        swissprot <- multipartToMatrix(row["swissprot"], swissprot_subfields)
        insert_NetAffx_multipart_field(conn, "swissprot", swissprot,
                                       mrna_insres, probeset_id, verbose)

        ## Extract and insert the "unigene" data
        unigene <- multipartToMatrix(row["unigene"], unigene_subfields, min.nsubfields=2)
        insert_NetAffx_multipart_field(conn, "unigene", unigene,
                                       mrna_insres, probeset_id, verbose)

        ## Extract and insert the "GO" data
        GO_biological_process <- multipartToMatrix(row["GO_biological_process"],
                                                   GO_biological_process_subfields)
        insert_NetAffx_multipart_field(conn, "GO_biological_process", GO_biological_process,
                                       gene_insres, probeset_id, verbose)
        GO_cellular_component <- multipartToMatrix(row["GO_cellular_component"],
                                                   GO_biological_process_subfields)
        insert_NetAffx_multipart_field(conn, "GO_cellular_component", GO_cellular_component,
                                       gene_insres, probeset_id, verbose)
        GO_molecular_function <- multipartToMatrix(row["GO_molecular_function"],
                                                   GO_biological_process_subfields)
        insert_NetAffx_multipart_field(conn, "GO_molecular_function", GO_molecular_function,
                                       gene_insres, probeset_id, verbose)

        ## Extract and insert the "pathway" data
        pathway <- multipartToMatrix(row["pathway"], pathway_subfields)
        insert_NetAffx_multipart_field(conn, "pathway", pathway,
                                       mrna_insres, probeset_id, verbose)

        next # that's all for now

        protein_domains <- multipartToMatrix(row["protein_domains"],
                                             protein_domains_subfields, min.nsubfields=3)
        insert_NetAffx_multipart_field(conn, "protein_domains", protein_domains,
                                       acc2id, new_accessions, probeset_id)

        protein_families <- multipartToMatrix(row["protein_families"], protein_families_subfields)
        insert_NetAffx_multipart_field(conn, "protein_families", protein_families,
                                       acc2id, new_accessions, probeset_id)

        #if (is(conn, "SQLiteConnection"))
        #    dbCommit(conn)
    }
}

### Typical use:
###   > library(RSQLite)
###   > transcript_csv_file <- "srcdata/HuEx-1_0-st-v2.na21.hg18.transcript.csv"
###   > dbImport_NetAffx_HuEx_transcript(transcript_csv_file, "test.sqlite", 200, TRUE)
###
dbImport_NetAffx_HuEx_transcript <- function(transcript_csv_file, db_file, nrows=-1, verbose=FALSE)
{
    ## Takes about 1 min to load file "HuEx-1_0-st-v2.na21.hg18.transcript.csv"
    ## (312368x17) into 'transcript_data' on gopher6.
    transcript_data <- read.table(transcript_csv_file, header=TRUE, sep=",", quote="\"",
                                  nrows=nrows, stringsAsFactors=FALSE)
    conn <- dbConnect(dbDriver("SQLite"), dbname=db_file)
    create_NetAffx_HuEx_transcript_tables(conn)
    insert_NetAffx_HuEx_transcript_data(conn, transcript_data, verbose)
    dbDisconnect(conn)
}

