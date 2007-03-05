### =========================================================================
### Utilities for converting a NetAffx HuEx transcript CSV file to SQL
### 
### ... work in progress ... (nothing from this file is exported)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### General purpose low-level SQL helper functions (should probably go
### somewhere else).
###

### 'conn' must be a DBIConnection object, a filename or a file object
.dbSendQuery <- function(conn, sql)
{
    if (is(conn, "DBIConnection"))
        dbSendQuery(conn, sql)
    else
        cat(sql, ";\n", file=conn, sep="", append=TRUE)
}

### 'vals' must be a character vector with NAs for the SQL NULL.
### Return a character vector of the same length as 'vals' (and with the
### same names if 'vals' has names).
toSQLValues <- function(vals, col2type)
{
    lv <- length(vals)
    if (lv == 0)
        stop("'vals' is empty")
    cols <- names(vals)
    if (is.null(cols)) {
        if (lv != length(col2type)) {
            cat("          vals  = ", vals, "\n", sep="|")
            cat("names(col2type) = ", names(col2type), "\n", sep="|")
            cat("      col2type  = ", col2type, "\n", sep="|")
            stop("when unamed, 'vals' must be of the same length as 'col2type'")
        }
    } else {
        col2type <- col2type[cols]
        if (any(is.na(col2type)))
            stop("no type found for some of the values in 'vals'")
    }
    ## From now, vals and col2type have the same length
    ## and the mapping between them is positional.
    for (i in 1:lv) {
        val <- vals[i]
        if (is.na(val)) {
            vals[i] <- "NULL"
            next
        }
        if (col2type[i] == "INTEGER") {
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
    .dbSendQuery(conn, sql)
}

### 
dbInsertRow <- function(conn, tablename, row, col2type)
{
    sqlvals <- toSQLValues(row, col2type)
    sql <- paste(sqlvals, collapse=", ")
    sql <- paste("VALUES (", sql, ")", sep="")
    cols <- names(sqlvals)
    if (!is.null(cols)) {
        cols <- paste(cols, collapse=", ")
        sql <- paste("(", cols, ") ", sql, sep="")
    }
    sql <- paste("INSERT INTO ", tablename, " ", sql, sep="")
    .dbSendQuery(conn, sql)
}

### Return NULL if 'row' is not found in 'tablename'.
### Note that it wouldn't make sense to pass a row such that 'row[unique_col]'
### is NA.
dbGetThisRow <- function(conn, tablename, unique_col, row, col2type)
{
    sqlval <- toSQLValues(row[unique_col], col2type[unique_col])
    sql <- paste("SELECT * FROM ", tablename, " WHERE ",
                 unique_col, "=", sqlval, " LIMIT 1", sep="")
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
gene_desc <- list(
    col2type=c(
        #`_gene_id`="INTEGER",       # internal id (PRIMARY KEY)
        entrez_gene_id="INTEGER",   # if can't be "INTEGER" then use _gene_id as prim key
        gene_symbol="TEXT",
        gene_title="TEXT",
        cytoband="TEXT"
    ),
    col2key=c(
        #`_gene_id`="PRIMARY KEY",
        # entrez_gene_id="UNIQUE",
        entrez_gene_id="PRIMARY KEY",
        gene_symbol="UNIQUE"
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
mrna_assignment_desc <- list(
    col2type=c(
        source_name="TEXT",
        description="TEXT",
        assignment_seqname="TEXT",
        assignment_score="INTEGER",
        assignment_coverage="INTEGER",
        direct_probes="INTEGER",
        possible_probes="INTEGER",
        xhyb="INTEGER",
        probeset_id="INTEGER",      # REFERENCES probeset(probeset_id)
        `_mrna_id`="INTEGER"        # REFERENCES mrna(_mrna_id)
    ),
    col2key=c(
        probeset_id="NOT NULL REFERENCES probeset(probeset_id)",
        `_mrna_id`="NOT NULL REFERENCES mrna(_mrna_id)"
    )
)

### The "swissprot" table.
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
        GO_evidence="TEXT",
        `_mrna_id`="INTEGER"        # REFERENCES mrna(_mrna_id)
    ),
    col2key=c(
        `_mrna_id`="NOT NULL REFERENCES mrna(_mrna_id)"
    )
)

### The "pathway" table.
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
protein_domains_desc <- list(
    col2type=c(
        source="TEXT",
        accession_or_domain_name="TEXT",
        domain_description="TEXT",
        `_mrna_id`="INTEGER"        # REFERENCES mrna(_mrna_id)
    ),
    col2key=c(
        `_mrna_id`="NOT NULL REFERENCES mrna(_mrna_id)"
    )
)

### The "protein_families" table.
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

### Global schema (12 tables).
NETAFFX_HUEX_TRANSCRIPT_DB_schema <- list(
    probeset=probeset_desc,
    gene=gene_desc,
    mrna=mrna_desc,
    mrna_assignment=mrna_assignment_desc,
    swissprot=swissprot_desc,
    unigene=unigene_desc,
    GO_biological_process=GO_biological_process_desc,
    GO_cellular_component=GO_biological_process_desc,
    GO_molecular_function=GO_biological_process_desc,
    pathway=pathway_desc,
    protein_domains=protein_domains_desc,
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
    }
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
            if (nsubfields < min.nsubfields || nsubfields > ncol)
                stop("bad number of subfields")
            if (nsubfields < ncol)
                length(rows[[i]]) <- ncol
        }
        data <- unlist(rows)
        data[data == "---"] <- NA
        mat <- matrix(data=data, ncol=ncol, byrow=TRUE)
    }
    dimnames(mat) <- list(mat[ , 1], subfields)
    mat
}

### Return the number of parts in the multipart field (should always be equal
### to the number of rows that are inserted in 'tablename').
insert_NetAffx_multipart_field <- function(conn, tablename, multipart_val, probeset_id)
{
    if (is.na(multipart_val))
        return(0)
    col2type <- NETAFFX_HUEX_TRANSCRIPT_DB_schema[[tablename]]$col2type
    vals <- strsplit(multipart_val, " /// ", fixed=TRUE)[[1]]
    for (val in vals) {
        row <- strsplit(val, " // ", fixed=TRUE)[[1]]
        row[row == "---"] <- NA
        if (tablename %in% c("unigene", "protein_domains")
         && length(row) == length(col2type) - 2)
            row <- c(row, NA)
        row <- c(row, probeset_id)
        #names(row) <- names(col2type)
        dbInsertRow(conn, tablename, row, col2type)
    }
    return(length(vals))
}

insert_NetAffx_HuEx_transcript_data <- function(conn, transcript)
{
    .mrna.id <- 0
    for (i in 1:nrow(transcript)) {
        if (is(conn, "SQLiteConnection"))
            dbBeginTransaction(conn)
        row <- unlist(transcript[i, ])
        probeset_id <- row[["probeset_id"]] # instead of [ ] to get rid of the name
        if (!is(conn, "DBIConnection"))
            .dbSendQuery(conn, paste("-- probeset_id ", probeset_id, sep=""))
        row[row == "---"] <- NA

        ## Extract simple fields
        probeset_row <- row[names(probeset_desc$col2type)]
        names(probeset_row) <- NULL # should make dbInsertRow() slightly faster
        dbInsertRow(conn, "probeset", probeset_row, probeset_desc$col2type)

        ## Extract multipart fields
        gene_assignment <- multipartToMatrix(row["gene_assignment"], gene_assignment_subfields)
        mrna_assignment <- multipartToMatrix(row["mrna_assignment"], mrna_assignment_subfields)
        for (i in nrow(mrna_assignment)) {
            accession <- mrna_assignment[i, "accession"]
            i2 <- which(gene_assignment[ , "accession"] %in% accession)
            if (length(i2) >= 2)
                stop("in CSV line for probeset_id=", probeset_id, ": ",
                     "\"gene_assignment\" has more than 1 gene linked to ", accession)
            if (length(i2) == 1) {
                gene_row <- gene_assignment[i2, names(gene_desc$col2key)]
                row0 <- dbGetThisRow(conn, "gene", "entrez_gene_id", gene_row, gene_desc$col2key)
                if (is.null(row0))
                    dbInsertRow(conn, "gene", gene_row, gene_desc$col2key)
                gene_assignment <- gene_assignment[-i2, ]
                entrez_gene_id <- gene_row["entrez_gene_id"]
            } else {
                entrez_gene_id <- NA
            }
            mrna_row <- c(NA, accession, entrez_gene_id)
            row0 <- dbGetThisRow(conn, "mrna", "accession", mrna_row, mrna_desc$col2type)
            if (is.null(row0)) {
                mrna_id0 <- .mrna.id <- .mrna.id + 1
                mrna_row[1] <- .mrna.id
                dbInsertRow(conn, "mrna", mrna_row, mrna_desc$col2key)
            } else {
                mrna_id0 <- row0[1]
            }
        }
        if (nrow(gene_assignment) != 0)
            stop("in CSV line for probeset_id=", probeset_id, ": ",
                 "\"gene_assignment\" has unlinked genes")

        #for (tablename in names(NETAFFX_HUEX_TRANSCRIPT_DB_schema)[-1]) {
        #    multipart_val <- row[tablename]
        #    insert_NetAffx_multipart_field(conn, tablename, multipart_val, probeset_id)
        #}
        if (is(conn, "SQLiteConnection"))
            dbCommit(conn)
    }
}

build_NetAffx_HuEx_transcript_SQL <- function(csv_file, sql_file)
{
    data <- read.table(csv_file, header=TRUE, sep=",", quote="\"", stringsAsFactors=FALSE)
    conn <- dbConnect(dbDriver("SQLite"), dbname=sql_file)
    create_NetAffx_HuEx_transcript_tables(conn)
    insert_NetAffx_HuEx_transcript_data(conn, data)
    dbDisconnect(conn)
}

