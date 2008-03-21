## This is the schema file
## for the Affy Gene ST arrays
## Started: March/08 - Benilton Carvalho

## Strategy:
## featureSet table contains only *required* fields from PGF
## sequence table contains info from probe.tab
## likely to have only pmfeature table.

BASE_A <- 1
BASE_C <- 2
BASE_G <- 3
BASE_T <- 4

createAffySTFeatureSetSql <- ('
create table featureSet (
    fsetid integer primary key,
    man_fsetid text,
    type text,
    start_atom integer,
    transcript_cluster_id integer,
    seqname text,
    strand integer,
    start integer,
    stop integer,
    total_probes integer,
    gene_assignment text,
    mrna_assignment text,
    swissprot text,
    unigene text,
    GO_biological_process text,
    GO_cellular_component text,
    GO_molecular_function text,
    pathway text,
    protein_domains text,
    crosshyb_type integer,
    category text
)
')

createAffySTFeatureSql <- ('
create table %s (
    fid integer,
    fsetid integer not null references "featureSet" ("fsetid"),
    pbase integer,
    tbase integer,
    atom integer,
    x integer,
    y integer,
    gc_count integer,
    strand integer,
    PRIMARY KEY (atom, fid))
')

createAffySTPm_MmSql <- ('
create table %s (
    pm_fid integer primary key references "pmfeature" ("fid"),
    mm_fid integer references "mmfeature" ("fid"))
')


##  info here could actually on on the *feature tables but on ST arrays,
##  probe ID (fid) are not necessarily unique to a probeset (fsetid) and
##  many of the things would be duplicated (but it is still duplicated
##  here with transcript_cluster and trancript+probeid do not form a key
##  =(

createAffySTSequenceSql <- ('
create table sequence (
    fid integer,
    tstrand integer,
    interrogation_position integer,
    transcript_cluster integer,
    seqname text,
    start integer,
    stop integer,
    seq text,
    category text)
')
