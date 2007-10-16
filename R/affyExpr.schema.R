BASE_A <- 1
BASE_C <- 2
BASE_G <- 3
BASE_T <- 4

createAffyExprFeatureSetSql <- ('
create table featureSet (
    fsetid integer primary key,
    man_fsetid text,
    alignment text,
    gene_symbol text,
    chrom text,
    ensembl text,
    strand integer)
')

createAffyExprFeatureSql <- ('
create table %s (
    fid integer primary key,
    fsetid integer not null references "featureSet" ("fsetid"),
    pbase integer,
    tbase integer,
    atom integer,
    x integer,
    y integer)
')

createAffyExprPm_MmSql <- ('
create table %s (
    pm_fid integer primary key references "pmfeature" ("fid"),
    mm_fid integer references "mmfeature" ("fid"))
')

createAffyExprSequenceSql <- ('
create table sequence (
    fid integer primary key,
    tstrand integer,
    interrogation_position integer,
    seq text)
')
