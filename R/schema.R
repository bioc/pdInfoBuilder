SENSE <- 0
ANTISENSE <- 1
ALLELE_A <- 0
ALLELE_B <- 1

setPageSizeSql <- ('
pragma page_size = 8192;
')

## dbsnp_rs_id could be integer if we strip the leading 'rs' Also, chrom could
## be integer and we could have a separate mapping tabel to the cromosome label
## (so that names are more meaningful).
##
## BC: Seth, I think the below should be something like "createSnpFeatureSetSql"
##     b/c later we want to have pdInfo's for TIling, Exon arrays... thoughts?

createFeatureSetSql <- ('
create table featureSet (
    fsetid integer primary key not null,
    man_fsetid text,
    affy_snp_id integer,
    dbsnp_rs_id text,
    chrom text,
    physical_pos integer,
    strand integer,
    allele_a text,
    allele_b text,
    UNIQUE("fsetid"))
')


createFeatureSql <- ('
create table %s (
    fid integer not null,
    strand integer,
    allele integer,
    fsetid integer not null references "featureSet" ("fsetid"),
    pos integer,
    x integer,
    y integer,
    UNIQUE("fid"))
')

createPm_MmSql <- ('
create table pm_mm (
    pm_fid integer primary key references "pmfeature" ("fid"),
    mm_fid integer references "mmfeature" ("fid"))
')

createSequenceSql <- ('
create table sequence (
    fid integer primary key,
    offset integer,
    tstrand text,
    tallele text,
    ispm integer,
    seq text)
')

## BC: We'll also need a table for the control probes
##     Examples of control probes are those with AFFX-12345
