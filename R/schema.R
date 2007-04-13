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

# NB -- VC added cytoband

createSnpFeatureSetSql <- ('
create table featureSet (
    fsetid integer primary key,
    man_fsetid text,
    affy_snp_id integer,
    dbsnp_rs_id text,
    chrom text,
    physical_pos integer,
    strand integer,
    cytoband text,
    allele_a text,
    allele_b text,
    fragment_length integer)
')

createSnpFeatureSql <- ('
create table %s (
    fid integer primary key,
    strand integer,
    allele integer,
    fsetid integer not null references "featureSet" ("fsetid"),
    pos integer,
    x integer,
    y integer)
')

createSnpPm_MmSql <- ('
create table %s (
    pm_fid integer primary key references "pmfeature" ("fid"),
    mm_fid integer references "mmfeature" ("fid"))
')

createSnpSequenceSql <- ('
create table sequence (
    fid integer primary key,
    offset integer,
    tstrand text,
    tallele text,
    seq text)
')
