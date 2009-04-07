SENSE <- 0
ANTISENSE <- 1
ALLELE_A <- 0
ALLELE_B <- 1

## setPageSizeSql <- ('
## pragma page_size = 8192;
## ')

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
    gene_assoc text,
    fragment_length integer,
    dbsnp integer,
    cnv text)
')

createSnp6FeatureSetSql <- ('
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
    gene_assoc text,
    fragment_length integer,
    fragment_length2 integer,
    dbsnp integer,
    cnv text)
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


createCnvFeatureSetSql <- ('
create table featureSetCNV (
    fsetid integer primary key,
    man_fsetid text,
    chrom text,
    chrom_start integer,
    chrom_stop integer,
    strand integer,
    cytoband text,
    gene_assoc text,
    fragment_length integer,
    xpar integer,
    cnv text)
')

createCnv6FeatureSetSql <- ('
create table featureSetCNV (
    fsetid integer primary key,
    man_fsetid text,
    chrom text,
    chrom_start integer,
    chrom_stop integer,
    strand integer,
    cytoband text,
    gene_assoc text,
    fragment_length text,
    xpar integer,
    cnv text)
')

createCnvFeatureSql <- ('
create table %s (
    fid integer primary key,
    strand integer,
    fsetid integer not null references "featureSetCNV" ("fsetid"),
    x integer,
    y integer)
')

createCnvSequenceSql <- ('
create table sequenceCNV (
    fid integer primary key,
    offset integer,
    tstrand text,
    seq text)
')
