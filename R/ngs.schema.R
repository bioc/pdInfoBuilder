## incorporate blocks
## feature set table
createNgsFeatureSetSql <- ('
create table featureSet (
    fsetid integer primary key,
    man_fsetid text,
    alignment text,
    gene_symbol text,
    chrom text,
    ensembl text,
    strand integer,
    comment text)
')

## primary tables for pm and mm information
createNgsFeatureSql <- ('
create table %s (
    fid integer primary key,
    container text,
    unit_id integer,
    match_index integer,
    atom integer,
    probe_id text,
    x integer,
    y integer,
    position integer,
    fsetid integer not null references "featureSet" ("fsetid"))
')

##
createNgsPm_MmSql <- ('
create table %s (
    pm_fid integer primary key references "pmfeature" ("fid"),
    mm_fid integer references "mmfeature" ("fid"))
')

## create table sequence
createNgsSequenceSql <- ('
create table sequence (
    fid integer primary key,
    interrogation_position integer,
    seq text)
')

