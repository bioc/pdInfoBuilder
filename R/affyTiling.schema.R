createAffyTilingFeatureSetSql <- ('
create table featureSet (
    fsetid integer primary key,
    man_fsetid text,
    groupname text,
    version text,
    fullname text,
    name text)
')

createAffyTilingFeatureSql <- ('
create table %s (
    fid integer primary key,
    fsetid integer not null references "featureSet" ("fsetid"),
    strand integer,
    startpos integer,
    atom integer,
    x integer,
    y integer)
')

createAffyTilingPm_MmSql <- ('
create table %s (
    pm_fid integer primary key references "pmfeature" ("fid"),
    mm_fid integer references "mmfeature" ("fid"))
')

createAffyTilingSequenceSql <- ('
create table sequence (
    fid integer primary key,
    seq text)
')
