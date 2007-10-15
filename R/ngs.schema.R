createNgsFeatureSetSql <- ('
create table featureSet (
    fsetid integer primary key,
    man_fsetid text)
')

## NDF Contents in CAPITAL
## (x)   "PROBE_DESIGN_ID": not needed
## (f)   "CONTAINER": Feature Table
## (x)   "DESIGN_NOTE": not needed
## (x)   "SELECTION_CRITERIA": not needed
## (S)   "SEQ_ID": FeatureSet table (man_fsetid)
## (s)   "PROBE_SEQUENCE": Sequence table (seq)
## (x)   "MISMATCH": Will point to what table (pm/mm)Feature
## (f)   "MATCH_INDEX": pair PM/MM - Feature Table
## (f)   "FEATURE_ID": unit (4:9) - (unit_id)
## (x)   "ROW_NUM": not required
## (x)   "COL_NUM": not requied
## (x)   "PROBE_CLASS": internal / not required (fiducial/control/blabla)
## (f)   "PROBE_ID": Feature Table
## (f)   "POSITION": Feature Table
## (x)   "DESIGN_ID": not required
## (f)   "X": Feature Table
## (f)   "Y": Feature Table
## (?)   "DMD"

createNgsFeatureSql <- ('
create table %s (
    fid integer primary key,
    container text,
    unit_id integer,
    match_index integer,
    probe_id text,
    position integer,
    x integer,
    y integer,
    fsetid integer)
')

createNgsPm_MmSql <- ('
create table %s (
    pm_fid integer primary key references "pmfeature" ("fid"),
    mm_fid integer references "mmfeature" ("fid"))
')

createNgsSequenceSql <- ('
create table sequence (
    fid integer primary key,
    seq text)
')
