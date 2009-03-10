## Gene ST Arrays have lots in common with Exon ST arrays.
## Most of the code can be reused.
## The problem is that it is *not* uncommon for
##    'probes' (in Gene ST) to map to *multiple* probesets.
## This breaks the code used for the Exon arrays,
##  because the pmfeature table relies on the the fact that
##  one_probe -> one_probeset (ie., 'fid' was the PRIMARY KEY)

## The solution here is to reimplement the pmfeature and bgfeature
##  tables to allow this 1:N mapping from probes to probesets.
## The easiest, and suboptimal, solution is to remove 'fid' as
##   primary key of such tables.
## Instead, I prefer to add another table (f2fset) that will be
##  in charge of mapping features (probes) to featureSets (probesets)

#######################################################################
## SECTION A - Db Schema
#######################################################################
f2fsetSchema <- list(col2type=c(
                       fid="INTEGER",
                       fsetid="INTEGER",
                       atom="INTEGER"),
                     col2key=c(
                       fsetid="REFERENCES featureSet(fsetid)"
                       ))

geneStPmFeatureSchema <- list(col2type=c(
                                fid="INTEGER",
                                x="INTEGER",
                                y="INTEGER"
                                ),
                              col2key=c(
                                fid="PRIMARY KEY"
                                ))

