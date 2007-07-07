pdinfo <- "@PDINFONAME@"

# the command
# eg. library(pd.mapping50k.hind240)
# will have been run in the harness

test_crlmm <- function() {
# this stuff will work for the 4 chips we know about now
# 100k hind240, xba240; 250k sty, nsp
size = as.character(2*as.numeric(gsub(".*g(.*)k.*", "\\1", pdinfo)))
enz = gsub(".*k.(.*$)", "\\1", pdinfo)
enz = gsub("240", "", enz)
hapmPackName = paste("hapmap", size, "k", enz, sep="")
library(oligo)
library(hapmPackName, character.only=TRUE)
xxr = justSNPRMA(dir(system.file( "celFiles", package=hapmPackName), full=TRUE))
xxc = crlmm(xxr, correctionFile="corr.rda")
if (exists("xxc")) return(TRUE)
return(FALSE)
}

