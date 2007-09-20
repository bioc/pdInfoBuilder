
library("pdInfoBuilder")

targ = "../pd.mapping50k.hind240"

if (file.exists(targ)) stop(paste("cannot have folder", targ,
   "present if you want to build there"))
parms_pref = "../../parms_store/pd.mapping50k.hind240"  # should be in parms_store
anno_pref = "../../anno_store/pd.mapping50k.hind240"
pat = function(x) paste(parms_pref, x, sep="/")
pat2 = function(x) paste(anno_pref, x, sep="/")

# next 4 lines gunzip contents of anno_store, if needed
curd = getwd()
setwd(anno_pref)
try(system( paste(c("gunzip", dir(patt="gz$")), collapse=" ")))
setwd(curd)

# REMEMBER TO GUNZIP for NOW ... eventually use gzfiles?
cdfFile <- pat2("Mapping50K_Hind240.cdf")
csvAnno <- pat2("Mapping50K_Hind240.na22.annot.csv")
csvSeq <- pat2("Mapping50K_Hind_probe_tab")
spline <- pat("pd.mapping50k.hind240.spline.params.rda")
refd <- pat("pd.mapping50k.hind240Ref.rda")
crlmmInf <- pat("pd.mapping50k.hind240CrlmmInfo.rda")

pkg <- new("AffySNPPDInfoPkgSeed",
           version="0.3.4",
           author="Vince Carey", email="stvjc@channing.harvard.edu",
           biocViews="AnnotationData",
           genomebuild="NCBI Build 36",
           cdfFile=cdfFile, csvAnnoFile=csvAnno, csvSeqFile=csvSeq,
           splineParamFile=spline, crlmmInfoFile=crlmmInf,
           referenceDistFile=refd)


makePdInfoPackage(pkg, destDir="..")

# next 4 lines gzip contents of anno_store
curd = getwd()
setwd(anno_pref)
try(system( paste(c("gzip", dir()), collapse=" ")))
setwd(curd)
