
library("pdInfoBuilder")

targ = "../pd.genomewidesnp.6"

if (file.exists(targ)) stop(paste("cannot have folder", targ,
   "present if you want to build there"))
parms_pref = "../../parms_store/pd.genomewidesnp.6"  # should be in parms_store
anno_pref = "../../anno_store/pd.genomewidesnp.6"
pat = function(x) paste(parms_pref, x, sep="/")
pat2 = function(x) paste(anno_pref, x, sep="/")

# next 4 lines gunzip contents of anno_store, if needed
curd = getwd()
setwd(anno_pref)
try(system( paste(c("gunzip", dir(patt="gz$")), collapse=" ")))
setwd(curd)

# REMEMBER TO GUNZIP for NOW ... eventually use gzfiles?
cdfFile <- pat2("GenomeWideSNP_6.cdf")
csvAnno <- pat2("GenomeWideSNP_6.na24.annot.csv")
csvSeq <- pat2("GenomeWideSNP_6.probe_tab")
#../../anno_store/pd.genomewidesnp.6/GenomeWideSNP_6.probe_tab
csvAnnoCnv <- pat2("GenomeWideSNP_6.cn.na24.annot.csv")
csvSeqCnv <- pat2("GenomeWideSNP_6.CN_probe_tab")
spline <- pat("pd.genomewidesnp.6.spline.params.rda")
refd <- pat("pd.genomewidesnp.6Ref.rda")
crlmmInf <- pat("pd.genomewidesnp.6CrlmmInfo.rda")

pkg <- new("AffySNPCNVPDInfoPkgSeed",
           version="0.3.5",
           author="Vince Carey", email="stvjc@channing.harvard.edu",
           biocViews="AnnotationData",
           genomebuild="NCBI Build 36",
           cdfFile=cdfFile, csvAnnoFile=csvAnno, csvSeqFile=csvSeq,
           csvAnnoFileCnv=csvAnnoCnv, csvSeqFileCnv=csvSeqCnv,
           splineParamFile=spline,
           crlmmInfoFile=crlmmInf, referenceDistFile=refd)

makePdInfoPackage(pkg, destDir="..")

# next 4 lines gzip contents of anno_store
curd = getwd()
setwd(anno_pref)
try(system( paste(c("gzip", dir()), collapse=" ")))
setwd(curd)
