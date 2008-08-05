## modified by Matt Settles June 2,2008

## Nimblegen Header
## (s)	"software": 
## (s)	"version":
## (s)   "imagefile":
## (s)	"designfile":
## (s)	"designname":
##	(s)	"designid":
## (d)	"date":
## (i)	"border":
## (i)	"ul_x":
## (i)	"ul_y":
## (i)   "ur_x":
## (i)	"ur_y":
## (i)	"lr_x":
## (i)	"lr_y":
## (i)	"ll_x":
## (i)   "ll_y":
## (i)	"score":
## (i)	"qcscore":
##	(b)	"locallyaligned":
## (b)	"correctAstig":
## (?)	"Knots":
## (b)	"auto":

## readPairHeader works with xys files or pair files
"readPairHeader" <- 
function (file)
{
    firstfield <- scan(file, what = "", sep = "\t", quote = "\"", #"
        nlines = 100, flush = TRUE, quiet = TRUE, blank.lines.skip = FALSE,
        multi.line = FALSE, allowEscape = FALSE)
    NHeaderRecords <- grep("# software=", firstfield)
    txt <- scan(file, what = "", sep = "\t", quote = "\"", nlines = NHeaderRecords, #"
     quiet = TRUE, allowEscape = FALSE)
    substring(text=txt[1],first=3)
    txt <- strsplit(txt,split="=")
    txt <- data.frame(Value=sapply(txt,function(x) x[2]),row.names=sapply(txt,function(x) x[1]),stringsAsFactors=FALSE)
    out <- list(NHeaderRecords = NHeaderRecords, BeginRawData = NHeaderRecords)
    out$Version <- txt["version",]
    out$Date <- txt["date",]
    out$DesignName <- txt["designname",]
    out$DesignId <- txt["designid",]
    out
}#readNimbleGenHeader

# XYS Contents in CAPITAL - contains Header 
##	(x)	"X":
##	(x)	"Y":
##	(x)	"SIGNAL":
##	(x)	"COUNT":

# PAIR Contents in CAPITAL - contains Header
## (x)	"IMAGE_ID": 
## (x)	"GENE_EXPR_OPTION": 
## (x)	"SEQ_ID": 
## (x)   "PROBE_ID":
## (x) 	"POSITION": 
##	(x)	"X":
## (x)   "Y":
## (x)   "MATCH_INDEX": 
## (x)	"SEQ_URL": 
## (x)	"PM":
## (x)	"MM":

## NDF Contents in CAPITAL - no Header
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

## NGD Contents in CAPITAL - no Header
## (S)	"SEQ_ID": 
## (S) 	"COMMENT": 

## CALL Contents in CAPITAL - no Header
##	(x)	"IMAGE_ID":
## (x)	"SEQ_ID":
## (x)	"PROBE_PAIRS":
##	(x)	"FILTERED_PROBE_PAIRS":
## (x)	"PM_AVG":
##	(x)	"PM_CV":
##	(x)	"MM_AVG":
##	(x)	"MM_CV":
##	(x)	"DIFF_AVG":
##	(x)	"DIFF_CV":
##	(x)	"GENE_EXPR_OPTION":
##	(x)	"GENE_INFO":

## POS Contents in CAPITAL - no Header
##	(x)	"PROBE_ID":
##	(x)	"SEQ_ID":
##	(x)	"CHROMOSOME":
##	(x)	"POSITION":
##	(x)	"COUNT":
##	(x)	"LENGTH":
## (x)	"GC":

## GFF Contents - contains Header ##gff-version	3 9 columns - no column header
##	(x)	"seqid":
##	(x)	"source":
##	(x)	"type":
##	(x)	"start":
##	(x)	"end":
##	(x)	"score":
##	(x)	"strand":
##	(x)	"phase":
##	(x)	"attributes":
