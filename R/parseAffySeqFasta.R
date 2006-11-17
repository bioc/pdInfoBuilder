## DEBUG
if (FALSE) {

    fastad <- "/Users/seth/proj/SnpAnnBuilder-git/AffyMapping500k-data"
    fastaf <- file.path(fastad, "Mapping250K_Nsp_flanking_sequences_fasta")

    fname <- fastaf
    batch_size <- 10

    con <- file(fname, open="r")
    close(con)
}
## check netaffx to verify correctness of sequence parsing

## END DEBUG

parseAffySeqFasta <- function(con, n=10) {
    "Return a named list of sequence vectors named by featureSet ID

     con - connection to Affy sequence fasta file
       n - number of sequences to fetch
"
    if ((n %% 2) != 0)
      n <- n - 1
    n <- 2 * n ## we need 2n lines to get n records
    lines <- readLines(con, n=n)
    if (length(lines) == 0)
      return(NULL)
    if (length(lines) < n && (length(lines) %% 2) != 0)
      stop("unmatched read")
    fsetids <- sub("^>affx:", "", lines[seq(1, n, by=2)])
    fsetseqs <- strsplit(lines[seq(2, n, by=2)], "")
    names(fsetseqs) <- fsetids
    fsetseqs
}

encodeSeq <- function(s) {
    "Return an integer vector of length 3*length(s)
s should be a character vector of DNA bases (atcg).  Each element is
converted to a 3-tuple that encodes the letter as:
    A 1 0 0
    C 0 1 0
    T 0 0 1
    G 1 1 1
The return value is the concatenation of the encodings of each
base in s.
"
    ## A C T (G is all 1's)
    ##
    zeros <- rep(0L, 3)
    codes <- list(a=c(1L, 0L, 0L), c=c(0L, 1L, 0L), t=c(0L, 0L, 1L),
                  g=c(1L, 1L, 1L),
                  K=zeros, M=zeros, R=zeros, S=zeros, W=zeros, Y=zeros)
    tmp <- codes[tolower(s)]
    names(tmp) <- NULL
    do.call(c, tmp)
}

