pmmmBlockToMat <- function(b) {
    direction <- b$direction
    b <- do.call(rbind, lapply(b$groups, function(g) {
        do.call(cbind, lapply(g[c("indices", "x", "y", "expos")], t))
    }))
    nrb <- nrow(b)
    allele <- rep(c(ALLELE_A, ALLELE_B), nrb/2)
    if (direction == 3) {
        strand <- rep(c(SENSE, SENSE, ANTISENSE, ANTISENSE), nrb/4)
    } else if (direction == 2) {
        strand <- rep(ANTISENSE, nrb)
    } else if (direction == 1) {
        strand <- rep(SENSE, nrb)
    }
    b <- cbind(b, strand, strand)
    b <- cbind(b, allele, allele)

    colnames(b) <- rep(c("indices", "x", "y", "expos", "strand",
                         "allele"), each=2)
    b
}

groupDirectionToInt <- function(dir) {
    switch(dir,
        sense=SENSE,
        antisense=ANTISENSE,
        stop("unknown direction: ", dir))
}

readCdfUnitToMat <- function(u, verify.pmmm=TRUE) {
    ## Return a matrix when given a unit returned from
    ## affxparser::readCdf.  Columns are x, y, indices, atom, indexpos,
    ## ispm, strand, allele.
    ##
    ## XXX: we rely on the probes always being given in alleleA, alleleB
    ##      order.
    cols <- c("x", "y", "indices", "atom", "indexpos", "ispm")
    mat <- do.call(rbind, lapply(u$groups, function(g) {
        dir <- groupDirectionToInt(g$groupdirection)
         N <- length(g$x)
        cbind(do.call(cbind, g[cols]), strand=rep(dir, N),
              allele=rep(c(ALLELE_A, ALLELE_B), each=2, length=N))
    }))
    if (verify.pmmm) {
        ## verify pm/mm match
        ispm <- as.logical(mat[, "ispm"])
        stopifnot(all(mat[ispm, "atom"] == mat[!ispm, "atom"]))
    }
    mat
}

