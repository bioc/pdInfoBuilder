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
