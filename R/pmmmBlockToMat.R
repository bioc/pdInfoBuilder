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
    nGroups <- length(u$groups)

    ## BEGIN: Benilton's modification
    
    ## The 2 lines below added by Benilton to fix the bug on the alleles
    ## Assuming either 2 or 4 groups. Need to keep in mind that if:
    ## 2 groups: Group 1 is all ALLELE A and Group 2 is all ALLELE B
    ## 4 groups: Group 1 (Allele A Sense) Group 2 (Allele B Sense)
    ##           Group 3 (Allele A Antisense) Group 4 (Allele B Antisense)
    allele <- rep(c(ALLELE_A, ALLELE_B), nGroups/2)
    for (i in 1:nGroups) u$groups[[i]]$allele=as.integer(allele[i])

    ## The code below is the original by Seth (commented uot)
    ## mat <- do.call(rbind, lapply(u$groups, function(g) {
    ##     dir <- groupDirectionToInt(g$groupdirection)
    ##      N <- length(g$x)
    ##     cbind(do.call(cbind, g[cols]), strand=rep(dir, N),
    ##           allele=rep(c(ALLELE_A, ALLELE_B), each=2, length=N))
    ## }))

    ## The code below is the replacement for the above
    ## maybe there's a better way of doing this, that's
    ## why I left all the original code.
    mat <- do.call(rbind, lapply(u$groups, function(g) {
      dir <- groupDirectionToInt(g$groupdirection)
      N <- length(g$x)
      cbind(do.call(cbind, g[cols]), strand=rep(dir, N),
            allele=rep(g$allele, length=N))}))

    ## End of Benilton's modifications
    
    if (verify.pmmm) {
        ## verify pm/mm match
        ispm <- as.logical(mat[, "ispm"])
        stopifnot(all(mat[ispm, "atom"] == mat[!ispm, "atom"]))
    }
    mat
}


readCdfUnitToMat.affyExpr <- function(u, verify.pmmm=TRUE) {
  cols <- c("x", "y", "indices", "pbase", "tbase", "atom", "ispm")
  bases <- c("A", "C", "G", "T")
  u$groups[[1]][["pbase"]] <- match(u$groups[[1]]$pbase, bases)
  u$groups[[1]][["tbase"]] <- match(u$groups[[1]]$tbase, bases)
  mat <- do.call(cbind, u$groups[[1]][cols])
  if (verify.pmmm) {
    ## verify pm/mm match
    ispm <- as.logical(mat[, "ispm"])
    stopifnot(all(mat[ispm, "atom"] == mat[!ispm, "atom"]))
  }
  mat
}

readCdfUnitToMat.affyTiling2 <- function(u, nx, verify.pmmm=TRUE) {
  cols.pm <- c("pmx", "pmy", "strand", "startpos")
  cols.mm <- c("mmx", "mmy", "strand", "startpos")

  mat <- rbind(cbind(do.call(cbind, u[cols.pm]), atom=1:length(u[[2]]), ispm=1),
               cbind(do.call(cbind, u[cols.mm]), atom=1:length(u[[2]]), ispm=0))
  colnames(mat) <- c("x", "y", "strand", "startpos", "atom", "ispm")
  mat <- cbind(mat, indices=(mat[,"x"]+1+mat[,"y"]*nx))
  if (verify.pmmm) {
    ## verify pm/mm match
    ispm <- as.logical(mat[, "ispm"])
    stopifnot(all(mat[ispm, "atom"] == mat[!ispm, "atom"]))
  }
  mat
}

readCdfUnitToMat.affyTiling <- function(u, nx, verify.pmmm=TRUE) {
  cols.pm <- c("pmx", "pmy", "strand", "startpos")
  cols.mm <- c("mmx", "mmy", "strand", "startpos")

  if (u$seqInfo$mapping == "pmmm"){
    mat <- rbind(cbind(do.call(cbind, u[cols.pm]), atom=1:length(u[[2]]), ispm=1),
                 cbind(do.call(cbind, u[cols.mm]), atom=1:length(u[[2]]), ispm=0))
  } else if (u$seqInfo$mapping == "onlypm"){
    mat <- cbind(do.call(cbind, u[cols.pm]), atom=1:length(u[[2]]), ispm=1)
  } else {
    stop("u$seqInfo$mapping is something other than pmmm/onlypm")
  }
  
  colnames(mat) <- c("x", "y", "strand", "startpos", "atom", "ispm")
  mat <- cbind(mat, indices=(mat[,"x"]+1+mat[,"y"]*nx))
  mat <- as.data.frame(mat)

  pm.seq <- u[["probeseq"]]
  if (u$seqInfo$mapping == "pmmm"){
    mid.base <- substr(pm.seq, 13, 13)
    iA <- which(toupper(mid.base) == "A")
    iC <- which(toupper(mid.base) == "C")
    iG <- which(toupper(mid.base) == "G")
    iT <- which(toupper(mid.base) == "T")
    mid.base[iA] <- "T"
    mid.base[iC] <- "G"
    mid.base[iG] <- "C"
    mid.base[iT] <- "A"
    mm.seq <- paste(substr(pm.seq, 1, 12), mid.base,
                    substr(pm.seq, 14, 25), sep="")
    mat <- data.frame(mat, probeseq=c(pm.seq, mm.seq),
                      stringsAsFactors=FALSE)
  }else{
    mat <- data.frame(mat, probeseq=pm.seq, stringsAsFactors=FALSE)
    verify.pmmm <- FALSE
  }
  
  if (verify.pmmm) {
    ## verify pm/mm match
    ispm <- as.logical(mat[, "ispm"])
    stopifnot(all(mat[ispm, "atom"] == mat[!ispm, "atom"]))
  }
  mat
}
