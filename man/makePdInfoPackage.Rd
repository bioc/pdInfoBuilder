\name{makePdInfoPackage}
\alias{makePdInfoPackage}
\alias{makePdInfoPackage,AffySNPPDInfoPkgSeed-method}
\alias{makePdInfoPackage,AffySNPCNVPDInfoPkgSeed-method}
\alias{makePdInfoPackage,AffyExpressionPDInfoPkgSeed-method}
\alias{makePdInfoPackage,AffyTilingPDInfoPkgSeed-method}
\alias{makePdInfoPackage,AffySTPDInfoPkgSeed-method}
\alias{makePdInfoPackage,NgsExpressionPDInfoPkgSeed-method}
\alias{makePdInfoPackage,NgsTilingPDInfoPkgSeed-method}

\title{Create a Platform Design Info Package}
\description{
  This generic function create a platform design info package based on
  the parameters contained in \code{object} which will generally be an
  instance of a subclass of \code{PkgSeed}.  The result is a new
  directory on the filesystem containing the source for the generated
  pdInfo package.
}
\usage{
makePdInfoPackage(object, destDir, batch_size = 10000, quiet = FALSE, unlink = FALSE)
}

\arguments{
  \item{object}{See \code{showMethods("makePdInfoPackage")} to see
  available methods.}
  \item{destDir}{Path where the resulting pdInfo package source
  directory will be written.}
  \item{batch_size}{An integer controlling the size of batches
  processed when reading the flatfiles and loading the DB.  In
  general, larger
  values of \code{batch_size} will use more memory and less time
  (unless you exceed physical memory, in which case more time will be
  used as well).}
  \item{quiet}{A logical value.  When \code{TRUE}, diagnostic and
  status messages are not printed.}
  \item{unlink}{A logical value. If 'TRUE', and 'destDir' already contains a
          file or directory with the name 'pkgname', try to unlink
          (remove) it.}
}
\details{
        In general, creating the SQLite database will be a time and
        memory intensive task.
}

\value{
  This function is called for its side-effect of producing a pdInfo
  source package directory.
}

\author{Seth Falcon}

\examples{
cdfFile <- "Mapping250K_Nsp.cdf"
csvAnno <- "Mapping250K_Nsp_annot.csv"
csvSeq <- "Mapping250K_Nsp_probe_tab"

pkg <- new("AffySNPPDInfoPkgSeed", 
           version="0.1.5",
           author="A.U. Thor", email="au@thor.net",
           biocViews="AnnotationData",
           genomebuild="NCBI Build 35, May 2004",
           cdfFile=cdfFile, csvAnnoFile=csvAnno, csvSeqFile=csvSeq)

\dontrun{
makePdInfoPackage(pkg, destDir=".")
}
}

\keyword{methods}

