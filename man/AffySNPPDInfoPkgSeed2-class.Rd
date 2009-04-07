\name{AffySNPPDInfoPkgSeed2-class}
\Rdversion{1.1}
\docType{class}
\alias{AffySNPPDInfoPkgSeed2-class}
%% \alias{chipName,AffySNPPDInfoPkgSeed2-method}

\title{Class "AffySNPPDInfoPkgSeed2"}
\description{A generic annotation package builder for Affymetrix SNP
  arrays. This is a simplified version of the annotation package and
  \code{crlmm} will *not* work for them.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffySNPPDInfoPkgSeed2", cdfFile, csvAnnoFile, csvSeqFile, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{cdfFile}:}{CDF file for the design.}
    \item{\code{csvAnnoFile}:}{Affymetrix CSV Annotation file.}
    \item{\code{csvSeqFile}:}{Affymetrix Probe Sequence file.}

    \item{\code{chipName}:}{Name of the chip or platform}
    \item{\code{manufacturer}:}{chip/platform manufacturer}
    \item{\code{url}:}{chip URL}
    \item{\code{genomebuild}:}{The genome build this platform is based
      upon.}
    \item{\code{organism}:}{ organism for chip. }
    \item{\code{species}:}{ species for chip. }  
    \item{\code{version}:}{A character vector giving the version
      number of the package.}
    \item{\code{license}:}{The license of the package}
    \item{\code{author}:}{Author of the package}
    \item{\code{email}:}{An email address to use in the Maintainer field}
    \item{\code{biocViews}:}{Character data for the biocViews field of
    the DESCRIPTION file}
  }
}
\section{Methods}{
  \describe{
    \item{chipName}{\code{signature(object = "AffySNPPDInfoPkgSeed2")}: ... }
  }
}

\note{
  This is a simplified annotation package. CRLMM won't work for these
  objects.
  }

\examples{
showClass("AffySNPPDInfoPkgSeed2")
}
\keyword{classes}
