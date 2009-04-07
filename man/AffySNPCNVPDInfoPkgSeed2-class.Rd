\name{AffySNPCNVPDInfoPkgSeed2-class}
\Rdversion{1.1}
\docType{class}
\alias{AffySNPCNVPDInfoPkgSeed2-class}
%% \alias{chipName,AffySNPCNVPDInfoPkgSeed2-method}
%% \alias{makePdInfoPackage,AffySNPCNVPDInfoPkgSeed2-method}

\title{Class "AffySNPCNVPDInfoPkgSeed2"}
\description{A generic annotation package builder for Affymetrix SNP/CNV
  arrays. This is a simplified version of the annotation package and
  \code{crlmm} will *NOT* work for them.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffySNPCNVPDInfoPkgSeed2", csvAnnoFileCnv, csvSeqFileCnv, ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{cdfFile}:}{Path to the CDF file for this.}
    \item{\code{csvAnnoFile}:}{Path to the Affymetrix CSV annotation for
      the SNP probes.}
    \item{\code{csvSeqFile}:}{Path to the (SNP) probe sequence file.}
    \item{\code{csvAnnoFileCnv}:}{Path to the Affymetrix CSV annotation
      for the CNV probes.}
    \item{\code{csvSeqFileCnv}:}{Path to the (CNV) probe sequence file.}
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
    \item{chipName}{\code{signature(object = "AffySNPCNVPDInfoPkgSeed2")}: ... }
    \item{makePdInfoPackage}{\code{signature(object = "AffySNPCNVPDInfoPkgSeed2")}: ... }
	 }
}
\author{Benilton Carvalho}
\note{Annotation packages created from this package do *NOT* support \code{crlmm}.}
\examples{
showClass("AffySNPCNVPDInfoPkgSeed2")
}
\keyword{classes}
