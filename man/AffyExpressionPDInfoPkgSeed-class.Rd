\name{AffyExpressionPDInfoPkgSeed-class}
\docType{class}
\alias{AffyExpressionPDInfoPkgSeed-class}

\title{Class "AffyExpressionPDInfoPkgSeed"}
\description{PD Info Package Seed for Affymetrix Expression Arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffyExpressionPDInfoPkgSeed", cdfFile, csvAnnoFile, tabSeqFile, ...)}.
}
\section{Slots}{
	 \describe{
	   \item{\code{cdfFile}:}{CDF filename}
	   \item{\code{celFile}:}{CEL filename}
	   \item{\code{tabSeqFile}:}{TAB sequence file}

%% these are from PDInfoPkgSeed-class
    \item{\code{chipName}:}{Name of the chip or platform}
    \item{\code{manufacturer}:}{chip/platform manufacturer}
    \item{\code{url}:}{chip URL}
    \item{\code{genomebuild}:}{The genome build this platform is based
    upon.}
	\item{\code{organism}:}{ organism for chip. }
	\item{\code{species}:}{ species for chip. }  
%% these are from PkgSeed-class
    \item{\code{name}:}{The name of the package.}
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
    \item{chipName}{chipName}
    \item{getGeometry}{initialize}
    \item{makePdInfoPackage}{package creator}
	 }
}
\examples{
showClass("AffyExpressionPDInfoPkgSeed")
}
\keyword{classes}
