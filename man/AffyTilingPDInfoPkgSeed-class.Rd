\name{AffyTilingPDInfoPkgSeed-class}
\docType{class}
\alias{AffyTilingPDInfoPkgSeed-class}

\title{Class "AffyTilingPDInfoPkgSeed"}
\description{PD Info Package Seed for Affymetrix Tiling Arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffyTilingPDInfoPkgSeed", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{bpmapFile}:}{BPMAP File - provided by Affymetrix}
    \item{\code{celFile}:}{CEL File - provided by Affymetrix}
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
    \item{makePdInfoPackage}{\code{signature(object = "AffyTilingPDInfoPkgSeed")}: ... }
    \item{chipName}{\code{signature(object = "AffyTilingPDInfoPkgSeed")}: ... }
  }
}

\examples{
showClass("AffyTilingPDInfoPkgSeed")
}
\keyword{classes}
