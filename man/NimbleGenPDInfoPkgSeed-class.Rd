\name{NimbleGenPDInfoPkgSeed-class}
\docType{class}
\alias{NimbleGenPDInfoPkgSeed-class}


\title{Class "NimbleGenPDInfoPkgSeed"}
\description{PDInfo package Seed for all NimbleGen arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NimbleGenPDInfoPkgSeed", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{manufacturer}:}{Manufacturer = NimbleGen}
    \item{\code{chipName}:}{Name of the chip or platform}
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
    \item{chipName}{\code{signature(object = "NimbleGenPDInfoPkgSeed")}: ... }
    \item{getGeometry}{\code{signature(object = "NimbleGenPDInfoPkgSeed")}: ... }
	 }
}
\examples{
showClass("NimbleGenPDInfoPkgSeed")
}
\keyword{classes}
