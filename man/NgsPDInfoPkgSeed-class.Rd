\name{NgsPDInfoPkgSeed-class}
\docType{class}
\alias{NgsPDInfoPkgSeed-class}


\title{Class "NgsPDInfoPkgSeed"}
\description{PDInfo package Seed for all NimbleGen arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NgsPDInfoPkgSeed", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{manufacturer}:}{Manufacturer = NimbleGen}
    %% these are from PDInfoPkgSeed-class
    \item{\code{chipName}:}{Name of the chip or platform}
    %\item{\code{manufacturer}:}{chip/platform manufacturer}
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
\section{Extends}{
Class \code{"\linkS4class{PDInfoPkgSeed}"}, directly.
Class \code{"\linkS4class{PkgSeed}"}, by class "PDInfoPkgSeed", distance 2.
}
\section{Methods}{
  \describe{
    \item{chipName}{\code{signature(object = "NgsPDInfoPkgSeed")}: ... }
    \item{getGeometry}{\code{signature(object = "NgsPDInfoPkgSeed")}: ... }
	 }
}
\examples{
showClass("NgsPDInfoPkgSeed")
}
\keyword{classes}
