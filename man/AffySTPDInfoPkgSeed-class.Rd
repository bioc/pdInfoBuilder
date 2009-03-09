\name{AffySTPDInfoPkgSeed-class}
\docType{class}
\alias{AffySTPDInfoPkgSeed-class}
\alias{AffyExonPDInfoPkgSeed-class}
\alias{AffyGenePDInfoPkgSeed-class}

\title{Class "AffySTPDInfoPkgSeed" for the Sense Target gene-level array }
\description{	 container for parameters related to pdmapping package construction for ST type arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffySTPDInfoPkgSeed", pgfFile, clfFile, probeFile, transFile, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{pgfFile}:}{Object of class \code{"ScalarCharacter"} path to pgf }
    \item{\code{clfFile}:}{Object of class \code{"ScalarCharacter"} path to clf }
    \item{\code{probeFile}:}{Object of class \code{"ScalarCharacter"}, path to probe sequence file (Optional) }
    \item{\code{transFile}:}{Object of class \code{"ScalarCharacter"}, path to trans file (Optional) }
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
\section{Extends}{
Class \code{"\linkS4class{PDInfoPkgSeed}"}, directly.
Class \code{"\linkS4class{PkgSeed}"}, by class "PDInfoPkgSeed", distance 2.
}
\section{Methods}{
  \describe{
    \item{chipName}{\code{signature(object = "AffySTPDInfoPkgSeed")}: ... }
    \item{getGeometry}{\code{signature(object = "AffySTPDInfoPkgSeed")}: ... }
    \item{makePdInfoPackage}{\code{signature(object = "AffySTPDInfoPkgSeed")}: ... }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ B. Carvalho }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("AffySTPDInfoPkgSeed")
}
\keyword{classes}