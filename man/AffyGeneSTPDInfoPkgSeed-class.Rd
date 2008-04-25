\name{AffyGenePDInfoPkgSeed-class}
\docType{class}
\alias{AffyGenePDInfoPkgSeed-class}
\alias{chipName,AffyGenePDInfoPkgSeed-method}
\alias{initialize,AffyGenePDInfoPkgSeed-method}
\alias{makePdInfoPackage,AffyGenePDInfoPkgSeed-method}

\title{Class "AffyGenePDInfoPkgSeed" for the Sense Target gene-level array }
\description{	 container for parameters related to pdmapping package construction for ST gene level array}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffyGenePDInfoPkgSeed", pgfFile, clfFile, probeFile, transFile, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{pgfFile}:}{Object of class \code{"ScalarCharacter"} path to pgf }
    \item{\code{clfFile}:}{Object of class \code{"ScalarCharacter"} path to clf }
    \item{\code{probeFile}:}{Object of class \code{"ScalarCharacter"} ~~ }
    \item{\code{transFile}:}{Object of class \code{"ScalarCharacter"} ~~ }
    \item{\code{chipName}:}{Object of class \code{"character"} ~~ }
    \item{\code{manufacturer}:}{Object of class \code{"character"} ~~ }
    \item{\code{genomebuild}:}{Object of class \code{"character"} ~~ }
    \item{\code{pdInfoObjectName}:}{Object of class \code{"character"} ~~ }
    \item{\code{geometry}:}{Object of class \code{"integer"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{version}:}{Object of class \code{"character"} ~~ }
    \item{\code{license}:}{Object of class \code{"character"} ~~ }
    \item{\code{author}:}{Object of class \code{"character"} ~~ }
    \item{\code{email}:}{Object of class \code{"character"} ~~ }
    \item{\code{url}:}{Object of class \code{"character"} ~~ }
    \item{\code{biocViews}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{PDInfoPkgSeed}"}, directly.
Class \code{"\linkS4class{PkgSeed}"}, by class "PDInfoPkgSeed", distance 2.
}
\section{Methods}{
  \describe{
    \item{chipName}{\code{signature(object = "AffyGenePDInfoPkgSeed")}: ... }
    \item{initialize}{\code{signature(.Object = "AffyGenePDInfoPkgSeed")}: ... }
    \item{makePdInfoPackage}{\code{signature(object = "AffyGenePDInfoPkgSeed")}: ... }
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
showClass("AffyGenePDInfoPkgSeed")
}
\keyword{classes}
