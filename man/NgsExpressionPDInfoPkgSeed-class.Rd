\name{NgsExpressionPDInfoPkgSeed-class}
\docType{class}
\alias{NgsExpressionPDInfoPkgSeed-class}
\alias{initialize,NgsExpressionPDInfoPkgSeed-method}
\alias{makePdInfoPackage,NgsExpressionPDInfoPkgSeed-method}

\title{Class "NgsExpressionPDInfoPkgSeed"}
\description{PDInfo package Seed for NimbleGen Expression arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NgsExpressionPDInfoPkgSeed", ndfFile, pairFile, xysFile, ngdFile ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{ndfFile}:}{NDF (NimbleGen Design) file}
    \item{\code{pairFile}:}{PAIR (pair) file - used as template}
    \item{\code{xysFile}:}{XYS File - used as template}
    \item{\code{ngdFile}:}{NGD File - annotation data}
    \item{\code{chipName}:}{Design identifier}
    \item{\code{manufacturer}:}{Manufacturer: Affymetrix/NimbleGen}
    \item{\code{genomebuild}:}{Genome Build}
    \item{\code{pdInfoObjectName}:}{object name}
    \item{\code{geometry}:}{array dimensions}
    \item{\code{name}:}{name}
    \item{\code{version}:}{Version}
    \item{\code{license}:}{License}
    \item{\code{author}:}{Author}
    \item{\code{email}:}{Author's email}
    \item{\code{url}:}{website}
    \item{\code{biocViews}:}{biocViews}
      }
}
\section{Extends}{
Class \code{"\linkS4class{PDInfoPkgSeed}"}, directly.
Class \code{"\linkS4class{PkgSeed}"}, by class "PDInfoPkgSeed", distance 2.
}
\section{Methods}{
  \describe{
    \item{chipName}{\code{signature(object = "NgsPDInfoPkgSeed")}: ... }
    \item{initialize}{\code{signature(.Object = "NgsPDInfoPkgSeed")}: ... }
	 }
}
\examples{
showClass("NgsExpressionPDInfoPkgSeed")
}
\keyword{classes}
