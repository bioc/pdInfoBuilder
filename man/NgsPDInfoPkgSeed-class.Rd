\name{NgsPDInfoPkgSeed-class}
\docType{class}
\alias{NgsPDInfoPkgSeed-class}
\alias{initialize,NgsPDInfoPkgSeed-method}
\alias{NgsTilingPDInfoPkgSeed-class}
\alias{makePdInfoPackage,NgsTilingPDInfoPkgSeed-method}
\alias{NgsExpressionPDInfoPkgSeed-class}
\alias{makePdInfoPackage,NgsExpressionPDInfoPkgSeed-method}

\title{Class "NgsPDInfoPkgSeed"}
\description{PDInfo package Seed for NimbleGen arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NgsPDInfoPkgSeed", ndfFile, posFile, xysFile, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{ndfFile}:}{NDF (NimbleGen Design) file}
    \item{\code{posFile}:}{POS (Positions) file}
    \item{\code{xysFile}:}{XYS File - used as template}
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
showClass("NgsPDInfoPkgSeed")
}
\keyword{classes}
