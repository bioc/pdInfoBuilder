\name{NgsTilingPDInfoPkgSeed-class}
\docType{class}
\alias{NgsTilingPDInfoPkgSeed-class}
\alias{initialize,NgsTilingPDInfoPkgSeed-method}
\alias{makePdInfoPackage,NgsTilingPDInfoPkgSeed-method}

\title{Class "NgsTilingPDInfoPkgSeed"}
\description{PDInfo package Seed for NimbleGen Tiling arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NgsTilingPDInfoPkgSeed", ndfFile, posFile, pairFile, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{ndfFile}:}{NDF (NimbleGen Design) file}
    \item{\code{posFile}:}{POS (Positions) file}
    \item{\code{pairFile}:}{XYS File - used as template}
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
showClass("NgsTilingPDInfoPkgSeed")
}
\keyword{classes}
