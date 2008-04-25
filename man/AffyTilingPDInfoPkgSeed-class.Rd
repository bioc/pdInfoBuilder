\name{AffyTilingPDInfoPkgSeed-class}
\docType{class}
\alias{AffyTilingPDInfoPkgSeed-class}
\alias{initialize,AffyTilingPDInfoPkgSeed-method}
\alias{makePdInfoPackage,AffyTilingPDInfoPkgSeed-method}
\alias{chipName,AffyTilingPDInfoPkgSeed-method}

\title{Class "AffyTilingPDInfoPkgSeed"}
\description{PD Info Package Seed for Affymetrix Tiling Arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffyTilingPDInfoPkgSeed", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{bpmapFile}:}{BPMAP File - provided by Affymetrix}
    \item{\code{cifFile}:}{CIF File - provided by Affymetrix}
    \item{\code{chipName}:}{Automatically determined.}
    \item{\code{manufacturer}:}{Automatically determined.}
    \item{\code{genomebuild}:}{Genome Build.}
    \item{\code{pdInfoObjectName}:}{Automatically determined.}
    \item{\code{geometry}:}{Automatically determined.}
    \item{\code{name}:}{Automatically determined.}
    \item{\code{version}:}{Version}
    \item{\code{license}:}{License}
    \item{\code{author}:}{Author}
    \item{\code{email}:}{Email}
    \item{\code{url}:}{URL}
    \item{\code{biocViews}:}{biocViews}
  }
}
\section{Extends}{
Class \code{"\linkS4class{PDInfoPkgSeed}"}, directly.
Class \code{"\linkS4class{PkgSeed}"}, by class "PDInfoPkgSeed", distance 2.
}
\section{Methods}{
  \describe{
    \item{makePdInfoPackage}{\code{signature(object = "AffyTilingPDInfoPkgSeed")}: ... }
	 }
    \item{chipName}{\code{signature(object = "AffyTilingPDInfoPkgSeed")}: ... }
  }

\examples{
showClass("AffyTilingPDInfoPkgSeed")
}
\keyword{classes}
