\name{AffyExpressionPDInfoPkgSeed-class}
\docType{class}
\alias{AffyExpressionPDInfoPkgSeed-class}
\alias{chipName,AffyExpressionPDInfoPkgSeed-method}
\alias{initialize,AffyExpressionPDInfoPkgSeed-method}
\alias{makePdInfoPackage,AffyExpressionPDInfoPkgSeed-method}

\title{Class "AffyExpressionPDInfoPkgSeed"}
\description{PD Info Package Seed for Affymetrix Expression Arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffyExpressionPDInfoPkgSeed", cdfFile, csvAnnoFile, tabSeqFile, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{cdfFile}:}{CDF filename}
    \item{\code{csvAnnoFile}:}{CSV annotation file}
    \item{\code{tabSeqFile}:}{TAB sequence file}
    \item{\code{chipName}:}{NA}
    \item{\code{manufacturer}:}{Affymetrix}
    \item{\code{genomebuild}:}{Genome Build}
    \item{\code{pdInfoObjectName}:}{NA}
    \item{\code{geometry}:}{Array's geometry}
    \item{\code{name}:}{name}
    \item{\code{version}:}{version}
    \item{\code{license}:}{license}
    \item{\code{author}:}{author}
    \item{\code{email}:}{email}
    \item{\code{url}:}{url}
    \item{\code{biocViews}:}{biocViews}
  }
}
\section{Extends}{
Class \code{"\linkS4class{PDInfoPkgSeed}"}, directly.
Class \code{"\linkS4class{PkgSeed}"}, by class "PDInfoPkgSeed", distance 2.
}
\section{Methods}{
  \describe{
    \item{chipName}{chipName}
    \item{initialize}{initialize}
    \item{makePdInfoPackage}{package creator}
	 }
}
\examples{
showClass("AffyExpressionPDInfoPkgSeed")
}
\keyword{classes}
