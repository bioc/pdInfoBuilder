\name{AffySNPPDInfoPkgSeed-class}
\docType{class}
\alias{AffySNPPDInfoPkgSeed-class}

\title{Class "AffySNPPDInfoPkgSeed"}
\description{
  This class represents Platform Design (PD) packages for Affymetrix
  mapping (SNP chip) arrays.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffySNPPDInfoPkgSeed", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{cdfFile}:}{Path to the CDF file for this
    platform.  Example: \code{"files/Mapping250K_Nsp.cdf"}}

    \item{\code{csvAnnoFile}:}{Path to the Affymetrix CSV annotation
    file.  Example: \code{"files/Mapping250K_Nsp_annot.csv"}}

    \item{\code{csvSeqFile}:}{Path to the Affymetrix probe sequence
    file.  Example: \code{Mapping250K_Nsp_probe_tab}}

%% these are from PDInfoPkgSeed-class, bummer for the duplication :-(
    \item{\code{chipName}:}{Name of the chip or platform}
    \item{\code{manufacturer}:}{chip/platform manufacturer}
    \item{\code{genomebuild}:}{The genome build this platform is based
    upon.}
    \item{\code{pdInfoObjectName}:}{Name of the R object that will be
    made available in a generated pdInfo package.}
    \item{\code{geometry}:}{Used to describe the geometry of the platform. }

%% these are from PkgSeed-class, bummer for the duplication :-(
%% these are from PkgSeed-class, bummer for the duplication :-(
    \item{\code{name}:}{The name of the package.}
    \item{\code{version}:}{A character vector giving the version
    number of the package.}
    \item{\code{license}:}{The license of the package}
    \item{\code{author}:}{Author of the package}
    \item{\code{email}:}{An email address to use in the Maintainer field}
    \item{\code{url}:}{URL}
    \item{\code{biocViews}:}{Character data for the biocViews field of
    the DESCRIPTION file}

  }
}
\section{Extends}{
Class \code{"\linkS4class{PDInfoPkgSeed}"}, directly.
Class \code{"\linkS4class{PkgSeed}"}, by class "PDInfoPkgSeed", distance 2.
}
\section{Methods}{
See \code{showMethods(classes="AffySNPPDInfoPkgSeed")}
}

\author{Seth Falcon}
\examples{
showClass("AffySNPPDInfoPkgSeed")

cdfFile <- "Mapping250K_Nsp.cdf"
csvAnno <- "Mapping250K_Nsp_annot.csv"
csvSeq <- "Mapping250K_Nsp_probe_tab"

pkg <- new("AffySNPPDInfoPkgSeed", 
           version="0.1.5",
           author="A. U. Thor", email="au@thor.net",
           biocViews="AnnotationData",
           genomebuild="NCBI Build 35, May 2004",
           cdfFile=cdfFile, csvAnnoFile=csvAnno, csvSeqFile=csvSeq)

showMethods(classes=class(pkg))

}
\keyword{classes}
