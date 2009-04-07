\name{NgsTilingPDInfoPkgSeed-class}
\docType{class}
\alias{NgsTilingPDInfoPkgSeed-class}

\title{Class "NgsTilingPDInfoPkgSeed"}
\description{PDInfo package Seed for NimbleGen Tiling arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NgsTilingPDInfoPkgSeed", ndfFile, xysFile, pairFile, posFile ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{ndfFile}:}{NDF (NimbleGen Design) file}
    \item{\code{xysFile}:}{XYS File - used as template}
    \item{\code{posFile}:}{POS (Positions) file}
    %% these are from PDInfoPkgSeed-class
    \item{\code{chipName}:}{Name of the chip or platform}
    \item{\code{manufacturer}:}{chip/platform manufacturer}
    \item{\code{url}:}{chip URL}
    \item{\code{genomebuild}:}{The genome build this platform is based
      upon.}
    \item{\code{organism}:}{ organism for chip. }
    \item{\code{species}:}{ species for chip. }  
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
    \item{makePdInfoPackage}{\code{signature(.Object = "NgsTilingPDInfoPkgSeed")}: ... }
    \item{chipName}{\code{signature(object = "NimbleGenPDInfoPkgSeed")}: ... }
    \item{getGeometry}{\code{signature(object = "NimbleGenPDInfoPkgSeed")}: ... }
	 }
}
\examples{
showClass("NgsTilingPDInfoPkgSeed")
}
\keyword{classes}
