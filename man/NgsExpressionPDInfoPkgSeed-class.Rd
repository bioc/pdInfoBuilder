\name{NgsExpressionPDInfoPkgSeed-class}
\docType{class}
\alias{NgsExpressionPDInfoPkgSeed-class}

\title{Class "NgsExpressionPDInfoPkgSeed"}
\description{PDInfo package Seed for NimbleGen Expression arrays}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NgsExpressionPDInfoPkgSeed", ndfFile, pairFile, xysFile, ngdFile ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{ndfFile}:}{NDF (NimbleGen Design) file}
    \item{\code{xysFile}:}{XYS File - used as template}
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
    \item{makePdInfoPackage}{\code{signature(.Object = "NgsExpressionPDInfoPkgSeed")}: ... }
    \item{chipName}{\code{signature(object = "NimbleGenPDInfoPkgSeed")}: ... }
    \item{getGeometry}{\code{signature(.Object = "NimbleGenPDInfoPkgSeed")}: ... }
  }
}
\examples{
showClass("NgsExpressionPDInfoPkgSeed")
}
\keyword{classes}
