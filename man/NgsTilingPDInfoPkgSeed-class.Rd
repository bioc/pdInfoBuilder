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
    \item{\code{pairFile}:}{Pair File - used as template}
    \item{\code{posFile}:}{POS (Positions) file}
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
Class \code{"\linkS4class{NgsTilingPDInfoPkgSeed}"}, directly
Class \code{"\linkS4class{PDInfoPkgSeed}"}, distance 2.
Class \code{"\linkS4class{PkgSeed}"}, by class "PDInfoPkgSeed", distance 3.
}
\section{Methods}{
  \describe{
    \item{makePdInfoPackage}{\code{signature(.Object = "NgsTilingPDInfoPkgSeed")}: ... }
    \item{chipName}{\code{signature(object = "NgsPDInfoPkgSeed")}: ... }
    \item{getGeometry}{\code{signature(object = "NgsPDInfoPkgSeed")}: ... }
	 }
}
\section{Note}{
	Nimblegen provides pairFiles for each array and color as the standard, however this file type
	contains a lot of redundant information that can be found in the ndf file. The xysFile
	is a reduced pairFile that can be obtained from Nimblegen upon request and provides a significant
	amount of memory savings. Therefore only one of xysFile or pairFile must be specified.
}
\section{seealso}
See
\examples{
showClass("NgsTilingPDInfoPkgSeed")
}
\keyword{classes}
