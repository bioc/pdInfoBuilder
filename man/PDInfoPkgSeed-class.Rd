\name{PDInfoPkgSeed-class}
\docType{class}
\alias{PDInfoPkgSeed-class}

\title{Class "PDInfoPkgSeed"}
\description{
  This class represents platform design (PD) packages.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PDInfoPkgSeed", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{chipName}:}{Name of the chip or platform}
    \item{\code{manufacturer}:}{chip/platform manufacturer}
    \item{\code{genomebuild}:}{The genome build this platform is based
    upon.}
    \item{\code{pdInfoObjectName}:}{Name of the R object that will be
    made available in a generated pdInfo package.}
    \item{\code{geometry}:}{Used to describe the geometry of the platform. }

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
Class \code{"\linkS4class{PkgSeed}"}, directly.
}
\section{Methods}{
No methods defined with class "PDInfoPkgSeed" in the signature.
}

\author{Seth Falcon}
\examples{
showClass("PDInfoPkgSeed")
}
\keyword{classes}
