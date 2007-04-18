\name{PkgSeed-class}
\docType{class}
\alias{PkgSeed-class}

\title{Class "PkgSeed"}
\description{
  This class represents R package details.  It is intended to provide
  a base upon which package generating methods can be built.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PkgSeed", ...)}.
}
\section{Slots}{
	 \describe{
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
\section{Methods}{
No methods defined with class "PkgSeed" in the signature.
}

\author{Seth Falcon}
\examples{
showClass("PkgSeed")
}
\keyword{classes}
