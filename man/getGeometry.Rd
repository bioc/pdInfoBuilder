\name{chipName}
\alias{getGeometry}
\alias{getGeometry,AffySNPPDInfoPkgSeed-method}
\alias{getGeometry,NgsPDInfoPkgSeed-method}
\alias{getGeometry,NgsExpressionPDInfoPkgSeed-method}
\alias{getGeometry,NgsTilingPDInfoPkgSeed-method}

\title{Return the Chip/Platform geometry}
\description{
  This generic function returns the geometry for a chip/platform.
}
\usage{
getGeometry(object)
}

\arguments{
  \item{object}{See \code{showMethods("getGeometry")}, but generally
  \code{object} will be a subclass of \code{PkgSeed}.}
}
\details{
  The idea is that the input files can be used to determine the
  geometry for each platform.  For example, the method for
  \code{AffySNPPDInfoPkgSeed} objects reads the header of the CDF file
  to extract the geometry.
}
\value{
  A character vector of length one giving the "rows;cols"
}

\author{Matt Settles}

\keyword{methods}

