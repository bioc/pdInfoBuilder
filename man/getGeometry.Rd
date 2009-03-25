\name{getGeometry}
\alias{getGeometry}
\alias{getGeometry,AffySNPPDInfoPkgSeed-method}
\alias{getGeometry,AffySNPCNVPDInfoPkgSeed-method}
\alias{getGeometry,AffyExpressionPDInfoPkgSeed-method}
\alias{getGeometry,AffyTilingPDInfoPkgSeed-method}
\alias{getGeometry,AffySTPDInfoPkgSeed-method}
\alias{getGeometry,NimbleGenPDInfoPkgSeed-method}

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
 A list with two elements nrows and ncols
}

\author{Matt Settles}

\keyword{methods}

