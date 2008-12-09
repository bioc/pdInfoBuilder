\name{chipName}
\alias{chipName}
\alias{chipName,AffySNPPDInfoPkgSeed-method}
\alias{chipName,AffySNPCNVPDInfoPkgSeed-method}
\alias{chipName,AffyExpressionPDInfoPkgSeed-method}
\alias{chipName,AffyTilingPDInfoPkgSeed-method}
\alias{chipName,AffySTPDInfoPkgSeed-method}
\alias{chipName,NgsPDInfoPkgSeed-method}

\title{Return an Official Chip/Platform Name}
\description{
  This generic function returns an official or standard chip/platform name.
}
\usage{
chipName(object)
}

\arguments{
  \item{object}{See \code{showMethods("chipName")}, but generally
  \code{object} will be a subclass of \code{PkgSeed}.}
}
\details{
  The idea is that the input files can be used to determine a standard
  name for each platform.  For example, the method for
  \code{AffySNPPDInfoPkgSeed} objects reads the header of the CDF file
  to extract a name.
}
\value{
  A character vector of length one giving a standard name for the
  platform.
}

\author{Seth Falcon}

\keyword{methods}

