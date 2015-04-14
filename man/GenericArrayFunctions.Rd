\name{cdf2table}
\alias{cdf2table}
\alias{sequenceParser}

\title{
Helper functions to assist the creation of an annotation package for a generic array
}
\description{
Helper functions to assist the creation of an annotation package for a generic array. This includes converting CDF files into flat tables and parsing probe sequence files.
}
\usage{
cdf2table(cdfFile)
sequenceParser(seqFile)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdfFile}{name of the CDF file to be used}
  \item{seqFile}{name of the probe sequence file}
}
\details{
\code{cdf2table} will convert a CDF to a flat table.

\code{seqFile} will extract a flat table containing physical location and probe sequences.
}
\keyword{ manip }
