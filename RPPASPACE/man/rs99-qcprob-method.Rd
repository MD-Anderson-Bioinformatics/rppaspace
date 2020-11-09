% $Id: rs99-qcprob-method.Rd

\name{qcprob-method}
\alias{qcprob}
\alias{qcprob-method}
\alias{qcprob,ANY-method}
\alias{qcprob,NULL-method}
\docType{methods}
\title{Method \dQuote{qcprob}}
\description{
  \code{qcprob} is a generic function used to produce a quality control
  probability based on the input object. The method invokes particular
  \code{\link{methods}} which depend on the \code{\link{class}} of the
  first argument.
}
\usage{
\S4method{qcprob}{ANY}(object, \dots)
\S4method{qcprob}{NULL}(object, \dots)
}
\arguments{
  \item{object}{an object for which a QC probability is desired}
  \item{\dots}{additional arguments affecting the QC probability produced}
}
\value{
  The form of the value returned by \code{qcprob} depends on the
  class of its argument. See the documentation of the particular methods
  for details of what is produced by that method.

  If the object is \code{NULL}, \code{NA} is returned.
}
\author{
  P. Roebuck \email{paul_roebuck@comcast.net},
  James M. Melott \email{jmmelott@mdanderson.org}
}
\keyword{methods}

