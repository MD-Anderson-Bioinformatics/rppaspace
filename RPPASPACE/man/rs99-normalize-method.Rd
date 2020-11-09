% $Id: rs99-normalize-method.Rd

\name{normalize-method}
\alias{normalize}
\alias{normalize-method}
\alias{normalize,ANY-method}
\alias{normalize,NULL-method}
\docType{methods}
\title{Method \dQuote{normalize}}
\description{
  \code{normalize} is a generic function used to normalize the data based
  on the input object. The method invokes particular \code{\link{methods}}
  which depend on the \code{\link{class}} of the first argument.
}
\usage{
\S4method{normalize}{ANY}(object, \dots)
\S4method{normalize}{NULL}(object, \dots)
}
\arguments{
  \item{object}{an object to be normalized}
  \item{\dots}{additional arguments affecting the normalization process}
}
\value{
  The form of the value returned by \code{normalize} depends on the
  class of its argument. See the documentation of the particular methods
  for details of what is produced by that method.

  If the object is \code{NULL}, \code{NA} is returned.
}
\author{
  P. Roebuck \email{paul_roebuck@comcast.net}
}
\keyword{methods}

