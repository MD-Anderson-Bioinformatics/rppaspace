% $Id: rs99-paramString-method.Rd

\name{paramString-method}
\alias{paramString}
\alias{paramString-method}
\alias{paramString,ANY-method}
\docType{methods}
\title{Method \dQuote{paramString}}
\description{
  \code{paramString} is a generic function used to produce a string
  representation of an object's slots. The method invokes particular
  \code{\link{methods}} which depend on the \code{\link{class}} of the
  first argument.
}
\usage{
\S4method{paramString}{ANY}(object, \dots)
}
\arguments{
  \item{object}{an object for which a string representation is desired}
  \item{\dots}{additional arguments affecting the string representation
    produced}
}
\value{
  The form of the value returned by \code{paramString} depends on the
  class of its argument. See the documentation of the particular methods
  for details of what is produced by that method.
  \cr
  The \code{paramString} method returns a character vector, possibly
  empty but never \code{NULL}.
}
\section{Warning}{
  The \code{paramString} method should not be called by user except for
  informational purposes. The content and format of the returned string
  may vary between different versions of this package.
}
\author{
  P. Roebuck \email{paul_roebuck@comcast.net},
  James M. Melott \email{jmmelott@mdanderson.org}
}
\keyword{methods}
\keyword{internal}

