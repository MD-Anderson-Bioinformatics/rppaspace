% $Id: rs98-Directory-class.Rd

\name{Directory-class}
\alias{Directory-class}
\alias{Directory}
\alias{coerce,character,Directory-method}
\alias{coerce,Directory,character-method}
\alias{is.Directory}
\docType{class}
\title{Class \dQuote{Directory}}
\description{
  The Directory class represents a file system directory.
}
\usage{
Directory(path)
is.Directory(x)
\S4method{coerce}{character,Directory}(from, to, strict=TRUE)
\S4method{coerce}{Directory,character}(from, to, strict=TRUE)
}
\arguments{
  \item{path}{character string specifying a directory}
  \item{x}{object of class \code{Directory}}
  \item{from}{object of class \code{Directory} or character string specifying
    pathname of directory}
  \item{to}{object of class \code{Directory} or character string specifying
    pathname of directory}
  \item{strict}{logical scalar. If \code{TRUE}, the returned object must be
    strictly from the target class; otherwise, any simple extension of the
    target class will be returned, without further change.}
}
\section{Objects from the Class}{
  Although objects of the class can be created by a direct call to
  \link[methods]{new}, the preferred method is to use the \code{Directory}
  generator function.
}
\section{Slots}{
  \describe{
    \item{\code{path}:}{character string specifying a directory}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{
      \code{signature(from = "Directory", to = "character")}:\cr
      Coerce an object of class \code{Directory} to its character string
      pathname equivalent.}
    \item{coerce}{
      \code{signature(from = "character", to = "Directory")}:\cr
      Coerce a character string specifying directory pathname to an
      equivalent object of class \code{Directory}.}
  }
}
\value{
  The \code{Directory} generator returns an object of class \code{Directory}.

  The \code{is.Directory} method returns \code{TRUE} if its
  argument is an object of class \code{Directory}.
}
\note{
  The coercion methods should not be called explicitly; instead, use an
  explicit call to the \code{as} method.
}
\author{
  P. Roebuck \email{paul_roebuck@comcast.net},
  James M. Melott \email{jmmelott@mdanderson.org}
}
\seealso{
  \code{\link[methods]{as}}
}
\examples{
  pkgdir <- Directory(system.file(package="RPPASPACE"))
  pkgdir.path <- as(pkgdir, "character")
}
\keyword{classes}
\keyword{file}
\keyword{methods}

