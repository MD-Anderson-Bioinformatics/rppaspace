% $Id: rs20-RPPASpatialParams-class.Rd

\name{RPPASpatialParams-class}
\alias{RPPASpatialParams-class}
\alias{RPPASpatialParams}
\alias{is.RPPASpatialParams}
\alias{paramString,RPPASpatialParams-method}
\docType{class}
\title{Class \dQuote{RPPASpatialParams}}
\description{
  The RPPASpatialParams class is used to bundle the parameter set together
  that control how to perform spatial adjustment into a reusable object.
}
\usage{
RPPASpatialParams(cutoff=0.8,
                  k=100,
                  gamma=0.1,
                  plotSurface=FALSE)
is.RPPASpatialParams(x)
\S4method{paramString}{RPPASpatialParams}(object, slots, \dots)
}
\arguments{
  \item{cutoff}{numeric scalar used to identify the background cutoff
    with value in closed interval [0..1]. Default is \code{0.8}.}
  \item{k}{numeric scalar used as smoothing model argument.
    Default is \code{100}.}
  \item{gamma}{numeric scalar used as model parameter with value in
    closed interval [0..2]. Default is \code{0.1}.}
  \item{plotSurface}{logical scalar. If \code{TRUE}, plots surfaces.
    Default is \code{FALSE}.}
  \item{object}{object of class \code{RPPASpatialParams}}
  \item{x}{object of class \code{RPPASpatialParams}}
  \item{slots}{strings specifying \code{RPPASpatialParams} slotnames to display
    (for debugging)}
  \item{\dots}{extra arguments for generic routines}
}
\details{
  The \code{cutoff} argument passed to \code{quantile} is percentile
  of the background estimates used to define the noise region of slide.

  The \code{k} argument passed to \code{s} sets upper limit on
  degrees of freedom associated with smoothing.

  The \code{gamma} argument passed to \code{gam} provides a constant
  multiplier used to inflate model degrees of freedom in the
  \acronym{GCV} or \acronym{UBRE}/\acronym{AIC} score.
}
\section{Objects from the Class}{
  Although objects of the class can be created by a direct call to
  \link[methods]{new}, the preferred method is to use the
  \code{RPPASpatialParams} generator function.
}
\section{Slots}{
  \describe{
    \item{\code{cutoff}:}{numeric scalar; see arguments above}
    \item{\code{k}:}{numeric scalar; see arguments above}
    \item{\code{gamma}:}{numeric scalar; see arguments above}
    \item{\code{plotSurface}:}{logical scalar; see arguments above}
  }
}
\section{Methods}{
  \describe{
    \item{paramString(object)}{Returns string representation of object.}
  }
}
\value{
  The \code{RPPASpatialParams} generator returns an object of class
  \code{RPPASpatialParams}.

  The \code{is.RPPASpatialParams} method returns \code{TRUE} if its
  argument is an object of class \code{RPPASpatialParams}.

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
\seealso{
  \code{\link{spatialCorrection}}
}
\examples{
showClass("RPPASpatialParams")
spatialparams <- RPPASpatialParams(cutoff=0.8,
								   k=100,
								   gamma=0.1,
								   plotSurface=FALSE)
paramString(spatialparams)
}
\keyword{classes}

