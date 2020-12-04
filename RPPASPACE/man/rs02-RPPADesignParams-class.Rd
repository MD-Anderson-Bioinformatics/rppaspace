% $Id: rs02-RPPADesign-class.Rd

\name{RPPADesignParams-class}
\alias{RPPADesignParams-class}
\alias{RPPADesignParams}
\alias{is.RPPADesignParams}
\alias{paramString,RPPADesignParams-method}
\alias{plot,RPPADesignParams-method}
\alias{plot,RPPA-method}
\docType{class}
\title{Class \dQuote{RPPADesignParams}}
\description{
  The \code{RPPADesignParams} class is used to bundle the design parameter 
  set together for easier re-use.
}
\usage{

RPPADesignParams(center = FALSE, seriesToIgnore = NULL )

is.RPPADesignParams(x)
\S4method{paramString}{RPPADesignParams}(object, slots, \dots)
\S4method{plot}{RPPA}(x, measure, main, \dots)
}
\arguments{
  \item{center}{logical scalar. If \code{TRUE}, then dilution steps are
    centered around \code{0}.}
  \item{x}{object of class \code{RPPADesignParams} (or \code{RPPA} in plot method)}
  \item{seriesToIgnore}{object of class \code{list} or \code{NULL} where list members 
	are numeric Series.Id values to ignore when fitting the data to a curve.}
  \item{object}{object of class \code{RPPADesignParams} in \code{paramString} method}
  \item{slots}{strings specifying \code{RPPADesignParams} slotnames to display
    (for debugging)}
  \item{main}{overall title for plot}
  \item{measure}{character string specifying measure to plot}
  \item{\dots}{extra arguments for generic or plotting routines}
}
\details{
  Allows control of some specific controls for how RPPA slides are processed.
  
}
\section{Objects from the Class}{
  Although objects of these classes can be created by a direct call to
  \link[methods]{new}, the preferred method is to start with the
  \code{RPPADesignParams} generator, followed by the
  \code{RPPADesignFromParams} function to construct the final object
  (the \code{RPPADesign} generator is directly implemented in this way).
}
\section{Slots}{
  For \code{RPPADesignParams} class:

  \describe{
    \item{\code{center}:}{see corresponding argument above}
    \item{\code{seriesToIgnore}:}{see corresponding argument above}
  }
}
\section{Methods}{
  \describe{
    \item{paramString}{
      \code{signature(object = "RPPADesignParams")}:\cr
      Returns string representation of object.}
  }
}
\value{
  The \code{RPPADesignParams} generator returns an object of class
  \code{RPPADesignParams}.

  The \code{is.RPPADesignParams} method returns \code{TRUE} if its
  argument is an object of class \code{RPPADesignParams}.

  The \code{paramString} method returns a character vector, possibly
  empty but never \code{NULL}.
}
\section{Warning}{
  The \code{paramString} method should not be called by user except for
  informational purposes. The content and format of the returned string
  may vary between different versions of this package.
}
\author{
  Kevin R. Coombes \email{coombes.3@osu.edu},
  P. Roebuck \email{paul_roebuck@comcast.net},
  James M. Melott \email{jmmelott@mdanderson.org}
}
\seealso{
  \code{\linkS4class{RPPA}}
}
\examples{
showClass("RPPADesignParams")
designparams <- designparams <- RPPADesignParams(center=FALSE,
								 seriesToIgnore=list()
								 )
paramString(designparams)
}

\keyword{classes}
\keyword{regression}

