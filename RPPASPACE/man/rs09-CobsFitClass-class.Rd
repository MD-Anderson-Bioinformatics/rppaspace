% $Id: rs09-CobsFitClass-class.Rd

\name{CobsFitClass-class}
\alias{CobsFitClass-class}
\alias{fitSeries,CobsFitClass-method}
\alias{fitSlide,CobsFitClass-method}
\alias{fitted,CobsFitClass-method}
\alias{trimConc,CobsFitClass-method}
\docType{class}
\title{Class \dQuote{CobsFitClass}}
\description{
  The \code{CobsFitClass} class represents models that were fit with the
  nonparametric model.
}
\usage{
\S4method{fitSeries}{CobsFitClass}(object,
          diln,
          intensity,
          est.conc,
          method="nls",
          silent=TRUE,
          trace=FALSE,
          \dots)
\S4method{fitSlide}{CobsFitClass}(object,
         conc,
         intensity,
         \dots)
\S4method{fitted}{CobsFitClass}(object, conc, \dots)
\S4method{trimConc}{CobsFitClass}(object,
         conc,
         intensity,
         design,
         trimLevel,
         \dots)
}
\arguments{
  \item{object}{object of class \code{CobsFitClass}}
  \item{diln}{numeric vector of dilutions for series to be fit}
  \item{intensity}{numeric vector of observed intensities for series to be fit}
  \item{est.conc}{numeric estimated concentration for EC50 dilution}
  \item{method}{character string specifying regression method to use to fit the
    series}
  \item{silent}{logical scalar. If \code{TRUE}, report of error messages will
    be suppressed in \code{try(nlsmeth(...))}}
  \item{trace}{logical scalar. Used in \code{nls} method.}
  \item{conc}{numeric vector containing estimates of the log concentration
    for each dilution series}
  \item{design}{object of class \code{RPPADesignParams} describing options for processing
    the array}
  \item{trimLevel}{numeric scalar multiplied to Median Absolute Deviation \acronym{MAD}}
  \item{\dots}{extra arguments for generic routines}
}
\section{Objects from the Class}{
  Objects are created internally by calls to the methods
  \code{\link[=FitClass-class]{fitSlide}} or \code{\link{RPPAFit}}.
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{object of class \code{cobs} summarizing nonparametric
      fit}
    \item{\code{lambda}:}{numeric}
  }
}
\section{Extends}{
  Class \code{\linkS4class{FitClass}}, directly.
}
\section{Methods}{
  \describe{
    \item{fitSeries}{
      \code{signature(object = "CobsFitClass")}:\cr
      Finds the concentration for an individual dilution series given the
      curve fit for the slide.}
    \item{fitSlide}{
      \code{signature(object = "CobsFitClass")}:\cr
      Uses the concentration and intensity series for an entire slide to
      fit a curve for the slide of intensity = f(conc).}
    \item{fitted}{
      \code{signature(object = "CobsFitClass")}:\cr
      Extracts fitted values of the model.}
    \item{trimConc}{
      \code{signature(object = "CobsFitClass")}:\cr
      Returns concentration and intensity cutoffs for the model.}
  }
}
\value{
  The \code{fitted} method returns a numeric vector.
}
\author{
  P. Roebuck \email{paul_roebuck@comcast.net},
  James M. Melott \email{jmmelott@mdanderson.org}
}
\seealso{
  \code{\linkS4class{FitClass}}
}
\keyword{classes}
\keyword{models}

