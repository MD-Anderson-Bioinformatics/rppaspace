% $Id: rs08-LogisticFitClass-class.Rd

\name{LogisticFitClass-class}
\alias{LogisticFitClass-class}
\alias{coef,LogisticFitClass-method}
\alias{coefficients,LogisticFitClass-method}
\alias{fitSeries,LogisticFitClass-method}
\alias{fitSlide,LogisticFitClass-method}
\alias{fitted,LogisticFitClass-method}
\alias{trimConc,LogisticFitClass-method}
\docType{class}
\title{Class \dQuote{LogisticFitClass}}
\description{
  The \code{LogisticFitClass} class represents models that were fit with the
  logistic model.
}
\usage{
\S4method{coef}{LogisticFitClass}(object, \dots)
\S4method{coefficients}{LogisticFitClass}(object, \dots)
\S4method{fitSeries}{LogisticFitClass}(object,
          diln,
          intensity,
          est.conc,
          method="nls",
          silent=TRUE,
          trace=FALSE,
          \dots)
\S4method{fitSlide}{LogisticFitClass}(object,
         conc,
         intensity,
         \dots)
\S4method{fitted}{LogisticFitClass}(object,
       conc,
       \dots)
\S4method{trimConc}{LogisticFitClass}(object,
         conc,
         intensity,
         design,
         trimLevel,
         \dots)
} 
\arguments{
  \item{object}{object of class \code{LogisticFitClass}}
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
    \item{\code{coefficients}:}{numeric vector of length 3, representing
      alpha, beta, and gamma respectively.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{FitClass}}, directly.
}
\section{Methods}{
  \describe{
    \item{coef}{
      \code{signature(object = "LogisticFitClass")}:\cr
      Extracts model coefficients from objects returned by modeling functions.}
    \item{coefficients}{
      \code{signature(object = "LogisticFitClass")}:\cr
      An alias for \code{coef}}
    \item{fitSeries}{
      \code{signature(object = "LogisticFitClass")}:\cr
      Finds the concentration for an individual dilution series given the
      curve fit for the slide.}
    \item{fitSlide}{
      \code{signature(object = "LogisticFitClass")}:\cr
      Uses the concentration and intensity series for an entire slide to
      fit a curve for the slide of intensity = f(conc).}
    \item{fitted}{
      \code{signature(object = "LogisticFitClass")}:\cr
      Extracts fitted values of the model.}
    \item{trimConc}{
      \code{signature(object = "LogisticFitClass")}:\cr
      Returns concentration and intensity cutoffs for the model.}
  }
}
\value{
  The \code{coef} and \code{coefficients} methods return a named vector
  of length three with logistic curve coefficients.

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

