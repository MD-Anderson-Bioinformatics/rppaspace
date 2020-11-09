% $Id: rs06-RPPAFitParams-class.Rd

\name{RPPAFitParams-class}
\alias{RPPAFitParams-class}
\alias{RPPAFitParams}
\alias{RPPAFitFromParams}
\alias{RPPAFit}
\alias{is.RPPAFit}
\alias{is.RPPAFitParams}
\alias{paramString,RPPAFitParams-method}
\docType{class}
\title{Fitting Dilution Curves to Protein Lysate Arrays with Class \dQuote{RPPAFitParams}}
\description{
  The \code{RPPAFit} function fits an intensity response model to the
  dilution series in a reverse-phase protein array experiment. Individual
  sample concentrations are estimated by first matching individual sample
  dilution series to the overall logistic response for the slide and then 
  fitting a second time using the specified model, usually cobs.
  The \code{RPPAFitParams} class is a convenient place to wrap the parameters
  that control the model fit into a reusable object.
}
\usage{
RPPAFit(rppa,
        measure,
        model="logistic",
        xform=NULL,
        method=c("nls", "nlrob", "nlrq"),
        trim=2,
        ci=FALSE,
        ignoreNegative=TRUE,
        trace=FALSE,
        verbose=FALSE,
        veryVerbose=FALSE,
        warnLevel=0)

RPPAFitParams(measure,
              model="logistic",
              xform=NULL,
              method=c("nls", "nlrob", "nlrq"),
              trim=2,
              ci=FALSE,
              ignoreNegative=TRUE,
              trace=FALSE,
              verbose=FALSE,
              veryVerbose=FALSE,
              warnLevel=0)

RPPAFitFromParams(rppa,
                  fitparams,
                  progmethod=NULL)

is.RPPAFit(x)
is.RPPAFitParams(x)
\S4method{paramString}{RPPAFitParams}(object, slots, \dots)
}
\arguments{
  \item{rppa}{object of class \code{\link{RPPA}} containing the raw data
    to be fit}
  \item{fitparams}{object of the class \code{RPPAFitParams}, bundling
    together the following arguments.}
  \item{progmethod}{user defined function that will take a string telling which 
	portion	of the process is running and do with it as the function specifies.
	Default is a function that does nothing.}
  \item{measure}{character string identifying the column of the raw RPPA data
    that should be used to fit to the model.}
  \item{model}{character string specifying the model for the response curve
    fitted for the slide. Valid values are:
    \tabular{ll}{
      \code{"logistic"} \tab assumes a logistic shape for the curve \cr
      \code{"loess"}    \tab fits a loess curve to the response \cr
      \code{"cobs"}     \tab fits a b-spline curve to the slide with the
                             constraint that curve be strictly increasing \cr
    }}
  \item{xform}{optional function that takes a single input vector and
    returns a single output vector of the same length. The \code{measure}
    column is transformed using this function before fitting the model.}
  \item{method}{character string specifying the method for matching
    the individual dilution series to the response curve fitted for the slide.
    Valid values are:
    \tabular{ll}{
      \code{"nls"}   \tab uses the optimal fit based on nonlinear least
                          squares \cr
      \code{"nlrob"} \tab uses \code{nlrob} which is robust nls
                          from \pkg{robustbase} package \cr
      \code{"nlrq"}  \tab uses \code{nlrq} which is robust median regression
                          from \pkg{quantreg} package \cr
    }}
  \item{trim}{numeric or logical scalar specifying trim level for
    concentrations. If positive, concentrations will be trimmed to reflect
    min and max concentrations we can estimate given the background noise.
    If \code{TRUE}, the trim level defaults to \code{2}, which was originally
    the hardcoded value; otherwise, raw concentrations are returned without
    trimming.}
  \item{ci}{logical scalar. If \code{TRUE}, computes 90\% confidence intervals
    on the log concentration estimates.}
  \item{ignoreNegative}{logical scalar. If \code{TRUE}, converts negative values
    to \code{NA} before fitting the model.}
  \item{trace}{logical scalar passed to \link[stats]{nls} in the \code{method}
    portion of the routine}
  \item{verbose}{logical scalar. If \code{TRUE}, prints updates while fitting
    the data}
  \item{veryVerbose}{logical scalar. If \code{TRUE}, prints voluminous updates
    as each individual dilution series is fitted}
  \item{warnLevel}{integer scalar used to set the \code{warn} option before
    calling \code{method}. Since this is wrapped in a \code{try} function, it
    won't cause failure but will give us a chance to figure out which dilution
    series are failing. Setting \code{warnLevel} to two or greater may change
    the values returned.}
  \item{object}{object of class \code{RPPAFitParams}}
  \item{x}{object of class \code{RPPAFit} (or \code{RPPAFitParams})}
  \item{slots}{strings specifying \code{RPPAFitParams} slotnames to display
    (for debugging)}
  \item{\dots}{extra arguments for generic routines.}
}
\details{
  The basic mathematical model is given by \deqn{Y = f(X-\delta_i),}
  where \eqn{Y} is the observed intensity, \eqn{X} is the designed dilution
  step and \eqn{f} is the model for the protein response curve.
  By fitting a joint model, we assume that the response curve is the same for
  all dilution series on the array. The real point of the model, however, is
  to be able to draw inferences on the \eqn{\delta_i}, which represent the
  (log) concentration of the protein present in different dilution series.

  As the first step in fitting the model, we compute crude estimates of the
  individual \eqn{\delta_i} assuming a rough logistic shape for the protein
  response curve.

  Next, we fit an overall response curve for the slide \eqn{f} using the
  estimated concentrations and observed intensities \eqn{Y = f(\delta_i)}.
  The model for \eqn{f} is specified in the \eqn{model} parameter.

  Next, we update the estimates of the individual \eqn{\delta_i} using our
  improved fitted model \eqn{f} for the overall slide response curve. These
  individual series are matched to the overall slide response curve using the
  algorithm specified in \code{method}. The default method is \code{nls}, a
  least squares match-up, but we also offer robust alternatives which can do
  better.

  Finally, we re-estimate \eqn{f} using the improved estimates for
  \eqn{\delta_i}. We continue to iterate between \eqn{f} and \eqn{\delta_i}.
  We do this twice since that seems to give reasonable convergence.

  If the \code{ci} argument is \code{TRUE}, then the function also computes
  confidence intervals around the estimates of the log concentration.
  Since this step can be time-consuming, it is not performed by default.
  Moreover, confidence intervals can be computed after the main model is fit
  and evaluated, using the \code{\link{getConfidenceInterval}} function.
}
\section{Objects from the Class}{
  Although objects of the class can be created by a direct call to
  \link[methods]{new}, the preferred method is to use the \code{RPPAFitParams}
  function.
}
\section{Slots}{
  \describe{
    \item{\code{measure}:}{character; see arguments above}
    \item{\code{xform}:}{function or \code{NULL}; see arguments above}
    \item{\code{method}:}{character; see arguments above}
    \item{\code{ci}:}{logical scalar; see arguments above}
    \item{\code{ignoreNegative}:}{logical scalar; see arguments above}
    \item{\code{trace}:}{logical scalar; see arguments above}
    \item{\code{verbose}:}{logical scalar; see arguments above}
    \item{\code{veryVerbose}:}{logical scalar; see arguments above}
    \item{\code{warnLevel}:}{numeric; see arguments above}
    \item{\code{trim}:}{numeric; see arguments above}
    \item{\code{model}:}{character; see arguments above}
  }
}
\section{Methods}{
  \describe{
    \item{paramString}{
      \code{signature(object = "RPPAFitParams")}:\cr
      Returns string representation of object.}
  }
}
\value{
  The \code{RPPAFit} generator and \code{RPPAFitFromParams} function return
  an object of class \code{\link{RPPAFit}}.

  The \code{RPPAFitParams} generator returns an object of class
  \code{RPPAFitParams}.

  The \code{is.RPPAFit} method returns \code{TRUE} if its
  argument is an object of class \code{RPPAFit}.

  The \code{is.RPPAFitParams} method returns \code{TRUE} if its
  argument is an object of class \code{RPPAFitParams}.

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
  Kevin R. Coombes \email{coombes.3@osu.edu},
  James M. Melott \email{jmmelott@mdanderson.org}
}
\seealso{
  \code{\linkS4class{RPPAFit}},
  \code{\link{RPPAFit-class}},
  \code{\link{RPPA}},
  \code{\link{RPPADesignParams}}
}

\keyword{classes}
\keyword{models}
\keyword{nonlinear}
\keyword{regression}
\keyword{robust}

