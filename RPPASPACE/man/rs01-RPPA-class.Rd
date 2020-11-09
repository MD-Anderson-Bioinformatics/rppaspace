% $Id: rs01-RPPA-class.Rd

\name{RPPA-class}
\alias{RPPA-class}
\alias{RPPA}
\alias{is.RPPA}
\alias{dim,RPPA-method}
\alias{image,RPPA-method}
\alias{summary,RPPA-method}
\alias{plot,RPPA,ANY-method}
\alias{seriesNames}
\alias{seriesToUseToMakeCurve}
\docType{class}
\title{Class \dQuote{RPPA}}
\description{
  The RPPA class represents the raw quantification data from a
  reverse-phase protein array experiment.
}
\usage{
RPPA(file,
     path=".",
     slideNumber=NA,
     antibody=NULL,
     tracking=NULL,
     seriesToIgnore=NULL,
	 warningsFileName="warnings.txt"
	 )
is.RPPA(x)
\S4method{dim}{RPPA}(x)
\S4method{image}{RPPA}(x, measure="Net.Value", 
     main = .mkPlotTitle(measure,
                 x@antibody),
     colorbar=FALSE, col=terrain.colors(256), \dots)
\S4method{summary}{RPPA}(object, \dots)
seriesNames(rppa)
seriesToUseToMakeCurve(rppa)
}
\arguments{
  \item{file}{character string or connection specifying text file
    containing quantifications of a reverse-phase protein array experiment}
  \item{path}{character string specifying the path from the current
    directory to the file. The default value assumes the file is
    contained in the current directory. If \code{file} is a connection,
    this argument is ignored.}
  \item{antibody}{character string specifying antibody name. If missing,
    default value is filename (referenced by \code{file} argument) without
    extension.}
  \item{slideNumber}{integer containing the index of the slide currently
	being processed.}
  \item{warningsFileName}{character string holding the name of the file
	to which to write out warning messages generated during processing.}
  \item{tracking}{data.frame used to track the points data from a slide and 
    how they are used. (see section \sQuote{Tracking} below)}
  \item{seriesToIgnore}{Comma separated list of series names to ignore. 
    These series will not be used to calculate the curve used to fit data. 
	Names in list must match series names in sample file. 
	}
  \item{object}{object of class \code{RPPA}}
  \item{x}{object of class \code{RPPA}}
  \item{measure}{character string containing the name of the measurement column
    in \code{data} that should be displayed by the \code{image} method}
  \item{main}{character string used to title the image plot}
  \item{colorbar}{logical scalar that determines whether to include a
    color bar in the plot. If \code{TRUE}, the image cannot be used as
    one panel in a window with multiple plots. Default is \code{FALSE}.}
  \item{col}{graphics parameter used by \link[graphics]{image}. }
  \item{\dots}{extra arguments for generic or plotting routines}
  \item{rppa}{object of class \code{RPPA}}
}
\details{
  The data frame slot (\code{data}) in a valid RPPA object constructed
  from a quantification file using the \code{RPPA} generator function
  is guaranteed to contain at least 14 columns of information:
  \tabular{ll}{
    \code{Order} \tab Spot number order in file \cr
    \code{Main.Row} \tab logical location of spot on the array \cr
    \code{Main.Col} \tab logical location of spot on the array \cr
    \code{Sub.Row}  \tab logical location of spot on the array \cr
    \code{Sub.Col}  \tab logical location of spot on the array \cr
    \code{Series.Id}   \tab unique numeric identifier of sample spotted at location \cr
    \code{Spot.Type} \tab type of spot at location \cr
    \code{Dilution} \tab measurement representing background-corrected
                         mean intensity of the spot \cr
    \code{Net.Value} \tab measurement representing background-corrected
                         mean intensity of the spot \cr
    \code{Raw.Value} \tab measurement representing 
                         mean intensity of the spot \cr
    \code{Background.Value} \tab measurement representing mean background
                         intensity of the spot \cr
    \code{Spot.X.Position} \tab X location of spot on graphic image \cr
    \code{Spot.Y.Position} \tab Y location of spot on graphic image \cr
    \code{Original.Order} \tab Spot number order in original input file \cr
  } 
  Taken together, the four components (Main.Row, Sub.Row, Main.Col, Sub.Col) give 
  the logical location of aspot on an array. Additional columns may be included.
}
\section{Tracking}{
  An object for tracking how points in the slide are to be used in the process.
  The information comes from the sample file of the first slide that has a valid layout.
  The layout of all other slides are compared to this and skipped if they don't have an identical layout.
  
  Slots
  \tabular{lllllllllll}{
	\code{spotType}   \tab Spot.Type according to design file. \tab No default \cr
	\code{isNegCtrl}  \tab Is point a Negative Control Point. \tab Default: FALSE \tab TRUE if value of Spot.Type is negative control type (Blank, Buffer, NegCtrl) \cr
	\code{isPosCtrl}  \tab Is point a Positive Control Point. \tab Default: FALSE \tab TRUE if value of Spot.Type is positive control type (PosCtrl, or PosCtrl-Noise)\cr
	\code{isCtrl}     \tab Is point a Control Point.          \tab Default: FALSE \tab TRUE if value of Spot.Type is any control type (Blank, Buffer, NegCtrl, PosCtrl or PosCtrl-Noise)\cr
	\code{applySpatialCorrection}  \tab Apply spatial correction to point  \tab Default: TRUE  \tab FALSE if Spot.Type in design file is any control type, but TRUE if an Noise type. \cr
	\code{makePartOfCurve} \tab Should point be used to create curve to which to fit data? \tab Default: TRUE \tab FALSE if SpotType is control type or noise type (Blank, Buffer, Noise, NegCtrl, PosCtrl or PosCtrl-Noise) or in SeriesToIgnore parameter\cr
	\code{fitToCurve} \tab Should point be fit to curve? \tab Default: TRUE \tab FALSE if point is control in design but not a noise point (Blank, Buffer, NegCtrl, or PosCtrl) \cr
	\code{isNoise}    \tab Should point be used in noise calculations \tab Default: FALSE \tab TRUE if Noise Point (Noise or PosCtrl-Noise) \cr
	\code{isSample}   \tab Is point a sample point? \tab Default: TRUE \tab FALSE if not Sample Point in design file\cr
	\code{badPoint}   \tab Does this point have a value that was not used or caused problems in processing and whos output accuracy should be questioned? \tab Default: FALSE \tab Used as status indicator during run. \cr
	\code{dilution}   \tab Dilution value for this point. (Decimal value from slide). \cr
  }
}
\section{Objects from the Class}{
  Although objects of the class can be created by a direct call to
  \link[methods]{new}, the preferred method is to use the
  \code{RPPA} generator function.
}
\section{Slots}{
  \describe{
    \item{\code{data}}{data.frame containing the contents of a
      quantification file}
    \item{\code{file}}{character string specifying the name of the file that
      the data was loaded from}
    \item{\code{slideNumber}:}{integer containing the index of the slide currently
	  being processed.}
    \item{\code{antibody}}{character string specifying name of antibody}
    \item{\code{tracking}}{data.frame used to track the points data from a slide and 
      how they are used. (see section \sQuote{Tracking} below)}
    \item{\code{seriesToIgnore}}{NULL or Comma separated list of series names to ignore. 
      These series will not be used to calculate the curve used to fit data. 
	  Names in list must match series names in sample file.}
    \item{\code{warningsFileName}}{character string holding the name of the file
	  to which to write out warning messages generated during processing.}
  }
}
\section{Methods}{
  \describe{
    \item{dim}{
      \code{signature(x = "RPPA")}:\cr
      Returns the dimensions of the slide layout.}
    \item{image}{
      \code{signature(x = "RPPA")}:\cr
      Produces a "geographic" image of the measurement column named by
      the \code{measure} argument. The colors in the image represent the
      intensity of the measurement at each spot on the array, and the
      display locations match the row and column locations of the spot.
      Any measurement column can be displayed using this function. An
      optional color bar can be added, placed along the right edge.}
    \item{summary}{
      \code{signature(object = "RPPA")}:\cr
      Prints a summary of the underlying data frame.}
  }
}
\value{
  The \code{RPPA} generator returns an object of class \code{RPPA}.

  The \code{is.RPPA} method returns \code{TRUE} if its
  argument is an object of class \code{RPPA}.

  The \code{dim} method returns a numeric vector of length 4.

  The \code{image} method invisibly returns the \code{RPPA} object on
  which it was invoked.

  The \code{summary} method returns a summary of the underlying data frame.

  The \code{seriesNames} function returns a character vector containing
  the names of the unique (non-control) dilution series on the array.

  The \code{seriesToUseToMakeCurve} function returns a character vector containing
  the names of the unique (non-control) dilution series on the array that are
  used to create a curve to fit samples to.
}

\author{
  Kevin R. Coombes \email{coombes.3@osu.edu},
  P. Roebuck \email{paul_roebuck@comcast.net},
  James M. Melott \email{jmmelott@mdanderson.org}
}
\seealso{
  \code{\linkS4class{RPPADesignParams}},
  \code{\linkS4class{RPPAFit}}
}
\keyword{color}
\keyword{classes}
\keyword{file}
\keyword{hplot}

