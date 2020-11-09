% $Id: rs14-spatialCorrection.Rd

\name{spatialCorrection}
\alias{spatialCorrection}
\alias{spatialAdjustment}
\alias{spatialAdjustmentFromParams}
\title{Spatial Correction}
\description{
  This function estimates a smoothed surface from positive control spots on
  an RPPA slide. The surface is used to perform spatial corrections (i.e.,
  because of uneven hybridization) on the array.
  It is used before \code{RPPAFit}, one slide at a time.
}
\usage{
spatialAdjustmentFromParams(rppa,
                            spatialparams)
spatialAdjustment(rppa,
                  cutoff=0.8,
                  k=100,
                  gamma=0.1,
                  plotSurface=FALSE)
spatialCorrection(rppa,
                  measure=c("Net.Value", "Raw.Value"),
                  cutoff=0.8,
                  k=100,
                  gamma=0.1,
                  plotSurface=FALSE)
}
\arguments{
  \item{rppa}{object of class \code{RPPA}}
  \item{spatialparams}{object of class \code{RPPASpatialParams} containing
    parameters used to perform spatial adjustment}
  \item{measure}{character string specifying fit measure to smooth}
  \item{cutoff}{numeric scalar used to identify the background cutoff
    with value in range [0..1]}
  \item{k}{numeric scalar used as smoothing model argument.}
  \item{gamma}{numeric scalar used as model parameter with value in
   range [0..2]}
  \item{plotSurface}{logical scalar. If \code{TRUE}, plots surfaces.}
}
\details{
  The observed spot intensities are assumed to be a combination of true
  signal, background noise, and hybridization effects according to the
  following model:

       \deqn{Y_rc = Y * H_rc + B_rc}

  where \eqn{Y_rc} is the observed intensity, \eqn{Y} is the true signal,
  \eqn{H_rc} is the effect of hybridization, and \eqn{B_rc} is the
  background noise. The subscripts "r" and "c" refer to the physical row
  and column of the spot on the array. Background noise is estimated
  locally by the array software. The hybridization effect is estimated
  fitting a generalized additive model (\acronym{GAM}) to positive control
  spots printed uniformly across the array.

  The estimated surface is used to scale the intensities on the array.
  Each intensity is adjusted by the amount that is needed to make the
  positive control surface flat at the value of the median of the surface.
  This is done by dividing each spot by the estimated surface value and
  then multiplying by the median of the surface.

  Positive control spots that are expressed below the cutoff for the
  noise region are excluded from the computation of the surface.

  Sometimes, positive control spots are printed in a dilution series to
  avoid saturation problems with these spots. When this happens, the
  observed intensities are adjusted by the positive control surface that
  has the most similar expression level.

  The \code{cutoff} argument passed to \code{quantile} is percentile
  of the background estimates used to define the noise region of slide.

  The \code{k} argument passed to \code{s} sets upper limit on
  degrees of freedom associated with smoothing.

  The \code{gamma} argument passed to \code{gam} provides a constant
  multiplier used to inflate model degrees of freedom in the
  \acronym{GCV} or \acronym{UBRE}/\acronym{AIC} score.
}
\value{
  Returns modified \code{rppa} with an additional measurement column
  named after the \code{measure} with an \code{Adj.} prefix. For example,
  if the measure was \code{Net.Value}, the name of the adjusted column
  would be \code{Adj.Net.Value}.
}
\references{
  Neeley ES, Baggerly KA, Kornblau SM.\cr
  \emph{Surface Adjustment of Reverse Phase Protein Arrays Using Positive
  Control Spots}\cr
  Cancer Informatics (2012) 11: 77-86.\cr
  \url{http://www.ncbi.nlm.nih.gov/pubmed/22550399}
}
\author{
  P. Roebuck \email{paul_roebuck@comcast.net},
  E. Shannon Neeley \email{sneeley@stat.byu.edu},
  James M. Melott \email{jmmelott@mdanderson.org}
}
\seealso{
  \code{\linkS4class{RPPASpatialParams}},
  \code{\link[stats]{quantile}},
  \code{\link[mgcv]{gam}},
  \code{\link[mgcv]{s}},
  \code{\link[mgcv]{choose.k}}
}
\keyword{smooth}

