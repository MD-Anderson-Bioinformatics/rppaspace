% $Id: rs80-registerNormalizationMethod.Rd

\name{registerNormalizationMethod}
\alias{getRegisteredNormalizationMethod}
\alias{getRegisteredNormalizationMethodKeys}
\alias{getRegisteredNormalizationMethodLabel}
\alias{registerNormalizationMethod}
\title{Normalization Method Registration Methods}
\description{
  These routines represent the high-level access for normalization method
  registration, which enables data-driven access by other routines.
  This represents the initial implementation and may change in the future.
}
\usage{
getRegisteredNormalizationMethod(key)
getRegisteredNormalizationMethodLabel(key)
getRegisteredNormalizationMethodKeys()
registerNormalizationMethod(key, method, ui.label=names(key))
}
\arguments{
  \item{key}{character string representing a registered normalization method}
  \item{method}{function to invoke for normalization}
  \item{ui.label}{character string specifying label to display by UI}
}
\value{
  \code{getRegisteredNormalizationMethod} returns the \code{method} associated
  with \code{key}.

  \code{getRegisteredNormalizationMethodLabel} returns the \code{ui.label}
  associated with \code{key}.

  \code{getRegisteredNormalizationMethodKeys} returns vector of \code{key}s
  for all registered normalization methods.

  \code{registerNormalizationMethod} is invoked for its side effect, which is
  registering \code{method} and \code{ui.label} by association to \code{key}.
}
\author{
  P. Roebuck \email{paul_roebuck@comcast.net},
  James M. Melott \email{jmmelott@mdanderson.org}
}
\seealso{
  \code{\link{getRegisteredObject}},
  \code{\link{getRegisteredObjectKeys}},
  \code{\link{registerMethod}}
}
\examples{
  \dontrun{
## Not run due to lack of capability to unregister methods
## Create new normalization method
normalize.testNorm <- function(concs, bar) {
    return(normconcs <- concs - bar)
}

## Register normalization method to enable its use by package
registerNormalizationMethod("testNorm", normalize.testNorm, "Registered Test Normalization Class")

## Use it...
concs <- matrix(runif(500), nrow=10)
normalize(concs, method="testNorm", bar=0.005)

## Show all registered fit models
sapply(getRegisteredNormalizationMethodKeys(),
       function(key) {
           c(key = key,
             label=getRegisteredNormalizationMethodLabel(key))
       })
  }
}
\keyword{data}

