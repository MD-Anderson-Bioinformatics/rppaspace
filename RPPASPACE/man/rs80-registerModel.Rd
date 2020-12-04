% $Id: rs80-registerModel.Rd

\name{registerModel}
\alias{getRegisteredModel}
\alias{getRegisteredModelKeys}
\alias{getRegisteredModelLabel}
\alias{registerModel}
\title{Model Registration Methods}
\description{
  These routines represent the high-level access for model registration,
  which enables data-driven access by other routines.
  This represents the initial implementation and may change in the future.
}
\usage{
getRegisteredModel(key)
getRegisteredModelLabel(key)
getRegisteredModelKeys()
registerModel(key, classname, ui.label=names(key))
}
\arguments{
  \item{key}{character string representing a registered model}
  \item{classname}{character string specifying Model class name to register}
  \item{ui.label}{character string specifying label to display by UI}
}
\value{
  \code{getRegisteredModel} returns the \code{classname} associated with
  \code{key}.

  \code{getRegisteredModelLabel} returns the \code{ui.label} associated with
  \code{key}.

  \code{getRegisteredModelKeys} returns vector of \code{key}s for all
  registered models.

  \code{registerModel} is invoked for its side effect, which is registering
  \code{classname} and \code{ui.label} by association to \code{key}.
}
\author{
  P. Roebuck \email{paul_roebuck@comcast.net},
  James M. Melott \email{jmmelott@mdanderson.org}
}
\seealso{
  \code{\link{getRegisteredObject}},
  \code{\link{getRegisteredObjectKeys}},
  \code{\link{registerClassname}}
}
\examples{
## Create new (but nonfunctional) fit model
## Not run due to lack of capability to unregister class
  \dontrun{
setClass("TestFitClass",
         representation("FitClass",
                        testfit="character"),
         prototype(testfit="TestFitClass"))

## Register fit model to enable its use by package
registerModel("testfit", "TestFitClass", "Registered Test Fit Class")

## Show all registered fit models
sapply(getRegisteredModelKeys(),
       function(key) {
           c(model=getRegisteredModel(key),
             label=getRegisteredModelLabel(key))
       })
   }
}
\keyword{data}

