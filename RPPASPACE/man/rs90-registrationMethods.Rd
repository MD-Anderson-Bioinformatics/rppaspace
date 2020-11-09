% $Id: rs90-registrationMethods.Rd

\name{registrationMethods}
\alias{getRegisteredObject}
\alias{getRegisteredObjectKeys}
\alias{getRegisteredMethod}
\alias{getRegisteredMethodKeys}
\alias{registerClassname}
\alias{registerMethod}
\title{Generic Registration Methods}
\description{
  These routines represent the lower-level access for registration,
  which enables data-driven access by other routines. This represents
  the initial implementation and may change in the future.
}
\usage{
getRegisteredObject(key,
                    envir,
                    objtype=c("method", "classname"))
getRegisteredObjectKeys(envir)
getRegisteredMethodKeys(envir)
getRegisteredMethod(key, envir)
registerClassname(key, classname, \dots, envir)
registerMethod(key, method, \dots, envir)
}
\arguments{
  \item{key}{character string representing a registered object}
  \item{envir}{environment where registered objects are stored}
  \item{objtype}{character string specifying type of registered object}
  \item{method}{function to register}
  \item{classname}{character string specifying class name to register}
  \item{\dots}{extra arguments for generic routines}
}
\details{
  The routines here basically use the \code{envir} passed by the caller
  as a private hash table. The \code{key} argument is used to retrieve
  the value stored during registration.

  \code{getRegisteredMethodKeys} is an alias for \code{getRegisteredObjectKeys}.
  \code{getRegisteredMethod} invokes \code{getRegisteredObject} with
  \code{objtype} of "method".
}
\value{
  \code{getRegisteredObject} returns list associated with \code{key}.
  List components will include named arguments provided upon registration.

  \code{getRegisteredObjectKeys} returns vector of \code{key}s for all
  registered objects.

  \code{registerClassname} is invoked for its side effect, which is registering
  \code{classname} (and optionally other named arguments) by association to
  \code{key}.

  \code{registerMethod} is invoked for its side effect, which is registering
  \code{method} (and optionally other named arguments) by association to
  \code{key}.
}
\author{
  P. Roebuck \email{paul_roebuck@comcast.net},
  James M. Melott \email{jmmelott@mdanderson.org}
}
\seealso{
  The \code{envir} argument of \code{\link[base]{eval}},
  \code{\link[base]{get}}, and \code{\link[base]{exists}}.

  \code{\link[base]{assign}},
  \code{\link[base]{get}}
}
\keyword{data}
\keyword{internal}

