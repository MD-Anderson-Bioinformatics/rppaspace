###
### $Id: registration.R
### Generic registration routines
###

#:KRC: This whole thing still smells as though it is rebuilding
# an inferior version of S4 derived classes

##
## Private Methods
##

##-----------------------------------------------------------------------------
.validate.classname <- function(classname) {
  #:KRC: user-hostile.
    if (!is.character(classname)) {
        stop(sprintf("argument %s must be character",
                     sQuote("classname")))
    } else if (!(length(classname) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("classname")))
    } else if (!nzchar(classname)) {
        stop(sprintf("argument %s must not be empty string",
                     sQuote("classname")))
    }
}


##-----------------------------------------------------------------------------
.validate.envir <- function(envir) {
  # lists and environments are inter-convertible, and yuo can look things
  # up ineither one. Why require an environment, insted of a successful
  # call to (as(, "environment") ?
    if (!is.environment(envir)) {
        stop(sprintf("argument %s must be environment",
                     sQuote("envir")))
    }
}


##-----------------------------------------------------------------------------
.validate.key <- function(key) {
  #:KRC: user-hostile.
    if (!is.character(key)) {
        stop(sprintf("argument %s must be character",
                     sQuote("key")))
    } else if (!(length(key) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("key")))
    } else if (!nzchar(key)) {
        stop(sprintf("argument %s must not be empty string",
                     sQuote("key")))
    }
}


##-----------------------------------------------------------------------------
.validate.method <- function(method) {
    if (!is.function(method)) {
        stop(sprintf("argument %s must be function",
                     sQuote("method")))
    }
}


##
## Package Methods
##

##-----------------------------------------------------------------------------
## Returns list associated with key. List will contain named arguments provided
## upon registration.
getRegisteredObject <- function(key,
                                envir,
                                objtype=c("method",
                                          "classname")) {
    ## Check arguments
    .validate.key(key)
    .validate.envir(envir)
    objtype <- match.arg(objtype)

    ## Begin processsing
    if (!exists(key, envir=envir)) {
        stop(sprintf("no registered %s associated with argument %s (%s)",
                     objtype,
                     sQuote(key),
                     key))
    }

    return(get(key, envir=envir))
}


##-----------------------------------------------------------------------------
## Returns vector containing "keys" for all registered objects.
getRegisteredObjectKeys <- function(envir) {
    ## Check arguments
    .validate.envir(envir)

    ## Begin processsing
    return(keys <- objects(envir=envir))
}

getRegisteredMethodKeys <- getRegisteredObjectKeys


##-----------------------------------------------------------------------------
## Returns list associated with key. List will contain method as well as any
## other named arguments provided upon registration. 
getRegisteredMethod <- function(key,
                                envir) {
    getRegisteredObject(key, envir, "method")
}


##-----------------------------------------------------------------------------
## Registers classname (and optionally other named arguments) by association to
## key.
registerClassname <- function(key,
                              classname,
                              ...,
                              envir) {
    ## Check arguments
    .validate.key(key)
    .validate.classname(classname)
    .validate.envir(envir)

    ## Begin processsing
    dots <- list(...)
    registerable <- c(list(classname=classname),
                      namedArgs <- dots[names(dots) != ""])
    assign(key, registerable, envir=envir)
}


##-----------------------------------------------------------------------------
## Registers method (and optionally other named arguments) by association to
## key.
registerMethod <- function(key,
                           method,
                           ...,
                           envir) {
    ## Check arguments
    .validate.key(key)
    .validate.method(method)
    .validate.envir(envir)

    ## Begin processsing
    dots <- list(...)
    registerable <- c(list(method=method),
                      namedArgs <- dots[names(dots) != ""])
    assign(key, registerable, envir=envir)
}

