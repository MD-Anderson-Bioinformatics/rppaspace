###
### $Id: Directory-class.R
###

#:KRC: This entire class is completely unnecessary.

##
## Classes
##

##=============================================================================
setClass("Directory",
         representation(path="character"))
setClassUnion("OptionalDirectory", c("Directory", "NULL"))


##-----------------------------------------------------------------------------
## Invoked by validObject() method.
validDirectory <- function(object) {

    #cat("validating", class(object), "object", "\n")
    msg <- NULL

    ## Validate path slot
    {
        path <- object@path

        ## Ensure path exists (and appropriate filesystem object)
        if (!file.exists(path)) {
            msg <- c(msg, sprintf("path %s does not exist",
                                  dQuote(path)))
        } else if (!file.info(path)$isdir) {
            msg <- c(msg, sprintf("path %s is not directory",
                                  dQuote(path)))
        }
    }

    ## Pass or fail?
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
}

setValidity("Directory", validDirectory)


##-----------------------------------------------------------------------------
is.Directory <- function(x) {
    is(x, "Directory")
}


##-----------------------------------------------------------------------------
## Generator method
Directory <- function(path) {
    ## Check arguments
    if (!is.character(path)) {
        stop(sprintf("argument %s must be character",
                     sQuote("path")))
    } else if (!(length(path) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("path")))
    }

    ## Create new class
    new("Directory",
        path=path)
}


##
## Generic Methods (Coercion)
##

##-----------------------------------------------------------------------------
## Coercion method
setAs("Directory", "character",
      function(from) {
          from@path
      })


##-----------------------------------------------------------------------------
## Coercion method
setAs("character", "Directory",
      function(from) {
          Directory(from)
      })


##
## Public Routines
##

##-----------------------------------------------------------------------------
## :TODO: make generic...
pathname <- function(object) {
    if (!is.Directory(object)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("object"), "Directory"))

    }

    as(object, "character")
}

