###
### $Id: NULL-methods.R
###

requireNamespace("methods")
options(warn=1)


##
## Methods
##


##-----------------------------------------------------------------------------
mkNumericHandlerMethod <- function(methodName) {
    stopifnot(is.character(methodName) && length(methodName) == 1)

    ##-------------------------------------------------------------------------
    setMethod(methodName,
        signature(object="NULL"),
        function(object) {
            NA_real_
        })
}


methodNames <- c("normalize", "qcprob")
sapply(methodNames, mkNumericHandlerMethod)
rm(mkNumericHandlerMethod)

