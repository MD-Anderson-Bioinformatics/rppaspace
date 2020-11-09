###
### $Id: rs4-normalize.R
### Used after RPPAFit to remove sample effect
###


##
## Module Variables
##
## :TODO: Migrate this to .onLoad since it's singleton-like
.NormEnv <- new.env(hash=TRUE)                   # Private environment
attr(.NormEnv, "name") <- "RPPASPACENormalizationMethods"


##=============================================================================
setClass("RPPANormalizationParams",
         representation(name="character",        # ui.name of norm method
                        method="character",      # key for norm method
                        arglist="OptionalList")) # optional named arguments


##-----------------------------------------------------------------------------
## Invoked by validObject() method.
validRPPANormalizationParams <- function(object) {

    #cat("validating", class(object), "object", "\n")
    msg <- NULL

    ## Validate arglist slot (if it exists)
    {
        arglist <- object@arglist
        if (!is.null(arglist)) {
            nm <- names(arglist)
            if (is.null(nm) || "" %in% nm) {
                msg <- c(msg, "all list components must be named")
            }
        }
    }

    ## Pass or fail?
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
}

setValidity("RPPANormalizationParams", validRPPANormalizationParams)


##-----------------------------------------------------------------------------
is.RPPANormalizationParams <- function(x) {
    is(x, "RPPANormalizationParams")
}


##-----------------------------------------------------------------------------
## Generates an RPPANormalizationParams object.
RPPANormalizationParams <- function(method=getRegisteredNormalizationMethodKeys(),
                                    arglist=NULL) {
    ## Check arguments
    method <- match.arg(method)
    if (!is.null(arglist)) stopifnot(is.list(arglist))

    ## Create new class
    new("RPPANormalizationParams",
        name=getRegisteredNormalizationMethodLabel(method),
        method=method,
        arglist=arglist)
}


##-----------------------------------------------------------------------------
## Returns a string representation of this instance. The content and format of
## the returned string may vary between versions. Returned string may be
## empty, but never null.
setMethod("paramString", signature(object="RPPANormalizationParams"),
          function(object,
                   slots=slotNames(object),
                   ...) {
    ## Check arguments
    stopifnot(is.character(slots) && length(slots) >= 1)

    ## :TODO: Implementation currently ignores the 'slots' argument
    ## and returns string containing parameters from various slots
    ## as though:
    ##     slotsToDisplay <- c("name", "method", "arglist")
    ##     paramString(np, slotsToDisplay)
    ##
    arglist <- paste(lapply(names(object@arglist),
                            function(nam, lst) {
                                val <- if (is.character(lst[[nam]])) {
                                           dQuote(lst[[nam]])
                                       } else {
                                           lst[[nam]]
                                       }
                                paste(nam, val, sep="=")
                            },
                            lst=object@arglist),
                     collapse=", ")
    paste(paste("name:", shQuote(object@name)), "\n",
          paste("method:", shQuote(object@method)), "\n",
          paste("arglist:", shQuote(arglist)), "\n",
          sep="")
})


##
## Private Methods
##

##-----------------------------------------------------------------------------
## Returns private environment for storing registered normalization methods
normenv <- function() {
    .NormEnv
}


##-----------------------------------------------------------------------------
registerPkgNormalizationMethods <- function() {
    registerNormalizationMethod("medpolish",
                                normalize.medpolish,
                                "Median Polish")
    registerNormalizationMethod("median",
                                normalize.median,
                                "Median")
    registerNormalizationMethod("house",
                                normalize.house,
                                "Housekeeping")
    registerNormalizationMethod("vs",
                                normalize.vs,
                                "Variable Slope")
    registerNormalizationMethod("none",
                                normalize.none,
                                "Do Not Normalize")
}

##-----------------------------------------------------------------------------
## Normalization method. Does not normalize.

normalize.none <- function(concs, ...) { 
    ## Do nothing, just return what was passed in.
    ## Present because RPPASPACE requires a normalization method to be declared
    message("No normalization of concentrations across slides will be performed.")
    stopifnot(is.matrix(concs) || is.data.frame(concs))

    normconcs <- concs
    ## Store processing info in "normalization" attribute
#    attr(normconcs, "normalization") <- c(list(method=method,
    attr(normconcs, "normalization") <- c(list(method="none",
                                               calc.medians=FALSE,
                                               sweep.cols=FALSE),
                                          as.list(new.env(hash=TRUE)),
                                          attr(normconcs, "normalization"))

    normconcs
}

##-----------------------------------------------------------------------------
## Normalization method. Fits additive model using Tukey's median polish.
normalize.medpolish <- function(concs, ...) {
    stopifnot(is.matrix(concs) || is.data.frame(concs))
    print("Medpolish normalization of concentrations across slides will be performed.")

    ## Median polish to normalize sample, slide effects
    ##   where:
    ##     row       - sample correction
    ##     residuals - polished concentrations
    ##
    pol <- medpolish(concs, trace.iter=FALSE, na.rm=TRUE)
    conc.medpol <- cbind(pol$row, pol$residuals)
    colnames(conc.medpol)[1] <- "Correction"

    normconcs <- conc.medpol
}


##-----------------------------------------------------------------------------
## Normalization method. Sample median is subtracted from each sample.
normalize.median <- function(concs, rowMedian, ...) {
    stopifnot(is.matrix(concs) || is.data.frame(concs))
    stopifnot(is.numeric(rowMedian))
	print("Median normalization of concentrations across slides will be performed.")
    normconcs <- sweep(concs, 1, rowMedian, FUN="-")
}


##-----------------------------------------------------------------------------
## Normalization method. Median of set of housekeeping antibodies is subtracted
## from each sample.
normalize.house <- function(concs,
                            antibodies,
                            ...) {
    stopifnot(is.matrix(concs) || is.data.frame(concs))
    stopifnot(is.character(antibodies) && length(antibodies) >= 1)
	print("House normalization of concentrations across slides will be performed.")

    if (!all(antibodies %in% colnames(concs))) {
        missingNames <- antibodies[!antibodies %in% colnames(concs)]
        stop(sprintf(ngettext(length(missingNames),
                              "argument %s specifies invalid %s column name: %s",
                              "argument %s specifies invalid %s column names: %s"),
                     sQuote("antibodies"),
                     sQuote("concs"),
                     paste(missingNames, collapse=", ")))
    }

    houseMedian <- apply(as.matrix(concs[, antibodies]),
                         1,
                         median,
                         na.rm=TRUE)
    normconcs <- sweep(concs, 1, houseMedian, FUN="-")
    ## Store method-specific info in "normalization" attribute
    attr(normconcs, "normalization") <- list(antibodies=antibodies,
                                             houseMedian=houseMedian)

    normconcs
}


##-----------------------------------------------------------------------------
## Normalization method (variable slope). Sample median is subtracted from
## each sample after applying multiplicative gamma.
normalize.vs <- function(concs, rowMedian,
                         ...) {
    stopifnot(is.matrix(concs) || is.data.frame(concs))
    stopifnot(is.numeric(rowMedian))

	print("Variable slope normalization of concentrations across slides will be performed.")

    ##-------------------------------------------------------------------------
    ## Estimates the multiplicative gamma terms from variable slope
    ## normalization. It takes as input the data matrix (with samples in
    ## the rows and antibodies in the columns). It is assumed that this matrix
    ## has already had the column median swept out from its columns.
    ## It outputs estimates of the gammas (multiplicative protein effects).
    estimateGamma <- function(Xhat) {
        stopifnot(is.matrix(Xhat) || is.data.frame(Xhat))

        nCol <- ncol(Xhat)
        gamma <- matrix(0, nrow=nCol, ncol=nCol)
        means <- apply(Xhat, 2, mean)

        for (i in seq(1, nCol-1)) {
            for (j in seq(i+1, nCol)) {
                r <- cor(Xhat[, i], Xhat[, j], use="complete.obs")
                a <- Xhat[, i]
                n <- length(a)
                tt <- r * sqrt((n-2) / (1-r^2))
                chk <- pt(tt, n-2, lower.tail=FALSE)
                if (chk < 0.05) {
                    eig <- eigen(var(cbind(Xhat[, i], Xhat[, j]), na.rm=TRUE))
                    tmp <- (-1) * eig$vectors[1, 2] / eig$vectors[2, 2]
                    gamma[i, j] <- tmp
                }
            }
        }
        gamma[gamma <= 0] <- 1
        upper <- upper.tri(gamma)
        ind.upper <- which(upper, arr.ind=TRUE)

        design <- matrix(0, ncol=nCol, nrow=nrow(ind.upper))
        for (i in seq_len(nrow(ind.upper))) {
            design[i, ind.upper[i, 1]] <- -1
            design[i, ind.upper[i, 2]] <- 1
        }

        loggamma <- log(gamma[upper])

        newrow <- rep((1 / nCol), nCol)
        nonsingular <- rbind(newrow, design)
        lestimateMean <- qr.solve(nonsingular, c(0, loggamma))

        estimate1 <- exp(lestimateMean)
    }

    if (ncol(concs) == 1) {
        print("Only one column.  Gammas not calculated.")
        gamma <- NA
        normconcs <- sweep(concs, 1, rowMedian, FUN="-")
    } else {
        gamma <- estimateGamma(concs)
        temp <- sweep(concs, 2, gamma, FUN="/")
        normconcs <- sweep(temp, 1, rowMedian, FUN="-")
    }
    ## Store method-specific info in "normalization" attribute
    attr(normconcs, "normalization") <- list(gamma=gamma)

    normconcs
}


##
## Public Methods
##

##-----------------------------------------------------------------------------
## Returns normalization method associated with key for invocation.
getRegisteredNormalizationMethod <- function(key) {
    method <- getRegisteredMethod(key, envir=normenv())$method
}


##-----------------------------------------------------------------------------
## Returns label associated with key for display by user interface.
getRegisteredNormalizationMethodLabel <- function(key) {
    ui.label <- getRegisteredMethod(key, envir=normenv())$ui.label
}


##-----------------------------------------------------------------------------
## Returns vector containing "keys" for all registered normalization methods.
getRegisteredNormalizationMethodKeys <- function() {
    keys <- getRegisteredMethodKeys(envir=normenv())
    if (length(keys) == 0) {
        stop("no registered normalization methods exist")
    }

    keys
}


##-----------------------------------------------------------------------------
## Registers specific normalization method for use by normalize() method.
registerNormalizationMethod <- function(key,
                                        method,
                                        ui.label=names(key)) {
    if (is.null(ui.label)) {
        ui.label <- key
    }
    ui.label <- as.character(ui.label)[1]

    registerMethod(key, method, ui.label=ui.label, envir=normenv())
}


##-----------------------------------------------------------------------------
## Method has two required input values:
##   1) the data matrix with samples in the rows and antibodies in the columns.
##   2) the name of the method of sample loading normalization. This argument
##      may be augmented with user-provided normalization methods.
##      Package-provided values are:
##
##      median - the sample median (row median) is subtracted from each sample
##      house  - housekeeping normalization. The median of a housekeeping
##               antibody or set of housekeeping antibodies are used. The
##               names of the antibodies to be used must be supplied as
##               a named argument to this method.
##      vs     - variable slope normalization. Here the sample median
##               is used along with a multiplicative gamma.
##      none   - do not normalize
##
setClassUnion("MatrixLike", c("matrix", "data.frame"))
setMethod("normalize", signature(object="MatrixLike"),
          function(object,
                   method=getRegisteredNormalizationMethodKeys(),
                   calc.medians=TRUE,
                   sweep.cols=calc.medians,
                   ...) {
    ## Check arguments
    method <- match.arg(method)

    if (method == "none") {
        calc.medians=FALSE
        sweep.cols=calc.medians
    }

    stopifnot(is.logical(calc.medians) && length(calc.medians) == 1)
    stopifnot(is.logical(sweep.cols) && length(sweep.cols) == 1)

    ## Call original function...
    rppaNormalize(object,
                  method=method,
                  calc.medians=calc.medians,
                  sweep.cols=sweep.cols,
                  ...)
})


##-----------------------------------------------------------------------------
## Performs normalization for sample loading after quantification.
rppaNormalize <- function(concs,
                          method,
                          calc.medians,
                          sweep.cols,
                          ...) {
    ## Check arguments
    if (!(is.matrix(concs) || is.data.frame(concs))) {
        stop(sprintf("argument %s must be matrix-like",
                     sQuote("concs")))
    }
    stopifnot(is.character(method) && length(method) == 1)
    stopifnot(is.logical(calc.medians) && length(calc.medians) == 1)
    stopifnot(is.logical(sweep.cols) && length(sweep.cols) == 1)
    if (!(method %in% getRegisteredNormalizationMethodKeys())) {
        stop(sprintf("argument %s must be registered normalization method",
                     sQuote("method")))
    }
    if (sweep.cols && !calc.medians) {
        stop(sprintf("argument %s must be TRUE if argument %s is TRUE",
                     sQuote("calc.medians"),
                     sQuote("sweep.cols")))
    }
    if (method == "vs" && !sweep.cols) {
        ## normalize.vs() method precondition
        stop(sprintf("argument %s must be TRUE if argument %s is %s",
                     sQuote("sweep.cols"),
                     sQuote("method"),
                     dQuote("vs")))
    }

    ## Begin processsing

    ## Create environment for method to use
    env <- new.env(hash=TRUE)
    normMethod <- getRegisteredNormalizationMethod(method)
    environment(normMethod) <- env

    ## Calculate medians?
    if (calc.medians) {
        env$rowMedian <- apply(concs, 1, median, na.rm=TRUE)
        env$colMedian <- apply(concs, 2, median, na.rm=TRUE)
    } else {
        env$rowMedian <- as.numeric(NA)
        env$colMedian <- as.numeric(NA)
    }

    ## Sweep columns first?
    if (sweep.cols) {
        concs <- sweep(concs, 2, env$colMedian, FUN="-")
    }

    ## Invoke method
    normconcs <- normMethod(concs, ...)

    ## Store processing info in "normalization" attribute
    attr(normconcs, "normalization") <- c(list(method=method,
                                               calc.medians=calc.medians,
                                               sweep.cols=sweep.cols),
                                          as.list(env),
                                          attr(normconcs, "normalization"))

    normconcs
}

