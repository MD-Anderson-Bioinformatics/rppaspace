###
### $Id: rs2-design.R
###


##=============================================================================
setClassUnion("OptionalFilename", c("character", "NULL"))
setClassUnion("OptionalList", c("list", "NULL"))
setClass("RPPADesignParams",
         representation(center="logical",
                        seriesToIgnore="OptionalList"
						))


##-----------------------------------------------------------------------------
is.RPPADesignParams <- function(x) {
    is(x, "RPPADesignParams")
}


##-----------------------------------------------------------------------------
RPPADesignParams <- function(
                             center=FALSE,
                             seriesToIgnore=NULL
							 ) {

	## Check arguments
 
    ## Convert numeric argument to logical counterpart
    if (is.numeric(center)) {
        center <- as.logical(center)
    }
    if (!is.logical(center)) {
        stop(sprintf("argument %s must be logical",
                     sQuote("center")))
    } else if (!(length(center) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("center")))
    }

    ## Create new class
    new("RPPADesignParams",
        center=center,
        seriesToIgnore=seriesToIgnore
		)
}

##-----------------------------------------------------------------------------
## Returns a string representation of this instance. The content and format of
## the returned string may vary between versions. Returned string may be
## empty, but never null.
setMethod("paramString", signature(object="RPPADesignParams"),
          function(object,
                   slots=slotNames(object),
                   ...) {
    ## Check arguments
    stopifnot(is.character(slots) && length(slots) >= 1)

    ## :TODO: Implementation currently ignores the 'slots' argument
    ## and returns string containing parameters from various slots
    ## as though:
    ##     slotsToDisplay <- c("center", "seriesToIgnore", "dilutionsInSeries")
    ##     paramString(dp, slotsToDisplay)
    ##
    paste(paste("center:", object@center), "\n",
          paste("seriesToIgnore:", shQuote(object@seriesToIgnore)), "\n",
          #paste("dilutionsInSeries:", shQuote(object@dilutionsInSeries)), "\n",
          sep="")
})


##-----------------------------------------------------------------------------
## Plot the series in an RPPA under a given design layout to see if the series
## makes sense under this layout.
## :TBD: Is this signature backwards?
setMethod("plot", signature(x="RPPA"),
          function(x,
                   measure="Net.Value",
                   main="",
                   ...) {
    ## Check arguments
    if (!is.character(measure)) {
        stop(sprintf("argument %s must be character",
                     sQuote("measure")))
    } else if (!(length(measure) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("measure")))
    } else if (!(measure %in% colnames(x@data))) {
        stop(sprintf("invalid measure %s",
                     sQuote(measure)))
    }

    if (!is.character(main)) {
        stop(sprintf("argument %s must be character",
                     sQuote("main")))
    } else if (!(length(main) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("main")))
    }

    ## Begin processing
    vert <- x@data[, measure]
    horz <- x@data$Steps

    if (!nzchar(main)) {
        main <- .mkPlotTitle(paste(measure, "Intensity vs. Dilution Step"),
                             x@antibody)
    }
    par(mfrow=c(1, 1))  # Avoid existing partitions of graphic device

    steps <- x@data$Steps[x@data$Spot.Type %in% spottype.sample]

    plot(c(min(steps), max(steps)),
         c(min(vert), max(vert)),
         main=main,
         sub=paste("File:", x@file),
         type="n",
         xlab="Dilution Step",
         ylab="Intensity")

	series <- x@data$Series.Id
    s <- seriesNames(x) # Strip out control spots
    bow <- rainbow(length(s))

    for (i in seq_along(s)) {
        lines(x=horz[series == s[i]],
              y=vert[series == s[i]],
              col=bow[i],
              type="b")
    }
})