###
### $Id: rs8-prefitqc.R
###


##=============================================================================
setClass("RPPAPreFitQC",
         representation("VIRTUAL"))


##=============================================================================
setClass("DS5RPPAPreFitQC",
         contains="RPPAPreFitQC",
         representation(slopediff="numeric",
                        cvs="numeric",
                        drdiffs="numeric",
                        slopes="numeric",
                        percentgood="numeric",
                        adjusted="logical",
                        antibody="character"),
         prototype(adjusted=FALSE))


##-----------------------------------------------------------------------------
## Returns TRUE if class of argument is subclass of RPPAPreFitQC.
is.RPPAPreFitQC <- function(x) {
    extends(class(x), "RPPAPreFitQC")
}


##-----------------------------------------------------------------------------
## Generates subclass instance of an RPPAPreFitQC object (factory-style)
RPPAPreFitQC <- function(rppa,
                         useAdjusted=FALSE) {
    ## Check arguments
    if (!is.RPPA(rppa)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("rppa"), "RPPA"))
    } else {
        measures <- c("Net.Value", "Raw.Value")
        if (as.logical(useAdjusted)) {
            measures <- paste("Adj", measures, sep=".")
        }
        reqdNames <- measures
        if (!(all(reqdNames %in% colnames(rppa@data)))) {
            missingNames <- reqdNames[!reqdNames %in% colnames(rppa@data)]
            stop(sprintf(ngettext(length(missingNames),
                                  "argument %s missing required column: %s",
                                  "argument %s missing required columns: %s"),
                         sQuote("rppa"),
                         paste(missingNames, collapse=", ")))
        }
    }

	## Requires columns
	reqdNames <- c("Spot.Type", "Dilution")
	if (!(all(reqdNames %in% colnames(rppa@data)))) {
		missingNames <- reqdNames[!reqdNames %in% colnames(rppa@data)]
		stop(sprintf(ngettext(length(missingNames),
							  "argument %s missing required column: %s",
							  "argument %s missing required columns: %s"),
					 sQuote("rppa"), paste(missingNames, collapse=", ")))
	}

	if (!(any(rppa@data$Spot.Type %in% spottype.positivecontrol))) {
		stop("contains no positive controls")
	}

	Spot.Type <- NA
    ## Begin processing
    ndilutions <- with(subset(rppa@data,
                    Spot.Type %in% spottype.positivecontrol),
                    length(unique(Dilution)))
					
##TODO: Only supports slides with 5 dilution steps
    switch(EXPR=as.character(ndilutions),
           "5"=DS5RPPAPreFitQC(rppa, ndilutions, measures, useAdjusted),
           stop(sprintf("unsupported slide layout - ndilutions = %d",
                        ndilutions)))
}


##-----------------------------------------------------------------------------
## Generates a DS5RPPAPreFitQC object
DS5RPPAPreFitQC <- function(rppa,
							numDilutionsInSeries,
                            measures,
                            useAdjusted=FALSE) {

    ## Check arguments
	stopifnot(numDilutionsInSeries == 5)
    stopifnot(is.RPPA(rppa))
    stopifnot(is.character(measures) && length(measures) == 2)
    stopifnot(is.logical(useAdjusted) && length(useAdjusted) == 1)

    tracking <- rppa@tracking

    ##-------------------------------------------------------------------------
    ## Calculate coefficient of variance
    CV <- function(thedata) {
        stopifnot(is.numeric(thedata))

        log.measure <- log(thedata)
        sd(log.measure, na.rm=TRUE) / mean(log.measure, na.rm=TRUE)
    }


    ##-------------------------------------------------------------------------
    percentGood <- function(df) {
        stopifnot(is.data.frame(df))

        mean.total <- df[[measure.total]]
        mean.net <- df[[measure.net]]
        bsr <- (mean.total - mean.net) / mean.net
        100 * (sum(bsr < 1) / length(bsr))
    }


    ##-------------------------------------------------------------------------
    slopeDiff <- function(df) {
        stopifnot(is.data.frame(df))

        #Transform dilutions for Positive controls 
        #to corresponding 0:(max # dilutions -1) values for use in model
        dilutions<-df$Dilution
        uniqueDilutions <- sort(unique(dilutions)) 
        steps <- match(dilutions, uniqueDilutions) - 1 

        mean.net <- df[[measure.net]]
        model <- lm(log(mean.net) ~ steps)
        slope <- coef(model)["steps"]
        ## See how far the observed slopes are from that of perfect dilution.
        perfect.slope <- 0.6931
		
        abs(slope - perfect.slope)
    }

    ## Begin processing
    if (!requireNamespace("timeDate")) {
        stop(sprintf("%s package required for pre-fit QC in the %s method",
                     sQuote("timeDate"),
                     sQuote("DS5RPPAPreFitQC")))
    }

    layout <- rppa@data[,!(names(rppa@data) %in% c("Series.Id"))]

	#Incorporates filtering series to be ignored by only using sample points that are also to be made part of the curve
    samples.df <- rppa@data[rownames(tracking[tracking$isSample & tracking$makePartOfCurve,]),]
	
	posCtrlRows <- rownames(tracking[tracking$isPosCtrl | tracking$isNoise,])
    positives_or_noise.df <- rppa@data[posCtrlRows,]
    rownames(samples.df) <- 1:nrow(samples.df)

    dilutions <- sort(unique(positives_or_noise.df$Dilution), decreasing=TRUE)
    ndilutions <- length(dilutions)  # number of dilution series per slide

    stopifnot(ndilutions == numDilutionsInSeries)
	
	strengths <- sort(unique(rppa@data$Steps), decreasing=TRUE)

    ## Set measure names to use
    measure.net   <- measures[1]
    measure.total <- measures[2]

    ## Make list containing 'Sub.Row' value pairs for each PC dilution
    ## These are typically listed as PosCtrl points
    plocats <- lapply(dilutions,
        function(dilution, layout.df) {
            pc.tf <- with(layout.df,
                ((Spot.Type %in% spottype.positivecontrol | Spot.Type %in% spottype.noise) &
                    Dilution == dilution))
            unique(layout.df[pc.tf, "Sub.Row"])
        },
        layout.df=rppa@data)
    names(plocats) <- dilutions

    positive.locats <- lapply(plocats,
        function(locat, pos.df) {
            locat.tf <- with(pos.df, Sub.Row %in% locat)
            pos.df[locat.tf, measure.net]
        },
        pos.df=positives_or_noise.df[, c("Sub.Row", measure.net)])

    ## For each dilution step of positive controls, ...
    slopes <- NULL
    cvs <- NULL

    for (positives in positive.locats) {
        ## Calculate slope for measure (good if near zero)
        y <- positives
        x <- seq_along(y)    # seq_len(nrows(positives_or_noise.df)/ndilutions) = 1:96
        model <- lm(y~x)
        slope <- abs(coef(model)["x"])
        slopes <- c(slopes, slope)

        ## Calculate coefficient of variances for measure
        cv <- CV(positives)
	
        cvs <- c(cvs, cv)
    }

    names(cvs) <- strengths

    ## Calculate slope of step series for measure of positive controls
    slopediff <- slopeDiff(positives_or_noise.df)

    ## Calculate dynamic range difference for measure of each dilution step
    ## of samples
    drdiffs <- NULL
    for (dilution in dilutions) {
        samp.tf <- with(samples.df,
            (Spot.Type %in% c("Sample") &
                Dilution == dilution))
        samples.mean.net <- samples.df[samp.tf, measure.net]
        drdiff <- diff(range(samples.mean.net))
        drdiffs <- c(drdiffs, drdiff)
    }

    names(drdiffs) <- strengths

    ## Calculate percent good
    percentgood <- percentGood(samples.df)

	printQCParts <- TRUE
	if (printQCParts) {
		z = (
			3.013 
			-(0.9585 * slopediff) 
			-(21.51  * cvs[1])
			-(43.06 * cvs[2]) 
			-(19.29 * cvs[4]) 
			-(0.01574 * slopes[5]) 
			+(0.00003885 * drdiffs[2]) 
			-(0.00004131 * drdiffs[4])
			+(0.01271 * percentgood))

		qc_score = exp(z) / (1 + exp(z))

#		print(paste(3.013, -(0.9585 * slopediff), -(21.51  * cvs[1]), 
#			-(43.06 * cvs[2]), -(19.29 * cvs[4]), -(0.01574 * slopes[5]), 
#			(0.00003885 * drdiffs[2]), -(0.00004131 * drdiffs[4]), 
#			(0.01271 * percentgood), "=", z , "==>", qc_score
#			))
	}
	
    ## Create new class
    new("DS5RPPAPreFitQC",
        antibody=rppa@antibody,
        adjusted=useAdjusted,
        slopediff=slopediff,
        cvs=cvs,
        drdiffs=drdiffs,
        slopes=slopes,
        percentgood=percentgood)
}


##-----------------------------------------------------------------------------
setMethod("qcprob", signature(object="RPPAPreFitQC"),
          function(object,
                   ...) {
    stop(sprintf("%s method must be implemented by any subclass of %s",
                 sQuote("qcprob"),
                 sQuote("RPPAPreFitQC")))
})


##-----------------------------------------------------------------------------
setMethod("qcprob", signature(object="DS5RPPAPreFitQC"),
          function(object,
                   ...) {
    ##-------------------------------------------------------------------------
    pred.model <- function(x) {
        slopediff   <- x@slopediff
        cv1         <- x@cvs[1]
        cv2         <- x@cvs[2]
        cv8         <- x@cvs[4]
        step16slope <- x@slopes[5]
        drdiff2     <- x@drdiffs[2]
        drdiff8     <- x@drdiffs[4]
        percentgood <- x@percentgood

        z <- 3.013 -
            (0.9585     * slopediff) -
            (21.51      * cv1) -
            (43.06      * cv2) -
            (19.29      * cv8) -
            (0.01574    * step16slope) +
            (0.00003885 * drdiff2) -
            (0.00004131 * drdiff8) +
            (0.01271    * percentgood)
        as.numeric(z)  ## Strip attributes
    }


    ## Begin processing
    z <- pred.model(object)
    exp(z) / (1 + exp(z)) # probability of good slide
})


##-----------------------------------------------------------------------------
setMethod("summary", signature(object="RPPAPreFitQC"),
          function(object,
                   ...) {
    stop(sprintf("%s method must be implemented by any subclass of %s",
                 sQuote("qcprob"),
                 sQuote("RPPAPreFitQC")))
})


##-----------------------------------------------------------------------------
setMethod("summary", signature(object="DS5RPPAPreFitQC"),
          function(object,
                   ...) {
    cat(sprintf("A %s object representing QC measures for antibody %s",
                class(object),
                object@antibody), "\n")
    cat("Measures spatially adjusted:", object@adjusted, "\n")
    cat("Difference from perfect slope:", object@slopediff, "\n")
    cat("Coefficients of variances:", "\n")
    cat("\t", paste(object@cvs, collapse=", "), "\n")
    cat("Difference in ranges:", "\n")
    cat("\t", paste(round(object@drdiffs, 2), collapse=", "), "\n")
    cat("Step slopes:", "\n")
    cat("\t", paste(round(object@slopes, 5), collapse=", "), "\n")
    cat("Percent good:", object@percentgood, "\n")
})

