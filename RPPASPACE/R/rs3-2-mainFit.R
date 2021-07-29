###
### $Id: rs3-2-mainFit.R
###


##
## Module Variables
##
## :TODO: Migrate this to .onLoad since it's singleton-like
.ModelEnv <- new.env(hash=TRUE)
attr(.ModelEnv, "name") <- "RPPASPACEModels"


##
## Private Methods
##

##-----------------------------------------------------------------------------
## Returns private environment for storing registered models
modelenv <- function() {
    return(.ModelEnv)
}


##-----------------------------------------------------------------------------
registerPkgFitModels <- function() {
    registerModel("logistic", "LogisticFitClass", "Logistic")
    registerModel("cobs", "CobsFitClass", "Monotone Increasing B-spline")
    registerModel("loess", "LoessFitClass", "Loess")
}


##
## Public Methods
##

##-----------------------------------------------------------------------------
## REQUIRED INPUTS:
## rppa    - an RPPA object
## measure - name of column in rppa@data to view as intensity
##
##-----------------------------------------------------------------------------
## All optional inputs are handled through RPPAFitParams
##
## OPTIONAL INPUTS CONTROLLING THE ALGORITHM:
## model   - which statistical model to fit
## method  - how to fit the alpha and beta values.
## xform   - function used to transform the values before fitting
## trim    - if 0, concentrations will not be trimmed / bounded; otherwise,
##           trim the concentrations using this value as the trim level.
## ci      - if true, calculate a 90% confidence interval for the LC50
##           rppaspace values
## ignoreNegative - if true, then values below zero are converted to NA before
##           estimating alpha and beta. Only matters for the 'quantiles' method.
##
## OPTIONAL INPUTS CONTROLLING VERBOSITY:
## trace   - passed to nls in the 'bayesian' portion
## verbose - if true, have the function tell you what it is doing
## veryVerbose - if true, have the function overwhelm you telling you what it
##           is doing
## warnLevel - used to set the 'warn' option before calling 'rlm'.
##           since this is wrapped in a 'try', it won't cause failure but will
##           give us a chance to figure out which dilution series failed.
##           Setting warnLevel to two or greater may change computed results.


##-----------------------------------------------------------------------------
## Refactored this to improve maintainability; everything now happens once.
## Retain this generator for backwards compatibility.
RPPAFit <- function(rppa,
                    measure,
                    model="logistic",
                    xform=NULL,
                    method=c("nls", "nlrob", "nlrq"),
                    trim=2,
                    ci=FALSE,
                    ignoreNegative=TRUE,
                    trace=FALSE,
                    verbose=FALSE,
                    veryVerbose=FALSE,
					warnLevel=0,
                    residualsrotation=as.integer(0)) {
    ## Check arguments
    ## [1] 'rppa' is checked in the call to RPPAFitFromParams
    ## [2] Everything else is checked in RPPAFitParams
    params <- RPPAFitParams(measure,
                            model,
                            xform,
                            method,
                            trim,
                            ci,
                            ignoreNegative,
                            trace,
                            verbose,
                            veryVerbose,
							warnLevel,
                            residualsrotation)
    RPPAFitFromParams(rppa, params)
}


##-----------------------------------------------------------------------------
## Collect parameters in one place and make sure they are reasonably sensible
RPPAFitParams <- function(measure,
                          model="logistic",
                          xform=NULL,
                          method=c("nls", "nlrob", "nlrq"),
                          trim=2,
                          ci=FALSE,
                          ignoreNegative=TRUE,
                          trace=FALSE,
                          verbose=FALSE,
                          veryVerbose=FALSE,
                          warnLevel=0,
						  residualsrotation=as.integer(0)) {
    ## Check arguments
  
    ## Start with the critical parameter, 'model', since it tells us which
    ## statistical model we are supposed to fit. At the time the curve is
    ## fit, 'model' has to be the name of a registered FitClass. You could
    ## conceivably define your own class and register it after setting up the
    ## fit parameters. Because of that fact, we defer some checking until
    ## later (in RPPAFitFromParams).
    if (!is.character(model)) {
        stop(sprintf("argument %s must be a character string",
                     sQuote("model")))
    }
    model <- model[1]

    ## This parameter is also critical, and we cannot complete checking it
    ## until later, when we know whether it refers to an actual column in
    ## the RPPA object being fit.
    if (!is.character(measure)) {
        stop(sprintf("argument %s must be character",
                     sQuote("measure")))
    } else if (!(length(measure) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("measure")))
    }

    if (!is.null(xform)) {
        if (!is.function(xform)) {
            stop(sprintf("argument %s must be function, if specified",
                         sQuote("xform")))
        }
    }

    ## Remaining parameters are less important, so we mostly just transform
    ## them into the correct type and length.
    method <- match.arg(method)
    ci <- as.logical(ci)[1]
    ignoreNegative <- as.logical(ignoreNegative)[1]
    trace <- as.logical(trace)[1]
    verbose <- as.logical(verbose)[1]
    veryVerbose <- as.logical(veryVerbose)[1]
    warnLevel <- as.integer(warnLevel)[1]

    if (!(is.numeric(trim) || is.logical(trim))) {
        stop(sprintf("argument %s must be numeric (or logical)",
                     sQuote("trim")))
    }
    trim <- if (is.logical(trim) && isTRUE(trim[1])) {
                ## Substitute numeric equivalent for original behavior
                formals(RPPAFitParams)$trim
            } else {
                as.numeric(trim)[1]
            }
    if (is.na(trim) || trim < 0) {
        stop(sprintf("argument %s must be a non-negative quantity",
                     sQuote("trim")))
    }

    ## Create new class
    new("RPPAFitParams",
        measure=measure,
        xform=xform,
        method=method,
        ci=ci,
        ignoreNegative=ignoreNegative,
        trace=trace,
        verbose=verbose,
        veryVerbose=veryVerbose,
        warnLevel=warnLevel,
        trim=trim,
        model=model)
}


##-----------------------------------------------------------------------------
## Returns a string representation of this instance. The content and format of
## the returned string may vary between versions. Returned string may be
## empty, but never null.
setMethod("paramString", signature(object="RPPAFitParams"),
          function(object,
                   slots=slotNames(object),
                   ...) {
    ## :TODO: Implementation currently ignores the 'slots' argument and
    ## returns string containing parameters from various slots as though:
    ##
    ##     slotsToDisplay <- c("measure", "model", "method", "trim",
    ##                         "ci", "ignoreNegative", "warnLevel")
    ##     paramString(fp, slotsToDisplay)
    ##
    paste(paste("measure:", shQuote(object@measure)), "\n",
          paste("model:", shQuote(object@model)), "\n",
          paste("method:", shQuote(object@method)), "\n",
          paste("trim:", object@trim), "\n",
          paste("ci:", object@ci), "\n",
          paste("ignoreNegative:", object@ignoreNegative), "\n",
          paste("warnLevel:", object@warnLevel), "\n",
          sep="")
})


##-----------------------------------------------------------------------------
## Compute a multiple of the logit transform.
## in practice, epsilon never changes.
.calcLogitz <- function(data,   #intensity
                        alpha,  #minimum sample value on slide of points used to make curve (after values less than 0 have been set to 0)
                        beta,   #max_range of the values on slide of points used to make curve (max - min -> 0:max_range)
                        gamma,  #hard-coded to 1 in calling method
                        epsilon=0.001) {
    z <- (data - alpha) / beta  #(individual intensity - lowest intensity) / range of intensities  (converts to range 0 : 1) 
    ## Set points at min and max of points used to make the dilution curve to epsilon and 1-epsilon respectively
    z[z < epsilon] <- epsilon
    z[z > (1-epsilon)] <- 1-epsilon
    temp <- log(z / (1-z)) / gamma
}


##-----------------------------------------------------------------------------
## Computes crude estimates of the parameters (i.e., alpha, beta, gamma, and
## the EC50s) for non-ignored sample points so routine can be initialized.
.firstPass <- function(yval,
                       rppa,
                       ignoreNegative,
                       epsilon=1e-4) {

    ## 'yval' is a vector of intensity values selected from the slide file.
    ## in practice, 'epsilon' never changes

    ## Check arguments
    stopifnot(is.numeric(yval))
    stopifnot(is.logical(ignoreNegative) && length(ignoreNegative) == 1)
    stopifnot(is.numeric(epsilon))

    ## Begin processing
    if (ignoreNegative) {
        yval[yval < 0] <- NA
    }

    ## 'temp' is a vector of intensity values of non-ignored sample points selected from the slide file.
    temp <- yval[rppa@tracking$makePartOfCurve]

    ## Compute a robust estimate for alpha and beta
    ## Tried to use robust 'rnls' from 'sfsmisc' to fit alpha and beta below
    ## but found it was still too sensitive to outliers.
    ## This is a crude way to get alpha and beta but it seems to be robust.

    ## HERE IS THE PROBLEM!
    ## Quantiles forces a double trim in conjunction with calcLogitz
    ## and epsilon above, artificially shrinking alpha and beta
    ## Try to match version 0.12
    # lBot <- quantile(temp, probs=0.05, na.rm=TRUE)
    # lTop <- quantile(temp, probs=0.95, na.rm=TRUE)

    lBot <- min(temp, na.rm=TRUE)
    lTop <- max(temp, na.rm=TRUE)

	#Only non-ignored sample points are used to calculate the alpha, beta, and gamma,
	#but the Logitz curve for starting values are calculated for all samples
    lindata <- .calcLogitz(yval, lBot, lTop-lBot, 1)
	
    sampleseries <- seriesNames(rppa) # limits to series sample points only
	nonignoredSeries <- seriesToUseToMakeCurve(rppa)

    ##-------------------------------------------------------------------------
    .estimateLogisticSlope <- function(ser,
                                       allSteps,
                                       seriesAll,
                                       ld) {

        ## Get the items in this dilution series
        items <- seriesAll == ser
        ## Get the log2 steps for this series
        steps <- allSteps[items]
        ## Get the transformed intensity values
        x <- ld[items]
        ## Only use valid values
        y <- x[!is.na(x) & !is.infinite(x)]

        if (all(is.na(y))) {
            NA
        } else {
            ## ychange / xchange
			#TODO: Possible issue here if max or min are not 
			#the max step point and min step points respectively. (JMM)
            (max(y) - min(y)) / (max(steps) - min(steps))
        }
    }

    ## Initial estimate of the logistic slope across the dilution steps.
    ## There should be a more comprehensible way to write this:
    dat <- unlist(lapply(nonignoredSeries,
                         .estimateLogisticSlope,
						 rppa@data$Steps,
                         rppa@data$Series.Id,
                         lindata))
	
    ## Force a nonzero slope for starting estimate
    minslope <- 0.001

    ## Filter out zero slopes from blank sample
    dat <- dat[dat > minslope]

    ## Use mean, not median. We may have many blank samples on a slide,
    ## in which case median slope will be 0
	## Gamma calculated only using nonignored series
    gamma <- mean(dat, trim=0.01)
    gamma <- max(gamma, minslope) # force nonzero slope for starting estimate

    ## For each (all) samples, estimate the offset
    passer <- rep(NA, length(sampleseries))
    layout <- rppa@data
    names(passer) <- sampleseries
    seriesdata <- lindata[rppa@tracking$fitToCurve]

    seriessteps <- layout$Steps[rppa@tracking$fitToCurve]
	
    for (this in sampleseries) {
        items <- rppa@data$Series.Id[rppa@tracking$fitToCurve] == this
        passer[names(passer) == this] <- median(seriesdata[items] / gamma - seriessteps[items],
                               na.rm=TRUE)
    }

    list(lBot=lBot,
         lTop=lTop,
         gamma=gamma,
         passer=passer)
}


##-----------------------------------------------------------------------------
## Returns model associated with key for invocation.
getRegisteredModel <- function(key) {
    return(classname <- getRegisteredObject(key,
                                            envir=modelenv(),
                                            "classname")$classname)
}


##-----------------------------------------------------------------------------
## Returns label associated with key for display by user interface.
getRegisteredModelLabel <- function(key) {
    return(ui.label <- getRegisteredObject(key,
                                           envir=modelenv(),
                                           "classname")$ui.label)
}


##-----------------------------------------------------------------------------
## Returns vector containing "keys" for all registered models.
getRegisteredModelKeys <- function() {
    keys <- getRegisteredObjectKeys(envir=modelenv())
    if (length(keys) == 0) {
        stop("no registered models exist")
    }

    return(keys)
}


##-----------------------------------------------------------------------------
## Registers specific model for use by RPPASPACE package.
## Register names of what code will consider "valid" models.
registerModel <- function(key,
                          classname,
                          ui.label=names(key)) {
    if (is.null(ui.label)) {
        ui.label <- key
    }
    ui.label <- as.character(ui.label)[1]
 
    ## Verify class
    tryCatch(new(classname),
             error=function(e) {
                 stop(sprintf("cannot create instance of classname %s",
                              sQuote(classname)))
             })

    if (!extends(classname, "FitClass")) {
       stop(sprintf("argument %s must be name of subclass of class %s",
                    sQuote("classname"),
                    sQuote("FitClass")))
    }

    registerClassname(key, classname, ui.label=ui.label, envir=modelenv())
}


##-----------------------------------------------------------------------------
RPPAFitFromParams <- function(rppa,
                              fitparams,
                              progmethod=NULL) {
  
    ## Check arguments
    if (!is.RPPA(rppa)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("rppa"), "RPPA"))
    }

    if (!is.RPPAFitParams(fitparams)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("fitparams"), "RPPAFitParams"))
    }

    if (!is.null(progmethod)) {
        if (!is.function(progmethod)) {
            stop(sprintf("argument %s must be function, if specified",
                         sQuote("progmethod")))
        }
    } else {
        ## Create a placeholder function
        progmethod <- function(phase) {}
    }

    ## :WORKAROUND: codetools (via R CMD check) complains unnecessarily about
    ## "no visible binding" as code below uses assign() rather than "<-".
    measure <- model <- xform <- method <- trim <- warnLevel <-
    ignoreNegative <- trace <- verbose <- veryVerbose <- ci <- NULL
	
	dilutionsInSeries <- length(unique(rppa@data$Dilution[rppa@data$Dilution > 0]))

    ## Create variables from 'fitparams' slots
    for (slotname in slotNames(fitparams)) {
        assign(slotname, slot(fitparams, slotname))
    }

    ## Need to make certain that the 'model' is a registered FitClass
    modelClass <- tryCatch(getRegisteredModel(model),
		error=function(e) {
			stop(sprintf("argument %s must be name of registered fit class",
				sQuote("model")))
		})

    ## Need to make sure that 'measure' refers to an actual data column
    if (missing("measure")) {
        stop("missing name of the measurement column to fit")
    }

    dn <- dimnames(rppa@data)[[2]]
    temp <- pmatch(measure, dn)
    if (is.na(temp)) {
        stop(paste("supply the name of a valid measurement column to fit:", measure, "is not valid.", sep=" " ))
    } else if (length(temp) > 1) {
        stop(sprintf("argument %s must identify unique column of argument %s",
                     sQuote("measure"),
                     sQuote("rppa")))
    }
    measure <- dn[temp]

    call <- match.call()
    intensity <- if (!is.null(xform)) {
                     xform(rppa@data[, measure])
                 } else {
                     rppa@data[, measure]
                 }

	silent <- warnLevel < 0

    ## Perform the first pass to initialize the estimates
    progmethod("firstpass")
    first <- .firstPass(intensity, rppa, ignoreNegative)
    passer    <- first$passer

    gamma.est <- first$gamma

    if (verbose) {
        cat(paste("Completed first pass. Parameters:",
                  paste("\t", "lBot =", first$lBot),
                  paste("\t", "lTop =", first$lTop),
                  paste("\t", "G =", first$gamma),
                  "\n",
                  sep="\n"))
        flush.console()
    }

    ## Put our current guess at the x and y values into vectors
    curveyval <- intensity[rppa@tracking$makePartOfCurve]
	fityval <- intensity[rppa@tracking$fitToCurve]
    ## Create new class
    fc <- new(modelClass)

    ## Do a two pass estimation, first using rough conc. estimates,
    ## then using better ones
    curvesteps <- rppa@data$Steps[rppa@tracking$makePartOfCurve]
    fitsteps <- rppa@data$Steps[rppa@tracking$fitToCurve]

    series <- seriesNames(rppa)
    curvesamplenames <- rppa@data$Series.Id[rppa@tracking$makePartOfCurve]
    fitsamplenames <- rppa@data$Series.Id[rppa@tracking$fitToCurve]
	noise.calc <- as.numeric(NA)

    for (pass in seq_len(2)) {

        pass.name <- if (pass == 1) "coarse" else "fine"
        curvexval <- if (pass == 1) {
                    curvesteps + passer[curvesamplenames]
                } else {
                    curvesteps + pass2[curvesamplenames]
                }

        fitxval <- if (pass == 1) {
                    fitsteps + passer[fitsamplenames]
                } else {
                    fitsteps + pass2[fitsamplenames]
                }

        ## Fit a response curve for the slide of the form yval = f(xval)
        progmethod(sprintf("%s fit slide", pass.name))

        fc <- fitSlide(fc,
                       conc=curvexval,
                       intensity=curveyval,
                       method=method)

        ## Conditional on the response curve fit for the slide
        ## Perform a separate fit of EC50 values for each dilution series.

        pass2 <- rep(NA, length(series))
        names(pass2) <- series

        allResid <- matrix(NA, nrow=length(series), ncol=dilutionsInSeries)
        rownames(allResid) <- series

        ss.ratio <- pass2
        warn2  <- rep("", length(series))
        names(warn2) <- series
        series.len <- length(series)
        report.incr <- as.integer(5)  
        i.this <- as.integer(1)

        progmethod(sprintf("%s fit series", pass.name))

		#All sample series, not just unignored ones
        for (this in series) {
            items <- fitsamplenames == this

            ## Report in 5 percent increments rather than iterations (speed)
            percent.this <- as.integer((i.this / series.len) * 100)
            if (!(percent.this %% report.incr)) {
                ## Skip 0 and 100 when reporting
                if (percent.this %% as.integer(100)) {
                    ## Notify progress update
                    progmethod(sprintf("%s fit series (%d%%)",
                                       pass.name,
                                       percent.this))
                }
            }
			
			#Fit one sample series of points to the curve generated from nonignored series
            fs <- fitSeries(fc,
                            diln=fitsteps[items],
                            intensity=fityval[items],
                            est.conc=passer[names(passer) == this],
                            method=method,
                            silent=silent,
                            trace=trace)
            pass2[names(pass2) == this] <- fs$est.conc
            warn2[names(warn2) == this] <- fs$warn

            for (i in 1:dilutionsInSeries) {
                allResid[rownames(allResid)==this, i] <- signif(fs$resids[i], 7)
            }

            ## Compute R^2 as sum(r[i]^2) / sum((y[i]-mean(y))^2),
            ## the fraction of variance explained for this series
            resids <- signif(fs$resids, 7)
            sse <- signif(sum(resids*resids),7)
            sst <- signif(var(fityval[items]) * (length(fityval[items])-1), 7)
            rsquared <- 1 - sse/sst
            ss.ratio[names(ss.ratio) == this] <- signif(rsquared, 7)
        }

        progmethod(sprintf("%s fit series complete", pass.name))

        if (verbose) {
            cat("Finished estimating EC50 values. Coefficients:", "\n")
            print(summary(pass2))
            cat("SS Ratio:", "\n")
            print(summary(ss.ratio))
            flush.console()
        }
		
		##--------------------- NOISE --------------------------
		## If there are any noise points defined,
		## then fit those points to the curves calculated above

		if (any(rppa@tracking$isNoise == TRUE)) {
			noise.calc <- .calculateNoise(intensity, rppa, dilutionsInSeries, fc, ignoreNegative, method, silent, trace, verbose)
		}
    }

    ## Create new class
    result <- new("RPPAFit",
                  call=call,
                  rppa=rppa,
                  measure=measure,
                  method=method,
                  trimset=c(lo.intensity=-100000,
                            hi.intensity=100000,
                            lo.conc=-1000,
                            hi.conc=1000),
                  model=fc,
                  concentrations=pass2,
                  lower=pass2,
                  upper=pass2,
                  intensities=signif(fitted(fc, pass2), 7),
                  ss.ratio=ss.ratio,
                  conf.width=0,
				  noise=noise.calc,
                  warn=warn2,
				  residualsrotation=fitparams@residualsrotation,
                  version=packageDescription("RPPASPACE", fields="Version"))
	  if (trim > 0) {
        if (verbose) {
            cat("Trimming concentrations...", "\n")
            flush.console()
        }

        progmethod("trim")
        tc <- tryCatch({
            trimConc(fc,
                conc=signif(fitted(result, "X"),7),
                intensity=intensity,
				steps=rppa@data$Steps[rppa@data$Spot.Type %in% spottype.sample],
                trimLevel=trim,
				antibody=rppa@antibody
				)
            },
                error=function(e) {
                    message(conditionMessage(e))
                    list(x = NaN, fm = NaN, iter = NaN)
            })

        if (any(is.nan(c(tc$lo.conc, tc$hi.conc)))) {
            warning("concentration values zero'd out due to NaNs",
                    immediate.=TRUE)
            series.conc <- rep(0, length(result@concentrations))
            names(series.conc) <- names(result@concentrations)
        } else {
            series.conc <- result@concentrations
            series.conc[series.conc < tc$lo.conc] <- tc$lo.conc
            series.conc[series.conc > tc$hi.conc] <- tc$hi.conc
        }

        minunique <- 5  ## :TBD: Magic# for undetermined limit
        nunique <- length(unique(sort(series.conc)))
        if (!(nunique > minunique)) {
            warning(sprintf("trim level %s is too high: #unique conc. values=%d",
                            as.character(trim),
                            nunique),
                    immediate.=TRUE)
        }

        result@concentrations <- series.conc
        result@trimset <- unlist(tc)
    }

    ## Should confidence intervals for the estimates be computed?
    if (ci) {
        if (verbose) {
            cat("Computing confidence intervals...", "\n")
            flush.console()
        }
        result <- getConfidenceInterval(result, progmethod=progmethod)
    }

    result
}

## Calculate noise on slide using positive control points
## specified as noise or posctrl-noise points in the slide file.
.calculateNoise <- function(intensity, rppa, dilutionsInSeries, fc, 
		ignoreNegative, method, silent, trace, verbose) {

	isNoise = rppa@tracking$isNoise
	noise.allSeries <- as.character(rppa@data$Series.Id[isNoise])
	noise.series <- unique(as.character(rppa@data$Series.Id[isNoise]))
	noise.samplenames <- as.character(rppa@data$Series.Id[isNoise])
	noise.steps <- rppa@data$Steps[isNoise]

	nval <- intensity[isNoise]
	
	if (ignoreNegative) {
		nval[nval < 0] <- NA
	}

	## Compute a robust estimate for alpha and beta
	noise.lBot <- min(nval, na.rm=TRUE)
	noise.lTop <- max(nval, na.rm=TRUE)

	noise.lindata <- .calcLogitz(nval, noise.lBot, noise.lTop-noise.lBot, 1)

	logisticSlope <- rep(NA, length(noise.series))
	names(logisticSlope) <- noise.series
	
	#Estimate logistic slope for noise series
	for (noise.ser in noise.allSeries) {
		noise.items <- noise.allSeries == noise.ser
		noise.tempSteps <- noise.steps[noise.items]
		noise.x <- noise.lindata[noise.items]
		noise.y <- noise.x[!is.na(noise.x) & !is.infinite(noise.x)]
	
		if (!all(is.na(noise.y))) {
            ## ychange / xchange
            logisticSlope[noise.ser] <- (max(noise.y) - min(noise.y)) / (max(noise.tempSteps) - min(noise.tempSteps))
        }
	}

	## Force a nonzero slope for starting estimate
    minslope <- 0.001

    ## Filter out zero slopes from blank sample
    logisticSlope <- logisticSlope[logisticSlope > minslope]

    ## Use mean, not median. We may have many blank samples on a slide,
    ## in which case median slope will be 0.  
	## This shouldn't apply to control points, but better to be consistent with the rest of the slide.
    ngamma <- mean(logisticSlope, trim=0.01)
    ngamma <- max(ngamma, minslope) # force nonzero slope for starting estimate

    ## For each sample, estimate the offset
    noise.passer <- rep(NA, length(noise.series))
    names(noise.passer) <- noise.series

    for (noise.ser in noise.allSeries) {
		noise.items <- noise.allSeries == noise.ser
        noise.passer[noise.ser] <- median(noise.lindata[noise.items] / ngamma - noise.steps[noise.items],
                               na.rm=TRUE)
    }

	noise.pass2 <- rep(NA, length(noise.series))
	names(noise.pass2) <- noise.series

	noise.allResid <- matrix(NA, nrow=length(noise.series), ncol=dilutionsInSeries)
	rownames(noise.allResid) <- noise.series
	
	noise.warn2 <- rep("", length(noise.series))

	## Calculate fit for Noise points
	for (noise.ser in noise.series) {
		noise.items <- noise.samplenames == noise.ser

		noise.fs <- fitSeries(fc,
						diln=noise.steps[noise.items],
						intensity=nval[noise.items],
						est.conc=noise.passer[noise.ser],
						method=method,
						silent=silent,
						trace=trace)
		noise.pass2[noise.ser] <- noise.fs$est.conc
		noise.warn2[noise.ser] <- noise.fs$warn
	}
	
	if (verbose) {
		cat("Finished estimating noise EC50 values. Coefficients:", "\n")
		print(summary(noise.pass2))
		cat("warnings:", "\n")
		print(paste(noise.warn2, collapse="; "))
		flush.console()
	}
	return(noise.pass2)
}


##-----------------------------------------------------------------------------
## Extracts model coefficients from objects returned by modeling functions.
setMethod("coef", signature(object="RPPAFit"),
          function(object,
                   ...) {
    callGeneric(object@model)
})


setMethod("coefficients", signature(object="RPPAFit"),
          function(object,
                   ...) {
    coef(object, ...)
})

