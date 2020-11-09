###
### $Id: rs3-0-dilutionFitMethods.R
###


##=============================================================================
setClass("FitClass",
         representation("VIRTUAL"))

##=============================================================================
setClass("LogisticFitClass",
         contains="FitClass",
         representation(coefficients="numeric"),
         prototype=prototype(coefficients=c(alpha=0, beta=0, gamma=0)))

##=============================================================================
setOldClass("cobs")
setClass("CobsFitClass",
         contains="FitClass",
         representation(model="cobs",
                        lambda="numeric"),
         prototype=prototype(lambda=0))

##=============================================================================
setOldClass("loess")
setClass("LoessFitClass",
         contains="FitClass",
         representation(model="loess"))


##-----------------------------------------------------------------------------
## Returns TRUE if class of argument is subclass of FitClass.
is.FitClass <- function(x) {
    extends(class(x), "FitClass")
}


####################################################################
## GENERIC METHODS FOR FitClass: Typically throw an error since they
## must be implemented by derived classes.

##-----------------------------------------------------------------------------
## Finds the concentration for an individual dilution series given the
## curve fit for the slide
##
## Inputs
## dilutions and intensities for a single dilution series
## est.conc = starting estimated concentration for dilution = 0
##
## Outputs
## est.conc = estimated concentration for dilution = 0
##
setMethod("fitSeries", signature(object="FitClass"),
          function(object,
                   diln,
                   intensity,
                   est.conc,
                   method="nls",
                   silent=TRUE,
                   trace=FALSE,
                   ...) {
    stop(sprintf("%s method must be implemented by any subclass of %s",
                 sQuote("fitSeries"),
                 sQuote("FitClass")))
})


##-----------------------------------------------------------------------------
## Use the conc and intensity series for an entire slide to
## fit a curve for the slide of intensity = f(conc)
setMethod("fitSlide", signature(object="FitClass"),
          function(object,
                   conc,
                   intensity,
                   ...) {
    stop(sprintf("%s method must be implemented by any subclass of %s",
                 sQuote("fitSlide"),
                 sQuote("FitClass")))
})


##-----------------------------------------------------------------------------
setMethod("fitted", signature(object="FitClass"),
          function(object,
                   conc,
                   ...) {
    stop(sprintf("%s method must be implemented by any subclass of %s",
                 sQuote("fitted"),
                 sQuote("FitClass")))
})


##-----------------------------------------------------------------------------
## Returns concentration and intensity cutoffs for the model
setMethod("trimConc", signature(object="FitClass"),
          function(object,
                   conc,
                   intensity,
				   steps,
                   trimLevel,
				   antibody, 
                   ...) {
    stop(sprintf("%s method must be implemented by any subclass of %s",
                 sQuote("trimConc"),
                 sQuote("FitClass")))
})


##-----------------------------------------------------------------------------
## Extracts model coefficients from objects returned by modeling functions.
## N.B.: Should be overridden by classes that have coefficients!
setMethod("coef", signature(object="FitClass"),
          function(object,
                   ...) {
    NULL
})


###################################################################
## Utility functions used to implement methods for derived classes
## :KRC: Should these be used for the FitClass method so we can
## both document them and use them if we decide to develop new
## and improved fitting algorithms in the future?


##-----------------------------------------------------------------------------
.slide.model <- function(conc) {
    ## Check arguments
    stopifnot(is.numeric(conc))

    ## Begin processing

    ## :TODO: Come up with another way to do this that doesn't involve
    ## writing into user workspace.
#    obj <- get(".RPPA.fit.model", envir=.GlobalEnv)
    obj <- get(".RPPA.fit.model", envir=RPPASPACE_Temp_Env)
    fitted(obj, conc)
#    fitted(.RPPA.fit.model, conc)
}


##-----------------------------------------------------------------------------
## Fit the dilution series to the model for the slide
.series.fit <- function(object,
                        diln,   #dilution step values (0, -1, -2, -3, -4) for 5 point dilution series
                        intensity,
                        est.conc,
                        method=c("nls", "nlrob", "nlrq"),
                        silent=TRUE,
                        trace=FALSE,
                        ...) {
    ## Check arguments
    stopifnot(is.FitClass(object))
    stopifnot(is.numeric(diln))
    stopifnot(is.numeric(intensity))
    stopifnot(is.numeric(est.conc))
    stopifnot(is.logical(silent) && length(silent) == 1)
    stopifnot(is.logical(trace) && length(trace) == 1)
    method <- match.arg(method)

    ## Begin processing

    ## Ensure necessary packages available
    if (method == "nlrob") {
        if (!requireNamespace("robustbase")) {
            stop(sprintf("%s package required for %s method",
                         sQuote("robustbase"),
                         sQuote(method)))
        }
    } else if (method == "nlrq") {
        if (!requireNamespace("quantreg")) {
            stop(sprintf("%s package required for %s method",
                         sQuote("quantreg"),
                         sQuote(method)))
        }
    }

    ## Define regression method
    nlsmeth <- switch(EXPR=method,
                      nls=stats::nls,
                      nlrob=robustbase::nlrob,
                      nlrq=function(...) {
                          params <- quantreg::nlrq.control(maxiter=10,
                                                           eps=1e-02)
                          quantreg::nlrq(control=params, ...)
                      },
                      stop(sprintf("unrecognized regression method %s",
                                   sQuote(method))))

    ## Function .slide.model references object back here for curve model
    ## :NOTE: Writing into global environment is considered rude.
#    assign(".RPPA.fit.model", object, envir=.GlobalEnv)
    assign(".RPPA.fit.model", object, envir=RPPASPACE_Temp_Env)
	#.RPPA.fit.model <- object

    tmp <- try({
					#nls doesn't work correctly if values fit model exactly
					#nlrob and nlrq both involve randomness and without setting
					#the seed ahead of time, might give different values each run. (JMM)
                    #nlsmeth(Y ~ RPPASPACE:::.slide.model(Steps+X),

                    nlsmeth(Y ~ .slide.model(Steps+X),
                            data=data.frame(Y=intensity,
                                            Steps=diln),
                            start=list(X=est.conc),
                            trace=trace)
               },
               silent=silent)

    if (inherits(tmp, "try-error")) {
        warn <- "unavoidable nls/rlm error"
        ## :TBD: Should 'est.conc' be set to something different on error?
        resids <- 0
        if (!silent) {
            warning(warn)
        }
    } else {
        ## Model fitting succeeded, so we can continue
        warn <- ""
        est.conc <- coef(tmp)
        resids <- residuals(tmp)
    }

    list(warn=warn,
         est.conc=est.conc,
         resids=resids)
}


##-----------------------------------------------------------------------------
## Returns estimate of background noise.
.est.bg.noise <- function(object,
                          conc,
                          intensity,
                          trimLevel) {
    ## Check arguments
    stopifnot(is.FitClass(object))
    stopifnot(is.numeric(conc))
    stopifnot(is.numeric(intensity))
    stopifnot(is.numeric(trimLevel) && length(trimLevel) == 1)

    ## Begin processing

    ## Trim the concentration estimates to bound lower and upper
    ## concentration estimates at the limits of what can be detected
    ## given our background noise.
    r <- fitted(object, conc) - intensity  # residuals
    sMad <- mad(r, na.rm=TRUE)

    ## Use trim level * (median absolute deviation of residuals)
    ## as estimate for background noise

    trim <- trimLevel * sMad
}


##-----------------------------------------------------------------------------
.generic.trim <- function(object,
                          conc,
                          intensity,
						  steps,
                          trimLevel,
						  antibody, 
                          ...) {
    ## Check arguments
    stopifnot(is.FitClass(object))
    stopifnot(is.numeric(conc))
    stopifnot(is.numeric(intensity))
    stopifnot(is.numeric(trimLevel) && length(trimLevel) == 1)

    ## Begin processing
    trim <- .est.bg.noise(object, conc, intensity, trimLevel)

    ## Determine high and low intensities
    lBot <- quantile(intensity, probs=0.01, na.rm=TRUE)
    lTop <- quantile(intensity, probs=0.99, na.rm=TRUE)
    lo.intensity <- lBot + trim
    ## In practice, we rarely see the response "top out".
    ## Do not trim at the top end.
    # hi.intensity <- lTop - trim
    hi.intensity <- max(intensity)

    ## Determine high and low concentrations

    ## Search fitted model to find conc corresponding to lo.intensity
    print(paste("Antibody ", antibody, " min conc=", min(conc, na.rm=TRUE), ", max conc=", max(conc, na.rm=TRUE), sep=""))
	lo.conc  <- tryCatch({
		bisection.search(
			min(conc, na.rm=TRUE),
			max(conc, na.rm=TRUE),
			function(x, object) {
				fitted(object, x) - lo.intensity
			},
			f.extra=object,
			tol=0.1)$x
		},
			error=function(e) {
				msg <- paste("Antibody (", antibody, ") Unable to calculate low concentration value. Setting value to NA", sep="")
				warning(msg)
				message(msg)
				NA
	})

    ## Adjust min allowable conc to point at undiluted spot
	if (!is.na(lo.conc)) {
		max.step <- max(steps)
		lo.conc <- lo.conc - max.step
	}
	
	
	hi.conc  <- tryCatch({
		bisection.search(
			min(conc, na.rm=TRUE),
			max(conc, na.rm=TRUE),
			function(x, object) {
				fitted(object, x) - hi.intensity
			},
			f.extra=object,
			tol=0.1)$x
		},
			error=function(e) {
				msg <- paste("Antibody (", antibody, ") Unable to calculate high concentration value. Setting value to NA", sep="")
				warning(msg)
				message(msg)
				NA
	})

    ## Adjust max allowable conc to point at most dilute spot
    if (!is.na(hi.conc)) {
		min.step <- min(steps)
		hi.conc <- hi.conc - min.step
	}
	
    list(lo.intensity=lo.intensity,
         hi.intensity=hi.intensity,
         lo.conc=lo.conc,
         hi.conc=hi.conc,
         level=trimLevel)
}


##
## Loess model class
##

##-----------------------------------------------------------------------------
setMethod("fitSlide", signature(object="LoessFitClass"),
          function(object,
                   conc,
                   intensity,
                   ...) {
    model <- loess(intensity ~ conc)

    ## Create new class
    new("LoessFitClass",
        model=model)
})


##-----------------------------------------------------------------------------
setMethod("fitted", signature(object="LoessFitClass"),
          function(object,
                   conc,
                   ...) {
    model <- object@model

    ## loess will not interpolate beyond the initial fitted conc. range
    lo <- min(model$x)
    conc <- pmax(min(model$x), conc)
    conc <- pmin(max(model$x), conc)
    conc.pred <- conc
    conc.pred[is.na(conc)] <- lo

    intensity <- predict(model, data.frame(conc=conc.pred))
    intensity[is.na(conc)] <- NA

    intensity
})


##-----------------------------------------------------------------------------
setMethod("fitSeries", signature(object="LoessFitClass"),
          function(object,
                   diln,
                   intensity,
                   est.conc,
                   method="nls",
                   silent=TRUE,
                   trace=FALSE,
                   ...) {
    .series.fit(object, diln, intensity, est.conc, method, silent, trace)
})


##-----------------------------------------------------------------------------
## Trim level default based on trying various cutoff levels on multiple slides.
setMethod("trimConc", signature(object="LoessFitClass"),
          function(object,
                   conc,
                   intensity,
				   steps,
                   trimLevel=2,  # arbitrary based on experimentation
				   antibody,
                   ...) {
    .generic.trim(object, conc, intensity, steps, trimLevel, antibody, ...)
})


##
## Cobs model class
##

##-----------------------------------------------------------------------------
setMethod("fitSlide", signature(object="CobsFitClass"),
          function(object,
                   conc,
                   intensity,
                   ...) {
    if (!requireNamespace("cobs")) {
        stop(sprintf("%s package required for %s method",
                     sQuote("cobs"),
                     sQuote("fitSlide")))
    }

    model <- cobs(conc,
                  intensity,
                  constraint="increase",
                  nknots=20,
                  lambda=object@lambda,
                  degree=2,
                  tau=0.5,
                  print.warn=FALSE,
                  print.mesg=FALSE,
				  maxiter=300
				  )

	## Create new class
    new("CobsFitClass",
        model=model,
        lambda=model$lambda)
})


##-----------------------------------------------------------------------------
.predict.spline <- function(xvec,
                            aknot,
                            acoef) {
    ## Check arguments
    stopifnot(is.numeric(xvec))
    stopifnot(is.numeric(aknot))
    stopifnot(is.numeric(acoef))

    ## Begin processing
    aknot1 <- aknot[1]
    aknotn <- aknot[length(aknot)]
    aknotnew <- c(rep(aknot1, 2),
                  aknot,
                  rep(aknotn, 2))

    adj <- 1e-8
    xvec[xvec < (aknot1 + adj)] <- aknot1 + adj
    xvec[xvec > (aknotn - adj)] <- aknotn - adj

    a <- splines::spline.des(aknotnew, xvec, ord=3)
    fvalvec <- (a$design) %*% acoef

    return(as.vector(fvalvec))
}


##-----------------------------------------------------------------------------
setMethod("fitted", signature(object="CobsFitClass"),
          function(object,
                   conc,
                   ...) {
    model <- object@model

    ## Predict missing values at min intensity
    lo <- min(model$x)
    conc.pred <- conc
    conc.pred[is.na(conc)] <- lo

    ## :TODO: Add argument to enable Jianhua's code, or remove it
    intensity <- if (TRUE) {
		 ## predict.cobs is irritating
		 ## It returns predicted values after sorting on the input
		 ## vector. So we get intensity ~ sort(fit)
		 ## We need to undo this sort to find predictions using the
		 ## original concentration ordering
		 n <- length(conc.pred)
		 if (n > 1) {

			 ## Undo sort on fit
			 o <- sort.list(conc.pred,
							method="quick",
							na.last=NA)
			 u <- rep(NA, n)
			 u[o] <- seq_along(conc.pred)
			 cobs.intensity <- predict(model, conc.pred[o])[, "fit"]
			 cobs.intensity[u]
		 } else {
			 ## Only one data point
			 predict(model, conc.pred)[, "fit"]
		 }
	 } else {
		 if (!requireNamespace("splines")) {
			 stop(sprintf("%s package required for %s method",
						  sQuote("splines"),
						  sQuote("fitted")))
		 }

		 ## The above sort and unsort process is yucky and a bit
		 ## slow. Jianhua did not use the cobs predict method and
		 ## instead evaluates the spline directly. Unfortunately,
		 ## there seems to be a bug in .predict.spline where it does
		 ## not have the correct number of coefficients sometimes.

		 .predict.spline(conc.pred,
						 model$knots,
						 model$coef)
	 }

    intensity[is.na(conc)] <- NA
	
    intensity
})


##-----------------------------------------------------------------------------
setMethod("fitSeries", signature(object="CobsFitClass"),
          function(object,
                   diln,
                   intensity,
                   est.conc,
                   method="nls",
                   silent=TRUE,
                   trace=FALSE,
                   ...) {
    .series.fit(object, diln, intensity, est.conc, method, silent, trace)
})


##-----------------------------------------------------------------------------
## Trim level default based on trying various cutoff levels on multiple slides.
setMethod("trimConc", signature(object="CobsFitClass"),
          function(object,
                   conc,
                   intensity,
				   steps,
                   trimLevel=2,  # arbitrary based on experimentation
				   antibody,
                   ...) {
    .generic.trim(object, conc, intensity, steps, trimLevel, antibody, ...)
})


##
## Logistic model class
##

##-----------------------------------------------------------------------------
## N.B.: rnls does not work with local functions
.ea <- function(x) {
    exp(x) / (1 + exp(x))
}


##-----------------------------------------------------------------------------
.coef.quantile.est <- function(intensity) {
    ## Check arguments
    stopifnot(is.numeric(intensity))

    ## Begin processing
    lBot <- quantile(intensity, probs=0.05, na.rm=TRUE)
    lTop <- quantile(intensity, probs=0.95, na.rm=TRUE)
    p.alpha <- lBot
    p.beta  <- lTop - lBot
    p.gamma <- log(2)  # Assume linear response on log2 scale as first guess
    list(alpha=p.alpha,
         beta=p.beta,
         gamma=p.gamma)
}


##-----------------------------------------------------------------------------
setMethod("fitSlide", signature(object="LogisticFitClass"),
          function(object,
                   conc,
                   intensity,
                   ...) {
    cf <- as.list(coef(object))
	
	magicOffset <- 5000

    if (cf$gamma == 0) {
        ## Initialize coefficients
        cf <- .coef.quantile.est(intensity)
    }

    nls.model <- try(nls(yval ~ log(alpha + beta*.ea(gamma*xval) + magicOffset),
                         data=data.frame(xval=conc,
                                         yval=log(intensity + magicOffset)),
                         start=list(alpha=cf$alpha,
                                    beta=cf$beta,
                                    gamma=cf$gamma),
                         control=nls.control(maxiter=100),
                         na.action="na.omit"))
    if (inherits(nls.model, "try-error")) {
        warning("unable to perform first pass overall slide fit. trying quantiles.")
        ## Crude (but robust) way to get alpha and beta
        cf <- .coef.quantile.est(intensity)
    } else {
        cf <- coef(nls.model)
        ## :TBD: Why is this done? Seemingly creates unnecessary list...
        cf <- list(alpha=cf["alpha"],
                   beta= cf["beta"],
                   gamma=cf["gamma"])
    }

    ## Create new class
    new("LogisticFitClass",
        coefficients=unlist(cf))
})


##-----------------------------------------------------------------------------
setMethod("fitted", signature(object="LogisticFitClass"),
          function(object,
                   conc,
                   ...) {
    cf <- as.list(object@coefficients)
    cf$alpha + cf$beta * .ea(cf$gamma * conc)
})


##-----------------------------------------------------------------------------
setMethod("fitSeries", signature(object="LogisticFitClass"),
          function(object,
                   diln,
                   intensity,
                   est.conc,
                   method="nls",
                   silent=TRUE,
                   trace=FALSE,
                   ...) {
    .series.fit(object, diln, intensity, est.conc, method, silent, trace)
})


##-----------------------------------------------------------------------------
## Trim level default based on trying various cutoff levels on multiple slides.
setMethod("trimConc", signature(object="LogisticFitClass"),
          function(object,
                   conc,
                   intensity,
				   steps,
                   trimLevel=2,  # arbitrary based on experimentation
				   antibody,
                   ...) {
    cf <- as.list(object@coefficients)
    noise <- .est.bg.noise(object, conc, intensity, trimLevel)
    trim <- noise / cf$beta

    if (trim <= 0 || trim >= 1) {
        warning(sprintf("Antibody (%s) trimConc: trim should be in interval (0, 1): trim=%s",
						antibody,
                        trim),
                immediate.=TRUE)
    }

    ## Determine high and low intensities
    lo.intensity <- cf$alpha + cf$beta * trim
    hi.intensity <- cf$alpha + cf$beta - cf$beta * trim

    ## Determine high and low concentrations
    max.step <- max(steps)
    min.step <- min(steps)

    if (!requireNamespace("boot")) {
        stop(sprintf("%s package required for %s method",
                     sQuote("boot"),
                     sQuote("trimConc")))
    }

    lo.logit <- tryCatch(boot::logit(trim),
                         error=function(cond) {
                             warning(sprintf("Antibody (%s) logit: %s: p=%f, odds=%f",
											 antibody, 
                                             conditionMessage(cond),
                                             p <- trim,
                                             p / (1 - p)),
                                     immediate.=TRUE)
                             NaN
                         })

    hi.logit <- tryCatch(boot::logit(1-trim),
                         error=function(cond) {
                             warning(sprintf("Antibody (%s) logit: %s: p=%f, odds=%f",
											 antibody, 
                                             conditionMessage(cond),
                                             p <- 1-trim,
                                             p / (1 - p)),
                                     immediate.=TRUE)
                             NaN
                         })

    lo.conc <- lo.logit / cf$gamma - max.step
    hi.conc <- hi.logit / cf$gamma - min.step

    list(lo.intensity=lo.intensity,
         hi.intensity=hi.intensity,
         lo.conc=lo.conc,
         hi.conc=hi.conc,
         level=trimLevel)
})


##-----------------------------------------------------------------------------
## Extracts model coefficients from LogisticFitClass.
setMethod("coef", signature(object="LogisticFitClass"),
          function(object,
                   ...) {
    object@coefficients
})

