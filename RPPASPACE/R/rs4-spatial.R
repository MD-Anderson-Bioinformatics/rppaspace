###
### $Id: rs4-spatial.R
###


##=============================================================================
setClass("RPPASpatialParams",
         representation(cutoff="numeric",
                        k="numeric",
                        gamma="numeric",
                        plotSurface="logical"))


##-----------------------------------------------------------------------------
is.RPPASpatialParams <- function(x) {
    is(x, "RPPASpatialParams")
}


##-----------------------------------------------------------------------------
RPPASpatialParams <- function(cutoff=0.8,
                              k=100,
                              gamma=0.1,
                              plotSurface=FALSE) {
    ## Check arguments
    stopifnot(is.numeric(cutoff) && length(cutoff) == 1)
    stopifnot(is.numeric(k) && length(k) == 1)
    stopifnot(is.numeric(gamma) && length(gamma) == 1)
    stopifnot(is.logical(plotSurface) && length(plotSurface) == 1)

    ## Create new class
    new("RPPASpatialParams",
        cutoff=cutoff,
        k=k,
        gamma=gamma,
        plotSurface=plotSurface)
}


##-----------------------------------------------------------------------------
## Invoked by validObject() method.
validSpatialParams <- function(object) {

    #cat("validating", class(object), "object", "\n")
    msg <- NULL

    ## Validate cutoff slot
    {
        cutoff <- object@cutoff
        min.cutoff <- 0
        max.cutoff <- 1

        if (!(cutoff >= min.cutoff && cutoff <= max.cutoff)) {
            msg <- c(msg, sprintf("cutoff must be in interval [%d, %d]",
                                  min.cutoff, max.cutoff))
        }
    }

    ## Validate k slot
    {
        k <- object@k                   # [passed directly to mgcv::s()]
        min.k <- 2.0

        ## Ensure k value acceptable [passed directly to mgcv::s()] 
        if (!(is.finite(k))) {
            msg <- c(msg, "k must be finite")
        } else if (!(k >= min.k)) {
            msg <- c(msg, sprintf("k must be greater than or equal %d",
                                  min.k))
        }
    }

    ## Validate gamma slot
    {
        gamma <- object@gamma           # [passed directly to mgcv::gam()]
        min.gamma <- 0
        max.gamma <- 2

        if (!(is.finite(gamma))) {
            msg <- c(msg, "gamma must be finite")
        } else if (!(gamma >= min.gamma && gamma <= max.gamma)) {
            msg <- c(msg, sprintf("gamma must be in interval [%d, %d]",
                                  min.gamma, max.gamma))
        }
    }

    ## Validate plotSurface slot
    if (!is.logical(object@plotSurface)) {
        msg <- c(msg, "plotSurface must be logical")
    }

    ## Pass or fail?
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
}

setValidity("RPPASpatialParams", validSpatialParams)


##-----------------------------------------------------------------------------
## Returns a string representation of this instance. The content and format of
## the returned string may vary between versions. Returned string may be
## empty, but never null.
setMethod("paramString", signature(object="RPPASpatialParams"),
          function(object,
                   slots=slotNames(object),
                   ...) {
    ## Check arguments
    stopifnot(is.character(slots) && length(slots) >= 1)

    ## :TODO: Implementation currently ignores the 'slots' argument
    ## and returns string containing parameters from various slots.
    ## as though:
    ##     slotsToDisplay <- c("cutoff", "k", "gamma", "plotSurface")
    ##     paramString(sp, slotsToDisplay)
    ##
    paste(paste("cutoff:", object@cutoff), "\n",
          paste("k:", object@k), "\n",
          paste("gamma:", object@gamma), "\n",
          paste("plotSurface:", object@plotSurface), "\n",
          sep="")
})


##-----------------------------------------------------------------------------
## Computes the background cutoff of a slide.
.computeBackgroundCutoff <- function(mydata, tracking, measure, cutoff) {
    ## Check arguments
    stopifnot(is.data.frame(mydata))
    stopifnot(is.character(measure) && length(measure) == 1)
    stopifnot(is.numeric(cutoff) && length(cutoff) == 1)

    ## Begin processing
    ## Find all negative controls
    negcon <- tracking$isNegCtrl

    ## Identify noise region
    bg <- if (any(negcon)) {
              mydata[[measure]][negcon]
          } else {
              ## Use background to compute noise region
              mydata$Raw.Value - mydata$Net.Value
          }

    ## Compute background cutoff using the quantile of 'cutoff' argument
    bgCut <- quantile(bg, probs=cutoff)
	
    ## If computed background cutoff too low, use a larger quartile
    if (bgCut <= 100) {  
        bgCut <- quantile(bg, probs=0.99)
        if (bgCut <= 100) {  
            bgCut <- max(bg[-which.max(bg)])
			warning(paste("Final bgCut value was less than threshhold of 100.  Value=", bgCut, sep=" "))
        }
    }
    return(bgCut)
}


##-----------------------------------------------------------------------------
spatialCorrection <- function(rppa,
                              measure=c("Net.Value", "Raw.Value"),
                              cutoff=0.8,
                              k=100,
                              gamma=0.1,
                              plotSurface=FALSE) {
    ## Check arguments
    if (!is.RPPA(rppa)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("rppa"), "RPPA"))
    }

	reqdNames <- c("Spot.Type", "Dilution")

	if (!(all(reqdNames %in% colnames(rppa@data)))) {
		missingNames <- reqdNames[!reqdNames %in% colnames(rppa@data)]
		stop(sprintf(ngettext(length(missingNames),
							  "argument %s missing required column: %s",
							  "argument %s missing required columns: %s"),
					 sQuote("rppa"), paste(missingNames, collapse=", ")))
	}

	# Verify slide is the same length as the tracking information
	expectedRowCount <- TRUE
	if (nrow(rppa@tracking) != nrow(rppa@data)) {
		msg <- paste( "The number of rows (", nrow(rppa@data), ") in slide ", rppa@antibody, 
			" do not match the number of rows (", nrow(rppa@tracking),
			") in the first slide. All slides in a run are required to have the same number of entries.", sep="")
		warning(msg)
		expectedRowCount <- FALSE
	}

	if (!(any(rppa@data$Spot.Type %in% spottype.positivecontrol))) {
		stop("contains no positive controls")
	}

    measure <- match.arg(measure)
    if (!(measure %in% colnames(rppa@data))) {
        stop(sprintf("argument %s missing column for measure %s",
                     sQuote("rppa"), sQuote(measure)))
    }

    adjMeasure <- paste("Adj", measure, sep=".")
    if (adjMeasure %in% colnames(rppa@data)) {
        ## Don't allow doing this again as the merge will screw up...
        warning(sprintf("argument %s has already been spatially corrected for measure %s",
                        sQuote("rppa"), sQuote(measure)))
        return(rppa)
    }

    if (!is.numeric(cutoff)) {
        stop(sprintf("argument %s must be numeric",
                     sQuote("cutoff")))
    } else if (!(length(cutoff) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("cutoff")))
    } else if (!(cutoff >= 0 && cutoff <= 1)) {
        stop(sprintf("argument %s must be in interval [%d, %d]",
                     sQuote("cutoff"), 0, 1))
    }

    ## Valid value is >= 2. [passed directly to mgcv::s()]
    if (!is.numeric(k)) {
        stop(sprintf("argument %s must be numeric",
                     sQuote("k")))
    } else if (!(length(k) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("k")))
    } else if (!is.finite(k)) {
        stop(sprintf("argument %s must be finite",
                     sQuote("k")))
    }

    ## Valid range is probably [0..2] [passed directly to mgcv::gam()]
    if (!is.numeric(gamma)) {
        stop(sprintf("argument %s must be numeric",
                     sQuote("gamma")))
    } else if (!(length(gamma) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("gamma")))
    } else if (!(is.finite(gamma) && gamma > 0)) {
        stop(sprintf("argument %s must be a positive finite quantity",
                     sQuote("gamma")))
    }

    if (!is.logical(plotSurface)) {
        stop(sprintf("argument %s must be logical",
                     sQuote("plotSurface")))
    } else if (!(length(plotSurface) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("plotSurface")))
    }

	if (!expectedRowCount) {
	} else {
		## Begin processing
		if (!requireNamespace("mgcv")) {
			stop(sprintf("%s package required for fitting the GAM in the %s method",
						 sQuote("mgcv"),
						 sQuote("spatialCorrection")))
		}

		## Set up the row and column variables
		tracking <- rppa@tracking
		mydata <- rppa@data

		mydata$Row <- with(mydata, (Main.Row-1)*max(Sub.Row) + Sub.Row)
		mydata$Col <- with(mydata, (Main.Col-1)*max(Sub.Col) + Sub.Col)

		## Create data frame with row/column indices used for predicting surface
		pd <- data.frame(Row=mydata$Row, Col=mydata$Col)

		## Compute background cutoff
		bgCut <- .computeBackgroundCutoff(mydata, tracking, measure, cutoff)
		
		## Remove positive controls less than computed background cutoff 
		## by replacing measure values with NAs
		poscon <- tracking$isPosCtrl
		is.na(mydata[poscon, measure]) <- mydata[poscon, measure] < bgCut
		positives <- mydata[poscon, ]

		## :NOTE: This code assumes that only one positive control
		## series exists on a slide. If multiple exist, they would
		## be blended together by this method

		## Determine positive control dilutions
		dilutions <- sort(unique(positives$Dilution), decreasing=TRUE)
		ndilutions <- length(dilutions)

		## Create surface names
		surfaces <- sapply(dilutions,
						   function(dilution) {
							   paste("surface", dilution, sep="")
						   })

						   ## Fits a generalized additive model to estimate a surface
		## from positive controls
		fmla <- switch(EXPR=measure,
					   Net.Value = Net.Value ~ s(Row, Col, bs="ts", k=adjK),
					   Raw.Value = Raw.Value ~ s(Row, Col, bs="ts", k=adjK))

		surfacesPlottable <- TRUE
		for (dilution in dilutions) {
			pcsub <- positives[positives$Dilution == dilution, ]

			## Make choice of k robust in case that the number of
			## available spots is less than k (YH).
			adjK <- k
			spotCount <- sum(!is.na(pcsub[, measure]))
			if (spotCount < k) {
				adjK <- round(spotCount / 3)  ## arbitrary magic number
			}

			if(adjK < 3) {
				warning(paste("Insufficient number of positive control points to create surface for dilution", dilution, "on slide", rppa@antibody))
				surfacesPlottable <- FALSE
			}
			else
			{
				b1 <- mgcv::gam(fmla, data=pcsub,gamma=gamma)

				surface <- paste("surface", dilution, sep="")
				assign(surface, signif(mgcv::predict.gam(b1, newdata=pd), 7))

				remove(b1)
			}

		}

		## Plot the different surfaces
		if (plotSurface) {
			if (!surfacesPlottable) {
				warning("Unable to plot surfaces for ", rppa@antibody)
			} else {
				imageRPPA <- getMethod("image", class(rppa))
				temprppa <- rppa
				par(ask=TRUE)
				for (surface in surfaces) {
					try (
					{
						main <- .mkPlotTitle(measure,sprintf("%s [%s]", rppa@antibody, surface))
						print(paste(measure,sprintf("%s [%s]", rppa@antibody, surface)))
						temprppa@data[, surface] <- eval(as.name(surface))
						imageRPPA(temprppa, colorbar=TRUE, measure=surface, main=main)
					}
					)
				}
				par(ask=FALSE)
			}
		}

		## Constrain the surfaces so they do not cross
		if (ndilutions > 1) {
			if (surfacesPlottable) {
				for (i in seq(2, ndilutions)) {
					surface <- surfaces[i]
					prev.surface <- surfaces[i-1]
					s2 <- eval(as.name(surface))
					s1 <- eval(as.name(prev.surface))
					s2[s1 < s2] <- s1[s1 < s2] - 0.5
					assign(surface, s2)
				}
			}
		}
	}
    ##-------------------------------------------------------------------------
    ## For each spot on the array, finds the positive control surface with the
    ## most similar expression level.
    which.surface <- function(meas,
                              surfs) {
        stopifnot(is.numeric(meas) && length(meas) == 1)
        stopifnot(is.numeric(surfs) && length(surfs) >= 1)

        n <- length(surfs)
        x.surf <- NA     # index of most similar surface
        for (i in seq_len(n-1)) {
            if ((surfs[i] == meas) || (surfs[i] > meas && meas > surfs[i+1])) {
                x.surf <- i
                break
            }
		}

		if (is.na(x.surf)) {
			if (meas >= surfs[1]) {
				x.surf <- 0
			} else if (meas <= surfs[n]) {
				x.surf <- n
			}
		}
		
        return(x.surf)
    }


    ##-------------------------------------------------------------------------
    ## Returns the linear interpolation when the level of spot expression falls
    ## between two surfaces; otherwise, 0.
    getp <- function(x.surf,
                     meas,
                     surfs) {

        stopifnot(is.numeric(x.surf) && length(x.surf) == 1)
        stopifnot(is.numeric(meas) && length(meas) == 1)
        stopifnot(is.numeric(surfs) && length(surfs) >= 1)

        n <- length(surfs)

        stopifnot(x.surf >= 0 && x.surf <= n)
        p <- if (x.surf != 0 && x.surf != n) {
                 (surfs[x.surf] - meas) / (surfs[x.surf] - surfs[x.surf+1])
             } else {
                 0
             }
    }


    ##-------------------------------------------------------------------------
    ## Finds the surface or linear interpolation of two surfaces that will be
    ## used to perform the scaling. Returns the overall adjustment.
    getadj <- function(x.surf,
                       p,
                       adjs) {
        stopifnot(is.numeric(x.surf) && length(x.surf) == 1)
        stopifnot(is.numeric(p) && length(p) == 1)
        stopifnot(is.numeric(adjs) && length(adjs) >= 1)

        n <- length(adjs)
        stopifnot(x.surf >= 0 && x.surf <= n)
        adj <- if (x.surf == 0) {
                   adjs[1]
               } else if (x.surf == n) {
                   adjs[n]
               } else {
                   adjs[x.surf]*(1-p) + adjs[x.surf+1]*p
               }
    }


    ##-------------------------------------------------------------------------
    ## Scales measurement using values of surface.
    scaleBySurface <- function(x,
                               surface,
                               n=1) {
        stopifnot(is.numeric(x))
        stopifnot(is.character(surface))
        stopifnot(is.numeric(n))

        s1 <- eval.parent(as.name(surface), n=n)
        stopifnot(is.numeric(s1))
        (x / s1) * median(s1)
    }


    adjustment <- if (ndilutions > 1) {
      ## Organize input data for "which.surface" function
      input.mat <- as.matrix(rppa@data[, measure])

      for (surface in surfaces) {
          s1 <- eval(as.name(surface))
          input.mat <- cbind(input.mat, s1)
      }
      dimnames(input.mat) <- list(NULL,
                                  c("measure", surfaces))
      ## Find closest positive control surface and fraction
      ## between two surfaces
      place <- apply(input.mat,
                     1,
                     function(x) {
                         which.surface(meas=x[1],
                                       surfs=x[-1])
                     })
      values.mat <- cbind(place, input.mat)

      rownames(values.mat) <- NULL
      p <- apply(values.mat,
                 1,
                 function(x) {
                     getp(x.surf=x[1],
                          meas=x[2],
                          surfs=x[-(1:2)])
                 })

      ## Perform scaling on each positive control surface
      x <- rppa@data[, measure]
      adjs <- sapply(surfaces,
                     scaleBySurface,
                     x=x,
                     n=3)
      dimnames(adjs) <- list(NULL,
                             paste("adj", dilutions, sep=""))

      ## Now retrieve appropriate adjustment based on the
      ## closest positive control surface and the fraction
      ## between two surfaces as computed by "getadj" function.
      adjust.mat <- cbind(place, p, adjs)
      rownames(adjust.mat) <- NULL

      adjusted <- apply(adjust.mat,
            1,
            function(x) {
                getadj(x.surf=x[1],
                       p=x[2],
                       adjs=x[-(1:2)])
            })
    } else {
        x <- rppa@data[, measure]
        as.numeric(scaleBySurface(x, surfaces[1]))
    }
	
    rppa@data[[adjMeasure]] <- adjustment

    return(rppa)
}


##-----------------------------------------------------------------------------
spatialAdjustment <- function(rppa,
                              cutoff=0.8,
                              k=100,
                              gamma=0.1,
                              plotSurface=FALSE) {
    ## Check arguments
    if (!is.RPPA(rppa)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("rppa"), "RPPA"))
    }

    ## Begin processing
    measures <- eval(formals(spatialCorrection)$measure)
    tf.measures <- measures %in% colnames(rppa@data)

    if (all(!tf.measures)) {
        message("cannot perform spatial adjustment")
        stop(sprintf("argument %s missing all measure columns: %s",
                     sQuote("rppa"), paste(measures, collapse=", ")))
    }

    for (i in seq_along(measures)) {
        measure <- measures[i]

        if (tf.measures[i]) {
            rppa <- spatialCorrection(rppa,
                                      measure=measure,
                                      cutoff=cutoff,
                                      k=k,
                                      gamma=gamma,
                                      plotSurface=plotSurface)
        } else {
            message(sprintf("cannot perform spatial adjustment for measure %s",
                            sQuote(measure)))
            warning(sprintf("argument %s missing measure column: %s",
                            sQuote("rppa"), measure))
        }
    }

	if(any(rppa@data$Adj.Net.Value < -65536 | rppa@data$Adj.Net.Value > 65536)) {
		warning(paste("Antibody ", rppa@antibody, ": Spatial Adjustments cause values outside range of [-65536, 65536].  Values will be capped at those levels.", sep=""))
		
		rppa@data$Adj.Net.Value[rppa@data$Adj.Net.Value < -65536] <- -65536
		rppa@data$Adj.Net.Value[rppa@data$Adj.Net.Value > 65536] <- 65536
	}
	
    return(rppa)
}


##-----------------------------------------------------------------------------
spatialAdjustmentFromParams <- function(rppa,
                                        spatialparams) {
    ## Check arguments
    if (!is.RPPASpatialParams(spatialparams)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("spatialparams"), "RPPASpatialParams"))
    }

    ## Begin processing
    spatialAdjustment(rppa,
                      cutoff=spatialparams@cutoff,
                      k=spatialparams@k,
                      gamma=spatialparams@gamma,
                      plotSurface=spatialparams@plotSurface)
}

