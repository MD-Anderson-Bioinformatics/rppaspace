###
### $Id: rs6-rppaSetSummary.R
### Summarize fit processing results
###


##=============================================================================
setClassUnion("OptionalNumeric", c("numeric", "NULL"))
setClassUnion("OptionalMatrix", c("matrix", "NULL"))

setClass("RPPASetSummary",
         representation(raw="matrix",              ## raw concentrations
                        ss="matrix",               ## sum squares ratio
                        norm="matrix",             ## normalized concentrations
                        probs="OptionalNumeric",   ## probability good slide
                        completed="matrix",        ## what worked/failed
						noise="OptionalMatrix",    ## noise values, if calculated
                        design="RPPADesignParams", ## design parameters for all slides
                        onlynormqcgood="logical",  ## filter norm'd by qc score
                        version="character"),      ## package version
         prototype(probs=NULL))


##-----------------------------------------------------------------------------
is.RPPASetSummary <- function(x) {
    is(x, "RPPASetSummary")
}


##-----------------------------------------------------------------------------
## Returns a slot in the array of fits as a simple matrix view.
.fitSlot <- function(rppaset,
                     slotname) {
    ## Check arguments
    stopifnot(is.RPPASet(rppaset))
    stopifnot(is.character(slotname) && length(slotname) == 1)

    rppafits.tf <- rppaset@completed[, 'fit']
    rppafits <- rppaset@fits[rppafits.tf]

    if (!(slotname %in% slotNames(rppafits[[1]]))) {
        stop(sprintf("invalid slotname %s",
                     sQuote(slotname)))
    }

    ## Begin processing
    sapply(rppafits,
           slot,
           name=slotname)
}


##-----------------------------------------------------------------------------
## Create an RPPASetSummary object
RPPASetSummary <- function(rppaset,
                           onlynormqcgood=ran.prefitqc(rppaset)) {
    ## Check arguments
    if (!is.RPPASet(rppaset)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("rppaset"), "RPPASet"))
    }
    if (!is.logical(onlynormqcgood)) {
        stop(sprintf("argument %s must be logical",
                     sQuote("onlynormqcgood")))
    } else if (!(length(onlynormqcgood) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("onlynormqcgood")))
    }

    ## Begin processing
    conc.raw <- .fitSlot(rppaset, "concentrations")
    conc.ss  <- .fitSlot(rppaset, "ss.ratio")
	#noise <- .fitSlot(rppaset, "noise")

    ## Generate probabilities (goodness) for each processed slide (if any)
    prefitqcs.tf <- rppaset@completed[, "prefitqc"]
    probs <- if (!all(is.na(prefitqcs.tf))) {
                 prefitqcs <- rppaset@prefitqcs[prefitqcs.tf]
                 sapply(prefitqcs, qcprob)
             } else {
                 NULL
             }

    ## Normalize the concentrations
    norm.tf <- if (onlynormqcgood) {
                   ## Remove "bad" slides...
                   local({
                       probs.tmp <- probs
                       probs.tmp[is.na(probs.tmp)] <- 0
                       good.cutoff <- 0.8
                       probs.tmp >= good.cutoff
                   })
               } else {
                   rep(TRUE, ncol(conc.raw))
               }
    normparams <- rppaset@normparams
    normalizeArgs <- c(list(), normparams@arglist)
    normalizeArgs$object <- conc.raw[, norm.tf, drop=FALSE]
    normalizeArgs$method <- normparams@method

    conc.norm <- do.call(normalize, normalizeArgs)

    ## To allow reorder on write
    rppafits.tf <- rppaset@completed[, "fit"]
    rppa <- rppaset@rppas[rppafits.tf][[1]]

	if (is.null(rppa)) {
		msg <- "Skipping noise value calculations as no valid slides were found to use for calculating noise values."
		print(msg)
		warning(msg)
		noise <- NULL
	} else {
		# Calculate number of noise dilutions in slide
		noiseDilutions <- rppa@tracking[rppa@tracking$isNoise,]$dilution
		numNoiseDilutions <- length(unique(noiseDilutions))

		if (numNoiseDilutions > 0) {
		
			##Create the matrix for noise information output
			numNoiseSeries <- sum(rppa@tracking$isNoise, na.rm=TRUE) / numNoiseDilutions
		
			mat <- matrix(as.numeric(NA), nrow=numNoiseSeries + 2, ncol=length(rppaset@fits))

			rownames(mat) <- c(sort(unique(as.character(rppa@data$Series.Id[rppa@tracking$isNoise==TRUE]))), "mean", "sd")
			colnames(mat) <- names(rppaset@fits)

			for (i in 1:length(rppaset@fits)) {
				currentFitName <- names(rppaset@fits)[i]
				currentFitNoise <- rppaset@fits[[i]]@noise
				
				for (j in 1:length(currentFitNoise)) {
					currentRowName <- names(currentFitNoise)[[j]]
					mat[currentRowName, currentFitName] <- currentFitNoise[[j]]
				}
				tempPoints <- mat[0:numNoiseSeries,currentFitName]
				tempPoints <- tempPoints[!is.na(tempPoints)]
				
				##Set value for mean of noise points for slide
				if (length(tempPoints) == 0) {
					mat["mean", currentFitName] <- NA
				}
				else {
					mat["mean", currentFitName] <- mean(mat[0:numNoiseSeries,currentFitName], na.rm=TRUE)
				}
				
				##Set value for standard deviation of noise points for slide
				if (length(tempPoints) <= 1) {
					mat["sd", currentFitName] <- NA
				}
				else {
					mat["sd", currentFitName] <- sd(mat[0:numNoiseSeries,currentFitName], na.rm=TRUE)
				}
			}
			noise <- mat
		}
		else {
			noise <- NULL
		}
	}

	firstsample.tf <- rppa@tracking$isSample & !duplicated(rppa@data$Series.Id)
	rownums <- rppa@data$Series.Id
	locations <- as.integer(rownums[firstsample.tf])[]

	attr(conc.raw,  "locations") <- locations
	attr(conc.ss,   "locations") <- locations
	attr(conc.norm, "locations") <- locations
	#attr(noise, "locations") <- locations
	
    ## Create new class
    new("RPPASetSummary",
        raw=conc.raw,
        ss=conc.ss,
        norm=conc.norm,
        probs=probs,
        completed=rppaset@completed,
		noise=noise,
        design=rppaset@design,
        onlynormqcgood=onlynormqcgood,
        version=packageDescription("RPPASPACE", fields="Version"))
}


##-----------------------------------------------------------------------------
## Provide a convenience function to save fit summary results as CSV/TSV files
setMethod("write.summary", signature(object="RPPASetSummary"),
          function(object,
                   path,
                   prefix="rppaspace",
                   ...) {
    ## Check arguments
    if (!is.character(path)) {
        stop(sprintf("argument %s must be character",
                     sQuote("path")))
    } else if (!(length(path) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("path")))
    } else if (!dir.exists(path)) {
        stop(sprintf("directory %s does not exist",
                     dQuote(path)))
    } else if (!dir.writable(path)) {
        stop(sprintf("directory %s is not writable",
                     dQuote(path)))
    }

    if (!is.character(prefix)) {
        stop(sprintf("argument %s must be character",
                     sQuote("prefix")))
    } else if (!(length(prefix) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("prefix")))
    }


    ##-------------------------------------------------------------------------
    ## Allows concentrations to be written back in different order.
    get_concs_ordered_for_write <- function(object,
                                            slotname) {
        stopifnot(is.RPPASetSummary(object))
        stopifnot(is.character(slotname) && length(slotname) == 1)

        if (!(slotname %in% slotNames(object))) {
            stop(sprintf("invalid slotname %s",
                         sQuote(slotname)))
        }

        concs <- slot(object, slotname)
		tempLocs <- attr(concs, "locations", exact=TRUE)

        if (!is.null(tempLocs)) {
			locations <- as.integer(tempLocs)
            concs[order(locations), ]
        } else {
            concs
        }
    }


    ##-------------------------------------------------------------------------
    ## Create informative information to the filename
    mknormtag <- function(object) {
        stopifnot(is.RPPASetSummary(object))

        attrs <- attr(object@norm, "normalization", exact=TRUE)

        normMethod <- attrs$method
        qctag <- if (object@onlynormqcgood) {
                     "qc"
                 }
        if (normMethod == "none") {
            normMethod
        } else {
            if (normMethod == "medpolish") {
                if (is.null(qctag)) {
                    normMethod
                } else {
                    paste(normMethod, qctag, sep="-")
                }
            } else {
                sweeptag <- if (attrs$sweep.cols) "rowcol" else "col"
                if (is.null(qctag)) {
                    paste(normMethod, sweeptag, sep="-")
                } else {
                    paste(normMethod, sweeptag, qctag, sep="-")
                }
            }
        }
    }


    ## Begin processing
    ## Write file for raw concentrations
	print("Writing raw concentrations file")
    filename <- sprintf("%s_conc_raw.csv", prefix)
    conc.raw <- signif(get_concs_ordered_for_write(object, "raw"),7)
    write.csv(conc.raw, file=file.path(path, .portableFilename(filename)))

    ## Write file for R^2 statistics
	print("Writing R^2 statistics file")
    filename <- sprintf("%s_ss_ratio.csv", prefix)
    conc.ss <- signif(get_concs_ordered_for_write(object, "ss"),7)
    write.csv(conc.ss, file=file.path(path, .portableFilename(filename)))

    ## Write file for normalized concentrations
	print("Writing normalized concentrations file")
	filename <- sprintf("%s_conc_norm_%s.csv", prefix, mknormtag(object))
	conc.norm <- signif(get_concs_ordered_for_write(object, "norm"),7)
	write.csv(conc.norm, file=file.path(path, .portableFilename(filename)))

	#Number of columns in combined qc output
	#If any new qc statistics are created, follow a similar pattern to code below 
	#for prefit qc and noise qc values to make sure they are written as their own file 
	#and that the relevant statistics are included in the combined qc output file.
	qc_column_count <- 0

    ## If positive control points set up to be used to calculate noise 
	noise_qc_calculated <- length(dim(object@noise)) == 2
	
    if (noise_qc_calculated) {
		qc_column_count <- qc_column_count + 4
		filename <- sprintf("%s_noise.csv", prefix)
		noise.out <- signif(slot(object, "noise"),7)
		print("Writing noise qc output file.")
		write.csv(noise.out, file=file.path(path, .portableFilename(filename)))
	}

    ## If QC processing was requested...
	prefit_qc_calculated <- !is.null(object@probs)
    if (prefit_qc_calculated) {
		qc_column_count <- qc_column_count + 1
        ## Write file for QC probabilities
        filename <- sprintf("%s_prefit_qc.csv", prefix)
        probs.df <- data.frame("Filename"=names(object@probs),
                               "Probabilities"=signif(object@probs,7),
                               row.names=seq_along(object@probs),
							   stringsAsFactors = FALSE
							   )
		print("Writing prefit qc output file.")
        write.csv(probs.df, file=file.path(path, .portableFilename(filename)))
    }
	

	## If QC processing was requested or positive control points set up
	## to be used to calculate noise, then write combined output file.
	if (qc_column_count > 0) {
		filename <- sprintf("%s_combined_qc.csv", prefix)
		antibodies <- colnames(object@raw)
		combined_qc_data <- data.frame(matrix(as.numeric(NA), ncol=qc_column_count, nrow=length(antibodies)), row.names = antibodies, stringsAsFactors = FALSE)
		
		current_column <- 1
		if (prefit_qc_calculated) {
			qc_col_names <- "prefit_qc"
			combined_qc_data[probs.df$Filename, current_column] <- probs.df$Probabilities
			current_column <- current_column + 1
		} else {
			print("No Pre-fit QC values calculated. Pre-Fit QC Will not be part of combined QC output.")
		}

		#If noise was calculated
		if (noise_qc_calculated) {
			noise_to_combine <- t(object@noise[rownames(object@noise) %in% c("mean","sd"),])
			
			if (prefit_qc_calculated) {
				qc_col_names <- c(qc_col_names, "noise_sd", "noise_mean", "noise_cv", "noise_n")

			} else {
				qc_col_names <- c("noise_sd", "noise_mean", "noise_cv", "noise_n")
			}
			
			#Set noise_sd column for calculated sds
			combined_qc_data[rownames(noise_to_combine), current_column] <- signif(noise_to_combine[rownames(noise_to_combine),"sd"],7)
			
			#Set noise_mean column for calculated means
			combined_qc_data[rownames(noise_to_combine), current_column + 1] <- signif(noise_to_combine[rownames(noise_to_combine),"mean"], 7)
				
			#Set noise_cv column for calculated sd/mean
			combined_qc_data[rownames(noise_to_combine), current_column + 2] <- 
				signif(noise_to_combine[rownames(noise_to_combine),"sd"] / noise_to_combine[rownames(noise_to_combine),"mean"], 7)
				
			#Set noise_n column for number of dilution series to use as noise 
			combined_qc_data[rownames(noise_to_combine), current_column + 3] <- length(rownames(object@noise)) - 2

			current_column <- current_column + 4
		} else {
			print("No Noise QC values calculated. Noise QC Will not be part of combined QC output.")
		}
		colnames(combined_qc_data) <- qc_col_names
		
		print("Writing combined qc output file.")
		write.csv(combined_qc_data, file=file.path(path, .portableFilename(filename)))
	}

    ## Write file for stage completion summary
	print("Writing summary file.")
    filename <- sprintf("%s_summary.tsv", prefix)
    write.table(object@completed,
                file=file.path(path, .portableFilename(filename)),
                sep='\t',
                col.names=NA)

    invisible(NULL)
})

