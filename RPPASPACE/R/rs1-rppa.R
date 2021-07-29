###
### $Id: rs1-rppa.R
###


##=============================================================================
setClassUnion("OptionalList", c("list", "NULL"))
setClass("RPPA",
         representation(data="data.frame",
                        file="character",
						slideNumber="integer",
                        antibody="character",
                        tracking="data.frame",
                        seriesToIgnore="OptionalList",
						warningsFileName="character"
                        ))


##-----------------------------------------------------------------------------
is.RPPA <- function(x) {
    is(x, "RPPA")
}


##-----------------------------------------------------------------------------
## Generates an RPPA object from a quantification file.
RPPA <- function(file,
                 path=".",
				 slideNumber=NA,
                 antibody=NULL,
                 tracking=NULL,
                 seriesToIgnore=NULL,
				 warningsFileName="warnings.txt"
                 ) {
				 
	warnings <- 0
    ## Check arguments
    if (is.character(file)) {
        if (!(length(file) == 1)) {
            stop(sprintf("argument %s must be of length 1",
                         sQuote("file")))
        } else if (!nzchar(file)) {
            stop(sprintf("argument %s must not be empty string",
                         sQuote("file")))
        }

        path_or_url <- if (.isAbsolutePathname(file) || .hasScheme(file)) {
                           file
                       } else {
                           if (!is.character(path)) {
                               stop(sprintf("argument %s must be character",
                                            sQuote("path")))
                           } else if (!(length(path) == 1)) {
                               stop(sprintf("argument %s must be of length 1",
                                            sQuote("path")))
                           }

                           if (.hasScheme(path)) {
                               paste(path, file, sep="/")
                           } else {
                               file.path(path, file)
                           }
                       }

        is_url <- .hasScheme(path_or_url)
        if (!is_url) {
            if (!file.exists(path_or_url)) {
                stop(sprintf("file %s does not exist",
                             dQuote(path_or_url)))
            }
        }

        ## Convert to connection object
        file <- if (is_url) {
                    url(path_or_url, "r")
                } else {
                    file(path_or_url, "r")
                }
        on.exit(close(file))
    }
    filename <- basename(summary(file)$description)

    if (!is.null(antibody)) {
        if (!is.character(antibody)) {
            stop(sprintf("argument %s must be character",
                         sQuote("antibody")))
        } else if (!(length(antibody) == 1)) {
            stop(sprintf("argument %s must be of length 1",
                         sQuote("antibody")))
        } else if (!nzchar(antibody)) {
            stop(sprintf("argument %s must not be empty string",
                         sQuote("antibody")))
        }
    } else {
        ## Use filename without extension as default value
        txt.re <- "\\.[tT][xX][tT]$"
        basename <- sub(txt.re, "", filename)
        antibody <- sub("[[:space:]]+$", "", basename)
    }
	
    ## Read quantification file
    quant.df <- readQuantification(file)

	#Add steps based on dilution values
	steps <- rep(NA, length(quant.df$Dilution))
	steps[quant.df$Dilution > 0] <- log(quant.df$Dilution[quant.df$Dilution > 0]/100)/log(2)
	steps[quant.df$Dilution == 0] <- 0
	quant.df$Steps <- steps
	##TODO:Centering needs to be handled appropriately - see original design class
	
    rownames(quant.df) <- paste( quant.df$Main.Row, 
            quant.df$Main.Col,  
            quant.df$Sub.Row, 
            quant.df$Sub.Col, sep = "_" )

	badSlide <- FALSE
	if (is.null(tracking)) {
		# Check that first slide object makes sense 
		series <- sort(quant.df[quant.df$Series.Id > 0 & quant.df$Dilution > 0,]$Series.Id)
		uniquesSeries <- unique(series)
		uniqueCharSeries <- as.character(uniquesSeries)
		dilutions <- quant.df[quant.df$Series.Id > 0 & quant.df$Dilution > 0,]$Dilution
		uniquePositiveDilutions <- unique(dilutions)
		
		# Test that all positive dilutions have entries for all series (and vice versa).
		for (dilution in uniquePositiveDilutions) {
			dilutionSeries <- sort(quant.df[quant.df$Series.Id > 0 & quant.df$Dilution == dilution,]$Series.Id)

			matchSeries <- all.equal(as.character(dilutionSeries), uniqueCharSeries)
			if ( matchSeries != "TRUE") {
				badSlide <- TRUE
				msg <- paste("Slide ", slideNumber, " (",  antibody, ") For dilution value ", dilution, " not all series are present. ", paste(matchSeries, collapse="; "), sep="")
				warning(msg)
				write(msg, warningsFileName, append=TRUE)
				warnings <- warnings + 1
			}
		}
		
		# Slide passes initial tests, so use it to create tracking object
		if (badSlide == FALSE) {
			print(paste("Creating tracking object using slide ", antibody), sep="")
			tracking <- createTracking(quant.df, antibody)
		}
	} else {
		# Verify slides are the same length
		if (nrow(tracking) != nrow(quant.df)) {
			msg <- paste("Slide ", slideNumber, " (",  antibody,  ") The number of rows (", nrow(quant.df), 
				") do not match the number of rows (", nrow(tracking),
				") in the first valid slide.", sep="")
			write(msg, warningsFileName, append=TRUE)
			warning(msg)
			warnings <- warnings + 1
		}

		# If tracking already exists then this is not the first slide
		# Verify spot types in tracking match ones in this slide
		spotTypesMatch <- as.character(all.equal(tracking$spotType, quant.df$Spot.Type))
				
		if (spotTypesMatch != "TRUE") {
			splitMatch <- unlist(strsplit(spotTypesMatch, " "))
			if (splitMatch[[3]] == "mismatch" || splitMatch[[3]] == "mismatches") {
				msg <- paste("Slide ", slideNumber, " (",  antibody, ") Spot types in slide do not match those in first valid slide.", paste(spotTypesMatch, collapse="; "), sep="")
				write(msg, warningsFileName, append=TRUE)
				warning(msg)
				warnings <- warnings + 1
			}
		}
		
		# Verify dilutions in tracking match ones in this slide
		dilutionsMatch <- all.equal(as.character(tracking$dilution), as.character(quant.df$Dilution))
				
		if (dilutionsMatch != "TRUE") {
			splitMatch <- unlist(strsplit(dilutionsMatch, " "))
			if (splitMatch[[3]] == "mismatch" || splitMatch[[3]] == "mismatches") {
				msg <- paste( "Slide ", slideNumber, " (",  antibody, ") Dilutions do not match those in first valid slide. ", paste(dilutionsMatch, collapse="; "), sep="")
				write(msg, warningsFileName, append=TRUE)
				warning(msg)
				warnings <- warnings + 1
			}
		}
	}

	if (badSlide == FALSE) {
		if (is.null(tracking) || (!is.data.frame(tracking))) {
		    stop(sprintf("Bad slide: unable to create tracking object."))
		}

		#Modify tracking for series to be ignored
		if (!is.null(seriesToIgnore) & length(seriesToIgnore) > 0) {
			ignoreItems <- subset(quant.df, quant.df$Series.Id %in% seriesToIgnore)		
			if (length(unique(seriesToIgnore)) > length(unique(ignoreItems$Series.Id))) {
				stop("The seriesToIgnore parameter contains series that do not exist in the Series.Id of the first valid slide.")
			}
			ignoreRows <- rownames(ignoreItems)
			tracking[ignoreRows,]$makePartOfCurve <- FALSE
		}
	}

	obj <- if ( warnings > 0) {
		NULL
	} else {
		## Create new class
		new("RPPA",
			data=quant.df,
			file=filename,
			antibody=antibody,
			tracking=tracking,
			seriesToIgnore=seriesToIgnore
			)
	}
}


##-----------------------------------------------------------------------------
setMethod("dim", signature(x="RPPA"),
          function(x) {
    .dimOfLayout(x@data)
})


##-----------------------------------------------------------------------------
setMethod("summary", signature(object="RPPA"),
          function(object,
                   ...) {
    cat(sprintf("An %s object loaded from file %s",
                class(object), dQuote(object@file)), "\n")
    cat("antibody:", object@antibody, "\n")
    cat("\n")
    print(dim(object))
    cat("\n")
    unneededColnames <- c(.locationColnames(), "Series.Id")
    summarizable <- !colnames(object@data) %in% unneededColnames
    print(summary(object@data[summarizable]))
})


##-----------------------------------------------------------------------------
setMethod("image", signature(x="RPPA"),
          function(x,
                   rot=as.integer(0),
                   measure="Net.Value",
                   main=.mkPlotTitle(measure, x@antibody),
                   colorbar=FALSE,
                   col=terrain.colors(256),
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

    # the whole thing should just use "try(as.logical())" and, if it
    # works, live with the result....
    if (is.numeric(colorbar)) {
        colorbar <- as.logical(colorbar)
    }

    if (!is.logical(colorbar)) {
        stop(sprintf("argument %s must be logical",
                     sQuote("colorbar")))
    } else if (!(length(colorbar) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("colorbar")))
    }
    ## Begin processing
    data.df <- x@data

    dim.rppa <- dim(x)
    my <- dim.rppa["Main.Row"] * dim.rppa["Sub.Row"]
    mx <- dim.rppa["Main.Col"] * dim.rppa["Sub.Col"]
    yspot <- 1 + my - (max(data.df$Sub.Row)*(data.df$Main.Row-1) + data.df$Sub.Row)
    xspot <- max(data.df$Sub.Col)*(data.df$Main.Col-1) + data.df$Sub.Col

    geo <- tapply(data.df[, measure],
                  list(xspot, yspot),
                  mean)

	if (rot == as.integer(90)) 
	{
		geo <- rotateMatrixClockwise90Degrees(geo, 1)
		temp <- my
		my <- mx
		mx <- temp
		temp <- yspot
		yspot <- xspot
		xspot <- temp
	}
	if (rot == as.integer(180)) 
	{
		geo <- rotateMatrixClockwise90Degrees(geo, 2)
	}
	if (rot == as.integer(270)) 
	{
		geo <- rotateMatrixClockwise90Degrees(geo, 3)
		temp <- my
		my <- mx
		mx <- temp
		temp <- yspot
		yspot <- xspot
		xspot <- temp
	}

    if (colorbar) {
        ## Get the size of the plotting region in relative units
        startPlt <- par()$plt

        ## We're only going to partition things on the x-axis, so only
        ## the first 2 coordinates are of interest. Define the boundaries
        ## for the two panels so that a 10/1 width ratio is attained.
        imagePlt    <- startPlt
        colorbarPlt <- startPlt
        startWidth  <- startPlt[2] - startPlt[1]

        imagePlt[2]    <- startPlt[1] + (10/12)*startWidth
        colorbarPlt[1] <- startPlt[2] - ( 1/12)*startWidth

        ## Draw the colorbar
        ## :TODO: Figure out how to set margins so it works in small windows...
        par(plt=colorbarPlt)
		
		minGeo <- min(geo, na.rm=TRUE)
		maxGeo <- min(geo, na.rm=TRUE)
		meanGeo <- mean(geo, na.rm=TRUE)
		if (minGeo >= meanGeo) { minGeo <- minGeo - 0.0000001 }
		if (maxGeo <= meanGeo) { maxGeo <- maxGeo + 0.0000001 }
		
		#Plot the colorbar (legend part of the plot) on the right side of the image box
        image(x=1,
              y=seq(minGeo, maxGeo, length=256),
              z=matrix(seq_len(256), nrow=1),
              col=col,
              xaxt="n",
              xlab="",
              yaxt="n",
              ylab="")
			  
        axis(4) # Put labeling at right
        box()

        ## Set things up to draw main image and revert back for next figure
        par(plt=imagePlt, new=TRUE)
        on.exit(par(plt=startPlt))
    }

	#Plot the actual image box
    image(x=seq_len(mx),
          y=seq_len(my),
          z=geo,
          col=col,
          main=main,
          sub=paste("File:", x@file),
          xaxt="n", #Supress plotting of x-axis
          xlab="",
          yaxt="n", #Supress plotting of y-axis
          ylab="",
          ...)

	at.x <- seq(from=dim.rppa["Sub.Col"],
                to=dim.rppa["Sub.Col"]*dim.rppa["Main.Col"],
                by=dim.rppa["Sub.Col"])
    at.y <- seq(from=dim.rppa["Sub.Row"],
                to=dim.rppa["Sub.Row"]*dim.rppa["Main.Row"],
                by=dim.rppa["Sub.Row"])
    axis(1, at=at.x)
    axis(2, at=at.y)

    abline(h=(0.5 + seq(0, my, length=1+dim.rppa["Main.Row"])))
    abline(v=(0.5 + seq(0, mx, length=1+dim.rppa["Main.Col"])))

	bv <- c(0:11)
	vlines <- c(5+11*bv,10+11*bv,11+11*bv)
	bh <- c(0:4)
	hlines <- c(5+bh*11,6+bh*11,11+bh*11)
	
	abline(v=(0.5 + vlines), col="green")
	abline(h=(0.5 + hlines), col="green")

    invisible(x)
})

##-----------------------------------------------------------------------------
seriesNames <- function(rppa) {
    ## Check arguments
    if (!is.RPPA(rppa)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("rppa"), "RPPA"))
    }

    series <- rppa@data$Series.Id[rppa@tracking$isSample]
    unique(series)
}

##-----------------------------------------------------------------------------
seriesToUseToMakeCurve <- function(rppa) {
    ## Check arguments
    if (!is.RPPA(rppa)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("rppa"), "RPPA"))
    }

    series <- unique(rppa@data$Series.Id[rppa@tracking$isSample])
	seriesToUse <- unique(rppa@data$Series.Id[rppa@tracking$makePartOfCurve])
    intersect(series, seriesToUse)
}