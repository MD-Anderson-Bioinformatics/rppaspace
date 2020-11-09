##-----------------------------------------------------------------------------
## Create a Tracking object to determine how points in dilution
## series should be handled within RPPASPACE.
createTracking <- function(slideData, antibody) {
		
    numRows <- length(rownames(slideData))
    NAS = rep(NA, numRows)        # Set of values all NA
    FALSES <- rep(FALSE, numRows) # Set of values all FALSE
    TRUES <- rep(TRUE, numRows)   # Set of values all TRUE

    ## tracking is a dataframe used to track the points data from a slide and how they are used
    ## Sample points are loaded with much of the information coming from the sample file for a single slide
    ## Additional information is either loaded from a design file or from object parameters passed in the R Script
    ## This design-related information is used to modify all slides that are part of a run of the software
    ## These sources of information determine which points from the same are used for the following:
    ## 1) Negative controls used to filter positive control points
    ## 2) Positive control points used to perform spatial adjustments on sample points
    ## 3) Sample points that are used to create a curve to represent the data for one slide/antibody
    ## 4) Sample points to fit to the curve created in step 3 above and provided as textual and graphical output
    ## 5) Positive Control points treated as noise values and fit to the curve and provided as textual output 
    ## 6) List of any points that caused issues during the above processing which have output values that should be suspect (badPoint)

    tracking <- data.frame(
		
        spotType = slideData$Spot.Type,

        # Is point a negative control point
        # Default FALSE
        # Set to TRUE if value of Spot.Type in design file is negative control type
        isNegCtrl = FALSES,

        # Is point a positive control point
        # Default FALSE
        # Set to TRUE if value of Spot.Type in design file is positive control type
        isPosCtrl = FALSES,

        # Is point a control point
        # Default FALSE
        # Set to TRUE if value of Spot.Type in design file is any control type
        isCtrl = FALSES,

        # Apply spatial correction to point
        # Default TRUE
        # Set to FALSE if value of Spot.Type in design file is any control type
        # Set to TRUE if value of Spot.Type in design file is any noise type
        applySpatialCorrection = TRUES,

        # Should point be used to create curve to which to fit data?
        # Default TRUE
        # Set to False if Spot.Type is control type or noise type in design file
        # Set to FALSE if in SeriesToIgnore
        makePartOfCurve = TRUES,

        # Should point be fit to curve?
        # Default TRUE
        # Set to False if point is control in design but not a noise point
        fitToCurve = TRUES,

        # Should point be used in noise calculations
        # Default FALSE
        # Set to TRUE if Noise Point in Design file
        isNoise = FALSES,

        # Is point a sample point?
        # Default TRUE
        # Set to False if Sample in design file
        isSample = TRUES,

        # Does this point have a value that was not used or caused problems in processing and whose output accuracy should be questioned?
        # Default FALSE
        # Set to TRUE if not a valid numeric value
        # Set to TRUE if becomes positive control point with negative value
        # Set to TRUE if causes error when creating curve
        # Set to TRUE if causes error when fitting to curve 
        # Set to TRUE if causes error when calculating noise
        badPoint=FALSES,
		
		# Dilution value for this spot on the slide
		dilution = slideData$Dilution,
		
        row.names=as.character(rownames(slideData)),
        stringsAsFactors = FALSE
    )

    #If point is a negative control point
    negnames <- which(slideData$Spot.Type %in% spottype.negativecontrol)
    if (length(negnames) > 0) {
        tracking[negnames,]$isNegCtrl <- TRUE
        tracking[negnames,]$isCtrl <- TRUE
        tracking[negnames,]$applySpatialCorrection <- FALSE
        tracking[negnames,]$isSample <- FALSE
        tracking[negnames,]$makePartOfCurve <- FALSE
        tracking[negnames,]$fitToCurve <- FALSE
    }
    #If point is a positive control point
    posnames <- which(slideData$Spot.Type %in% spottype.positivecontrol)
    if (length(posnames) > 0) {
        tracking[posnames,]$isPosCtrl <- TRUE
        tracking[posnames,]$isCtrl <- TRUE
        tracking[posnames,]$applySpatialCorrection <- FALSE
        tracking[posnames,]$isSample <- FALSE
        tracking[posnames,]$makePartOfCurve <- FALSE
        tracking[posnames,]$fitToCurve <- FALSE
    } 
	
	#If either positive or negative controls are not provided, don't do spatial corrections
	if (length(negnames) == 0 || length(posnames) == 0)  {
		tracking$applySpatialCorrection = FALSES
	}

    #If point is a noise point
    noisenames <- which(slideData$Spot.Type %in% spottype.noise)
    if (length(noisenames) > 0) {
        tracking[noisenames,]$applySpatialCorrection <- TRUE
        tracking[noisenames,]$makePartOfCurve <- FALSE
        tracking[noisenames,]$fitToCurve <- TRUE
        tracking[noisenames,]$isNoise <- TRUE
        tracking[noisenames,]$isSample <- FALSE
    }

	tracking
}