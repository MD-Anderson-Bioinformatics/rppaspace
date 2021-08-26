###
### $Id: rs5-rppaSet.R
### Fit a set of slides with a common layout
###


##=============================================================================
setClassUnion("OptionalRPPASpatialParams", c("RPPASpatialParams", "NULL"))

setClass("RPPASet",
         representation(call="call",                 ## function invocation
						design="RPPADesignParams",
                        rppas="array",               ## vector of RPPAs
                        spatialparams="OptionalRPPASpatialParams",
                        prefitqcs="array",           ## vector of QC values
                        fitparams="RPPAFitParams",
                        fits="array",                ## set of fits
                        completed="matrix",          ## what worked/failed
                        normparams="RPPANormalizationParams",
                        version="character",
						residualsrotation="integer",
						warningsFileName="character",
						errorsFileName="character"						
						))        ## package version


##-----------------------------------------------------------------------------
is.RPPASet <- function(x) {
    is(x, "RPPASet")
}

##-----------------------------------------------------------------------------
## Method to rotate a matrix 90 degrees, by default clockwise.
rotateMatrix  <- function(x, clockwise = T) 
{
	if (clockwise) 
	{ 
		t( apply(x, 2, rev))
	} 
	else 
	{
		apply( t(x),2, rev)
	} 
}

##-----------------------------------------------------------------------------
## Method to rotate a matrix numRots increments x 90 degrees clockwise.
rotateMatrixClockwise90Degrees <- function(m, numRots = c(0,1,2,3)) 
{
	#numRots
	# 0: 0 degree rotation, return original matrix 
	# 1: 90 degree rotation, rotate matrix clockwise once
	# 2: 180 degree rotation, rotate matrix clockwise twice
	# 3: 270 degree rotation, rotate matrix counterclockwise once
	return (switch(numRots+1, m, rotateMatrix(m), rotateMatrix(rotateMatrix(m)), rotateMatrix(m,FALSE)))
}

##-----------------------------------------------------------------------------
## Create the fit graphs and save them as PNG files

createResidualsPNG <- function (
	topPlotFileName,
	outputfilename,
	measure,
	antibody,
	m,
	majorxpoints = as.integer(NA),
	majorypoints = as.integer(NA),
	plotcolors = NA,
	missingvaluecolor = "#AAAAAA",
	gridcolor = "black",
	colorlegend = FALSE,
	cellsize = 3,
	rotation
	) 
{
	#---------------------------------------------
	#Setup for matrix colors
	#---------------------------------------------
	#Set up where breaks between colors should be.
	#10 color breaks with 11 color values over range from .4 to 1.
	
	#Set default Red-Yellow-Green color range if no color range was specified
	if (anyNA(plotcolors)) 
	{
		#Colors to use for rectangles in plot 
		#Go from Red (.4) to Yellow (.7) to Green (1.0)
		RYG <- c("#A50026", #Red
				 "#D73027",
				 "#F46D43",
				 "#FDAE61",
				 "#FEE066", 
				 "#FFFF88", #Yellow
				 "#D9EF88",
				 "#A6D96A",
				 "#66BD63",
				 "#1A9850",
				 "#006837"  #Green
				 ) 
				 
		plotcolors <- RYG 	
	}
	
	numcolors <- length(plotcolors)
	mincolorvalue <- .39999999999
	rangeofvalues <- 1 - mincolorvalue  
	breakpoints <- mincolorvalue + (0:(numcolors-1))/(numcolors/rangeofvalues)
	
	if (colorlegend) 
	{
		nacolorbarheight <- 12
		nacolordivider <- 10
		colorbarxsize <- 40
		colorbarysize <- as.integer(256/numcolors) * numcolors 
		if (colorbarysize / numcolors > 12)
		{
			colorbarysize <- 12 * numcolors 
		}
		colorbartotalheight <- colorbarysize + nacolordivider + nacolorbarheight
	}

	#Apply requested rotations to data values and switch grid labels as needed
	if (rotation == as.integer(90)) 
	{
		m <- rotateMatrixClockwise90Degrees(m, 1)
		
		temp <- majorxpoints
		majorxpoints <- majorypoints
		majorypoints <- temp
	} 
	else if (rotation == as.integer(180)) 
	{
		m <- rotateMatrixClockwise90Degrees(m, 2)
	} 
	else if (rotation == as.integer(270))
	{
		m <- rotateMatrixClockwise90Degrees(m, 3)
		
		temp <- majorxpoints
		majorxpoints <- majorypoints
		majorypoints <- temp
	}

	values <- as.vector(m)
	xpoints <- ncol(m)
	ypoints <- nrow(m)

	# Make matrices to plot series of rectangles based on number of x and y elements
	# times the sizes of the cells. This calculates the left (mxl) and bottom (myb)
	# coordinates of each rectangle to plot.
	mxl <- matrix(rep(0:(xpoints-1) * cellsize, ypoints), ncol=xpoints, nrow=ypoints, byrow=TRUE)
	myb <- matrix(rep(0:(ypoints-1) * cellsize, xpoints), ncol=xpoints, nrow=ypoints, byrow=FALSE)
	
	#Set index values for all colors in data matrix for slide RR2 values
	colindex <- as.integer(cut(as.array(m), breaks=c( unlist(breakpoints), 1)))
	#Color all blocks according to color scale adjusted from .4 to 1
	gcolors <- plotcolors[colindex]

	#Color all NA (missing) values gray.
	gcolors[is.na(gcolors)] <- missingvaluecolor

	#Base x grid labels off x major grid sections
	#gridxlabels <- (0:majorxpoints) * as.integer(xpoints/majorxpoints)
	gridxlabels <- (0:as.integer(xpoints/majorxpoints)) * majorxpoints
	
	#Base y grid labels off y major grid sections
	#If y grid labels are defined, attempt to use those, otherwise base them off y major grid sections
	#gridylabels <- ((majorypoints:0) * as.integer(ypoints/majorypoints))
	gridylabels <- ((as.integer(ypoints/majorypoints):0) * majorypoints)

	#Calculate how much to adjust grid locations to move each one a 2 pixels from axis per grid line > 0
	gridxadjust <- (0:(length(gridxlabels)-1)) * 2
	#if (min(gridxlabels) > 0) gridxadjust <- gridxadjust + 1
	gridyadjust <- ((length(gridylabels)-1):0) * 2
	#if (min(gridylabels) > 0) gridyadjust <- gridyadjust + 1

	#Set the pixel locations on grid for axis tick marks
	gridxlocs <- as.integer(gridxlabels * cellsize + gridxadjust - cellsize/2) + 1
	gridylocs <- as.integer(gridylabels * cellsize + gridyadjust - cellsize/2) + 1
	
	#Adjust left and top lines to be at 0 instead 1
	gridxlocs[1] <- -1
	gridylocs[length(gridylocs)] <- gridylocs[length(gridylocs)] 
	
	#If the last x grid location is not included in the generated locations, then add it.
	if ((xpoints %% majorxpoints != 0) && (length(gridxlocs) < xpoints + 1))
	{
		gridxlocs <- append(gridxlocs, as.integer(xpoints* cellsize + max(gridxadjust) + 1))
	}
	#If the last y grid location is not included in the generated locations, then add it.
	if ((ypoints %% majorypoints != 0) && (length(gridylocs) < ypoints + 1))
	{
		gridylocs <- append(gridylocs, as.integer(ypoints * cellsize + max(gridyadjust)))
	}

	#Adjust x coordinates of boxes by 2 pixels for each grid mark before it
	for(labloc in 1:length(gridxlabels))
	{
		mxl[,(1:xpoints)[(1:xpoints) > gridxlabels[labloc]]] <- mxl[,(1:xpoints)[(1:xpoints) > gridxlabels[labloc]]] + 2
	}
	
	#Adjust x coordinates of boxes by 2 pixels for each grid mark before it
	for(labloc in 1:length(gridylabels))
	{
		myb[(1:ypoints)[(1:ypoints) > gridylabels[labloc]],] <- myb[(1:ypoints)[(1:ypoints) > gridylabels[labloc]],] + 2
	}

	#Calculate plot dimensions of RSS portion of plot based on size of rectangles and gridlines to be drawn
	pdlx <- min(mxl) 
	pdrx <- max(mxl) + cellsize
	pdby <- min(myb)
	pdty <- max(myb) + cellsize
	
	#Total plot size in pixels
	rssxsize <- pdrx - pdlx
	rssysize <- pdty - pdby 
	
	#bottom, left, top, right
	top_image_height <- 320
	text_height_top <- 20
	text_height_bottom <- 20
	margins_all <- 30
	mar_top <- margins_all + text_height_top + top_image_height
	mar_left <- margins_all + cellsize
	mar_bottom <- margins_all
	mar_right <- margins_all - 5

	pngxsize <- as.integer(rssxsize + mar_left + mar_right)
	pngysize <- as.integer(rssysize + mar_top + mar_bottom + text_height_bottom)

	#Adjust width to include color legend if it has been included.
	if (colorlegend)
	{
		pngxsize <- pngxsize + colorbarxsize + mar_right
		pngysize <- as.integer(max(rssysize, colorbarysize) + text_height_bottom + mar_top + mar_bottom)
	}

	#If width is less than 640 pixels, then increase it to 640 to match with width of other plot.
	if (pngxsize < 640)
	{
		expandx <- 640 - pngxsize
		pngxsize <- 640
		mar_left <- mar_left + as.integer((expandx + 1) / 2)
		mar_right <- mar_right + as.integer(expandx / 2)
	}

	png(filename=outputfilename, bg="white", width = pngxsize, height = pngysize, units="px", type="cairo-png")
	plot.new()

	viewport(x=unit(mar_left, "native"), y=unit(mar_bottom, "native"))

	grid.raster(readPNG(topPlotFileName), x = unit(0, "native"),  y = unit(0, "native"), just=c("left","top"))
	
	grid.rect(
		x=unit(mxl+mar_left, "native"), 
		y=unit(myb+mar_top, "native"), 
		width=unit(cellsize, "native"), 
		height=unit(cellsize, "native"),
		gp=gpar(fill=gcolors, col=gcolors, lwd=0, alpha=1)
	)
	
	label_pos <- 1

	for(x_loc in (gridxlocs+mar_left))
	{
		grid.lines(x = unit(x_loc, "native"),
			  y = unit(c(min(gridylocs) + mar_top, max(gridylocs) + mar_top + 5 ), "native"),
			  arrow = NULL, name = NULL, gp=gpar(col=gridcolor,lwd=unit(2,"native"),lex=1), draw = TRUE)

		if (!is.na(gridxlabels[label_pos]))
		{
			grid.text(gridxlabels[label_pos], x = unit(x_loc, "native"), y = unit( max(gridylocs) + mar_top + 10, "native"),
				just = "top", hjust = NULL, vjust = NULL, rot = 0,
				check.overlap = FALSE, default.units = "native",
				name = NULL, gp = gpar(color="black", cex=1, fontfamily="Times"), draw = TRUE)
		}
		label_pos <- label_pos + 1
	}	
	
	label_pos <- 1
	
	for(y_loc in (gridylocs + mar_top))
	{
		grid.lines(y = unit(y_loc, "native"),
			x = unit(c(min(gridxlocs) + mar_left - 5, max(gridxlocs) + mar_left), "native"),
			gp=gpar(col=gridcolor, lwd=unit(2,"native"),lex=1), 
			draw = TRUE
		)

		if (!is.na(gridylabels[label_pos]))
		{
			grid.text(gridylabels[label_pos], x = unit(mar_left - max(c(8,cellsize/2 + 8)), "native"), y = unit(y_loc, "native"),
				just = "right", rot = 0,
				check.overlap = FALSE, default.units = "native",
				gp = gpar(color="black", cex=1, fontfamily="Times"), 
				draw = TRUE
			)
		}
		
		label_pos <- label_pos + 1
	}
	
	if (colorlegend) 
	{
		cgcellheight <- colorbarysize / numcolors
		cgcellwidth <- colorbarxsize %/% 2

		cgleft <- pngxsize - mar_right - colorbarxsize
		cgtop <- mar_top 
		cgcelltops <- cgtop + cgcellheight * (0:(numcolors-1))

		cggridmarklabels <- c("1.0", "0.9", "0.8", "0.7", "0.6", "0.5", "0.4")
		cgtextlocsarea <- max(cgcelltops) - min(cgcelltops) + 0 * cgcellheight
		cgtextdist <- cgtextlocsarea / (length(cggridmarklabels)-1)
		cgtextylocs <- cgtop + (0:(length(cggridmarklabels)-1)) * cgtextdist + 1

		#Draw color scale grid
		grid.rect(
			x=unit(cgleft, "native"), 
			y=unit(cgcelltops, "native"), 
			width=unit(cgcellwidth, "native"), 
			height=unit(cgcellheight, "native"),
			just=c("left","bottom"),
			gp=gpar(fill=rev(plotcolors),lwd=0, alpha=1)
		)
		
		#Draw color scale area border
		grid.rect(
			x=unit(cgleft, "native"), 
			y=unit(cgtop, "native"), 
			width=unit(cgcellwidth, "native"), 
			height=unit(cgcellheight * numcolors, "native"),
			just=c("left","bottom"),
			gp=gpar(fill=NA, col="black", lwd=unit(2,"native"), alpha=1)
		)
		
		#Place color scale ticks and labels
		for (yloc in cgtextylocs)
		{
			grid.lines(
				x = unit(c(cgleft + colorbarxsize/2,cgleft + colorbarxsize/2 + 5), "native"), 
				y = unit(c(yloc,yloc), "native"),
				gp=gpar(col="black", lwd=unit(2,"native"),lex=1),
				draw = TRUE
			)
		}		
		
		grid.text(cggridmarklabels, 
			x = unit(cgleft + colorbarxsize/2 + 8, "native"), 
			y = unit(cgtextylocs, "native"),
			just = c("left","centre"), 
			check.overlap = FALSE, 
			default.units = "native",
			gp = gpar(color="black", cex=1, fontfamily="Times"),
			draw = TRUE
		)
		  
		#Add missing color legend elements
		grid.rect(
			x = unit(cgleft, "native"), 
			y = unit(max(cgcelltops) + cgcellheight + nacolordivider, "native"), 
			width = unit(cgcellwidth, "native"), 
			height = unit(nacolorbarheight, "native"),
			just = c("left","bottom"),
			gp = gpar(fill=missingvaluecolor, col="black", lwd=unit(2,"native"), alpha=1)
		)
		
		grid.lines(
			x = unit(c(cgleft + cgcellwidth, cgleft + cgcellwidth + 3), "native"), 
			y = unit(c(max(cgcelltops) + cgcellheight + nacolordivider + nacolorbarheight/2, 
				max(cgcelltops) + cgcellheight + nacolordivider + nacolorbarheight/2), "native"), 
			gp = gpar(col="black", lwd=unit(2,"native"),lex=1), 
			draw = TRUE
		)

		grid.text("NA", 
			x = unit(cgleft + cgcellwidth + 5, "native"), 
			y = unit(max(cgcelltops) + cgcellheight + nacolordivider + nacolorbarheight/2, "native"),
			just = c("left","centre"),
			check.overlap = FALSE, 
			default.units = "native",
			gp = gpar(color="black", cex=1, fontfamily="Times"), 
			draw = TRUE
		)
	}

	#Plot title text
	grid.text(
		paste(measure, ":", antibody),
		x = unit(as.integer(pngxsize / 2), "native"), 
		y = unit(as.integer(mar_top - text_height_top), "native"), 
		just = c("center","bottom"),
		default.units = "native",
		gp = gpar(color="black",  cex=1.2, fontface = "bold", fontfamily="Times"),
		draw = TRUE)

	#Footer text (file name)
	grid.text(
		paste("File:", outputfilename),
		x = unit(as.integer(pngxsize / 2), "native"), 
		y = unit(as.integer(pngysize - mar_bottom / 2), "native"), 
		just = c("center","center"),
		default.units = "native",
		gp = gpar(color="black", cex=1, fontfamily="Times"),
		draw = TRUE)
	
	dev.off()
	
	return(c(x = pngxsize, y = pngysize))
}

##-----------------------------------------------------------------------------
## Create the fit graphs and save them as PNG files
.createFitGraphs <- function(rppaset,
                             path,
                             prefix,
							 residualsrotation,
							 majorXDivisions = rppaset@design@majorXDivisions,
							 majorYDivisions = rppaset@design@majorYDivisions
							 ) {
    ## Check arguments
    stopifnot(is.RPPASet(rppaset))
    stopifnot(is.character(path)   && length(path) == 1)
    stopifnot(is.character(prefix) && length(prefix) == 1)

    ## Begin processing
    saved.par <- par(no.readonly=TRUE)
    on.exit(par(saved.par))
    fitdev <- dev.cur()

    ## Use red/yellow/green palette for residual plots.
    ## From RColorBrewer palette RdYlGn
    RYG <- c("#A50026",
             "#D73027",
             "#F46D43",
             "#FDAE61",
             "#FEE08B",
             "#FFFFBF",
             "#D9EF8B",
             "#A6D96A",
             "#66BD63",
             "#1A9850",
             "#006837")
			 
	plotcolors <- RYG
    fitxform <- rppaset@fitparams@xform
    antibodies <- rownames(rppaset@fits)
    for (i in seq_along(antibodies)) {
        antibody <- antibodies[i]
        rppafit <- rppaset@fits[[i]]

		print(paste("Creating fit images for slide ", i, " (", antibody, ")"))
		#print(paste("Creating fit image 1, top part, for slide ", i, " (", antibody, ")"))
        if (!is.RPPAFit(rppafit)) {
			msg <- paste("Slide ", i, " (", antibody, ") will not have fit graphs due to missing or invalid RPPAFit object.", sep="")
            warning(msg)
			write(msg, warningsFileName, append=TRUE)
            next
        }
		measure <- "ResidualsR2"
        main <- .mkPlotTitle(rppafit@measure, antibody)

		#-----------------------------------------------------------------
        ##
        ## First pair of plots
        ##
		par(bg="white", mfrow=c(1, 1))
        ## Plot sigmoid curve graph
		#Write it to a temporary file to be added to the later plot.
        tryCatch(plot(rppafit,
                      main=main,
                      xform=fitxform,
                      xlim=c(-15, 15)),
                 error=function(e) {
					msg <- paste("Slide ", i, " (", antibody, ") unable to plot fit curve due to following issue:", conditionMessage(e), sep="")
					write(msg, warningsFileName, append=TRUE)
                    message(sprintf("cannot plot sigmoid curve for %s", antibody))
                    warning(conditionMessage(e), immediate.=TRUE)
                 })

        topPlotFileName <- sprintf("temp_%s_%s_1.png", prefix, antibody)
        outputfilename <- sprintf("%s_%s_1.png", prefix, antibody)
        dev.copy(png,
                 file.path(path, .portableFilename(topPlotFileName)),
                 width=640,
                 height=320)
		dev.off()

		#print(paste("Creating fit image 1, bottom part, for slide ", i, " (", antibody, ")"))
		rppa <- rppafit@rppa
		residualData <- residuals(rppafit, "r2")

		## Begin processing
		dim.rppa <- dim(rppa)
		
		mx <- max(rppa@data$Main.Col) * max(rppa@data$Sub.Col)
		my <- max(rppa@data$Main.Row) * max(rppa@data$Sub.Row)
			
		m <- matrix(residualData, ncol=mx, nrow=my, byrow=TRUE)

		if (is.na(majorXDivisions))
		{
			majorxpoints <- dim.rppa["Main.Col"]
		}
		else
		{
			majorxpoints <- majorXDivisions
		}
		
		if (majorxpoints <= 1)
		{
			warning(paste("Invalid majorXDivisions setting:", majorxpoints, ". Defaulting to value of 10."))
			majorxpoints <- 10
		}
		
		if (is.na(majorYDivisions))
		{
			majorypoints <- dim.rppa["Main.Row"]
		}
		else
		{
			majorypoints <- majorYDivisions
		}
		
		if (majorypoints <= 1)
		{
			warning(paste("Invalid majorYDivisions setting:", majorypoints, "Defaulting to value of 10."))
			majorypoints <- 10
		}
		
		residualsPlotDimensions <- 
			createResidualsPNG(
				file.path(path, .portableFilename(topPlotFileName)),
				outputfilename,
				measure,
				antibody,
				m,
				majorxpoints,
				majorypoints,
				plotcolors,
				missingvaluecolor = "#AAAAAA",
				gridcolor = "black",
				colorlegend = TRUE,
				cellsize = 3,
				rotation = residualsrotation
			)

        dev.off()
		file.remove(topPlotFileName)
        dev.set(fitdev)

		#-----------------------------------------------------------------
        ##
        ## Second pair of plots
        ##
		#print(paste("Creating fit image 2 for slide ", i, " (", antibody, ")"))

        par(bg="white", mfrow=c(2, 1))

        ## Plot residuals graph
        tryCatch(plot(rppafit,
                      main=main,
                      type="resid",
                      xform=fitxform,
                      xlim=c(-15, 15)),
                 error=function(e) {
					msg <- paste("Slide ", i, " (", antibody, ") unable to plot residuals graph due to following issue:", conditionMessage(e), sep="")
					write(msg, warningsFileName, append=TRUE)
					message(sprintf("cannot plot residuals graph for %s", antibody))
                    warning(conditionMessage(e), immediate.=TRUE)
                 })

        ## Plot steps graph
        tryCatch(plot(rppafit,
                      main=main,
                      type="steps",
                      xform=fitxform,
                      xlim=c(-15, 15)),
                 error=function(e) {
					msg <- paste("Slide ", i, " (", antibody, ") unable to plot cannot plot steps graph due to following issue:", conditionMessage(e), sep="")
					write(msg, warningsFileName, append=TRUE)
					message(sprintf("cannot plot steps graph for %s", antibody))
                    warning(conditionMessage(e), immediate.=TRUE)
                 })

        filename <- sprintf("%s_%s_2.png", prefix, antibody)
        dev.copy(png,
                 file.path(path, .portableFilename(filename)),
                 width=640,
                 height=640)
        dev.off()
        dev.set(fitdev)
    }
	
	return(c(plot_1=residualsPlotDimensions, plot_2 = c(x=640,y=640)))
}


##-----------------------------------------------------------------------------
## Examine system() return code to determine if execution of the shell failed.
.execShellFailed <- if (getRversion() <= "2.12") {
                        function(rc) { rc == 32512 }
                    } else {
                        function(rc) { rc == 127 }
                    }


##-----------------------------------------------------------------------------
## Merge output graph png files with source slide image file, save it as JPG file
.mergeGraphsAndImage <- function(antibody,
                                 prefix,
                                 outputdir,
                                 slideImageName,
								 slideImageRotation = 0,
								 missingSlide = FALSE,
								 #dimensions of two graphs already written to disk
								 graphDimensions = c(plot_1 = c(x=0, y=0), plot_2 = c(x=0, y=0))  
								 ) {
	fitdev <- dev.cur()
	
	## Check arguments
    stopifnot(is.character(antibody)  && length(antibody) == 1)
    stopifnot(is.character(prefix)    && length(prefix) == 1)
    stopifnot(is.character(outputdir) && length(outputdir) == 1)
    stopifnot(is.character(slideImageName) && length(slideImageName) == 1)
	
	rc <- FALSE
	
	tryCatch({
		if (!(slideImageRotation %in% c(0, 90, 180, 270))) {
			msg <- paste("Invalid slideImageRotation value =", 
				slideImageRotation, 
				". Acceptable values are 0, 90, 180, and 270. ",
				"Rotation request ignored.  Default rotation of 0 used.",
				sep="")
			write(msg, warningsFileName, append=TRUE)
			message(msg)
			slideImageRotation = 0
		}

		## Begin processing
		filename <- sprintf("%s_%s_1.png", prefix, antibody)
		pg1 <- file.path(outputdir, .portableFilename(filename))

		filename <- sprintf("%s_%s_2.png", prefix, antibody)
		pg2 <- file.path(outputdir, .portableFilename(filename))

		filename <- sprintf(".%s.temp.png", antibody)
		tempImageName <- file.path(outputdir, .portableFilename(filename))
		
		filename <- sprintf("%s.png", antibody)
		outputImageName <- file.path(outputdir, .portableFilename(filename))


		#The imager package fails to display the pngs if trying to append a 16 bit monochrome tif converted to color and color pngs
		#If we save the slide image as a png or jpg and then reload it, the merge works fine.
		#Converting before rotating or scaling takes much less time than saving and loading after add.color
		lenSlideImageName = nchar(slideImageName)

		if (missingSlide == FALSE && (substr(slideImageName, lenSlideImageName - 3, lenSlideImageName) %in% c(".tif", "tiff"))) {
			img <- readTIFF(slideImageName, native=TRUE, convert=TRUE)
			#writeJPEG(img, target=tempImageName, quality=1)
			writePNG(img, target=tempImageName)
			slideImage <- load.image(tempImageName)
			file.remove(tempImageName)
		} else {
			#Load the slide Image file
			slideImage <- load.image(slideImageName)
		}

		#Get image heights and widths
		slideWidth <- dim(slideImage)[1]
		slideHeight <- dim(slideImage)[2]

		minDim <- min(slideWidth, slideHeight)
		if (minDim > 640)
		{
			slideImage <- imresize(slideImage, scale = 640 / minDim, interpolation = 6)
			#Get new image heights and widths
			slideWidth <- dim(slideImage)[1]
			slideHeight <- dim(slideImage)[2]
		}

		#Rotate slide image if requested
		if (slideImageRotation != 0) {
			slideImage <- imrotate(slideImage, slideImageRotation)
		}

		if (dim(slideImage)[4] == 1) {  #Monochrome Image
			#Color depth will be lost if was initially 16 bit grayscale, being switched to 8 bit rgb.
			slideImage <- add.color(slideImage)
		}

		#Get image heights and widths
		slideWidth <- dim(slideImage)[1]
		slideHeight <- dim(slideImage)[2]

		save.image(slideImage, tempImageName)
		rm(slideImage)

		png1x <- graphDimensions["plot_1.x"]
		png1y <- graphDimensions["plot_1.y"]
		png2x <- graphDimensions["plot_2.x"]
		png2y <- graphDimensions["plot_2.y"]

		#Determine which slide image dimension is greater, width or height. 
		
		if (slideHeight > slideWidth)
		{
			orientation <- 1 # alternate layout: graphs left, slide right
			combinedGraphsWidth <- max(png1x, png2x)
			combinedGraphsHeight <- png1y + png2y
			combinedWidth <- combinedGraphsWidth + slideWidth
			combinedHeight <- max(combinedGraphsHeight, slideHeight)
		}
		else
		{
			orientation <- 0 # default layout: graphs top, slide bottom
			combinedGraphsWidth <- png1x + png2x
			combinedGraphsHeight <- max(png1y, png2y)
			combinedWidth <- max(combinedGraphsWidth, slideWidth)
			combinedHeight <- combinedGraphsHeight + slideHeight
		}
		print(paste("Creating combined output graphic", filename))
		png(filename=outputImageName, bg="white", width = combinedWidth, height = combinedHeight, units="px", type="cairo-png")
		plot.new()
		
		if (orientation == 1)
		{
			grid.raster(readPNG(pg1), x = unit(0, "native"), y = unit(0, "native"), just=c("left","top"), interpolate = FALSE,  
				width = unit(png1x, "native"), height = unit(-png1y, "native"))
				
			grid.raster(readPNG(pg2), x = unit(0, "native"), y = unit(png1y, "native"), just=c("left","top"), interpolate = FALSE, 
				width = unit(png2x, "native"), height = unit(-png2y, "native"))
				
			grid.raster(readPNG(tempImageName), x = unit(combinedGraphsWidth, "native"),  y = unit(0, "native"), interpolate = FALSE, 
				just=c("left","top"), width = unit(slideWidth, "native"), height = unit(-slideHeight, "native"))
		}
		else
		{
			grid.raster(readPNG(pg1), x = unit(0, "native"), y = unit(0, "native"), just=c("left","top"), interpolate = FALSE, 
				width = unit(png1x, "native"), height = unit(-png1y, "native"))
				
			grid.raster(readPNG(pg2), x = unit(png1x, "native"), y = unit(0, "native"), just=c("left","top"), interpolate = FALSE,  
				width = unit(png2x, "native"), height = unit(-png2y, "native"))
			
			grid.raster(readPNG(tempImageName), x = unit(0, "native"), y = unit(combinedGraphsHeight, "native"), just=c("left","top"), 
				interpolate=FALSE, width = unit(slideWidth, "native"), height = unit(-slideHeight, "native"))
		}
		dev.off()
		
		file.remove(tempImageName)

		rc <- TRUE
	},	
	error=function(cond) {
		print("Error merging graphs and slide image")
		message("###stacktrace###")
		dump.frames()
		message("<<<ERROR>>>", cond)
	})
	dev.set(fitdev)

    return(rc)
}


##-----------------------------------------------------------------------------
.plotProbabilityOfGoodSlide <- function(qcprobs,
                                        good.cutoff=0.8) {
    stopifnot(is.numeric(qcprobs)     && length(qcprobs) >= 1)
    stopifnot(is.numeric(good.cutoff) && length(good.cutoff) == 1)
    stopifnot(.isProbability(qcprobs))
    stopifnot(.isProbability(good.cutoff))

    nslides <- length(qcprobs)
    rating <- rep("poor", len=nslides)
    rating[which(qcprobs >= good.cutoff)] <- "good"

    rating.fac <- ordered(rating, levels=c("poor", "good"))
    col.qcprobs <- c("red", "green")
    stopifnot(nlevels(rating.fac) == length(col.qcprobs))

    plot(qcprobs,
         las=1,
         main="Predicted Slide Quality",
         sub=sprintf("#Good = %d, #Poor = %d",
                     ngood <- sum(rating == "good"),
                     npoor <- nslides - ngood),
         type='n',
         xaxt='n',
         xlim=c(1, nslides),
         yaxt='n',
         ylab="Probability of Good Slide",
         ylim=0:1)
    mtext(side=3, paste("cutoff =", good.cutoff))
    axis(1, at=seq_len(nslides))
    axis(2, at=seq(0, 1, by=0.1), las=1)
    rect(xleft=1,
         ybottom=c(0, good.cutoff),
         xright=nslides,
         ytop=c(good.cutoff, 1),
         col=c('lightpink', 'lightgreen'),
         border=NA)
    text(x=(nslides+1)/2,
         y=c(good.cutoff/2, 1-((1-good.cutoff)/2)),
         labels=toupper(levels(rating.fac)),
         cex=2)
    abline(h=good.cutoff)
    points(qcprobs,
           bg=col.qcprobs[rating.fac],
           col='black',
           pch=21)
}


##-----------------------------------------------------------------------------
## Return TRUE if PreFit QC was performed.
ran.prefitqc <- function(rppaset) {
    ## Check arguments
    stopifnot(is.RPPASet(rppaset))

    ## Begin processing
    prefitqcs.tf <- rppaset@completed[, "prefitqc"]
    return(!all(is.na(prefitqcs.tf)))
}


##-----------------------------------------------------------------------------
setMethod("normalize", signature(object="RPPASet"),
          function(object, ...) {
    concs <- .fitSlot(object, "concentrations")
    normparams <- object@normparams
    arglist <- c(list(concs,
                      method=normparams@method),
                 normparams@arglist,
                 ...)

    ## Assemble matrix of concentrations from all fits in object
    do.call(callGeneric, arglist)
})


##-----------------------------------------------------------------------------
## See R FAQ (8.1 How should I write summary methods?)
setMethod("summary", signature(object="RPPASet"),
          function(object,
                   onlynormqcgood=ran.prefitqc(object),
                   ...) {
    dots <- list(...)

    RPPASetSummary(object,
                   onlynormqcgood)
})


##-----------------------------------------------------------------------------
## Provide a convenience function to save fit results to disk
setMethod("write.summary", signature(object="RPPASet"),
          function(object,
                   path,
                   prefix="rppaspace",
                   graphs=TRUE,
				   createcombinedoutputimage=FALSE,
                   imagedir=NULL,
                   onlynormqcgood=ran.prefitqc(object),
				   imageextension=".tif",
				   imagerotation=as.integer(0),
				   residualsrotation=as.integer(0),
				   majorXDivisions=object@design@majorXDivisions,
				   majorYDivisions=object@design@majorYDivisions,
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

    if (is.numeric(graphs)) {
        graphs <- as.logical(graphs)
    }

    if (!is.logical(graphs)) {
        stop(sprintf("argument %s must be logical",
                     sQuote("graphs")))
    } else if (!(length(graphs) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("graphs")))
    }

    if (is.numeric(onlynormqcgood)) {
        onlynormqcgood <- as.logical(onlynormqcgood)
    }

    if (!is.logical(onlynormqcgood)) {
        stop(sprintf("argument %s must be logical",
                     sQuote("onlynormqcgood")))
    } else if (!(length(onlynormqcgood) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("onlynormqcgood")))
    }

    ## Begin processing

    ## Make sure at least one fit exists
    rppafits <- object@fits
    if (all(sapply(rppafits, is.null))) {
        stop("cannot summarize as no fits were stored")
    }

	graphDimensions <- c(plot_1 = c(x=0, y=0), plot_2 = c(x=0, y=0))
    ## Graph fits, if requested
    if (graphs) {
		tryCatch({
			#t <- proc.time()
			#print("Saving fit graphs")
			tm <- proc.time()

			## Save fit graphs
			dev.new(title="Fit Graphs")
			graphDimensions <- .createFitGraphs(object, path, prefix, residualsrotation, majorXDivisions, majorYDivisions)
			#t <- proc.time() - tm
			#print(paste("Fit graph output time:",t[3]))
		},
		error=function(cond) {
			print("Error creating Fit Graphs")
			message("###stacktrace###")
			dump.frames()
			message("<<<ERROR>>>", cond)
		})
		
	}

	if (createcombinedoutputimage){
	    #t <- proc.time()

		pkgimgdir <- system.file("images", package="RPPASPACE")
		print(paste("imagedir:",imagedir))
        if (is.null(imagedir)) {
            ## Assume the tif images are in a sibling directory named "tif"
			print("null imagedir")
            imagedir <- normalizePath(file.path(dirname(path), "tif"))
            if (!dir.exists(imagedir)) {
                ## As last resort, use package directory for missing image
                message(sprintf("image directory unspecified and sibling directory %s does not exist",
                                dQuote(imagedir)))
                imagedir <- pkgimgdir
            }
        }

        if (!is.character(imagedir)) {
            stop(sprintf("argument %s must be character",
                         sQuote("imagedir")))
        } else if (!(length(imagedir) == 1)) {
            stop(sprintf("argument %s must be of length 1",
                         sQuote("imagedir")))
        } else if (!dir.exists(imagedir)) {
            stop(sprintf("directory %s does not exist",
                         dQuote(imagedir)))
        }
		
		if(!(imageextension %in% c(".tif", ".png", ".bmp", ".gif", ".jpg"))) {
			stop(paste("Supported image extensions for slides include .tif, .tiff, .png, .bmp, .jpg.",
			"Only one slide image file type can be used in a single run of slides."))
		}

        ## Merge output graphs with source tiff file for each antibody
        imgfiles <- {
			txtfiles <- sapply(rppafits,
							   function(fit) {
								   if (!is.null(fit) && is.RPPAFit(fit)) {
									   fit@rppa@file
								   } else {
									   NA
								   }
							   })
			txt.re <- "\\.[tT][xX][tT]$"
			sub(txt.re, imageextension, txtfiles)
		}

        ## For each antibody...
        antibodies <- names(rppafits)

		#When using the imager package, do not do graphics merge in parallel or
		#it will probably either take much longer or completely lock up the computer. 
		i<- 0
		
		foreach (i = seq_along(antibodies)) %do% {
			antibody <- antibodies[i]

			rppafit <- rppafits[[i]]
			#cat(c("Testing ", antibody, "\n" ))

			if (!is.null(rppafit) && is.RPPAFit(rppafit)) {
				#print(paste("Merging graphs and image for", antibody))
				flush.console()

				## If no corresponding image exists, substitute "missing" image
				imgfile <- file.path(imagedir, imgfiles[antibody])
				missingSlide <- FALSE

				if (is.na(imgfile) || !file.exists(imgfile)) {
					#imgfile <- file.path(pkgimgdir, "missing_slide.jpg")
					imgfile <- file.path(pkgimgdir, paste("missing_slide", ".jpg", sep=""))
					if (is.na(imgfile) || !file.exists(imgfile)) {
						print("No slide image was provided and the missing_slide.jpg image file can not be found!")
					}
					rotation = imagerotation
					missingSlide <- TRUE
				} else {
					rotation = imagerotation
				}

				## Create merged image
				.mergeGraphsAndImage(antibody, prefix, path, imgfile, rotation, missingSlide, graphDimensions)
				#print(paste("Finished image processing for ", antibody, sep=""))
			}
			else {
				msg <- paste("Slide ", i, " (", antibody, ") will not be have images created due to missing or invalid RPPAFit object.", sep="")
				message(msg)
				write(msg, warningsFileName, append=TRUE)
			}
		}
	    #t <- proc.time() - tm
		#print(paste("Combined output image file creation time:",t[3]))

    }

    ## Write CSV files
    callGeneric(summary(object,
                        onlynormqcgood=onlynormqcgood),
                path,
                prefix)
})


##-----------------------------------------------------------------------------
## Create an RPPA set from a directory of slides.
RPPASet <- function(path,
                    designparams,
                    fitparams,
                    spatialparams=NULL,
                    normparams,
                    doprefitqc=FALSE,
                    parallelClusterSize,
					residualsrotation=as.integer(0),
					warningsFileName="warnings.txt",
					printTimings=TRUE
					) {
	ptm <- proc.time()
	if (printTimings) { print("Starting RPPASet") }
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
    }

    if (!is.RPPADesignParams(designparams)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("designparams"), "RPPADesignParams"))
    }

    if (!is.RPPAFitParams(fitparams)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("fitparams"), "RPPAFitParams"))
    }

    if (!is.null(spatialparams)) {
        if (!is.RPPASpatialParams(spatialparams)) {
            stop(sprintf("argument %s must be object of class %s",
                         sQuote("spatialparams"), "RPPASpatialParams"))
        }
    }

    if (!is.RPPANormalizationParams(normparams)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("normparams"), "RPPANormalizationParams"))
    }

    if (!is.logical(doprefitqc)) {
        stop(sprintf("argument %s must be logical",
                     sQuote("doprefitqc")))
    } else if (!(length(doprefitqc) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("doprefitqc")))
    } else if (is.na(doprefitqc)) {
        doprefitqc <- FALSE
		msg <- sprintf("argument %s converted from NA to FALSE", dQuote(doprefitqc))
        warning(msg, immediate.=FALSE)
		write(msg, warningsFileName, append=TRUE)
    }

	if (printTimings) { 
		t <- proc.time() - ptm
		print(paste("Variables checked:",t[3]))
	}
	
    ##-------------------------------------------------------------------------
    ## Returns the names of all TXT files in directory argument.
    ## :TBD: Should this get the list of slides from a file ('proteinAssay.tsv'
    ## or 'targets.txt') instead of assuming all .txt files are slides?
    getQuantificationFilenames <- function(path) {

		if (printTimings) { 
			t <- proc.time() - ptm
			print(paste("Getting file names:",t[3]))
		}
		
        stopifnot(is.character(path) && length(path) == 1)

        ## Assumes all .txt files in the directory are slides
        txt.re <- "\\.[tT][xX][tT]$"
        txtfiles <- list.files(path=path, pattern=txt.re)
        ## If SuperCurveGUI's input and output directories refer to the same
        ## path, then its settings file in TEXT format could be present...
        settingsfile.tf <- txtfiles %in% "rs-settings.txt"
        txtfiles[!settingsfile.tf]
    }

	if (printTimings) { 
		t <- proc.time() - ptm
		print(paste("Begin processing:",t[3]))
	}
	
    ## Begin processing
    call <- match.call()

    ## Get filenames of slides to process
    slideFilenames <- getQuantificationFilenames(path)
    if (length(slideFilenames) == 0) {
        stop(sprintf("no quantification files found in directory %s",
                     dQuote(path)))
    }
	
	if (printTimings) { 
		t <- proc.time() - ptm
		print(paste("Getting antibody info:",t[3]))
	}
    ab.list <- vector("list", length(slideFilenames))

    ## Fill in missing values with generated defaults
    x.which <- which(sapply(ab.list, is.null))
    txt.re <- "\\.[tT][xX][tT]$"
    for (x in x.which) {
        ab.list[[x]] <- sub(txt.re, "", slideFilenames[x])
    }
	
	if (printTimings) { 
		t <- proc.time() - ptm
		print(paste("Making antibody names unique:",t[3]))
	}
	
    ## Ensure antibody names are unique
    antibodies <- make.unique(abnames <- unlist(ab.list, use.names=FALSE))
    if (!identical(antibodies, abnames)) {
		msg <- "Slide names were not unique - This should not be possible with current slide input method." 
        warning(msg)
		write(msg, warningsFileName, append=TRUE)
    }
    remove(abnames)

    ## Tracking success/failure of each step
    input.tf    <- logical(length(slideFilenames))
    prefitqc.tf <- logical(length(slideFilenames))
    spatial.tf  <- logical(length(slideFilenames))
    fits.tf     <- logical(length(slideFilenames))

    ## Load slides to process
	message(paste("Reading slides from ", path, ".", sep=""))
	
	if (printTimings) { 
		tm <- proc.time()
		t <- proc.time() - ptm
		#print(paste("Reading first slide to determine layout:",t[3]))
	}
	tracking <- NULL
	

	if (printTimings) { 
		tm <- proc.time()
		t <- proc.time() - ptm
		print(paste("Reading slides:",t[3]))
	
	}
	firstGoodSlide <- FALSE
	trackingSource <- NA
	goodSlideCount <- 0
    rppas <- array(list(), length(slideFilenames), list(antibodies))
    for (i in seq_along(slideFilenames)) {

        slideFilename <- slideFilenames[i]
        antibody <- antibodies[i]

        print(paste("Reading", slideFilename))
        flush.console()

		if (printTimings) { 
			t <- proc.time() - ptm
			print(paste("Adding slide :", slideFilename, ":" ,t[3]))
		
		}
        rppa <- tryCatch({
            RPPA(slideFilename,
                 path=path,
                 antibody=antibody,
				 slideNumber=i,
                 tracking=tracking,
                 seriesToIgnore=designparams@seriesToIgnore,
				 warningsFileName
				)
        },
        error=function(e) {
			msg <- paste("Error creating RPPA object for slide ", slideFilename, ":", conditionMessage(e), sep="") 
            message(msg)
			write(msg, warningsFileName, append=TRUE)
            NULL
        })

		if (!is.null(rppa) && is.RPPA(rppa) && is.null(tracking)) {
			firstGoodSlide <- TRUE
			#print("Setting generic tracking object to one constructed using first valid slide." )
			tracking <- rppa@tracking
			trackingSource <- antibody
		} else {
			firstGoodSlide <- FALSE
		}
	
		## Update only on success
        if (is.RPPA(rppa)) {
			rppas[[i]] <- rppa
			input.tf[i] <- TRUE
			goodSlideCount <- goodSlideCount + 1
        } else {
			rppas[[i]] <- NULL
			input.tf[i] <- FALSE
			
			if (!is.null(tracking) ) {
				msg <- paste("Slide ", i, " (", antibody, ") will not be processed due to layout mismatches with first valid slide, ", trackingSource, ".", sep="")
			}
			if (goodSlideCount == 0) {
				msg <- paste("Slide ", i, " (", antibody, ") will not be processed due to format errors.", sep="")
			}
			message(msg)
			warning(msg)
			write(msg, warningsFileName, append=TRUE)
		}

        ## Plot the first possible slide as a quick design check
        if (firstGoodSlide)
        {
			print("Plotting design check graph")
			dev.new(title="Design Check")
			plot(rppa,
			   fitparams@measure)
			firstGoodSlide <- FALSE
        }
    }
	
	if (printTimings) { 
		t <- proc.time() - tm
		print(paste("Slide Reading Time:",t[3]))
	}

    ## This will trigger if no RPPAs exist.
    if (length(rppa) == 0) {
        stop("no slides can be processed")
    }

	if (printTimings) { 
		t <- proc.time() - ptm
		print(paste("Performing pre-fit qc:",t[3]))
	}	
	
    ## Perform pre-fit QC, if enabled
    prefitqcs <- array(list(), length(slideFilenames), list(antibodies))

	doQC <- FALSE
    if (doprefitqc) {
		doQC <- TRUE
		
		#Checking number of control points to see if enough available to do qc or spatial adjustments
		numPosCtrlPoints <- sum(tracking$isPosCtrl, na.rm=TRUE)
		numNegCtrlPoints <- sum(tracking$isNegCtrl, na.rm=TRUE)
		
		#Must have some negative controls
		#Must have at least 3 positive controls per unique dilution point value other than 0 to make a surface, but requiring at least 3
		if (numNegCtrlPoints < 1) { 
			msg <- paste("Quality control calculations requested, but not enough NegCtrl points provided in first valid slide ", trackingSource ,". QC will not be performed.", sep="")
			warning(msg) 
			write(msg, warningsFileName, append=TRUE)
			doQC <- FALSE
		}
		if (numPosCtrlPoints < 3) { 
			msg <- paste("Quality control calculations requested, but not enough PosCtrl or PosCtrl-Noise points provided in first valid slide ", trackingSource ,". QC will not be performed.", sep="")
			warning(msg) 
			write(msg, warningsFileName, append=TRUE)
			doQC <- FALSE
		}
	}

    if (doQC) {
        for (i in seq_along(slideFilenames)) {
            antibody <- antibodies[i]
			
			if (printTimings) { 
				t <- proc.time() - ptm
				print(paste("Slide:", antibody, ":",t[3]))
			}
			
            rppa <- rppas[[i]]
            if (!is.null(rppa)) {
				print(paste("Quality checking slide ", i, " : ", antibody, ".", sep=""))
				flush.console()
				
                prefitqc <- tryCatch({
					#RPPAPreFitQC(rppa, designparams@dilutionsInSeries)
					RPPAPreFitQC(rppa)
				},
				error=function(e) {
					traceback()
					warning(conditionMessage(e))
					write(conditionMessage(e), warningsFileName, append=TRUE)
					NULL
				})
                ## Update only on success
                if (is.RPPAPreFitQC(prefitqc)) {
                    prefitqcs[[i]] <- prefitqc
                    prefitqc.tf[i] <- TRUE
                }
            } else {
				msg <- paste("Slide ", i, " (", antibody, ") will not be assessed for quality control because of missing RPPA object.", sep="")
                warning(msg)
				message(msg)
				write(msg, warningsFileName, append=TRUE)
            }
        }
		
		if (printTimings) { 
			t <- proc.time() - ptm
			print(paste("Plotting goodness of slide values:",t[3]))
		}
		
        ## Plot 'goodness of slide' values
        dev.new(title="Predicted Slide Quality Plot")
        tryCatch({
                     qcprobs <- sapply(prefitqcs, qcprob)
                     qcprobs[is.na(qcprobs)] <- 0
                     .plotProbabilityOfGoodSlide(qcprobs)
                 },
                 error=function(e) {
					msg <- paste("Cannot plot quality control graph because of following error: ", conditionMessage(e), sep="")
					write(msg, warningsFileName, append=TRUE)
                    message(msg)
                    warning(conditionMessage(e), immediate.=TRUE)
                 })
    } else {
        prefitqc.tf <- rep(NA, length(prefitqc.tf))
    }

	if (printTimings) { 
		t <- proc.time() - tm
		print(paste("QC Time:",t[3])) 
	}
	
    ##-------------------------------------------------------------------------
    ## Determines if spatial adjustment processing is warranted
    shouldPerformSpatialAdj <- function(spatialparams, fitparams, tracking, warningsFileName) {
        stopifnot(is.RPPAFitParams(fitparams))
		
		if(is.null(spatialparams)){
			FALSE
		}
		else
		{
			measures <- eval(formals(spatialCorrection)$measure)
			numPosCtrlPoints <- sum(tracking$isPosCtrl, na.rm=TRUE)
			numNegCtrlPoints <- sum(tracking$isNegCtrl, na.rm=TRUE)
			
			#Must have some negative controls
			#Must have at least 3 positive controls per unique dilution point value other than 0 to make a surface, but requiring at least 3
			if (numNegCtrlPoints < 1) { 
				msg <- paste("Spatial adjustments requested, but not enough NegCtrl points provided in first valid slide, ", trackingSource, ".", sep="")
				write(msg, warningsFileName, append=TRUE)
				warning(msg) 
			}
			if (numPosCtrlPoints < 3) {
				msg <- paste("Spatial adjustments requested, but not enough PosCtrl or PosCtrl-Noise points provided in first valid slide, ", trackingSource, ".", sep="")
				write(msg, warningsFileName, append=TRUE)
				warning(msg) 
			}
			
			is.RPPASpatialParams(spatialparams) && (fitparams@measure %in% measures) && (numNegCtrlPoints > 0) && (numPosCtrlPoints > 3)
		}
    }

	if (printTimings) { 
		print("Checking if need to do spatial adjustment")
		tm <- proc.time()
	}
	
    ## Perform spatial adjustment, if enabled
    doSpatialAdj <- shouldPerformSpatialAdj(spatialparams, fitparams, tracking, warningsFileName)
    if (doSpatialAdj) {
	
		if (printTimings) { 
			t <- proc.time() - ptm
			print(paste("Doing spatial adjustment:",t[3]))
		}
		
        if (spatialparams@plotSurface) {
            ## Open new device if surface plots were requested
            dev.new(title="Surface Plots")
            message("*** watch for prompts to plot on R console ***")
        }

        message(paste("Performing spatial adjustments on slides."))
        flush.console()
		
		i <- 0
        if (spatialparams@plotSurface == TRUE || parallelClusterSize == 1) {
			#Parallel plotting of surfaces breaks, so never done in parallel if plotting surfaces.
            spatiallyAdjustedSlidesData <- foreach (i=seq_along(slideFilenames)) %do% {
                doSpatiallyAdjustSlide(antibody = antibodies[i], rppa = rppas[[i]], index=i, 
                spatialparams=spatialparams)
            }
        } else {
			spatiallyAdjustedSlidesData <- foreach (i=seq_along(slideFilenames), .export="doSpatiallyAdjustSlide") %dopar% {
#            spatiallyAdjustedSlidesData <- foreach (i=seq_along(slideFilenames)) %do% {
                doSpatiallyAdjustSlide(antibody = antibodies[i], rppa = rppas[[i]], index=i, 
                spatialparams=spatialparams)
            }
        }

		#Parallelizing the doSpatiallyAdjustSlide calls above causes error messages 
		#at the time to potentially overwrite each other in the warnings text file
		#To get past this problem, the messages are being added to the return list
		#and printed to the warnings file here.
		j <- 0
		foreach (j=seq_along(slideFilenames)) %do% {
	
			if (!is.null(spatiallyAdjustedSlidesData[[j]]$msg)) {
				msg <- paste(spatiallyAdjustedSlidesData[[j]]$msg, sep="")
				warning(msg)
				write(msg, warningsFileName, append=TRUE)
			}
		}

        rppas <- as.array(lapply(spatiallyAdjustedSlidesData, function(x){x$rppa}))

        spatial.tf <- as.logical(lapply(spatiallyAdjustedSlidesData, function(x){x$tf}))
        remove(spatiallyAdjustedSlidesData)
		
		if (printTimings) { 
			t <- proc.time() - ptm
			print(paste("Saving spatial adjustment to file. : ",t[3]))
		}
		
		## Save the spatial adjustments to a file
		i <- 0
		tempSpatial <- foreach (i=seq_along(slideFilenames)) %do% {
			
			if (!is.null(rppas[[i]])) {
				tempAdjMean <- rppas[[i]]@data$Adj.Net.Value
				if (is.null(tempAdjMean)) {
					rppas[[i]]@data$Adj.Net.Value <- rppas[[i]]@data$Net.Value
					msg <- paste("Slide ", i, " (", antibodies[i], 
						") will not be spatially adjusted due to issues with the values of the positive controls on the slide.",
						sep="")
					warning(msg)
					write(msg, warningsFileName, append=TRUE)
				}
				rppas[[i]]@data$Adj.Net.Value - rppas[[i]]@data$Net.Value
			} 
			#else {
			#	#Slides not properly input will have been assigned all values of NA
			#	as.double(0.0)
			#}
		}
		names(tempSpatial) <- antibodies
		
		# Write out the spatial adjustments done to each slide.
		print("Writing spatial adjustments output.")
		write.table(tempSpatial[spatial.tf], file=file.path("spatial_adjustments.tsv"), sep='\t')
		rm(tempAdjMean)
		rm(tempSpatial)
    } else {
		msg <- "Spatial adjustments will not be calculated."
		warning(msg)
		write(msg, warningsFileName, append=TRUE)
        spatial.tf <- rep(NA, length(spatial.tf))
    }

	if (printTimings) { 
		t <- proc.time() - tm
		print(paste("Spatial Adjustment Time:",t[3])) 
	}
	
    ## Create fits
    tm <- proc.time()

	if (printTimings) { 
		t <- proc.time() - ptm
		print(paste("Fitting slides:",t[3]))
	}

    adj.fitparams <- fitparams
    if (doSpatialAdj) {
        message("Performing fits using spatially adjusted measure.")
        adjMeasure <- paste("Adj", fitparams@measure, sep=".")
        adj.fitparams@measure <- adjMeasure
    } else {
		message("Performing fits using original measure in file.")
	}
    fits <- as.array(list(), length(slideFilenames), list(antibodies))
    fits.tf <- as.logical(FALSE, length(slideFilenames))

    slideNamesSeq <- seq_along(slideFilenames)

	tempfits <- foreach (i=slideNamesSeq, .export=c("dofitSlide")) %dopar% 
	{
		if (!is.null(rppas[[i]])) {
			tryCatch(
			{
				if (printTimings) { 
					t <- proc.time() - ptm
					print(paste("Fitting:", antibodies[i], ":" ,t[3]))
				}
				
				tempfit <- dofitSlide(antibody = antibodies[i], rppa = rppas[[i]], index=i, fitparams=adj.fitparams)

				if (is.null(tempfit) || !tempfit$tf){
					msg <- paste("Slide ", i, " (", antibodies[i], ") will not be fitted due to the dofitSlide method not producing a proper fit.", sep="")
					warning(msg)
					list( index=i, fit=NA, tf=FALSE, msg=msg)
				}
				else {
					tempfit
				}
			},
			error=function(e) {
				warning(conditionMessage(e))
				write(conditionMessage(e), warningsFileName, append=TRUE)
				warning(paste("Error fitting slide ", i, " ", antibody))
				msg <- paste("Slide ", i, " (", antibodies[i], ") will not be fitted due to some error in the dofitSlide method.", sep="")
				list( index=i, fit=NA, tf=FALSE, msg=msg)
			})
		}
		else {
			msg <- paste("Slide ", i, " (", antibodies[i], ") will not be fitted due to missing RPPA object.", sep="")
			warning(msg)
			list( index=i, fit=NA, tf=FALSE, msg=msg)
		}
    }

	if (printTimings) { 
		t <- proc.time() - ptm
		print(paste("Updating successful fit curves:",t[3]))
	}
	
    for ( fitentry in tempfits)
    {
		# Write warning messages created in fitting process to warning file and remove entry from tempfit object.  
		# Must be done in non-parallel portion of code to prevent seperate threads from overwriting warning messages.
		if (!is.null(fitentry$msg)) {
			msg <- paste(fitentry$msg, sep="")
			warning(msg)
			write(msg, warningsFileName, append=TRUE)
		}
        index <- fitentry[["index"]]
        fit <- fitentry[["fit"]]
        tf <- fitentry[["tf"]]
		fits[[antibodies[index]]] <- fit
		fits.tf[index] <- tf
    }
	
    remove(tempfits)

    fits <- as.array(fits)

    t <- proc.time() - tm
    print(paste("Curve Fitting Time:",t[3])) 

    ## Create matrix for tracking what succeeded/failed
    completed <- cbind(input.tf,
                       prefitqc.tf,
                       spatial.tf,
                       fits.tf)
    rownames(completed) <- slideFilenames
    colnames(completed) <- names(getStages()[1:4])

	if (printTimings) { 
		t <- proc.time() - ptm
		print(paste("Creating RPPASet object:",t[3]))
	}
	
    ## Create new class
    new("RPPASet",
        call=call,
		design=designparams,
        rppas=rppas,
        spatialparams=spatialparams,
        prefitqcs=prefitqcs,
        fits=fits,
        fitparams=fitparams,
        normparams=normparams,
        completed=completed,
        version=packageDescription("RPPASPACE", fields="Version"),
		residualsrotation=residualsrotation
		)
}

dofitSlide <- function(antibody, rppa, index, fitparams) {
        tf <- FALSE
        fit <- NULL
        print(paste("Fitting slide for ", antibody, ".", sep=""))
        flush.console()

        #rppa <- rppas[[i]]
        if (!is.null(rppa)) {
            fit <- 
                tryCatch({
                    RPPAFitFromParams(rppa,
                        fitparams=fitparams)
                    },
                    error=function(e) {
                        message(conditionMessage(e))
                        print(paste("Error in slide ", index, " ", antibody))
                        NULL
                    })

            ## Update only on success
            if (!is.null(fit) && is.RPPAFit(fit)) {

                #print(c("Setting tf = TRUE for index", index))
                tf <- TRUE
                fit.tf <- TRUE
            }
            else{
                #print(c("Error fitting slide ", index, " ", antibody))
                tf <- FALSE
                fit <- NULL
                fit.tf <- FALSE
            }
        } else {
			msg <- paste("no slide to fit for", antibody)
			warning(msg)
			write(msg, warningsFileName, append=TRUE)
        }
        flush.console()


        list(index=index, fit=fit, tf=tf)
}

doSpatiallyAdjustSlide <- function(antibody, rppa, index, spatialparams){
    tf <- FALSE
	msg <- NULL
    if (!is.null(rppa)) {
		print(paste("Spatially adjusting slide ", index, " : ", antibody, sep =""))
        rppaNew <- tryCatch({
            spatialAdjustmentFromParams(rppa, spatialparams)
        },
        error=function(e) {
			msg <- paste("Slide ", index, " (", antibody, ") will not be spatially adjusted due to the following problem: ", conditionMessage(e),
				"; If this message indicated a surface was not found, there were not at least three positive control points left to form ", 
				"a spatial adjustment surface for the specified dilution level after eliminating the positive control points below calculations ",
				"based on the slide's negative control values.", sep="")
			warning(msg)
			#write(msg, warningsFileName, append=TRUE)
            traceback()
            NULL
        })

        if (!is.null(rppaNew) && is.RPPA(rppaNew)) {

            ncols.list <- lapply(c(rppaNew, rppa),
                function(x) {
                     ncol(df <- slot(x, "data"))
                })

            if (!do.call("identical", ncols.list)) {
                ## Update only on modification
                rppa <- rppaNew
                tf <- TRUE
            }

            remove(ncols.list)
        }
    } 
    else {
		msg <- paste("Slide ", index, " (", antibody, ") will not be spatially adjusted because of missing RPPA object.", sep="")
		message(msg)
    }

    list(rppa = rppa, tf = tf, msg = msg)
}

doQcDataPrefit <- function(antibody, rppa, designparams, prefitqcOld){
    tf <- FALSE
    prefitqc <- prefitqcOld
	msg <- ""

    if (!is.null(rppa)) {
        prefitqcNew <- tryCatch({
            #RPPAPreFitQC(rppa, designparams@dilutionsInSeries)
            RPPAPreFitQC(rppa)
        },
        error=function(e) {
			msg <- paste("Error prefitting qc data for antibody ", antibody, ":", conditionMessage(e))
            message(msg)
            traceback()
            NULL
        })

        ## Update only on success
        if (is.RPPAPreFitQC(prefitqcNew)) {
            prefitqc <- prefitqcNew
            tf <- TRUE
        }
        else {
            prefitqc <- prefitqcOld
        }
    } 
    else {
		msg <- paste("No slide to quality check for", antibody)
		warning(msg)
    }

    list(prefitqc = prefitqc, tf = tf, msg=msg)
}
