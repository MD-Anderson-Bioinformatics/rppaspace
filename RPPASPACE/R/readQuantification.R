###
### $Id: readQuantification.R
###

##-----------------------------------------------------------------------------
## Reads slide input files in standard format and returns data frame containing
## the desired information
readQuantification <- function(conn) {
    ## Check arguments
    if (!inherits(conn, "connection")) {
        stop(sprintf("argument %s must be connection",
                     sQuote("conn")))
    } else if (!isOpen(conn, "r")) {
        stop(sprintf("connection %s not open for read",
                     dQuote(summary(conn)$description)))
    }

	#TODO: Add check for type of object for dimension size and valid values
#	if (!typeOf(expectedDimensions) 
#	expectedDimensions
	
    ## Begin processing
	reqdNames <- c("Order", 
		"Main.Row", "Main.Col", "Sub.Row", "Sub.Col", 
		"Series.Id", "Spot.Type", "Dilution", 
		"Net.Value", "Raw.Value", "Background.Value", 
		"Spot.X.Position", "Spot.Y.Position", "Original.Order")

    ## Read data from file
	slideData <- tryCatch(
		read.table(conn, quote="", header=TRUE, sep="\t", stringsAsFactors=FALSE ),
		error=function(e) {
			badformat <- "more columns than column names"
			if (conditionMessage(e) == badformat) {
				stop(sprintf('Data does not follow standard format: %s', badformat), call.=FALSE)
			} else {
				stop(e)
			}
		}
	)

	slideColNames <- colnames(slideData)
	
	if (!(all(reqdNames %in% slideColNames))) {
		missingNames <- reqdNames[!reqdNames %in% colnames(slideData)]
		stop(sprintf(ngettext(length(missingNames),
			"argument %s missing required column: %s",
			"argument %s missing required columns: %s"),
			sQuote("rppa"),
			paste(missingNames, collapse=", ")))
	}
	
	extraNames = list()
	if (length(reqdNames) < length(slideData)) {
		for(colName in slideColNames[!slideColNames %in% reqdNames]) {
			if(!startsWith(colName, "Original.")) {
				extraNames <- c(extraNames, colName)
				slideData[[colName]] <- NULL
			}
		}
		if (length(extraNames) > 0) {
			message("Non-standard columns in input must start with text 'Original.'")
			message(paste(
				"Extra columns detected and removed from input: ",
				paste(extraNames, collapse=", ")
			))
		}
	}

	slideData[, "Order"] 	<- as.integer(slideData[, "Order"])
	slideData[, "Main.Row"] <- as.integer(slideData[, "Main.Row"])
	slideData[, "Main.Col"] <- as.integer(slideData[, "Main.Col"])
	slideData[, "Sub.Row"] 	<- as.integer(slideData[, "Sub.Row"])
	slideData[, "Sub.Col"] 	<- as.integer(slideData[, "Sub.Col"])
	slideData[, "Series.Id"] <- as.integer(slideData[, "Series.Id"])
	
	slideData[, "Spot.Type"] <- as.character(slideData[, "Spot.Type"])
	slideData[, "Dilution"] <- as.numeric(slideData[, "Dilution"])

	slideData[, "Net.Value"] <- as.numeric(slideData[, "Net.Value"])
	slideData[, "Raw.Value"] <- as.numeric(slideData[, "Raw.Value"])
	slideData[, "Background.Value"] <- as.numeric(slideData[, "Background.Value"])
	slideData[, "Spot.X.Position"] <- as.numeric(slideData[, "Spot.X.Position"])
	slideData[, "Spot.Y.Position"] <- as.numeric(slideData[, "Spot.Y.Position"])
	
	slideData[, "Original.Order"] <- as.integer(slideData[, "Original.Order"])
	
	spotTypesInFile <- unique(slideData$Spot.Type)
	invalidSpotTypesInFile <- setdiff(spotTypesInFile, spottype.all)
	if (length(invalidSpotTypesInFile) > 0) {
		#Check for correct spot types in wrong case and correct
		correctableSpots <- intersect(tolower(invalidSpotTypesInFile), spottype.lower)
		if (length(correctableSpots) > 0) {
			slideData[, "Spot.Type.Lower"] <- tolower(slideData[, "Spot.Type"])
			
			#If the lower case version of spot type matches but not already in 
			# expected case then change to expected case
			for(spotType in spottype.all) {
				lowerSpotType <- tolower(spotType)
				if(lowerSpotType %in% correctableSpots)
				{
					slideData[slideData$"Spot.Type.Lower" == lowerSpotType,]$Spot.Type <- spotType
				}
			}
		
			#Remove temporary column
			slideData$"Spot.Type.Lower" <- NULL
		}
	}
	spotTypesInFile <- unique(slideData$Spot.Type)
	invalidSpotTypesInFile <- setdiff(spotTypesInFile, spottype.all)
	
	if (length(invalidSpotTypesInFile) > 0) {
		msg <- paste("Invalid spot types found in file: ",  invalidSpotTypesInFile, sep="")
		stop(msg)
	}	

    slideData
}
