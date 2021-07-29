#Converts MD Anderson ArrayPro formatted files into the Standard format used by RPPASPACE

#List of columns in RPPASPACE standard quantification files
outNames <- c("Order", 
	"Main.Row", "Main.Col", "Sub.Row", "Sub.Col", 
	"Series.Id", "Spot.Type", "Dilution", 
	"Net.Value", "Raw.Value", "Background.Value", 
	"Spot.X.Position", "Spot.Y.Position", "Original.Order")

#List of column names in ArrayPro files
reqdNames <- c("Main Row", "Main Col", "Sub Row", "Sub Col", 
	"Gene ID", "Slide", "Control", "Dilution", 
	"Net intensity (mean) {A}", "Raw intensity (mean) {A}", "Background (mean) {A}", 
	"Net intensity (med) {A}", "Raw intensity (med) {A}", "Background (med) {A}",
	"Net intensity (sum) {A}", "Raw intensity (sum) {A}", "Background (sum) {A}",
	"Raw intensity (cv) {A}", "Background (cv) {A}",
	"Spot position X {A}", "Spot position Y {A}",
	"Net intensity (np) {A}", "Background (np) {A}",
	"Area (sum) {A}", "Aspect {A}",
	"Radius Ratio {A}", "Ignored cells")

##-------------------------------------------------------------------------
## Returns the names of all TXT files in directory argument.
getQuantificationFilenames <- function(path) {
	stopifnot(is.character(path) && length(path) == 1)

	## Assumes all .txt files in the directory are slides
	txt.re <- "\\.[tT][xX][tT]$"
	txtfiles <- list.files(path=path, pattern=txt.re)
	## If SuperCurveGUI's input and output directories refer to the same
	## path, then its settings file in TEXT format could be present...
	settingsfile.tf <- txtfiles %in% "sc-settings.txt"
	txtfiles[!settingsfile.tf]
}


#Note: This assumes that multiple sets of input files are in a single directory (in this case, /arrapro_format)
# and sets in that directory each have a directory named Set### where ### is a three digit number, and that
# the actual text quantification files for that set are in a txt directory beneath the set directory.
# The output for each set will respectively be written to a standard_format/Set###/txt directory.

for (set in c(100:102 )) {   #Sets 100 through 102

	inDir <- paste("/arraypro_format/Set", set, "/txt/", sep="")
	outDir <- paste("/standard_format/Set", set, "/txt/", sep="")
	dir.create(outDir, recursive = TRUE)

	slideFilenames <- getQuantificationFilenames(inDir)

	for (i in seq_along(slideFilenames)) {

		slideFilename <- slideFilenames[i]
		## Read data from file
		message(paste("reading", paste(inDir, slideFilename, sep="")))
		slideData <- tryCatch(
			read.table(paste(inDir, slideFilename, sep=""), quote="", header=TRUE, sep="\t", stringsAsFactors=FALSE ),
			error=function(e) {
				badformat <- "more columns than column names"
				if (conditionMessage(e) == badformat) {
					stop(sprintf('Data does not follow standard format: %s', badformat), call.=FALSE)
				} else {
					stop(e)
				}
			}
		)
		#Create output slide in standard format
		slideOut <- data.frame(
			Order = slideData$X,
			Main.Row = slideData$"Main.Row",
			Main.Col = slideData$"Main.Col",
			Sub.Row = slideData$"Sub.Row",
			Sub.Col = slideData$"Sub.Col",
			Series.Id = slideData$"Gene.ID",
			Spot.Type = rep(as.character("Sample"), nrow(slideData)),
			Dilution = slideData$"Dilution",
			Net.Value = slideData$"Net.intensity..mean...A.",
			Raw.Value = slideData$"Raw.intensity..mean...A.",
			Background.Value = slideData$"Background..mean...A.",
			Spot.X.Position = slideData$"Spot.position.X..A.",
			Spot.Y.Position = slideData$"Spot.position.Y..A.",
			Original.Order = slideData$X,
			stringsAsFactors = FALSE
		)
		#Set positive control SpotType and Series.Id
		posCtrlPoints <- slideData$"Gene.ID" == "p-Ctrl"
		slideOut$Spot.Type[posCtrlPoints==TRUE] <- "PosCtrl"
		posCtrlSeriesIds <- list()
		for (i in 0:7)
		{
			posCtrlSeriesIds <- c(unlist(posCtrlSeriesIds), c(rep((1057+i*12):(1068+i*12), 5)))
		}
		
		slideOut$Series.Id[posCtrlPoints == TRUE] <- as.character(posCtrlSeriesIds)

		#Set negative control SpotType
		negCtrlPoints <- slideData$"Gene.ID" == "n-Ctrl"
		slideOut$Spot.Type[negCtrlPoints] <- "NegCtrl"
		slideOut$Series.Id[negCtrlPoints] <- "0"
		#browser()
		
		#Reorder to old style Superslide order so QC metric slopes will be calculated as expected
		## Logical dimensions of SuperSlide format
		nmainrow <- 4
		nmaincol <- 12
		nsubrow  <- 11
		nsubcol  <- 11

		nspot.mr <- nmaincol * nsubrow * nsubcol  # number of spots in main row

		## Ensure file is actually SuperSlide single subgrid layout
		dim.singlesubgrid <- as.integer(c(1,
										  1,
										  (nmainrow*nsubrow),
										  (nmaincol*nsubcol)))
		dim.data.df <- c(max(slideOut$Main.Row),
						 max(slideOut$Main.Col),
						 max(slideOut$Sub.Row),
						 max(slideOut$Sub.Col))
		if (!identical(dim.singlesubgrid, dim.data.df)) {
			stop(sprintf("dim of file %s (%s) does not match SuperSlide single subgrid (%s)",
						 dQuote(pathname),
						 paste(dim.data.df, collapse="x"),
						 paste(dim.singlesubgrid, collapse="x")))
		}

		slideOut$Order <- 1:5808
		
		print(paste("Writing:", paste(outDir, slideFilename, sep="")))
		write.table(slideOut, paste(outDir, slideFilename, sep=""), sep="\t", quote=FALSE, row.names=FALSE)
	}
}