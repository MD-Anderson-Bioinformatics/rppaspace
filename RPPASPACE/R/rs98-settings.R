###
### $Id: rs98-settings.R
###


##=============================================================================
setClassUnion("OptionalString", c("character", "NULL"))

setClass("RPPASPACESettings",
	representation(
		txtdir="Directory",
		imgdir="OptionalDirectory",
		outdir="Directory",
		designparams="RPPADesignParams",
		spatialparams="OptionalRPPASpatialParams",
		doprefitqc="logical",
		fitparams="RPPAFitParams",
		normparams="RPPANormalizationParams",
		onlynormqcgood="logical",
		seriesToIgnore="OptionalList",
		parallelClusterSize="integer", 
		createcombinedoutputimage="logical",
		imageextension="character", 
		imagerotation="integer",
		residualsrotation="integer",
		warningsFileName="character",
		errorsFileName="character"
	),
	prototype(
		parallelClusterSize=as.integer(1), 
		createcombinedoutputimage=as.logical(FALSE),
		imageextension=".tif", 
		imagerotation=as.integer(0),
		residualsrotation=as.integer(0),
		warningsFileName="warnings.txt",
		errorsFileName="errors.txt"
	)
)


##-----------------------------------------------------------------------------
## Invoked by validObject() method.
validRPPASPACESettings <- function(object) {

    #cat("validating", class(object), "object", "\n")
    msg <- NULL

    ## Validate txtdir slot
    {
        path <- object@txtdir@path

        ## Ensure directory contains TEXT files
        txt.re <- "\\.*[tT][xX][tT]$"
        txtfiles <- list.files(path, pattern=txt.re)
        if (length(txtfiles) == 0) {
            msg <- c(msg, paste("Slide directory (", path, ") contains no slide files with 'txt' extension.", sep = ""))
        }
    }

    ## Validate imgdir slot
    {
        if (!is.null(object@imgdir)) {
            path <- object@imgdir@path
			ext <- tolower(object@imageextension)
			image.re <- "\\.*[tT][iI][fF]{1,2}$"

            ## Issue warning if image directory does not contain files of designated type with names matching slide names.
			if (ext == ".jpg") {
				image.re <- "\\.*[jJ][pP][gG]{1,2}$"
			} else if (ext == ".png") {
				image.re <- "\\.*[pP][nN][gG]{1,2}$"
			} else if (ext == ".bmp") {
				image.re <- "\\.*[bB][mM][pP]{1,2}$"
			} else if (ext == ".gif") {
				image.re <- "\\.*[gG][iI][fF]{1,2}$"
			}
			
            imagefiles <- list.files(path, pattern=image.re)
            if (length(imagefiles) == 0) {
                warning(sprintf("image directory %s contains no %s files. Default missing image file graphic will be used for all output images.",
                                dQuote(path), object@imageextension))
            } else {
                ## :TODO: Do they correspond to ANY of the TEXT files?
            }
        }
    }

    ## Validate outdir slot
    {
        path <- object@outdir@path

        ## Ensure directory is writable
        if (!dir.writable(path)) {
            msg <- c(msg, "Output directory is not writable.")
        }
    }

    ## Validate onlynormqcgood slot against doprefitqc slot
    if (object@onlynormqcgood && !object@doprefitqc) {
        msg <- c(msg, "Cannot normalize only good slides unless qc performed")
    }

    if (object@parallelClusterSize < 1 || object@parallelClusterSize > detectCores(all.tests = FALSE, logical = TRUE)) {
        msg <- c(msg, paste(
            "parallelClusterSize must be 1 or greater and no greater than the number of cores on the system (",
            detectCores(all.tests = FALSE, logical = TRUE), ")"))
    }
	
	## Validate warningsFileName and errorsFileName (Names files to write warnings and errors to, respectively)
	if (is.null(object@warningsFileName)) {
		msg <- c(msg, "Must provide valid file name for warnings output file name.")
	}
	
	if (is.null(object@errorsFileName)) {
		msg <- c(msg, "Must provide valid file name for errors output file name.")
	}

    ## Pass or fail?
    if (is.null(msg)) {
        TRUE
    } else {
		if (!is.null(object@warningsFileName)) {
			write(msg, warningsFileName, append=TRUE)
		}
        msg
    }
}

setValidity("RPPASPACESettings", validRPPASPACESettings)


##-----------------------------------------------------------------------------
is.RPPASPACESettings <- function(x) {
    is(x, "RPPASPACESettings")
}


##-----------------------------------------------------------------------------
## Generator method
RPPASPACESettings <- function(txtdir,
							imgdir,
							outdir,
							designparams,
							fitparams,
							spatialparams=NULL,
							normparams,
							doprefitqc=FALSE,
							onlynormqcgood=doprefitqc,
							parallelClusterSize=as.integer(1), 
							createcombinedoutputimage=FALSE,
							imageextension=".tif", 
							imagerotation=as.integer(0),
							residualsrotation=as.integer(0),
							warningsFileName="warnings.txt",
							errorsFileName="errors.txt"
							) {

    ## Check arguments
    if (!is.character(txtdir)) {
        stop(sprintf("argument %s must be character",
                     sQuote("txtdir")))
    }

    if (!is.null(imgdir)) {
        if (!is.character(imgdir)) {
            stop(sprintf("argument %s must be character",
                         sQuote("imgdir")))
        }
    }

    if (!is.character(outdir)) {
        stop(sprintf("argument %s must be character",
                     sQuote("outdir")))
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
    }

    if (!is.logical(onlynormqcgood)) {
        stop(sprintf("argument %s must be logical",
                     sQuote("onlynormqcgood")))
    } else if (!(length(onlynormqcgood) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("onlynormqcgood")))
    }

    if (!is.logical(createcombinedoutputimage)) {
        stop(sprintf("argument %s must be logical",
                     sQuote("createcombinedoutputimage")))
    } else if (!(length(createcombinedoutputimage) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("createcombinedoutputimage")))
    }

    ## Create new class
    new("RPPASPACESettings",
        txtdir=as(txtdir, "Directory"),
        imgdir=if (!is.null(imgdir)) as(imgdir, "Directory") else NULL,
        outdir=as(outdir, "Directory"),
        designparams=designparams,
        spatialparams=spatialparams,
        doprefitqc=doprefitqc,
        fitparams=fitparams,
        normparams=normparams,
        onlynormqcgood=onlynormqcgood,
        parallelClusterSize=parallelClusterSize, 
		createcombinedoutputimage=createcombinedoutputimage,
		imageextension=imageextension, 
		imagerotation=imagerotation,
		residualsrotation=residualsrotation,
		warningsFileName=warningsFileName,
		errorsFileName=errorsFileName
	)
}


##-----------------------------------------------------------------------------
setMethod("write.summary", signature(object="RPPASPACESettings"),
          function(object,
                   path=as(object@outdir, "character"),
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

    ##---------------------------------------------------------------------
    makeFileHeader <- function(string) {
        stopifnot(is.character(string) && length(string) == 1)

        paste("###",
              paste("###", string),
              "###",
              "\n",
              sep="\n")
    }


    ## Begin processing
    version <- packageDescription("RPPASPACE", fields="Version")
    cat(makeFileHeader("RPPASPACE settings"),
        paramString(object),
        paste("rppaspace version:", version), "\n",
        "\n",  # blank line at EOF
        sep="",
        file=file.path(path, "rs-settings.txt"))

    invisible(NULL)
})


##-----------------------------------------------------------------------------
## Returns a string representation of this instance. The content and format of
## the returned string may vary between versions. Returned string may be
## empty, but never null.
setMethod("paramString", signature(object="RPPASPACESettings"),
          function(object,
                   designparams.slots,
                   fitparams.slots,
                   spatialparams.slots,
                   normparams.slots,
                   ...) {
    if (missing(designparams.slots)) {
        designparams.slots <- c("center",
                                "designfile")
    }

    if (missing(fitparams.slots)) {
        fitparams.slots <- c("measure",
                             "model",
                             "method",
                             "trim",
                             "ci",
                             "ignoreNegative",
                             "warnLevel")
    }

    if (missing(spatialparams.slots)) {
        spatialparams.slots <- c("cutoff",
                                 "k",
                                 "gamma",
                                 "plotSurface")
    }

    if (missing(normparams.slots)) {
        normparams.slots <- c("name",
                              "method",
                              "arglist")
    }


    ##---------------------------------------------------------------------
    indent <- function(params.text,
                       indention="  ") {
        paste(unlist(lapply(strsplit(params.text, '\n'),
                            function(x, indention) {
                                paste(indention, x)
                            },
                            indention)),
              collapse="\n")
    }


    ## Handle unspecified image directory
    imgdir <- if (!is.null(object@imgdir)) {
                  object@imgdir@path
              } else {
                  NULL
              }

    ## Handle parameters
    designparams  <- paramString(object@designparams, designparams.slots)
    spatialparams <- if (!is.null(object@spatialparams)) {
                         paramString(object@spatialparams, spatialparams.slots)
                     } else {
                         NULL
                     }
    fitparams     <- paramString(object@fitparams, fitparams.slots)
    normparams    <- paramString(object@normparams, normparams.slots)

    ## Create param string
    paste(sprintf("txtdir: %s\n", shQuote(object@txtdir@path)),
          sprintf("imgdir: %s\n", shQuote(imgdir)),
          sprintf("outdir: %s\n", shQuote(object@outdir@path)),
          sprintf("designparams:\n%s\n", indent(designparams)),
          if (!is.null(spatialparams)) {
              sprintf("spatialparams:\n%s\n", indent(spatialparams))
          } else {
              sprintf("dospatialadj: %s\n", FALSE)
          },
          if (!is.null(object@doprefitqc)) {
              sprintf("doprefitqc: %s\n", object@doprefitqc)
          },
          sprintf("fitparams:\n%s\n", indent(fitparams)),
          sprintf("normparams:\n%s\n", indent(normparams)),
          sprintf("onlynormqcgood: %s\n", object@onlynormqcgood),
		  sprintf("parallelClusterSize: %s\n", shQuote(object@parallelClusterSize)),
		  sprintf("createcombinedoutputimage: %s\n", object@createcombinedoutputimage),
		  sprintf("imageextension: %s\n", shQuote(object@imageextension)),
		  sprintf("imagerotation: %s\n", shQuote(object@imagerotation)),
		  sprintf("residualsrotation: %s\n", shQuote(object@residualsrotation)),
		  sprintf("warningsFileName: %s\n", shQuote(object@warningsFileName)),
		  sprintf("errorsFileName: %s\n", shQuote(object@errorsFileName)),
          sep="")
})


##-----------------------------------------------------------------------------
## Returns list of prerequisite packages based on requested processing.
getPrerequisitePackages <- function(settings) {
    ## Check arguments
    if (!is.RPPASPACESettings(settings)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("settings"), "RPPASPACESettings"))
    }

    ## Begin processing

    ## Get model-specific prerequisites
    model.prereqs <- switch(EXPR=settings@fitparams@model,
                            cobs="cobs",
                            logistic="boot")
    prerequisites <- model.prereqs

    ## Get fitmethod-specific prerequisites
    method.prereqs <- switch(EXPR=settings@fitparams@method,
                             nlrq="quantreg",
                             nlrob="robustbase")
    prerequisites <- c(prerequisites, method.prereqs)

    ## Get processing-specific prerequisites
    if (!is.null(settings@spatialparams)) {
        prerequisites <- c(prerequisites, "mgcv")
    }

    if (settings@doprefitqc) {
        prerequisites <- c(prerequisites, "timeDate")
    }

    prerequisites
}

