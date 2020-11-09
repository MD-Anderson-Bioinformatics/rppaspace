###
### $Id: rs99-fitcurve.R
###


##-----------------------------------------------------------------------------
fitCurveAndSummarizeFromSettings <- function(settings) {
    ## Check arguments
    if (!is.RPPASPACESettings(settings)) {
        stop(sprintf("argument %s must be object of class %s",
                     sQuote("settings"), "RPPASPACESettings"))
    }
    validObject(settings, complete=TRUE)  ## Invokes stop() if invalid

    ## Begin processing
    txtdir <- as(settings@txtdir, "character")
    outdir <- as(settings@outdir, "character")
    imgdir <- if (is.Directory(settings@imgdir)) {
                  as(settings@imgdir, "character")
              } else {
                  NULL
              }

    rppasetArgs <- list(path=txtdir,
                        designparams=settings@designparams,
                        fitparams=settings@fitparams,
                        spatialparams=settings@spatialparams,
                        normparams=settings@normparams,
                        doprefitqc=settings@doprefitqc,
                        parallelClusterSize=settings@parallelClusterSize,
						warningsFileName = warningsFileName
						)

    ## Set up code to work in parallel
    clusterSize <- settings@parallelClusterSize
    if (clusterSize > 1) {
        print(paste("Using", clusterSize, "cores" )) 
    }
    cluster <- makeCluster(clusterSize, outfile=file.path(outdir, "/parallel.log"))
    registerDoParallel(cluster)

    ## Perform analysis
    rppaset <- do.call(RPPASet, rppasetArgs)

    stopCluster(cluster)

    ## Save results (as rppaset takes forever to generate)
    rda.filename <- "rs-rppaset.RData"
    save(rppaset, file=file.path(outdir, rda.filename))

    ## Summarize the results
    temp.summary <- write.summary(rppaset,
                  path=outdir,
                  graphs=TRUE,
				  createoutputjpg=settings@createoutputjpg,
                  imagedir=imgdir,
                  onlynormqcgood=settings@onlynormqcgood,
				  imageextension=settings@imageextension,
				  imagerotation=settings@imagerotation,
				  residualsrotation=settings@residualsrotation
				  )


    return(temp.summary)
}

