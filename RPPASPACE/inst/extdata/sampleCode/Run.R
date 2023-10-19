#analysishome <- [Insert directory here]
analysishome <- "./"

## Pathnames (preferred layout)
txtdir <- file.path(analysishome, "txt" )
imgdir <- file.path(analysishome, "img" )
outdir <- file.path(analysishome, "out")
number_cpus_to_use <- 3

warningsFileName <- "warnings.txt"
errorsFileName <- "errors.txt"

dir.create(outdir)
setwd(outdir)

if (!dev.interactive()) {
    options(device="x11")
}
print("===========================================================")
print(paste("Starting RPPASPACE run in ", analysishome))

## Create settings
designparams <- RPPADesignParams(center=FALSE,
                                 seriesToIgnore=list(),
                                 majorXDivisions=as.integer(10),
                                 majorYDivisions=as.integer(10)
                                 )

spatialparams <- RPPASpatialParams(cutoff=0.8,
                                   k=100,
                                   gamma=0.1,
                                   plotSurface=FALSE)

fitparams <- RPPAFitParams(measure="Net.Value",
                           method="nls",
                           model="cobs",
                           trim=2,
                           ci=FALSE,
                           ignoreNegative=FALSE,
                           warnLevel=-1
                           )

normparams <- RPPANormalizationParams(method="none")

settings <- RPPASPACESettings(txtdir=txtdir,
                               imgdir=imgdir,
                               outdir=outdir,
                               designparams=designparams,
                               spatialparams=spatialparams,
                               doprefitqc=TRUE,
                               fitparams=fitparams,
                               normparams=normparams,
                               onlynormqcgood=FALSE,
                               imageextension=".jpg",
                               warningsFileName=warningsFileName,
                               parallelClusterSize=as.integer(number_cpus_to_use),
                               createcombinedoutputimage=FALSE)
write.summary(settings)

## Process slides
tryCatch({
    ptm <- proc.time()
    sapply(RPPASPACE:::getPrerequisitePackages(settings),
           function(pkgname) {
               do.call("library", list(package=pkgname))
           })
    fitCurveAndSummarizeFromSettings(settings)
    t <- proc.time() - ptm
    print(paste("Total time:",t[3]))
    print(sessionInfo())
},
error=function(cond) {
	message("###stacktrace###")
	dump.frames()
	invisible(sapply(names(last.dump),
					 function(acall) {
						 message(paste("   ", acall))
					 },
					 USE.NAMES=TRUE))
	message("<<<ERROR>>>", cond)
})
graphics.off()

setwd(analysishome)
