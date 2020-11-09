###
### $Id: AllGenerics.R
###


##
## S3 or non-generics converted to S4
##

## defined for: FitClass, LogisticFitClass, RPPAFit
if (!isGeneric("coef")) {
    setGeneric("coef",
               function(object, ...) standardGeneric("coef"))
}

## defined for: FitClass, LogisticFitClass, RPPAFit
if (!isGeneric("coefficients")) {
    setGeneric("coefficients",
               function(object, ...) standardGeneric("coefficients"))
}

## defined for: FitClass, LoessFitClass, CobsFitClass, LogisticFitClass, RPPAFit
if (!isGeneric("fitted")) {
    setGeneric("fitted",
               function(object, ...) standardGeneric("fitted"))
}

## defined for: RPPAFit
if (!isGeneric("hist")) {
    setGeneric("hist",
               function(x, ...) standardGeneric("hist"))
}

## defined for: RPPA, RPPAFit
if (getRversion() == "2.15.0") {
    ## R-2.15.0 was broken in respect to this generic. It was patched but damage done
    setGeneric("image")
} else {
    if (!isGeneric("image")) {
        setGeneric("image",
                   function(x, ...) standardGeneric("image"))
    }
}

## defined for: RPPAFit
if (!isGeneric("plot")) {
    setGeneric("plot",
               function(x, y, ...) standardGeneric("plot"))
}

## defined for: RPPAFit
if (!isGeneric("resid")) {
    setGeneric("resid",
               function(object, ...) standardGeneric("resid"))
}

## defined for: RPPAFit
if (!isGeneric("residuals")) {
    setGeneric("residuals",
               function(object, ...) standardGeneric("residuals"))
}

## defined for: RPPA, RPPAFit, RPPASet, DS5RPPAPreFitQC
if (!isGeneric("summary")) {
    setGeneric("summary",
               function(object, ...) standardGeneric("summary"))
}

##
## Borrowed generic definitions
##

## also defined by: affy
## defined for: RPPASet, matrix, data.frame
if (!isGeneric("normalize")) {
    setGeneric("normalize",
               function(object, ...) standardGeneric("normalize"))
}

##
## Brand new generics
##

## defined for: FitClass, LoessFitClass, CobsFitClass, LogisticFitClass
if (!isGeneric("fitSeries")) {
    setGeneric("fitSeries",
               function(object, ...) standardGeneric("fitSeries"))
}

## defined for: FitClass, LoessFitClass, CobsFitClass, LogisticFitClass
if (!isGeneric("fitSlide")) {
    setGeneric("fitSlide",
               function(object, ...) standardGeneric("fitSlide"))
}

## defined for: RPPADesignParams, RPPAFitParams, RPPASpatialParams,
##              RPPANormalizationParams, RPPASPACESettings
if (!isGeneric("paramString")) {
    setGeneric("paramString",
               function(object, ...) standardGeneric("paramString"))
}

## defined for: RPPAPreFitQC, DS5RPPAPreFitQC
if (!isGeneric("qcprob")) {
    setGeneric("qcprob",
               function(object, ...) standardGeneric("qcprob"))
}

## defined for: FitClass, LoessFitClass, CobsFitClass, LogisticFitClass
if (!isGeneric("trimConc")) {
    setGeneric("trimConc",
               function(object, ...) standardGeneric("trimConc"))
}

## defined for: RPPASet, RPPASetSummary, RPPASPACESettings
if (!isGeneric("write.summary")) {
    setGeneric("write.summary",
               function(object, ...) standardGeneric("write.summary"))
}

