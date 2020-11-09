###
### Package spell checking configuration file
###
### $Id: defaults.R 976 2015-07-27 01:10:08Z proebuck $
###


##
## PREREQ: Install GNU Aspell and appropriate language dictionary.
##

##
## Spellcheck the meat of the package as follows:
## 
## R> library(utils)
## R> pkgname <- "mypackage"           # replace value with name of package
## R> pkg.pathname <- system.file(package=pkgname)
## R> aspell_package_Rd_files(pkg.pathname)
## R> aspell_package_vignettes(pkg.pathname)
## R> aspell_package_R_files(pkg.pathname)
## R> aspell_package_C_files(pkg.pathname)
##
## To get just the list of questioned words on manpages, try this instead:
##
## R> asa <- aspell_package_Rd_files(pkg.pathname)
## R> unique(sort(asa$Original))
##


##
## Custom dictionaries
##
## Create new personal dictionary as follows:
##
## R> lang <- "en"                     # Use language code or locale name
## R> words <- c("foo", "bar", "baz")
## R> dict.basename <- sprintf("%s-%s", lang, deparse(substitute(words)))
## R> dict.filename <- sprintf("%s.rds", dict.basename)
## R> dict.pathname <- system.file(".aspell", dict.filename, package=pkgname)
## R> saveRDS(words, dict.pathname)
##
## Then add its basename to the "aspell_dicts_pkg" variable below.
##
## Add new words to existing dictionary as follows:
##
## R> words <- readRDS(dict.pathname)
## R> words <- c(words, "qux", "quux") # Add new more words..
## R> words <- sort(words)
## R> saveRDS(words, dict.pathname)
##
aspell_dicts_R   <- "en_stats"
aspell_dicts_pkg <- c("en-antibodies",
                      "en-developer",
                      "en-microarray",
                      "en-software",
                      "en-statistics")

dictionaries <- c(aspell_dicts_pkg,
                  aspell_dicts_R)

##
## Manpage sections we don't bother to spellcheck
##
drop <- c("\\author",
          "\\references",
          "\\source")

##
## Settings used by "aspell-utils" functions
##
Rd_files <- list(dictionaries=dictionaries,
                 drop=drop)
vignettes <- list(dictionaries=dictionaries)
R_files <- list(dictionaries=dictionaries)
C_files <- list(dictionaries=dictionaries)

