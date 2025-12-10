pkgname <- "UtilityFunctions"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('UtilityFunctions')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("calculate_h2_from_list")
### * calculate_h2_from_list

flush(stderr()); flush(stdout())

### Name: calculate_h2_from_list
### Title: Calculate Heritability from a List of ASReml Models
### Aliases: calculate_h2_from_list

### ** Examples

## Not run: 
##D # Assuming 'fm' is a list like fm <- list(fm0=model0, fm1=model1, ...)
##D heritabilities <- calculate_h2_from_list(fm, id = "Variety")
##D print(heritabilities)
## End(Not run)



cleanEx()
nameEx("plot.fa_asreml")
### * plot.fa_asreml

flush(stderr()); flush(stdout())

### Name: plot.fa_asreml
### Title: Master Visualization Suite for Factor Analytic Models
### Aliases: plot.fa_asreml

### ** Examples

## Not run: 
##D # Standard Selection View
##D plot(results, type = "fast", n_label = 5, highlight = c("CheckA", "CheckB"))
##D 
##D # Investigate GxE Drivers (Factor 2)
##D plot(results, type = "latent_reg", factor = 2)
##D 
##D # View Network Structure
##D plot(results, type = "biplot", factor = c(1, 2))
##D plot(results, type = "heatmap")
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
