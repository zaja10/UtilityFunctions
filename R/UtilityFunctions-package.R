#' UtilityFunctions: Factor Analytic and Design Tools
#'
#' A suite of tools for processing and visualizing Factor Analytic models
#' in plant breeding (ASReml-R).
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom dplyr %>% filter mutate select summarize
#' @importFrom rlang .data
#' @importFrom stats as.formula coef cov2cor density dnorm lm predict qnorm residuals sd var rbinom rnorm runif resid model.matrix setNames xtabs update
#' @importFrom utils head
#' @importFrom methods as
#' @importFrom graphics axis box grid image layout legend par plot points
#' @importFrom grDevices hcl.colors
NULL

# Register Global Variables to silence R CMD check
utils::globalVariables(c(
    ".", "Yield", "Genotype", "Trial", "Year", "Row", "Column",
    "MinT", "MaxT", "Precip", "ET0", # Weather
    "Var1", "Var2", "PlotValue", "Date", # Plotting/Envirotyping
    "RMSD", "OP", "Type", "Correlation", "Status", ".data"
))
