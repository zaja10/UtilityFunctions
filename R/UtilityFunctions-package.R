#' UtilityFunctions: Factor Analytic and Design Tools
#'
#' A suite of tools for processing and visualizing Factor Analytic models
#' in plant breeding (ASReml-R).
#'
#' ASReml-R is proprietary and therefore only suggested. Functions that fit or
#' interrogate \code{asreml} objects check for it at run time and fail with an
#' informative error when it is unavailable.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom rlang .data .env
#' @importFrom stats as.formula coef cor cov2cor dnorm lm na.omit predict qnorm
#' @importFrom stats quantile rbinom resid rnorm runif sd setNames terms update var
#' @importFrom utils head read.csv write.csv
#' @importFrom patchwork plot_annotation
NULL

# `%||%` is exported by base R only from 4.4.0; keep a private copy.
`%||%` <- function(x, y) if (is.null(x)) y else x
