#' ASReml Utility Wrapper
#'
#' @name utils_asreml
NULL

#' Check Availability
#' @export
check_asreml_availability <- function() {
    if (!requireNamespace("asreml", quietly = TRUE)) {
        stop("ASReml-R is required for this function.")
    }
    TRUE
}

#' Force Convergence (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Use \code{\link{mkConv}} instead. This function now wraps \code{mkConv} for backward compatibility.
#'
#' @param model asreml object.
#' @param max_iter Numeric limit.
#' @export
force_convergence <- function(model, max_iter = 10) {
    warning("force_convergence() is deprecated. Please use mkConv() instead.")
    mkConv(model)
}

#' Force Convergence by Variance Change
#'
#' Iteratively updates the model until the percentage change in variance components
#' is below a threshold (default 1%).
#'
#' @param mod asreml object.
#' @return Converged asreml object.
#' @export
mkConv <- function(mod) {
    if (!inherits(mod, "asreml")) stop("Not an asreml object")

    # Helper to safely extract percent change
    get_pct_chg <- function(m) {
        summ <- summary(m)
        if (!"varcomp" %in% names(summ)) {
            return(0)
        }
        vc <- summ$varcomp
        if (!"%ch" %in% names(vc)) {
            return(0)
        }
        return(vc[, "%ch"])
    }

    pctchg <- get_pct_chg(mod)
    tries <- 1

    # Loop while any change > 1% and max tries < 20
    while (any(abs(pctchg) > 1, na.rm = TRUE) && tries < 20) {
        tries <- tries + 1
        mod <- suppressWarnings(update(mod))
        pctchg <- get_pct_chg(mod)
    }

    if (tries >= 20) warning("mkConv reached max iterations (20) without full stabilization.")
    return(mod)
}
