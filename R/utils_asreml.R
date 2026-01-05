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

#' Force Convergence
#'
#' Repeatedly updates an ASReml model until convergence or max iterations.
#' @param model asreml object.
#' @param max_iter Numeric limit.
#' @export
force_convergence <- function(model, max_iter = 10) {
    if (!inherits(model, "asreml")) stop("Not an asreml object.")

    for (i in 1:max_iter) {
        if (model$converge) {
            message("Converged.")
            return(model)
        }
        message(paste("Update", i))
        model <- update(model)
    }
    warning("Did not converge.")
    model
}
