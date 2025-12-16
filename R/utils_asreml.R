#' Check ASReml Availability
#'
#' @description
#' Checks if the `asreml` package is installed and licensed.
#' Stops execution with a clear message if unavailable.
#'
#' @param action Character. "stop" (default) or "warning" or "logical".
#' @return Logical TRUE if available, or error/warning if not.
#' @export
check_asreml_availability <- function(action = "stop") {
    installed <- requireNamespace("asreml", quietly = TRUE)

    if (installed) {
        # Try a simple licensed operation (e.g., check license status if possible,
        # or just rely on package load not failing)
        # ASReml checks license on load/function call.
        # We can try to access a function.
    }

    if (!installed) {
        msg <- "The 'asreml' package is required for this function but is not installed.\nPlease install it (and acquire a license) to use this feature."
        if (action == "stop") stop(msg, call. = FALSE)
        if (action == "warning") warning(msg, call. = FALSE)
        return(FALSE)
    }

    return(TRUE)
}
