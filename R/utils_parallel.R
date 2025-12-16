#' Run Function in Parallel (OS Agnostic)
#'
#' @param X Vector or List to iterate over.
#' @param FUN Function to apply.
#' @param cores Number of cores.
#' @param ... Args for FUN.
#' @importFrom parallel mclapply makeCluster stopCluster parLapply detectCores
#' @export
run_parallel <- function(X, FUN, cores = NULL, ...) {
    if (is.null(cores)) cores <- max(1, parallel::detectCores(logical = FALSE) - 1)

    if (cores == 1) {
        return(lapply(X, FUN, ...))
    }

    if (.Platform$OS.type == "unix") {
        return(parallel::mclapply(X, FUN, ..., mc.cores = cores))
    } else {
        cl <- parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl))
        # Note: For Windows, you might need to export libraries inside FUN
        # or use clusterEvalQ if FUN relies on loaded packages.
        return(parallel::parLapply(cl, X, FUN, ...))
    }
}
