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

#' Bind Parallel FA Results
#'
#' Aggregates a list of results from \code{run_parallel} into a single dataframe.
#' Useful when processing 50 sites and needing a master summary table.
#'
#' @param result_list A list of objects (usually containing extraction results).
#' @param element Character. The specific element to extract (e.g., "fast" or "vaf").
#' @return A combined dataframe with a "Source" column.
#' @export
bind_fa_results <- function(result_list, element = "fast") {
    # Filter NULLs or Errors
    valid_results <- result_list[!sapply(result_list, is.null)]

    if (length(valid_results) == 0) {
        warning("No valid results to bind.")
        return(NULL)
    }

    combined <- do.call(rbind, lapply(names(valid_results), function(nm) {
        obj <- valid_results[[nm]]
        if (!element %in% names(obj)) {
            return(NULL)
        }

        df <- obj[[element]]
        # Add ID column if list is named
        if (!is.null(nm) && nm != "") df$Source_ID <- nm
        return(df)
    }))

    return(combined)
}
