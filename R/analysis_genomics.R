#' Predict Cross Utility (Single or Batch)
#'
#' Calculates the expected genetic variance (UC) for one or multiple crosses.
#'
#' @param parents Either:
#'    1. A numeric matrix (2 x M) for a single cross (Standard).
#'    2. A dataframe with columns "Parent1" and "Parent2" (Batch Mode).
#' @param ... Additional arguments.
#'    - If Batch Mode: must provide \code{markers} (Genotype Matrix) and \code{effects}, \code{map}.
#'    - If Single Mode: provide \code{effects}, \code{map}, \code{max_dist}.
#'
#' @return A numeric value (Single) or a Dataframe with UC values (Batch).
#' @export
predict_cross_utility <- function(parents, ...) {
    args <- list(...)

    # --- BATCH MODE (Dataframe) ---
    if (is.data.frame(parents)) {
        if (!all(c("Parent1", "Parent2") %in% names(parents))) {
            stop("Batch input must have 'Parent1' and 'Parent2' columns.")
        }
        # Unpack args for batch
        if (!"markers" %in% names(args)) stop("Batch mode requires 'markers' matrix.")
        markers <- args[["markers"]]
        effects <- args[["effects"]] # Renamed from marker_effects in call? Check consistency.
        map <- args[["map"]]

        # Handle naming inconsistency in optimize_mating_list call vs definition
        # optimize uses: predict_cross_utility(set, markers, marker_effects, map)
        # So args[[2]] is marker_effects
        if (is.null(effects)) effects <- args[[2]]
        if (is.null(map)) map <- args[[3]]

        # Pre-check existence
        all_parents <- unique(c(parents$Parent1, parents$Parent2))
        missing_p <- setdiff(all_parents, rownames(markers))
        if (length(missing_p) > 0) stop(paste("Parents missing from marker matrix:", head(missing_p)))

        # Loop prediction (Future optimization: Move loop to C++)
        # We use a simple apply here for clarity
        results <- parents
        results$UC <- NA_real_

        # Ensure sorted map for C++ safety
        if (is.unsorted(map)) {
            ord <- order(map)
            map <- map[ord]
            effects <- effects[ord]
            markers <- markers[, ord, drop = FALSE]
        }

        # Iterate
        # Note: This is the bottleneck. Ideally, we write 'calc_batch_variance_cpp'.
        uc_vals <- numeric(nrow(parents))

        for (i in 1:nrow(parents)) {
            p1_id <- parents$Parent1[i]
            p2_id <- parents$Parent2[i]

            # Extract rows
            p1_vec <- markers[p1_id, ]
            p2_vec <- markers[p2_id, ]

            # Combine into 2xM matrix expected by C++ wrapper logic
            # But wait, the C++ function takes vectors.
            # We can call C++ directly here for speed.
            uc_vals[i] <- calc_cross_variance_cpp(p1_vec, p2_vec, effects, map, 3.0)
        }

        results$UC <- uc_vals
        return(results)
    }

    # --- SINGLE MODE (Matrix) ---
    else {
        # Original logic (Legacy support)
        # Expects: parents (2xM), effects, map, max_dist
        effects <- args[[1]]
        map <- args[[2]]
        max_dist <- if (length(args) > 2) args[[3]] else 3.0

        # Validation
        if (nrow(parents) != 2) stop("Parents matrix must have 2 rows.")

        # Sort Safety
        if (is.unsorted(map)) {
            ord <- order(map)
            map <- map[ord]
            effects <- effects[ord]
            parents <- parents[, ord, drop = FALSE]
        }

        return(calc_cross_variance_cpp(parents[1, ], parents[2, ], effects, map, max_dist))
    }
}
