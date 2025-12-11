#' Predict Cross Utility (Usefulness Criterion)
#'
#' Calculates the Usefulness Criterion (UC) for potential crosses between a set of parents.
#'
#' @param parents Character vector of parent IDs to cross (all pairwise combinations will be generated).
#'        Alternatively, a dataframe with columns "Parent1" and "Parent2".
#' @param markers Matrix. Genotype matrix (n x m) with rownames matching parent IDs. Coding should be numeric (0, 1, 2) or (0, 2) for inbreds.
#' @param effects Named numeric vector of marker effects (BLUPs from a Training Set).
#' @param map Named numeric vector of marker positions in Morgans (NOT cM).
#' @param i_sel Numeric. Selection intensity (standardized). Default is 2.06 (Top 5%).
#'
#' @return A dataframe containing the columns Parent1, Parent2, Mu, Sigma, and UC (Usefulness Criterion).
#' @export
predict_cross_utility <- function(parents, markers, effects, map, i_sel = 2.06) {
    if (!is.matrix(markers)) stop("markers must be a matrix.")

    common_mk <- intersect(colnames(markers), names(effects))
    common_mk <- intersect(common_mk, names(map))

    if (length(common_mk) == 0) stop("No common markers between genotype matrix, effects, and map.")
    if (length(common_mk) < length(effects)) warning("Some markers with effects are missing in genotype/map. They will be ignored.")

    map <- map[common_mk]
    map <- map[order(map)]
    ordered_mk <- names(map)

    M <- markers[, ordered_mk, drop = FALSE]
    eff <- effects[ordered_mk]
    map_vec <- as.numeric(map)

    if (is.data.frame(parents)) {
        cross_df <- parents
        if (!all(c("Parent1", "Parent2") %in% names(cross_df))) stop("parents dataframe must have Parent1 and Parent2 columns.")
    } else {
        cross_df <- expand.grid(Parent1 = parents, Parent2 = parents, stringsAsFactors = FALSE)
        cross_df <- cross_df[cross_df$Parent1 <= cross_df$Parent2, ]
    }

    all_p <- unique(c(cross_df$Parent1, cross_df$Parent2))
    missing_p <- setdiff(all_p, rownames(M))
    if (length(missing_p) > 0) stop(paste("Parents not found in marker matrix:", paste(missing_p, collapse = ", ")))

    n_cross <- nrow(cross_df)
    results <- cross_df
    results$Mu <- NA
    results$Sigma <- NA
    results$UC <- NA

    mrk_eff_mat <- as.matrix(eff)
    gebv <- (M[all_p, ] %*% mrk_eff_mat)[, 1]

    cat(sprintf("Evaluating %d crosses...\n", n_cross))

    for (k in 1:n_cross) {
        p1_id <- cross_df$Parent1[k]
        p2_id <- cross_df$Parent2[k]

        results$Mu[k] <- (gebv[p1_id] + gebv[p2_id]) / 2

        if (p1_id == p2_id) {
            results$Sigma[k] <- 0
        } else {
            p1_vec <- M[p1_id, ]
            p2_vec <- M[p2_id, ]

            # Using C++ function
            var_prog <- calc_cross_variance_cpp(p1_vec, p2_vec, eff, map_vec)
            results$Sigma[k] <- sqrt(var_prog)
        }
    }

    results$UC <- results$Mu + (i_sel * results$Sigma)
    results <- results[order(results$UC, decreasing = TRUE), ]

    return(results)
}
