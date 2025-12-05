#' Population Genetics & Cross Prediction
#'
#' @description
#' Tools for predicting the value of crosses based on genomic data.

#' Predict Cross Utility (Usefulness Criterion)
#'
#' @description
#' Calculates the Usefulness Criterion (UC) for potential crosses between a set of parents.
#' \deqn{UC = \mu + (i \times \sigma_{prog})}
#' where \eqn{\mu} is the mid-parent value and \eqn{\sigma_{prog}} is the predicted gametic standard deviation.
#'
#' @param parents Character vector of parent IDs to cross (all pairwise combinations will be generated).
#'        Alternatively, a dataframe with columns "Parent1" and "Parent2".
#' @param markers Matrix. Genotype matrix (n x m) with rownames matching parent IDs. Coding should be numeric (0, 1, 2) or (0, 2) for inbreds.
#' @param effects Named numeric vector of marker effects (BLUPs from a Training Set).
#' @param map Named numeric vector of marker positions in Morgans (NOT cM).
#' @param i_sel Numeric. Selection intensity (standardized). Default is 2.06 (Top 5%).
#'
#' @return A dataframe containing:
#' \describe{
#'   \item{Parent1, Parent2}{The cross combination.}
#'   \item{Mu}{Mid-parent breeding value.}
#'   \item{Sigma}{Predicted standard deviation of the progeny.}
#'   \item{UC}{The Usefulness Criterion.}
#' }
#'
#' @details
#' This function uses a C++ backend (\code{calc_cross_variance_cpp}) to efficiently calculate the
#' expected variance of F2 progeny accounting for linkage disequilibrium between markers.
#'
#' @export
predict_cross_utility <- function(parents, markers, effects, map, i_sel = 2.06) {
    # 1. SETUP & VALIDATION
    if (!is.matrix(markers)) stop("markers must be a matrix.")

    # Ensure map and effects match markers
    common_mk <- intersect(colnames(markers), names(effects))
    common_mk <- intersect(common_mk, names(map))

    if (length(common_mk) == 0) stop("No common markers between genotype matrix, effects, and map.")
    if (length(common_mk) < length(effects)) warning("Some markers with effects are missing in genotype/map. They will be ignored.")

    # Subset and Align
    # It is CRITICAL that effects and map are sorted in the same order as markers columns
    # and map order (physical/genetic order) is respected for LD calc.

    # Sort map first (Genome order)
    map <- map[common_mk]
    map <- map[order(map)] # Assuming input map is comparable (e.g. cumulative Morgans or by Chr)
    # TODO: Handle multi-chromosomes. If map is just position, we need chr info to break LD.
    # For V1, we assume map is cumulative or user provides single chromosome?
    # Or we update C++ to handle chromosomes.
    # Simple fix: If distance > 0.5 Morgans (50cM), r = 0.5 (Unlinked).
    # The Haldane function handles this naturally: (1 - exp(-2d))/2. If d is large, r -> 0.5.
    # So Cumulative map works fine across chromosomes if gaps are large.

    ordered_mk <- names(map)

    M <- markers[, ordered_mk, drop = FALSE]
    eff <- effects[ordered_mk]
    map_vec <- as.numeric(map)

    # 2. GENERATE CROSS LIST
    if (is.data.frame(parents)) {
        cross_df <- parents
        if (!all(c("Parent1", "Parent2") %in% names(cross_df))) stop("parents dataframe must have Parent1 and Parent2 columns.")
    } else {
        # All pairwise
        cross_df <- expand.grid(Parent1 = parents, Parent2 = parents, stringsAsFactors = FALSE)
        # Remove selfs? Keeps selfs for checking.
        # Remove duplicates (A x B vs B x A)?
        cross_df <- cross_df[cross_df$Parent1 <= cross_df$Parent2, ]
    }

    # Check parents exist
    all_p <- unique(c(cross_df$Parent1, cross_df$Parent2))
    missing_p <- setdiff(all_p, rownames(M))
    if (length(missing_p) > 0) stop(paste("Parents not found in marker matrix:", paste(missing_p, collapse = ", ")))

    # 3. CALCULATE UTILITY
    n_cross <- nrow(cross_df)
    results <- cross_df
    results$Mu <- NA
    results$Sigma <- NA
    results$UC <- NA

    # Calculate GEBVs for parents (for Mu)
    # GEBV = M * eff
    # Using subsets
    mrk_eff_mat <- as.matrix(eff)
    gebv <- (M[all_p, ] %*% mrk_eff_mat)[, 1]

    cat(sprintf("Evaluating %d crosses...\n", n_cross))

    for (k in 1:n_cross) {
        p1_id <- cross_df$Parent1[k]
        p2_id <- cross_df$Parent2[k]

        # Mu
        results$Mu[k] <- (gebv[p1_id] + gebv[p2_id]) / 2

        # Sigma
        if (p1_id == p2_id) {
            results$Sigma[k] <- 0
        } else {
            p1_vec <- M[p1_id, ]
            p2_vec <- M[p2_id, ]

            # Call C++
            # We need to assume the package is loaded/compiled.
            # If pure source, this function might fail if not compiled.
            # We'll use tryCatch or assume.
            var_prog <- calc_cross_variance_cpp(p1_vec, p2_vec, eff, map_vec)
            results$Sigma[k] <- sqrt(var_prog)
        }
    }

    results$UC <- results$Mu + (i_sel * results$Sigma)
    results <- results[order(results$UC, decreasing = TRUE), ]

    return(results)
}
