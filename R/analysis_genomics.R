#' Predict Cross Utility (Genetic Variance)
#'
#' Calculates the expected genetic variance of a biparental cross given
#' parental genotypes and marker effects. Enforces map sorting for C++ engine safety.
#'
#' @param parents Matrix (2 x M) of parental genotypes (scored 0, 1, 2).
#' @param effects Vector (M) of marker effects.
#' @param map Vector (M) of genetic map positions (in Morgans).
#' @param max_dist Numeric. Maximum distance (Morgans) to consider linkage.
#'
#' @return A numeric value representing the expected variance.
#' @export
predict_cross_utility <- function(parents, effects, map, max_dist = 3.0) {
    # 1. Validation
    if (ncol(parents) != length(effects) || length(effects) != length(map)) {
        stop("Dimensions of parents, effects, and map must match.")
    }

    if (nrow(parents) != 2) {
        stop("Parents matrix must have exactly 2 rows.")
    }

    # 2. CRITICAL: Enforce Sort Order
    # The C++ engine uses a sliding window optimization that assumes
    # markers are sorted by position. If not sorted, the loop breaks prematurely.
    if (is.unsorted(map)) {
        ord <- order(map)
        map_sorted <- map[ord]
        eff_sorted <- effects[ord]
        # Reorder columns of parent matrix
        par_sorted <- parents[, ord, drop = FALSE]
    } else {
        map_sorted <- map
        eff_sorted <- effects
        par_sorted <- parents
    }

    # 3. Call C++ Engine
    p1 <- par_sorted[1, ]
    p2 <- par_sorted[2, ]

    var_g <- calc_cross_variance_cpp(p1, p2, eff_sorted, map_sorted, max_dist)

    return(var_g)
}
