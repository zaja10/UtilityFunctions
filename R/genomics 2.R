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
#' Optimize Mating List (Optimal Haplotype Stacking)
#'
#' Uses a Genetic Algorithm (Villiers et al., 2024) to select a set of crosses that maximizes
#' the total Usefulness Criterion (UC) of the progeny population, subject to budget constraints.
#'
#' @param parents Character vector of all available parent IDs.
#' @param marker_effects Named numeric vector of marker effects (BLUPs).
#' @param map Named numeric vector of marker positions (Morgans).
#' @param markers Matrix. Genotype matrix (n_parents x n_markers).
#' @param n_crosses Integer. Target number of crosses to select.
#' @param max_usage Integer. Maximum number of times a single parent can be used (constraint).
#' @param iterations Integer. Number of GA generations.
#' @param mutation_rate Numeric. Probability of mutating a cross in the set.
#'
#' @return A dataframe of selected crosses with their UC values.
#' @export
optimize_mating_list <- function(parents, marker_effects, map, markers, n_crosses = 100, max_usage = 10, iterations = 1000, mutation_rate = 0.1) {
    # Initial Population: Random Selection
    cross_pool <- expand.grid(Parent1 = parents, Parent2 = parents, stringsAsFactors = FALSE)
    cross_pool <- cross_pool[cross_pool$Parent1 < cross_pool$Parent2, ] # Unique pairs

    # Pre-calculate GEBVs for mean calc (faster than full UC for everything)
    # Actually, we need full UC for the fitness function.
    # Evaluating all 500k crosses is too slow.
    # Better strategy: Evaluate UC on the fly for the *selected* set.

    # Start with random set
    current_set <- cross_pool[sample(nrow(cross_pool), n_crosses), ]

    # Calculate initial fitness
    # Batch calculate UC for current set
    current_res <- predict_cross_utility(current_set, markers, marker_effects, map)
    current_fitness <- sum(current_res$UC)

    message(sprintf("Starting GA... Initial Fitness: %.2f", current_fitness))

    best_set <- current_res
    best_fitness <- current_fitness

    for (i in 1:iterations) {
        # Mutation: Swap X crosses
        n_mut <- max(1, round(n_crosses * mutation_rate))

        # Select candidates to remove (e.g. lowest UC in current set, or random)
        # Random for exploration
        remove_idx <- sample(n_crosses, n_mut)

        # Select candidates to add from pool (random)
        # Avoid duplicates in set
        new_candidates_pool <- cross_pool[!paste(cross_pool$Parent1, cross_pool$Parent2) %in%
            paste(best_set$Parent1, best_set$Parent2), ]

        if (nrow(new_candidates_pool) < n_mut) break # Pool exhausted

        add_candidates <- new_candidates_pool[sample(nrow(new_candidates_pool), n_mut), ]

        # Calculate UC for new candidates
        add_res <- predict_cross_utility(add_candidates, markers, marker_effects, map)

        # Proposed new set
        prop_set <- rbind(best_set[-remove_idx, ], add_res)

        # Check Constraints (Max Usage)
        usage <- table(c(prop_set$Parent1, prop_set$Parent2))
        if (any(usage > max_usage)) {
            # Invalid solution, skip
            next
        }

        # Calculate new fitness
        prop_fitness <- sum(prop_set$UC)

        # Selection
        if (prop_fitness > best_fitness) {
            best_set <- prop_set
            best_fitness <- prop_fitness
            if (i %% 100 == 0) message(sprintf("Iter %d: New Best Fitness: %.2f", i, best_fitness))
        }
    }

    best_set$Rank <- rank(-best_set$UC)
    best_set <- best_set[order(best_set$Rank), ]

    return(best_set)
}

