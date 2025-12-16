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
    current_res <- predict_cross_utility(current_set, markers = markers, effects = marker_effects, map = map)
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
        add_res <- predict_cross_utility(add_candidates, markers = markers, effects = marker_effects, map = map)

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
