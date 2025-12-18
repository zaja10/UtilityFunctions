#' Calculate Multi-Trait Stability Index (Enhanced)
#'
#' Computes a composite selection index across multiple traits, incorporating
#' both performance (OP) and stability (RMSD) for each trait.
#' Supports Independent Culling Levels (ICL).
#'
#' @param fa_objects_list A named list of objects of class \code{fa_model} (output from fit_fa_model).
#' @param weights Named numeric vector of economic weights (w) for each trait.
#' @param penalties Named numeric vector of stability penalties (lambda) for each trait.
#' @param thresholds Named numeric vector. Minimum acceptable raw value for a trait.
#'        E.g., c(Yield = 2.0, Protein = 10.5). Genotypes below this are culled.
#' @param cull_strict Logical. If TRUE, remove culled genotypes from the output.
#'        If FALSE (default), keep them but set their Index to -Inf and Status to "Culled".
#'
#' @return A dataframe with standardized scores per trait, raw values, and the final Index.
#' @export
calculate_mt_index <- function(fa_objects_list, weights, penalties,
                               thresholds = NULL, cull_strict = FALSE) {
    if (is.null(names(fa_objects_list))) stop("fa_objects_list must be named (Trait names).")

    traits <- names(fa_objects_list)
    common_genos <- NULL

    # 1. Extract and Standardize per trait
    res_list <- list()

    for (tr in traits) {
        obj <- fa_objects_list[[tr]]

        # Robustness: Check if this is a valid fa_model object
        if (is.null(obj$fast)) {
            warning(paste("Trait", tr, "does not have a valid $fast element. Skipping."))
            next
        }

        df <- obj$fast

        if (is.null(common_genos)) {
            common_genos <- df$Genotype
        } else {
            common_genos <- intersect(common_genos, df$Genotype)
        }

        # Store Raw OP for Gain Calculation (Critical for Step 4.2)
        df$OP_Raw <- df$OP

        # Standardize (Z-score)
        # Handle cases where SD is 0 to avoid NaNs
        sd_op <- sd(df$OP, na.rm = TRUE)
        sd_rmsd <- sd(df$RMSD, na.rm = TRUE)

        df$OP_Z <- if (!is.na(sd_op) && sd_op > 0) (df$OP - mean(df$OP, na.rm = TRUE)) / sd_op else 0
        df$RMSD_Z <- if (!is.na(sd_rmsd) && sd_rmsd > 0) (df$RMSD - mean(df$RMSD, na.rm = TRUE)) / sd_rmsd else 0

        res_list[[tr]] <- df
    }

    if (length(common_genos) == 0) stop("No common genotypes found across traits.")

    # 2. Identify Culled Genotypes (ICL)
    culled_genos <- c()
    if (!is.null(thresholds)) {
        for (tr in names(thresholds)) {
            if (!tr %in% names(res_list)) next

            limit <- thresholds[[tr]]
            sub_df <- res_list[[tr]]

            # Identify bad genotypes based on Raw OP
            # Assumes 'higher is better'. If trait is negative (e.g., Disease),
            # ensure your threshold/logic matches.
            bad <- sub_df$Genotype[sub_df$OP_Raw < limit]
            culled_genos <- unique(c(culled_genos, bad))
        }
    }

    # 3. Combine
    final_df <- data.frame(Genotype = common_genos, stringsAsFactors = FALSE)
    final_df$Index <- 0

    for (tr in traits) {
        if (!tr %in% names(res_list)) next

        w <- if (tr %in% names(weights)) weights[[tr]] else 1.0
        lam <- if (tr %in% names(penalties)) penalties[[tr]] else 0.0

        sub <- res_list[[tr]]
        match_idx <- match(final_df$Genotype, sub$Genotype)

        op_z <- sub$OP_Z[match_idx]
        rmsd_z <- sub$RMSD_Z[match_idx]
        op_raw <- sub$OP_Raw[match_idx]

        # Component Index: w * (OP - lambda * RMSD)
        # RMSD is "bad" (instability), so subtracting it improves the index.
        comp_val <- w * (op_z - (lam * rmsd_z))

        final_df[[paste0(tr, "_OP_Z")]] <- round(op_z, 2)
        final_df[[paste0(tr, "_RMSD_Z")]] <- round(rmsd_z, 2)
        final_df[[paste0(tr, "_Comp")]] <- comp_val
        final_df[[paste0(tr, "_raw")]] <- op_raw # Preserve raw for prediction

        final_df$Index <- final_df$Index + comp_val
    }

    # 4. Handle Culling Status
    final_df$Status <- "Selected"
    if (length(culled_genos) > 0) {
        is_culled <- final_df$Genotype %in% culled_genos
        final_df$Status[is_culled] <- "Culled"

        if (cull_strict) {
            final_df <- final_df[!is_culled, ]
        } else {
            # Penalize Index to push to bottom
            final_df$Index[is_culled] <- -Inf
        }
    }

    # Sort by Index (Descending)
    final_df <- final_df[order(final_df$Index, decreasing = TRUE), ]
    final_df$Rank <- 1:nrow(final_df)

    return(final_df)
}

#' Predict Response to Selection
#'
#' Calculates the expected response to selection (Realized Selection Differential)
#' based on the provided index.
#'
#' @param index_df Dataframe returned by \code{calculate_mt_index}.
#' @param select_pct Percentage of the population to select (default 10).
#'
#' @return A dataframe showing Population Mean, Selected Mean, Differential, and % Change per trait.
#' @export
predict_response <- function(index_df, select_pct = 10) {
    # Filter valid population (exclude culled lines)
    pop <- index_df[index_df$Status != "Culled", ]

    if (nrow(pop) == 0) stop("No valid genotypes in population.")

    # Calculate number to select
    n_sel <- max(1, ceiling(nrow(pop) * (select_pct / 100)))
    sel <- head(pop, n_sel)

    # Identify trait columns (ending in _raw)
    trait_cols <- grep("_raw$", names(pop), value = TRUE)
    if (length(trait_cols) == 0) stop("No '_raw' trait columns found. Did you use the enhanced calculate_mt_index?")

    res <- data.frame(Trait = sub("_raw", "", trait_cols), stringsAsFactors = FALSE)

    # Calculate means
    res$Pop_Mean <- colMeans(pop[, trait_cols, drop = FALSE], na.rm = TRUE)
    res$Sel_Mean <- colMeans(sel[, trait_cols, drop = FALSE], na.rm = TRUE)
    res$Differential <- res$Sel_Mean - res$Pop_Mean

    # Pct Change (Handle division by zero)
    res$Pct_Change <- ifelse(res$Pop_Mean == 0,
        NA,
        (res$Differential / abs(res$Pop_Mean)) * 100
    )

    # Formatting
    res$Pop_Mean <- round(res$Pop_Mean, 3)
    res$Sel_Mean <- round(res$Sel_Mean, 3)
    res$Differential <- round(res$Differential, 3)
    res$Pct_Change <- round(res$Pct_Change, 2)

    return(res)
}

#' Master Selection Pipeline
#'
#' A unified wrapper that integrates Multi-Trait Index calculation, Independent Culling,
#' and Response Prediction into a single workflow.
#'
#' @param fa_model_list A named list of \code{fa_model} objects.
#'        (This is the output of \code{fit_fa_model}, NOT the raw asreml object).
#' @param weights Named numeric vector of economic weights.
#' @param penalties Named numeric vector of stability penalties.
#' @param thresholds Named numeric vector for independent culling levels (optional).
#' @param check_geno Character string. Name of the Check variety to use as a benchmark (optional).
#' @param top_n Integer. Number of genotypes to select for prediction metrics (default 10).
#'
#' @return An object of class \code{mt_selection_results}.
#' @export
run_selection_pipeline <- function(fa_model_list, weights, penalties,
                                   thresholds = NULL, check_geno = NULL, top_n = 10) {
    if (is.null(names(fa_model_list))) stop("fa_model_list must be named (e.g., list(Yield=obj1, ...)).")

    message("--- Starting Multi-Trait Selection Pipeline ---")

    # 1. Calculate Index & Apply Culling
    message("1. Calculating Index and applying thresholds...")
    index_df <- calculate_mt_index(fa_model_list, weights, penalties,
        thresholds = thresholds, cull_strict = FALSE
    )

    # 2. Predict Response (Genetic Gain)
    message("2. Predicting Response to Selection...")
    total_valid <- sum(index_df$Status != "Culled")
    sel_pct <- if (total_valid > 0) (top_n / total_valid) * 100 else 10

    gain_df <- tryCatch(
        {
            predict_response(index_df, select_pct = sel_pct)
        },
        error = function(e) {
            warning("Could not calculate response: ", e$message)
            return(NULL)
        }
    )

    # 3. Construct Result Object
    res <- list(
        index = index_df,
        prediction = gain_df,
        params = list(
            weights = weights,
            penalties = penalties,
            thresholds = thresholds,
            check_geno = check_geno,
            top_n = top_n
        )
    )

    class(res) <- "mt_selection_results"

    message("--- Pipeline Complete. Use plot() to visualize. ---")
    return(res)
}

#' @export
print.mt_selection_results <- function(x, ...) {
    cat("Multi-Trait Selection Results\n")
    cat("=============================\n")
    cat("Traits: ", paste(names(x$params$weights), collapse = ", "), "\n")

    n_culled <- sum(x$index$Status == "Culled")
    cat("Culled Genotypes: ", n_culled, "\n")

    cat("\nTop", x$params$top_n, "Selections:\n")
    valid_top <- head(x$index[x$index$Status != "Culled", ], x$params$top_n)
    print(valid_top[, c("Genotype", "Rank", "Index")])

    if (!is.null(x$prediction)) {
        cat("\nPredicted Gain (Top ", x$params$top_n, "):\n", sep = "")
        print(x$prediction[, c("Trait", "Differential", "Pct_Change")])
    }
}
