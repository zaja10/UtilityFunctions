#' Calculate Multi-Trait Stability Index
#'
#' Computes a composite selection index across multiple traits, incorporating
#' both performance (OP) and stability (RMSD) for each trait.
#'
#' @param fa_objects_list A named list of \code{fa_asreml} objects (one per trait).
#' @param weights Named numeric vector of economic weights (w) for each trait.
#' @param penalties Named numeric vector of stability penalties (lambda) for each trait.
#'
#' @return A dataframe with standardized scores per trait and the final Index.
#' @export
calculate_mt_index <- function(fa_objects_list, weights, penalties) {
    if (is.null(names(fa_objects_list))) stop("fa_objects_list must be named (Trait names).")

    traits <- names(fa_objects_list)
    common_genos <- NULL

    # 1. Extract and Standardize per trait
    res_list <- list()

    for (tr in traits) {
        obj <- fa_objects_list[[tr]]
        if (is.null(obj$fast)) {
            # Try to calculate if missing
            obj <- tryCatch(calculate_op(obj), error = function(e) NULL)
        }

        if (is.null(obj) || is.null(obj$selection)) {
            warning(paste("Skipping trait", tr, "- no FAST indices available."))
            next
        }

        df <- obj$selection
        if (is.null(common_genos)) {
            common_genos <- df$Genotype
        } else {
            common_genos <- intersect(common_genos, df$Genotype)
        }

        # Standardize (Z-score)
        df$OP_Z <- (df$OP - mean(df$OP, na.rm = TRUE)) / sd(df$OP, na.rm = TRUE)
        df$RMSD_Z <- (df$RMSD - mean(df$RMSD, na.rm = TRUE)) / sd(df$RMSD, na.rm = TRUE)

        res_list[[tr]] <- df
    }

    if (length(common_genos) == 0) stop("No common genotypes found across traits.")

    # 2. Combine
    final_df <- data.frame(Genotype = common_genos)
    final_df$Index <- 0

    for (tr in traits) {
        w <- if (tr %in% names(weights)) weights[[tr]] else 1.0
        lam <- if (tr %in% names(penalties)) penalties[[tr]] else 0.0

        sub <- res_list[[tr]]
        match_idx <- match(final_df$Genotype, sub$Genotype)

        op_z <- sub$OP_Z[match_idx]
        rmsd_z <- sub$RMSD_Z[match_idx]

        # Component Index: w * (OP - lambda * RMSD)
        # Note: RMSD is "bad", so -RMSD is "good" stability
        # Formula: w * (OP_Z - lambda * RMSD_Z)

        comp_val <- w * (op_z - (lam * rmsd_z))

        final_df[[paste0(tr, "_OP_Z")]] <- round(op_z, 2)
        final_df[[paste0(tr, "_RMSD_Z")]] <- round(rmsd_z, 2)
        final_df[[paste0(tr, "_Comp")]] <- comp_val

        final_df$Index <- final_df$Index + comp_val
    }

    final_df <- final_df[order(final_df$Index, decreasing = TRUE), ]
    final_df$Rank <- 1:nrow(final_df)

    return(final_df)
}
