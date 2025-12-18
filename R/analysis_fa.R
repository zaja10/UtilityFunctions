#' Extract and Rotate Factor Analytic Model Parameters
#'
#' A universal parser for ASReml-R Factor Analytic (FA) and Reduced Rank (RR) models.
#' Implements the "Smith & Cullis" Factor Analytic Selection Tools (FAST) framework.
#'
#' @param model A fitted object of class \code{asreml}.
#' @param classify Character (Optional). The term to extract. E.g., \code{"fa(Site, 2)"}.
#'        If NULL, the function attempts to find the first FA/RR term in \code{model$call$random}.
#' @param psi_term Character (Optional). Specific variance term for RR models.
#' @param rotate Character or Logical. Rotation strategy. Options:
#'        \itemize{
#'          \item \code{"varimax"} or \code{TRUE}: Standard PCA/SVD rotation (Max Var).
#'          \item \code{"mean"}: Rotate Factor 1 to align with the average loading magnitude (Overall Performance).
#'          \item \code{FALSE}: Raw solution.
#'        }
#' @param annotation Optional dataframe or named list to annotate sites immediately.
#'        Passed to \code{\link{annotate_model}}.
#'
#' @return An object of class \code{fa_model} with 'loadings', 'scores', 'fast', and 'variance'.
#' @importFrom stats coef cov2cor sd
#' @import cli
#' @export
fit_fa_model <- function(model, classify = NULL, psi_term = NULL, rotate = "mean", annotation = NULL) {
    cli::cli_h1("Extracting FA/RR Model Parameters")

    if (is.logical(rotate)) {
        rotate <- if (rotate) "varimax" else "none"
    }

    if (!inherits(model, "asreml")) cli::cli_abort("Object must be of class {.cls asreml}.")

    # 1. Identify Model Structure -----------------------------------------------
    call_random <- paste(deparse(model$call$random), collapse = "")

    if (is.null(classify)) {
        # Auto-detect FA term
        match_fa <- regmatches(call_random, regexec("(fa|rr)\\(([^,]+),\\s*([0-9]+)\\)", call_random))
        if (length(match_fa[[1]]) < 4) {
            cli::cli_abort("Could not auto-detect simple 'fa(Site, k)' term in random formula. Please provide {.arg classify} argument.")
        }
        type <- match_fa[[1]][2]
        site_var <- match_fa[[1]][3]
        k <- as.numeric(match_fa[[1]][4])
        classify <- match_fa[[1]][1]
        cli::cli_alert_info("Auto-detected term: {.val {classify}} (k={k}, Site={site_var})")
    } else {
        clean_str <- gsub("\\s+", "", classify)
        k_match <- regmatches(clean_str, regexec("(fa|rr)\\([^,]+,\\s*([0-9]+)\\)", clean_str))
        if (length(k_match[[1]]) >= 3) {
            k <- as.numeric(k_match[[1]][3])
            type <- k_match[[1]][2]
            site_match <- regmatches(clean_str, regexec("(fa|rr)\\(([^,]+),", clean_str))
            site_var <- site_match[[1]][3]
        } else {
            cli::cli_abort("Could not parse 'k' from provided classify string: {.val {classify}}")
        }
    }

    # 2. Extract Terms ----------------------------------------------------------
    vp <- summary(model)$varcomp
    vp_names <- rownames(vp)

    # Loadings
    lambda_list <- lapply(1:k, function(i) {
        suffix <- paste0("!((fa|rr|comp)[_]?", i, ")$")
        rows <- grep(suffix, vp_names, ignore.case = TRUE)
        if (length(rows) == 0) {
            return(NULL)
        }

        vals <- vp[rows, "component"]
        full_names <- vp_names[rows]
        tokens <- strsplit(full_names, "!")
        sites <- sapply(tokens, function(x) x[length(x) - 1])
        data.frame(Group = sites, Value = vals, Factor = i)
    })

    lambda_list <- lambda_list[!sapply(lambda_list, is.null)]
    lambda_df <- do.call(rbind, lambda_list)

    if (is.null(lambda_df) || nrow(lambda_df) == 0) {
        cli::cli_abort("No Factor Loadings found matching detected k={k}. Check model convergence.")
    }

    lambda_mat <- xtabs(Value ~ Group + Factor, data = lambda_df)
    class(lambda_mat) <- "matrix"

    if (ncol(lambda_mat) < k) {
        pad <- matrix(0, nrow(lambda_mat), k - ncol(lambda_mat))
        lambda_mat <- cbind(lambda_mat, pad)
    }

    # Specific Variances
    target_sites <- rownames(lambda_mat)
    psi_vals <- numeric(length(target_sites))
    names(psi_vals) <- target_sites
    found_psi <- FALSE
    for (site in target_sites) {
        pat <- paste0("!", site, "!var$")
        match_idx <- grep(pat, vp_names, value = FALSE)
        if (length(match_idx) == 1) {
            psi_vals[site] <- vp[match_idx, "component"]
            found_psi <- TRUE
        } else {
            psi_vals[site] <- 0
        }
    }
    if (!found_psi && type == "fa") {
        cli::cli_warn("No specific variances (!var) found for FA model sites. Assuming 0 (Reduced Rank behavior).")
    }
    psi_df <- data.frame(Group = target_sites, Psi = psi_vals)

    # Scores
    coefs <- coef(model)$random
    coef_names <- rownames(coefs)
    scores_list <- vector("list", k)
    gen_col_name <- "Genotype"

    for (i in 1:k) {
        pat <- paste0("(Comp|Fac|fa|rr)[_]?", i, ".*")
        rows <- grep(pat, coef_names, ignore.case = TRUE)
        if (length(rows) > 0) {
            matches <- coef_names[rows]
            vals <- coefs[rows, 1]
            clean_ids <- sub(paste0(".*", "(Comp|Fac|fa|rr)[_]?", i, "(_|:)?"), "", matches)
            clean_ids <- sub("^:", "", clean_ids)
            scores_list[[i]] <- data.frame(Genotype = clean_ids, Value = vals, Factor = i)
        }
    }
    scores_df <- do.call(rbind, scores_list)
    has_scores <- !is.null(scores_df) && nrow(scores_df) > 0

    f_mat <- if (has_scores) xtabs(Value ~ Genotype + Factor, data = scores_df) else NULL
    if (!is.null(f_mat)) class(f_mat) <- "matrix"
    if (has_scores && ncol(f_mat) < k) {
        pad <- matrix(0, nrow(f_mat), k - ncol(f_mat))
        f_mat <- cbind(f_mat, pad)
    }

    # 3. Rotation Logic ---------------------------------------------------------
    lambda_rot <- lambda_mat
    f_rot <- f_mat
    var_exp <- rep(NA, k)

    if (rotate != "none" && nrow(lambda_mat) > 1 && ncol(lambda_mat) > 0) {
        # PCA (SVD)
        svd_res <- svd(lambda_mat)
        V <- svd_res$v
        lambda_pca <- lambda_mat %*% V
        f_pca <- if (has_scores) f_mat %*% V else NULL

        if (rotate == "mean" && k >= 1) {
            # Rotate to Mean strategy:
            # Align Factor 1 to be the "General Performance" dimension.
            # Simple heuristic: Ensure sum of F1 loadings is positive.
            # (For rigorous 'Rotate to Mean', we would rotate to align with unit vector,
            #  but SVD F1 is usually this if positive manifold exists).

            if (sum(lambda_pca[, 1]) < 0) {
                lambda_pca[, 1] <- -lambda_pca[, 1]
                V[, 1] <- -V[, 1]
                if (has_scores) f_pca[, 1] <- -f_pca[, 1]
            }
            lambda_rot <- lambda_pca
            f_rot <- f_pca
        } else {
            lambda_rot <- lambda_pca
            f_rot <- f_pca
        }

        # Recalculate Variance Explained
        tot_var <- sum(lambda_rot^2) # Trace(L'L)
        col_vars <- colSums(lambda_rot^2)
        var_exp <- (col_vars / tot_var) * 100
    }

    colnames(lambda_rot) <- paste0("Fac", 1:ncol(lambda_rot))
    if (has_scores) colnames(f_rot) <- paste0("Fac", 1:ncol(f_rot))

    # 4. Reconstruct & VAF ------------------------------------------------------
    G_est <- (lambda_rot %*% t(lambda_rot)) + diag(psi_df$Psi)
    C_est <- tryCatch(cov2cor(G_est), error = function(e) G_est)

    G_diag <- diag(G_est)
    vaf_df <- data.frame(Group = target_sites, stringsAsFactors = FALSE)
    total_vaf <- numeric(length(target_sites))

    for (i in 1:ncol(lambda_rot)) {
        v_fac <- lambda_rot[, i]^2
        pct <- ifelse(G_diag > 1e-9, (v_fac / G_diag) * 100, 0)
        vaf_df[[paste0("VAF_Fac", i)]] <- round(pct, 2)
        total_vaf <- total_vaf + pct
    }
    vaf_df$Total_VAF <- round(total_vaf, 2)

    # 5. FAST Indices Call ------------------------------------------------------
    temp_res <- list(loadings = list(rotated = lambda_rot), scores = list(rotated = f_rot))
    fast_df <- NULL
    if (has_scores) {
        fast_df <- calculate_fast_indices(temp_res, k = k)
    }

    res <- list(
        loadings = list(raw = lambda_mat, rotated = lambda_rot),
        scores = list(raw = f_mat, rotated = f_rot),
        var_comp = list(psi = psi_df, vaf = vaf_df),
        matrices = list(G = G_est, Cor = C_est),
        fast = fast_df,
        meta = list(k = k, group = site_var, genotype = gen_col_name, var_explained = var_exp)
    )
    class(res) <- "fa_model"

    if (!is.null(annotation)) {
        if (is.data.frame(annotation)) {
            res <- annotate_model(res, df = annotation)
        } else if (is.list(annotation)) {
            res <- do.call(annotate_model, c(list(object = res), annotation))
        }
    }

    cli::cli_alert_success("Extraction standard complete (k={k}, rotation={rotate}).")
    return(res)
}

#' Calculate FAST Indices (OP and RMSD)
#'
#' Calculates "Factor Analytic Selection Tools" indices from a FA Model.
#'
#' @param fa_object An object of class `fa_model` or a list with `loadings` and `scores`.
#' @param k Number of factors.
#' @return A dataframe with Genotype, OP, and RMSD.
#' @export
calculate_fast_indices <- function(fa_object, k = NULL) {
    L <- fa_object$loadings$rotated
    F_sc <- fa_object$scores$rotated

    if (is.null(L) || is.null(F_sc)) {
        return(NULL)
    }
    if (is.null(k)) k <- ncol(L)

    # OP (Overall Performance)
    # OP = f_g1 * mean(lambda_1)
    mean_load_1 <- mean(L[, 1])
    op <- F_sc[, 1] * mean_load_1

    # RMSD (Stability)
    rmsd <- rep(0, nrow(F_sc))

    if (k > 1 && ncol(L) > 1) {
        # Interaction Matrix
        L_int <- L[, 2:k, drop = FALSE]
        F_int <- F_sc[, 2:k, drop = FALSE]
        I_mat <- F_int %*% t(L_int)
        rmsd <- sqrt(rowMeans(I_mat^2))
    }

    df <- data.frame(Genotype = rownames(F_sc), OP = op, RMSD = rmsd)
    df <- df[order(df$OP, decreasing = TRUE), ]
    return(df)
}

#' Determine Interaction Classes (iClasses)
#'
#' Classifies genotypes based on their interaction scores (Significant vs NS).
#'
#' @param fa_object An FA model object.
#' @param alpha Significance level for specific interactions (approx).
#' @return Dataframe with iClass columns.
#' @export
calculate_i_classes <- function(fa_object, alpha = 0.05) {
    # TO BE IMPLEMENTED: Robust iClass determination
    # For now, simplistic sign-based classification on Factor 2
    F_sc <- fa_object$scores$rotated
    k <- ncol(F_sc)

    if (k < 2) {
        return(NULL)
    }

    # Simple classification on Factor 2 scores
    # + = Positive Interaction, - = Negative Interaction
    scores <- F_sc[, 2]
    class_labels <- ifelse(scores > 0, "Positive", "Negative")

    data.frame(Genotype = rownames(F_sc), iClass = class_labels)
}

#' Rank Genotypes using Selection Index
#'
#' Ranks genotypes based on a Smith-Cullis index: I = OP - b * RMSD.
#'
#' @param fast_df Dataframe containing OP and RMSD columns.
#' @param weight_stability Numeric. penalty 'b'.
#' @return The dataframe sorted by Index.
#' @export
rank_genotypes <- function(fast_df, weight_stability = 1) {
    if (inherits(fast_df, "fa_model")) fast_df <- fast_df$fast

    if (!all(c("OP", "RMSD") %in% names(fast_df))) {
        cli::cli_abort("Dataframe must contain 'OP' and 'RMSD' columns.")
    }

    fast_df$Index <- fast_df$OP - (weight_stability * fast_df$RMSD)
    fast_df <- fast_df[order(fast_df$Index, decreasing = TRUE), ]
    fast_df$Rank <- 1:nrow(fast_df)
    return(fast_df)
}
