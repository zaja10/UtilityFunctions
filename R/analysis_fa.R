#' Extract and Rotate Factor Analytic Model Parameters
#'
#' The primary engine of the toolkit. This function parses a fitted ASReml-R model object to
#' extract loadings, scores (BLUPs), and specific variances from a Factor Analytic (FA)
#' or Reduced Rank (RR) model term.
#'
#' It performs a Singular Value Decomposition (SVD) rotation to align the solution
#' with the Principal Component axis. This rotation is required to validly calculate
#' the **Factor Analytic Selection Tools (FAST)** indices (Overall Performance and Stability)
#' as proposed by Smith & Cullis (2018).
#'
#' @param model A fitted object of class \code{asreml}.
#' @param classify A character string defining the main latent model term.
#'        Must match the ASReml syntax used in the model call.
#'        \itemize{
#'          \item Standard FA: \code{"fa(Site, 2):Genotype"}
#'          \item Genomic FA: \code{"fa(Site, 2):vm(Genotype, grm)"}
#'          \item Reduced Rank: \code{"rr(Site, 2):Genotype"}
#'        }
#' @param psi_term A character string defining the specific variance term (Optional).
#'        Required only for **Reduced Rank** models where the diagonal variance is fitted
#'        separately (e.g., \code{"diag(Site):Genotype"}). If NULL for an RR model,
#'        specific variances are assumed to be zero.
#' @param rotate Logical. If \code{TRUE} (default), loadings and scores are rotated
#'        to the Principal Component solution via SVD. If \code{FALSE}, raw parameters
#'        are returned (warning: FAST indices may be invalid without rotation).
#'
#' Extract and Rotate Factor Analytic Model Parameters (Robust)
#'
#' A universal parser for ASReml-R Factor Analytic (FA) and Reduced Rank (RR) models.
#' Handles standard MET models `fa(Site, k):Genotype`, Multi-trait models `fa(Trait, k):Genotype`,
#' and Genomic models `fa(Site, k):vm(Genotype, G)`.
#'
#' @param model A fitted object of class \code{asreml}.
#' @param classify Character. The term to extract. E.g., \code{"fa(Site, 2):Genotype"}.
#' @param psi_term Character (Optional). For RR models, the term containing specific variances.
#'        E.g., \code{"diag(Site):Genotype"}. If NULL for RR, attempts automatic detection.
#' @param rotate Logical. Perform SVD rotation to PC solution? Default TRUE.
#'
#' @return An object of class \code{fa_model}.
#' @importFrom stats coef cov2cor sd terms
#' @import cli
#' @export
fit_fa_model <- function(model, classify, psi_term = NULL, rotate = TRUE) {
    cli::cli_h1("Extracting FA/RR Model Parameters")

    # 1. SETUP & PARSING --------------------------------------------------------
    if (!inherits(model, "asreml")) cli::cli_abort("Object must be of class {.cls asreml}.")

    # Robust extraction (Handle mocks that already have varcomp vs real objects needing summary)
    vc <- if (!is.null(model$varcomp)) model$varcomp else tryCatch(summary(model)$varcomp, error = function(e) NULL)

    if (is.null(vc)) cli::cli_abort("Could not extract variance components from model.")


    vc_names <- rownames(vc)
    coefs <- coef(model)$random
    coef_names <- rownames(coefs)

    # Initialize optional stats
    n_obs_long <- NULL
    rep_df <- NULL
    blues_df <- NULL

    # Robust Term Parsing
    clean_str <- gsub("\\s+", "", classify)

    # Split interaction: "fa(Site,2)" and "Genotype"
    # Logic: finding the part with "fa(" or "rr("
    parts <- strsplit(clean_str, ":")[[1]]
    latent_idx <- grep("(fa|rr)\\(", parts)

    if (length(latent_idx) == 0) cli::cli_abort("No 'fa()' or 'rr()' term found in {.arg classify}.")

    latent_part <- parts[latent_idx[1]] # e.g., "fa(Site,2)"
    gen_part_raw <- parts[-latent_idx[1]] # e.g., "vm(Genotype,G)"

    # Extract Group Name (Site/Trait) and k factors
    # Matches "fa(Name, k)"
    group_var <- sub("(fa|rr)\\(([^,]+),.*", "\\2", latent_part)
    k_str <- sub(".*,([0-9]+)\\).*", "\\1", latent_part)
    k <- as.numeric(k_str)

    if (is.na(k)) cli::cli_abort("Could not determine number of factors (k) from string.")

    # Extract Genotype Name (Handle vm(), ide(), etc.)
    # Removes outer function wrapper if present
    if (length(gen_part_raw) > 0) {
        if (grepl("\\(", gen_part_raw)) {
            gen_col_name <- sub("^[a-z]+\\(([^,)]+).*", "\\1", gen_part_raw)
        } else {
            gen_col_name <- gen_part_raw
        }
    } else {
        # Fallback if classify is just "fa(Site,2)" (rare but possible in simple structures)
        gen_col_name <- "Genotype"
    }

    is_rr <- grepl("rr\\(", clean_str)
    type_lbl <- if (is_rr) "Reduced Rank (RR)" else "Factor Analytic (FA)"

    cli::cli_alert_info("Model: {type_lbl} | Group: {.val {group_var}} | Genotype: {.val {gen_col_name}} | k={k}")

    # 2. IDENTIFY ASREML TERM PREFIX --------------------------------------------
    # ASReml outputs variance components with prefixes like "fa(Site,2)!Site!..."
    # We need to find the exact string ASReml uses.

    term_regex <- NULL

    if (!is_rr) {
        # Standard FA: Look for a variance term ending in !var
        # e.g., "fa(Site, 2)!Site!var"
        pat <- paste0(group_var, ".*!var$")
        candidates <- grep(pat, vc_names, value = TRUE)

        if (length(candidates) == 0) {
            # Try looser match
            candidates <- grep("!var$", vc_names, value = TRUE)
            # Filter for group_var
            candidates <- candidates[grepl(group_var, candidates)]
        }

        if (length(candidates) > 0) {
            # Extract prefix: "fa(Site, 2)" from "fa(Site, 2)!SiteA!var"
            # Everything before the first "!"
            term_regex <- strsplit(candidates[1], "!")[[1]][1]
        }
    } else {
        # RR: Look for loading term 1
        # e.g., "rr(Site, 2)!Site!fa_1"
        pat <- paste0(group_var, ".*!(fa|rr)[_]?1$")
        candidates <- grep(pat, vc_names, value = TRUE)
        if (length(candidates) > 0) {
            # Everything before the first "!"
            term_regex <- strsplit(candidates[1], "!")[[1]][1]
        }
    }

    if (is.null(term_regex)) cli::cli_abort("Could not locate model terms in varcomp output. Check spelling of {.arg classify}.")

    # Escape regex special chars for subsequent searches
    term_regex_safe <- gsub("\\(", "\\\\(", term_regex)
    term_regex_safe <- gsub("\\)", "\\\\)", term_regex_safe)

    # 3. EXTRACT LOADINGS (Lambda) ----------------------------------------------
    lambda_list <- vector("list", k)

    for (i in 1:k) {
        # Pattern: Prefix ! GroupLevel ! fa_i
        # Note: Sometimes ASReml uses "fa_1", sometimes "rr_1"
        pat <- paste0("^", term_regex_safe, "!.*!(fa|rr)[_]?", i, "$")
        rows <- grep(pat, vc_names)

        if (length(rows) == 0) next

        full_terms <- vc_names[rows]
        vals <- vc[rows, "component"]

        # Extract Group Levels (Site Names)
        # Remove prefix and suffix
        temp <- sub(paste0("^", term_regex_safe, "!"), "", full_terms)
        grps <- sub("!(fa|rr)[_]?[0-9]+$", "", temp)

        lambda_list[[i]] <- data.frame(
            Group = grps,
            Value = vals,
            Factor = i,
            stringsAsFactors = FALSE
        )
    }

    lambda_df <- do.call(rbind, lambda_list)
    if (is.null(lambda_df) || nrow(lambda_df) == 0) cli::cli_abort("No loadings extracted. Check model convergence.")

    # Matrix Conversion (Base R)
    lambda_mat <- xtabs(Value ~ Group + Factor, data = lambda_df)
    class(lambda_mat) <- "matrix" # strip xtabs class

    # Ensure all factors present
    if (ncol(lambda_mat) < k) {
        pad <- matrix(0, nrow(lambda_mat), k - ncol(lambda_mat))
        lambda_mat <- cbind(lambda_mat, pad)
    }

    # 4. EXTRACT SPECIFIC VARIANCES (Psi) ---------------------------------------
    psi_df <- NULL

    if (!is_rr) {
        # FA: Psi is internal
        pat <- paste0("^", term_regex_safe, "!.*!var$")
        rows <- grep(pat, vc_names)
        full_terms <- vc_names[rows]
        vals <- vc[rows, "component"]

        temp <- sub(paste0("^", term_regex_safe, "!"), "", full_terms)
        grps <- sub("!var$", "", temp)

        psi_df <- data.frame(Group = grps, Psi = vals, stringsAsFactors = FALSE)
    } else {
        # RR: Psi is external
        # Heuristic: Look for variance terms that match the Group names found in Lambda
        target_groups <- rownames(lambda_mat)

        # Identify candidate variance rows (those ending in !var or R!variance or just names)
        # Exclude the loadings we just found
        cand_rows <- vc_names[!grepl(term_regex_safe, vc_names)]

        found_grps <- c()
        found_vals <- c()

        # If user provided a specific term (e.g. "diag(Site)"), filter by it
        if (!is.null(psi_term)) {
            # Normalize user term for regex
            psi_clean <- gsub("\\(", "\\\\(", psi_term)
            psi_clean <- gsub("\\)", "\\\\)", psi_clean)
            cand_rows <- grep(psi_clean, cand_rows, value = TRUE)
        }

        # Match logic: The term must contain the group name
        for (g in target_groups) {
            # Regex: explicit boundary or component match
            # e.g. "diag(Site)!SiteA" or "SiteA!var"
            matches <- grep(paste0("!", g, "(!var|$)"), cand_rows, value = TRUE)

            if (length(matches) == 0) {
                # Try simple match (e.g. "at(Site, A):Genotype")
                matches <- grep(g, cand_rows, value = TRUE)
            }

            if (length(matches) > 0) {
                # Best guess: take the first match
                idx <- which(vc_names == matches[1])
                found_grps <- c(found_grps, g)
                found_vals <- c(found_vals, vc[idx, "component"])
            }
        }

        if (length(found_grps) > 0) {
            psi_df <- data.frame(Group = found_grps, Psi = found_vals, stringsAsFactors = FALSE)
        } else {
            cli::cli_warn("Could not auto-detect specific variances for RR model. Assuming 0.")
            psi_df <- data.frame(Group = rownames(lambda_mat), Psi = 0, stringsAsFactors = FALSE)
        }
    }

    # 5. SCORES -----------------------------------------------------------------
    scores_list <- vector("list", k)

    for (i in 1:k) {
        # ASReml scores: "Comp_1_Genotype" or "Fac_1_Genotype"
        p1 <- paste0("Comp[_]?", i, ".*", gen_col_name)
        p2 <- paste0("Fac[_]?", i, ".*", gen_col_name)

        rows <- unique(c(grep(p1, coef_names), grep(p2, coef_names)))

        if (length(rows) > 0) {
            full_terms <- coef_names[rows]
            vals <- coefs[rows, 1]

            # Clean Genotype Name
            # Remove "Comp_1_" prefix and "Genotype" prefix if repeated
            # Heuristic: Remove everything up to the Genotype column name and subsequent separator
            # .*: Greedy match ensures we get the last occurrence if repeated, but we anchor to gen_col_name
            genos <- sub(paste0(".*", gen_col_name, ".*[:_]"), "", full_terms)

            scores_list[[i]] <- data.frame(
                Genotype = genos,
                Value = vals,
                Factor = i,
                stringsAsFactors = FALSE
            )
        }
    }

    scores_df <- do.call(rbind, scores_list)
    has_scores <- FALSE
    f_mat <- matrix(0, 1, k)

    if (!is.null(scores_df) && nrow(scores_df) > 0) {
        f_mat <- xtabs(Value ~ Genotype + Factor, data = scores_df)
        f_mat <- matrix(f_mat,
            nrow = nrow(f_mat), ncol = ncol(f_mat),
            dimnames = dimnames(f_mat)
        )
        has_scores <- TRUE
    }

    # 6. ROTATION ---------------------------------------------------------------
    rot_mat <- diag(k)
    if (rotate && nrow(lambda_mat) > 1) {
        svd_res <- svd(lambda_mat)
        V <- svd_res$v
        lambda_rot <- lambda_mat %*% V

        # Sign convention: Force first element of Fac1 to be positive (or sum positive)
        # to ensure reproducibility
        if (sum(lambda_rot[, 1]) < 0) {
            lambda_rot[, 1] <- -lambda_rot[, 1]
            V[, 1] <- -V[, 1]
        }

        f_rot <- if (has_scores) f_mat %*% V else f_mat
        rot_mat <- V
    } else {
        lambda_rot <- lambda_mat
        f_rot <- f_mat
    }
    colnames(lambda_rot) <- colnames(f_rot) <- paste0("Fac", 1:k)

    # 7. RECONSTRUCTION ---------------------------------------------------------
    common_grps <- intersect(rownames(lambda_rot), psi_df$Group)
    if (length(common_grps) == 0) stop("Group mismatch between Lambda and Psi.")

    lam_ord <- lambda_rot[common_grps, , drop = FALSE]
    # Ensure Psi aligned
    psi_vec <- psi_df$Psi[match(common_grps, psi_df$Group)]
    psi_ord <- diag(psi_vec, nrow = length(psi_vec))

    G_est <- (lam_ord %*% t(lam_ord)) + psi_ord
    C_est <- tryCatch(cov2cor(G_est), error = function(e) G_est)

    G_diag <- diag(G_est)

    # Calculate VAF
    vaf_df <- data.frame(Group = common_grps, stringsAsFactors = FALSE)
    total_vaf <- numeric(length(common_grps))

    for (i in 1:k) {
        v_fac <- lam_ord[, i]^2
        pct_fac <- ifelse(G_diag > 1e-8, (v_fac / G_diag) * 100, 0)
        vaf_df[[paste0("VAF_Fac", i)]] <- round(pct_fac, 2)
        total_vaf <- total_vaf + pct_fac
    }
    vaf_df$Total_VAF <- round(total_vaf, 2)

    # 8. FAST INDICES -----------------------------------------------------------
    fast_df <- NULL
    if (has_scores) {
        OP <- f_rot[, 1] * mean(lambda_rot[, 1], na.rm = TRUE)

        if (k > 1) {
            # Vectorized RMSD
            # (Scores_k * Lambda_k') -> Deviation Matrix
            dev_mat <- f_rot[, 2:k, drop = FALSE] %*% t(lambda_rot[, 2:k, drop = FALSE])
            RMSD <- sqrt(rowMeans(dev_mat^2))
        } else {
            RMSD <- rep(0, nrow(f_rot))
        }

        fast_df <- data.frame(
            Genotype = rownames(f_rot),
            OP = OP,
            RMSD = RMSD,
            stringsAsFactors = FALSE
        )
        # Sort desc OP
        fast_df <- fast_df[order(fast_df$OP, decreasing = TRUE), ]
    }

    # 9. SITE BLUP RECONSTRUCTION (Regressed) -----------------------------------
    site_blups_long <- NULL
    if (has_scores) {
        # Prediction: G = F %*% L'
        reg_blups <- f_rot %*% t(lambda_rot)

        # Convert to Long Format manually
        # Create vectors for Genotype, Site, and Value
        genos_vec <- rep(rownames(reg_blups), times = ncol(reg_blups))
        sites_vec <- rep(colnames(reg_blups), each = nrow(reg_blups))
        vals_vec <- as.vector(reg_blups)

        site_blups_long <- data.frame(
            Genotype = genos_vec,
            Site = sites_vec,
            Pred_Value = vals_vec,
            stringsAsFactors = FALSE
        )

        # Merge with N_Obs if available
        if (!is.null(n_obs_long)) {
            site_blups_long <- merge(site_blups_long, n_obs_long,
                by = c("Genotype", "Site"),
                all.x = TRUE
            )
        }
    }

    # 8. OUTPUT AESTHETICS ------------------------------------------------------
    names(vaf_df)[1] <- group_var
    names(psi_df)[1] <- group_var

    out <- list(
        loadings = list(raw = lambda_mat, rotated = lambda_rot),
        scores = list(
            raw = if (has_scores) f_mat else NULL,
            rotated = f_rot
        ),
        var_comp = list(psi = psi_df, vaf = vaf_df),
        matrices = list(G = G_est, Cor = C_est),
        fast = fast_df,
        meta = list(k = k, type = type_lbl, group = group_var, genotype = gen_col_name)
    )
    class(out) <- "fa_model"
    return(out)
}

#' Print Method for FA Object
#' @export
print.fa_model <- function(x, ...) {
    cat("\n=== FACTOR ANALYTIC REPORT (ASRemlFAST) ===\n")
    cat(sprintf(" Model Type: %s\n", x$meta$type))
    cat(sprintf(" Rotation  : SVD (PC Solution)\n"))
    if (!is.null(x$fast)) {
        cat("\n--- Top 5 Genotypes (by OP) ---\n")
        print(head(x$fast[, c("Genotype", "OP", "RMSD")], 5))
    }
    invisible(x)
}

#' Factor Analytic Selection Tools (FAST)
#'
#' Functions to calculate Smith & Cullis (2018) selection indices from a Factor Analytic model.
#'
#' @name fast_tools
#' @return A \code{fast_selection} object or dataframes of indices.
NULL

#' Calculate Overall Performance (OP) and RMSD
#'
#' Calculates OP (Factor 1 performance) and RMSD (deviation/instability) and returns a selection object.
#'
#' @param model An object of class \code{fa_model}.
#' @return A \code{fast_selection} object.
#' @export
calculate_op <- function(model) {
    if (!inherits(model, "fa_model")) stop("Model must be of class fa_model")

    lam <- model$loadings$rotated
    sco <- model$scores$rotated

    if (is.null(sco)) stop("No genotype scores available in model.")

    k <- model$meta$k

    # OP = Factor 1 Score * Mean(Factor 1 Loading)
    mean_lam1 <- mean(lam[, 1], na.rm = TRUE)
    op <- sco[, 1] * mean_lam1

    # RMSD calculation
    if (k > 1) {
        effects_dev <- sco[, 2:k, drop = FALSE] %*% t(lam[, 2:k, drop = FALSE])
        rmsd <- apply(effects_dev, 1, function(x) sqrt(mean(x^2)))
    } else {
        rmsd <- rep(0, nrow(sco))
    }

    df <- data.frame(
        Genotype = rownames(sco),
        OP = op,
        RMSD = rmsd,
        stringsAsFactors = FALSE
    )
    df <- df[order(df$OP, decreasing = TRUE), ]

    obj <- list(
        selection = df,
        source_model = model
    )
    class(obj) <- "fast_selection"
    return(obj)
}

#' Calculate RMSD (Wrapper)
#'
#' @param model An object of class \code{fa_model}.
#' @return A \code{fast_selection} object (same as calculate_op).
#' @export
calculate_rmsd <- function(model) {
    return(calculate_op(model))
}

#' Calculate D-Optimality (Network Efficiency)
#'
#' Quantifies the information content of the trial network.
#'
#' @param object An object of class \code{fa_model}.
#' @return A list containing total_d and site_impact.
#' @export
calculate_d_optimality <- function(object) {
    lam <- object$loadings$rotated
    M_total <- t(lam) %*% lam
    total_det <- det(M_total)

    impact_scores <- numeric(nrow(lam))
    for (i in seq_len(nrow(lam))) {
        lam_min <- lam[-i, , drop = FALSE]
        impact_scores[i] <- (total_det - det(t(lam_min) %*% lam_min)) / total_det * 100
    }
    site_impact <- data.frame(Site = rownames(lam), Impact_Pct = round(impact_scores, 2))
    site_impact <- site_impact[order(site_impact$Impact_Pct, decreasing = TRUE), ]
    return(list(total_d = total_det, site_impact = site_impact))
}

#' Calculate Interaction Classes (iClasses)
#'
#' Implements Smith et al. (2021) to detect Specific Adaptation.
#'
#' @param object An object of class \code{fa_model}.
#' @param factor Integer. Which factor defines the interaction? (Default 2).
#' @param threshold Numeric. Loading magnitude for class assignment.
#'
#' @return A list with site classes and genetic effects.
#' @export
calculate_i_classes <- function(object, factor = 2, threshold = 0.1) {
    if (is.null(object$scores)) stop("Genotype scores (BLUPs) not found.")
    if (factor > object$meta$k) stop("Factor not found.")

    lam <- object$loadings$rotated[, factor]
    sco <- object$scores$rotated[, factor]
    lam1 <- object$loadings$rotated[, 1]
    sco1 <- object$scores$rotated[, 1]

    site_class <- rep("Neutral", length(lam))
    names(site_class) <- names(lam)
    site_class[lam > threshold] <- "Class_Pos"
    site_class[lam < -threshold] <- "Class_Neg"

    ml_pos <- mean(lam[site_class == "Class_Pos"], na.rm = TRUE)
    ml_neg <- mean(lam[site_class == "Class_Neg"], na.rm = TRUE)
    m_perf <- mean(lam1, na.rm = TRUE)

    pred_pos <- if (is.nan(ml_pos)) rep(NA, length(sco)) else (sco1 * m_perf) + (sco * ml_pos)
    pred_neg <- if (is.nan(ml_neg)) rep(NA, length(sco)) else (sco1 * m_perf) + (sco * ml_neg)

    gen_eff <- data.frame(
        Genotype = names(sco), Pred_Pos = pred_pos, Pred_Neg = pred_neg,
        Differential = pred_pos - pred_neg
    )
    gen_eff <- gen_eff[order(gen_eff$Differential, decreasing = TRUE, na.last = TRUE), ]

    return(list(
        site_classes = data.frame(Site = names(lam), Class = site_class),
        gen_effects = gen_eff, meta = list(factor = factor)
    ))
}

#' Calculate Selection Index
#'
#' Combines Overall Performance (OP) and specific Stability penalty (RMSD).
#'
#' @param object An object of class \code{fa_model}.
#' @param weight Numeric. Risk aversion weight applied to RMSD (default 1.0).
#' @return A dataframe sorted by the calculated Index.
#' @export
calculate_index <- function(object, weight = 1.0) {
    if (is.null(object$fast)) stop("No FAST indices found in object.")
    df <- object$fast
    df$Index <- df$OP - (weight * df$RMSD)
    df$Rank <- rank(-df$Index)
    df <- df[order(df$Index, decreasing = TRUE), ]
    return(df[, c("Rank", "Genotype", "OP", "RMSD", "Index")])
}

#' Plot FAST Selection
#'
#' @param x A \code{fast_selection} object.
#' @param type Character. "op_rmsd" (default).
#' @param ... Additional arguments.
#' @export
plot.fast_selection <- function(x, type = "op_rmsd", ...) {
    df <- x$selection
    if (type == "op_rmsd") {
        plot(df$RMSD, df$OP,
            xlab = "Stability (RMSD) [Lower is Stable]",
            ylab = "Overall Performance (OP) [Higher is Better]",
            main = "FAST Selection: OP vs RMSD",
            pch = 19, col = "steelblue", ...
        )
        grid()
        abline(h = mean(df$OP), col = "red", lty = 2)
        abline(v = mean(df$RMSD), col = "red", lty = 2)
    }
}

#' @export
fa.asreml <- function(...) {
    warning("Deprecated: 'fa.asreml' is deprecated. Use 'fit_fa_model()' instead.")
    fit_fa_model(...)
}
