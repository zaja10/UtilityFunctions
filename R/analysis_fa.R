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
#' @return An object of class \code{fa_model} containing:
#' \describe{
#'   \item{loadings}{A list containing \code{raw} and \code{rotated} site loading matrices.}
#'   \item{scores}{A list containing \code{raw}, \code{rotated} genotype scores, and \code{blups_in_met} (the Site-specific predicted values reconstructed from the latent factors).}
#'   \item{blues}{A dataframe of the centered BLUEs (Fixed Effects) for the Genotype, if fitted.}
#'   \item{data_stats}{A dataframe containing replication statistics (\code{n_obs}, \code{replicated}) derived from the original dataset.}
#'   \item{var_comp}{A list containing the extracted specific variances (\code{psi}) and a table of Variance Accounted For (\code{vaf}).}
#'   \item{matrices}{The estimated Genetic Covariance (\code{G}) and Genetic Correlation (\code{Cor}) matrices.}
#'   \item{fast}{A dataframe containing the FAST indices: Overall Performance (\code{OP}) and Stability (\code{RMSD}).}
#'   \item{meta}{Metadata regarding the model structure (k factors, model type).}
#' }
#'
#' @importFrom stats coef cov2cor sd
#' @export
fit_fa_model <- function(model, classify, psi_term = NULL, rotate = TRUE) {
    message("--- Starting FA/RR Extraction ---")

    # 1. SETUP ------------------------------------------------------------------
    vc <- summary(model)$varcomp
    vc_names <- rownames(vc)
    coefs <- coef(model)$random
    coef_names <- rownames(coefs)

    # Clean spaces
    clean_str <- gsub("\\s+", "", classify)

    # Robust Parsing of 'fa(Site,k):Genotype'
    # Use standard grep to find the fa/rr term
    term_labels <- attr(terms(as.formula(paste("~", classify))), "term.labels")
    latent_idx <- grep("(fa|rr)\\(", term_labels)
    if (length(latent_idx) == 0) {
        # Fallback to simple split if formula parsing fails
        parts <- strsplit(clean_str, ":")[[1]]
        latent_idx <- grep("(fa|rr)\\(", parts)
        if (length(latent_idx) == 0) stop("Could not identify 'fa()' or 'rr()' term in classify string.")
        latent_part <- parts[latent_idx[1]]
        gen_part_raw <- parts[-latent_idx[1]]
    } else {
        # Extract purely using regex on the input string to be safe
        # Assume format fa(Site, k)
        parts <- strsplit(clean_str, ":")[[1]]
        latent_part <- parts[grep("(fa|rr)\\(", parts)][1]
        gen_part_raw <- parts[grep("(fa|rr)\\(", parts, invert = TRUE)][1]
    }

    # Extract Site Name and K
    # Matches "fa(SiteName,2)" -> extracts "SiteName" and "2"
    site_col <- sub("(fa|rr)\\(([^,]+),.*", "\\2", latent_part)
    k <- as.numeric(sub(".*,([0-9]+)\\).*", "\\1", latent_part))

    # Extract Genotype Name
    # Handles "Genotype" or "vm(Genotype, ...)"
    if (grepl("\\(", gen_part_raw)) {
        gen_col_name <- sub("^[a-z]+\\(([^,)]+).*", "\\1", gen_part_raw)
    } else {
        gen_col_name <- gen_part_raw
    }

    is_rr <- grepl("rr\\(", clean_str)
    model_type <- if (is_rr) "Reduced Rank (RR)" else "Factor Analytic (FA)"

    message(sprintf("-> Type: %s | Site: '%s' | Factors: %d", model_type, site_col, k))

    # 1.1 DATA STATS (REPLICATION) ----------------------------------------------
    # Use Base R table() instead of dplyr grouping
    data_name <- as.character(model$call$data)
    rep_df <- NULL
    n_obs_long <- NULL

    if (exists(data_name)) {
        raw_df <- get(data_name)
        # Ensure columns exist
        if (all(c(site_col, gen_col_name) %in% names(raw_df))) {
            # Basic Aggregate: Count observations
            counts <- table(raw_df[[gen_col_name]], raw_df[[site_col]])

            # Convert to Long Format (Genotype, Site, n_obs)
            n_obs_long <- as.data.frame(counts)
            colnames(n_obs_long) <- c("Genotype", "Site", "n_obs")

            # Replicated status (Overall)
            total_n <- rowSums(counts)
            rep_df <- data.frame(
                Genotype = names(total_n),
                n_obs_total = as.numeric(total_n),
                replicated = as.numeric(total_n) > 1,
                stringsAsFactors = FALSE
            )
        }
    }

    # 1.2 BLUEs (Fixed Effects) -------------------------------------------------
    fixed_part <- coef(model)$fixed
    # Strict matching to avoid capturing interactions
    # e.g., if Genotype is "G", avoid matching "Site:G"
    # ASReml fixed effects usually named "Term_Level"
    blue_rows <- grep(paste0("^", gen_col_name, "(_|$)"), rownames(fixed_part))

    blues_df <- NULL
    if (length(blue_rows) > 0) {
        b_vals <- fixed_part[blue_rows, 1]
        v_fixed <- tryCatch(sqrt(model$vcoeff$fixed[blue_rows]), error = function(e) rep(NA, length(b_vals)))

        # Clean Genotype Names
        g_names <- sub(paste0(".*", gen_col_name, "(_|)?"), "", rownames(fixed_part)[blue_rows])

        blues_df <- data.frame(
            Genotype = g_names,
            BLUE = b_vals,
            BLUE_SE = v_fixed,
            stringsAsFactors = FALSE
        )

        # Center BLUEs
        blues_df$BLUE <- blues_df$BLUE - mean(blues_df$BLUE, na.rm = TRUE)
    }

    # 2. IDENTIFY PREFIX --------------------------------------------------------
    actual_term_prefix <- NULL
    if (!is_rr) {
        # Standard FA: Look for "!var"
        # Example: "fa(Site, 2)!Site_Level!var" or "fa(Site, 2)!var"
        psi_candidates <- grep(paste0(site_col, ".*!var$"), vc_names, value = TRUE)
        if (length(psi_candidates) == 0) stop(paste("Could not find !var for site:", site_col))
        actual_term_prefix <- strsplit(psi_candidates[1], "!")[[1]][1]
    } else {
        # RR: Look for loadings
        load_candidates <- grep(paste0(site_col, ".*!(rr|fa)[_]?1"), vc_names, value = TRUE)
        if (length(load_candidates) == 0) stop(paste("Could not find RR loadings for site:", site_col))
        actual_term_prefix <- strsplit(load_candidates[1], "!")[[1]][1]
    }
    # term_regex <- gsub("\\(", "\\\\(", actual_term_prefix) %>% gsub("\\)", "\\\\)", .)
    # No pipe!
    term_regex <- gsub("\\(", "\\\\(", actual_term_prefix)
    term_regex <- gsub("\\)", "\\\\)", term_regex)

    # 3. LOADINGS ---------------------------------------------------------------
    # Loop k factors
    lambda_list <- vector("list", k)

    for (i in 1:k) {
        pat <- paste0("^", term_regex, "!.*!(fa|rr)[_]?", i, "$")
        rows <- grep(pat, vc_names)
        if (length(rows) == 0) next

        # Extract Site Names using regex capture groups
        # Pattern: prefix!SITE_NAME!(fa|rr)_i
        # We replace the prefix and the suffix with empty string to get site
        # Note: This assumes the Site part is exactly between the first ! and the last !(fa|rr)

        full_terms <- vc_names[rows]
        vals <- vc[rows, "component"]

        # Robust extraction
        # Remove prefix
        temp <- sub(paste0("^", term_regex, "!"), "", full_terms)
        # Remove suffix
        sites <- sub(paste0("!(fa|rr)[_]?", i, "$"), "", temp)

        lambda_list[[i]] <- data.frame(
            Site = sites,
            Value = vals,
            Factor = i,
            stringsAsFactors = FALSE
        )
    }

    lambda_df <- do.call(rbind, lambda_list)
    if (is.null(lambda_df) || nrow(lambda_df) == 0) stop("No loadings found.")

    # Pivot to Matrix (Site x Factor)
    # Base R alternative to pivot_wider: xtabs
    lambda_mat <- xtabs(Value ~ Site + Factor, data = lambda_df)
    # Convert 'table' class to pure matrix
    lambda_mat <- matrix(lambda_mat,
        nrow = nrow(lambda_mat), ncol = ncol(lambda_mat),
        dimnames = dimnames(lambda_mat)
    )

    # Ensure correct column order
    if (ncol(lambda_mat) < k) {
        # Pad with 0 if some factors missing (rare)
        pad <- matrix(0, nrow(lambda_mat), k - ncol(lambda_mat))
        lambda_mat <- cbind(lambda_mat, pad)
    }

    # 4. PSI (SPECIFIC VARIANCES) -----------------------------------------------
    psi_df <- NULL
    if (!is_rr) {
        psi_rows <- grep(paste0("^", term_regex, "!.*!var$"), vc_names)
        full_terms <- vc_names[psi_rows]
        vals <- vc[psi_rows, "component"]

        # Extract Site
        temp <- sub(paste0("^", term_regex, "!"), "", full_terms)
        sites <- sub("!var$", "", temp)

        psi_df <- data.frame(Site = sites, Psi = vals, stringsAsFactors = FALSE)
    } else {
        if (!is.null(psi_term)) {
            message("-> Hunting for specific variances...")
            candidate_rows <- vc_names[!grepl(paste0("^", term_regex, "!"), vc_names)]
            sites_to_find <- rownames(lambda_mat)

            found_sites <- c()
            found_vals <- c()

            for (site in sites_to_find) {
                # Look for term ending in !Site or !Site!var or _Site
                # This is heuristic and might need tuning for complex diag terms
                match <- grep(site, candidate_rows, fixed = TRUE, value = TRUE)
                # Filter for var/component
                if (length(match) > 0) {
                    # Pick the one that looks like a variance
                    idx <- which(vc_names == match[1])
                    found_sites <- c(found_sites, site)
                    found_vals <- c(found_vals, vc[idx, "component"])
                }
            }

            if (length(found_sites) > 0) {
                psi_df <- data.frame(Site = found_sites, Psi = found_vals, stringsAsFactors = FALSE)
            } else {
                stop("Could not find specific variances for your sites.")
            }
        } else {
            psi_df <- data.frame(Site = rownames(lambda_mat), Psi = 0, stringsAsFactors = FALSE)
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
            # Heuristic: Remove everything up to the Genotype column name
            genos <- sub(paste0(".*", gen_col_name, "[:_]"), "", full_terms)

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
    common_sites <- intersect(rownames(lambda_rot), psi_df$Site)
    if (length(common_sites) == 0) stop("Site mismatch between Lambda and Psi.")

    lam_ord <- lambda_rot[common_sites, , drop = FALSE]
    # Ensure Psi aligned
    psi_vec <- psi_df$Psi[match(common_sites, psi_df$Site)]
    psi_ord <- diag(psi_vec)

    G_est <- (lam_ord %*% t(lam_ord)) + psi_ord
    C_est <- tryCatch(cov2cor(G_est), error = function(e) G_est)

    G_diag <- diag(G_est)

    # Calculate VAF
    vaf_df <- data.frame(Site = common_sites, stringsAsFactors = FALSE)
    total_vaf <- numeric(length(common_sites))

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

    out <- list(
        loadings = list(raw = lambda_mat, rotated = lambda_rot),
        scores = if (has_scores) list(raw = f_mat, rotated = f_rot, blups_in_met = site_blups_long) else NULL,
        blues = blues_df,
        data_stats = rep_df,
        var_comp = list(psi = psi_df, vaf = vaf_df),
        matrices = list(G = G_est, Cor = C_est),
        fast = fast_df,
        rotation_matrix = rot_mat,
        meta = list(k = k, type = model_type, classify = classify)
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
