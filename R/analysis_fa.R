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
#' @importFrom dplyr %>% mutate select arrange desc bind_rows distinct left_join
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom stringr str_extract
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom stats coef cov2cor sd
#' @export
fit_fa_model <- function(model, classify, psi_term = NULL, rotate = TRUE) {
    cat("--- Starting FA/RR Extraction ---\n")

    # 1. SETUP ------------------------------------------------------------------
    vc <- summary(model)$varcomp
    vc_names <- rownames(vc)
    coefs <- coef(model)$random
    coef_names <- rownames(coefs)

    clean_str <- gsub("\\s+", "", classify)
    parts <- strsplit(clean_str, ":")[[1]]
    latent_part <- parts[grep("(fa|rr)\\(", parts)]
    gen_part_raw <- parts[grep("(fa|rr)\\(", parts, invert = TRUE)]

    site_col <- sub("(fa|rr)\\(([^,]+),.*", "\\2", latent_part)
    k <- as.numeric(sub(".*,([0-9]+)\\).*", "\\1", latent_part))

    if (grepl("\\(", gen_part_raw)) {
        gen_col_name <- sub("^[a-z]+\\(([^,)]+).*", "\\1", gen_part_raw)
    } else {
        gen_col_name <- gen_part_raw
    }

    is_rr <- grepl("rr\\(", clean_str)
    model_type <- if (is_rr) "Reduced Rank (RR)" else "Factor Analytic (FA)"

    cat(sprintf("-> Type: %s | Site: '%s' | Factors: %d\n", model_type, site_col, k))

    # 1.1 DATA STATS (REPLICATION) ----------------------------------------------
    data_name <- as.character(model$call$data)
    rep_df <- NULL
    if (exists(data_name)) {
        raw_df <- get(data_name)
        # Ensure columns exist
        if (all(c(site_col, gen_col_name) %in% names(raw_df))) {
            # Basic Aggregate: Count observations
            counts <- table(raw_df[[gen_col_name]], raw_df[[site_col]])
            n_obs_long <- as.data.frame(counts)
            colnames(n_obs_long) <- c("Genotype", "Site", "n_obs")

            # Replicated status (Overall)
            total_n <- rowSums(counts)
            rep_df <- data.frame(Genotype = names(total_n), n_obs_total = as.numeric(total_n)) %>%
                mutate(replicated = n_obs_total > 1)
        }
    }

    # 1.2 BLUEs (Fixed Effects) -------------------------------------------------
    fixed_part <- coef(model)$fixed
    blue_rows <- grep(gen_col_name, rownames(fixed_part))
    blues_df <- NULL
    if (length(blue_rows) > 0) {
        # Extract
        b_vals <- fixed_part[blue_rows, 1]

        # Try to get SE if available
        v_fixed <- tryCatch(sqrt(model$vcoeff$fixed[blue_rows]), error = function(e) rep(NA, length(b_vals)))

        # Labels: Remove "Genotype_" prefix if present
        # Usually ASReml output is "Genotype_A" or "GenotypeA" depending on factor
        # We try to strip everything before the level name
        # Heuristic: Remove the column name
        g_names <- sub(paste0(".*", gen_col_name, "(_|)?"), "", rownames(fixed_part)[blue_rows])

        blues_df <- data.frame(Genotype = g_names, BLUE = b_vals, BLUE_SE = v_fixed)

        # Center BLUEs (as requested)
        blues_df$BLUE <- blues_df$BLUE - mean(blues_df$BLUE, na.rm = TRUE)
    }

    # 2. IDENTIFY PREFIX --------------------------------------------------------
    actual_term_prefix <- NULL
    if (!is_rr) {
        # Standard FA
        psi_candidates <- grep(paste0(site_col, ".*!var$"), vc_names, value = TRUE)
        if (length(psi_candidates) == 0) stop(paste("Could not find !var for site:", site_col))
        actual_term_prefix <- strsplit(psi_candidates[1], "!")[[1]][1]
    } else {
        # Reduced Rank
        load_candidates <- grep(paste0(site_col, ".*!(rr|fa)[_]?1"), vc_names, value = TRUE)
        if (length(load_candidates) == 0) stop(paste("Could not find RR loadings for site:", site_col))
        actual_term_prefix <- strsplit(load_candidates[1], "!")[[1]][1]
    }
    term_regex <- gsub("\\(", "\\\\(", actual_term_prefix) %>% gsub("\\)", "\\\\)", .)

    # 3. LOADINGS ---------------------------------------------------------------
    lambda_df <- data.frame()
    for (i in 1:k) {
        pat <- paste0("^", term_regex, "!.*!(fa|rr)[_]?", i, "$")
        rows <- grep(pat, vc_names)
        if (length(rows) == 0) next
        tmp <- data.frame(FullTerm = vc_names[rows], Value = vc[rows, "component"], Factor = i) %>%
            mutate(Site = str_extract(FullTerm, paste0("(?<=", term_regex, "!)(.*?)(?=!(fa|rr)[_]?", i, ")")))
        lambda_df <- bind_rows(lambda_df, tmp)
    }
    if (nrow(lambda_df) == 0) stop("No loadings found.")
    lambda_mat <- lambda_df %>%
        select(Site, Factor, Value) %>%
        pivot_wider(names_from = Factor, values_from = Value) %>%
        column_to_rownames("Site") %>%
        as.matrix()

    # 4. PSI (SPECIFIC VARIANCES) -----------------------------------------------
    psi_df <- NULL
    if (!is_rr) {
        psi_rows <- grep(paste0("^", term_regex, "!.*!var$"), vc_names)
        psi_df <- data.frame(FullTerm = vc_names[psi_rows], Psi = vc[psi_rows, "component"]) %>%
            mutate(Site = str_extract(FullTerm, paste0("(?<=", term_regex, "!)(.*?)(?=!var)")))
    } else {
        if (!is.null(psi_term)) {
            cat("-> Hunting for specific variances...\n")
            candidate_rows <- vc_names[!grepl(paste0("^", term_regex, "!"), vc_names)]
            sites_to_find <- rownames(lambda_mat)
            found_psi <- list()

            for (site in sites_to_find) {
                site_pat <- paste0("(_|!)", site, "(!var)?$")
                match <- grep(site_pat, candidate_rows, value = TRUE)

                if (length(match) > 0) {
                    idx <- which(vc_names == match[1])
                    found_psi[[site]] <- data.frame(Site = site, Psi = vc[idx, "component"])
                }
            }
            if (length(found_psi) > 0) {
                psi_df <- do.call(rbind, found_psi)
            } else {
                stop("Could not find specific variances for your sites.")
            }
        } else {
            psi_df <- data.frame(Site = rownames(lambda_mat), Psi = 0)
        }
    }

    # 5. SCORES -----------------------------------------------------------------
    scores_df <- data.frame()
    for (i in 1:k) {
        p1 <- paste0("Comp[_]?", i, ".*", gen_col_name)
        p2 <- paste0("Fac[_]?", i, ".*", gen_col_name)
        rows <- unique(c(grep(p1, coef_names), grep(p2, coef_names)))
        if (length(rows) > 0) {
            tmp <- data.frame(FullTerm = coef_names[rows], Value = coefs[rows, 1], Factor = i) %>%
                mutate(Genotype = sub(paste0(".*", gen_col_name, "[:_]"), "", FullTerm))
            scores_df <- bind_rows(scores_df, tmp)
        }
    }
    has_scores <- FALSE
    f_mat <- matrix(0, 1, k)
    if (nrow(scores_df) > 0) {
        f_mat <- scores_df %>%
            distinct(Genotype, Factor, .keep_all = T) %>%
            select(Genotype, Factor, Value) %>%
            pivot_wider(names_from = Factor, values_from = Value, values_fill = 0) %>%
            column_to_rownames("Genotype") %>%
            as.matrix()
        has_scores <- TRUE
    }

    # 6. ROTATION ---------------------------------------------------------------
    rot_mat <- diag(k)
    if (rotate && nrow(lambda_mat) > 1) {
        svd_res <- svd(lambda_mat)
        V <- svd_res$v
        lambda_rot <- lambda_mat %*% V
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
    psi_ord <- diag(psi_df$Psi[match(common_sites, psi_df$Site)])
    G_est <- (lam_ord %*% t(lam_ord)) + psi_ord
    C_est <- tryCatch(cov2cor(G_est), error = function(e) G_est)

    G_diag <- diag(G_est)
    vaf_list <- list(Site = common_sites)
    total_vaf <- numeric(length(common_sites))

    for (i in 1:k) {
        v_fac <- lam_ord[, i]^2
        pct_fac <- ifelse(G_diag > 1e-8, (v_fac / G_diag) * 100, 0)
        vaf_list[[paste0("VAF_Fac", i)]] <- round(pct_fac, 2)
        total_vaf <- total_vaf + pct_fac
    }
    vaf_list$Total_VAF <- round(total_vaf, 2)
    vaf_df <- as.data.frame(vaf_list)

    # 8. FAST INDICES -----------------------------------------------------------
    fast_df <- NULL
    if (has_scores) {
        OP <- f_rot[, 1] * mean(lambda_rot[, 1], na.rm = TRUE)
        RMSD <- if (k > 1) apply(f_rot[, 2:k, drop = F] %*% t(lambda_rot[, 2:k, drop = F]), 1, function(x) sqrt(mean(x^2))) else rep(0, nrow(f_rot))
        fast_df <- data.frame(Genotype = rownames(f_rot), OP = OP, RMSD = RMSD) %>% arrange(desc(OP))
    }

    # 9. SITE BLUP RECONSTRUCTION (Regressed) -----------------------------------
    site_blups_long <- NULL
    if (has_scores) {
        # Formula: G = F %*% L'
        # This gives the Genetic Value predicted by the Latent Factors for each Site
        reg_blups <- f_rot %*% t(lambda_rot)

        # Convert to Long Format
        site_blups_long <- as.data.frame(reg_blups) %>%
            tibble::rownames_to_column("Genotype") %>%
            tidyr::pivot_longer(-Genotype, names_to = "Site", values_to = "Pred_Value")

        # Merge with N_Obs if available
        if (exists("n_obs_long")) {
            # Be careful with Site names matching
            # lambda_rot rownames are Sites
            site_blups_long <- site_blups_long %>%
                dplyr::left_join(n_obs_long, by = c("Genotype", "Site"))
        }
    }

    out <- list(
        loadings = list(raw = lambda_mat, rotated = lambda_rot),
        scores = if (has_scores) list(raw = f_mat, rotated = f_rot, blups_in_met = site_blups_long) else NULL,
        blues = blues_df,
        data_stats = rep_df,
        var_comp = list(psi = psi_df, vaf = vaf_df),
        matrices = list(G = G_est, Cor = C_est),
        fast = fast_df, rotation_matrix = rot_mat, meta = list(k = k, type = model_type, classify = classify)
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
        RMSD = rmsd
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

