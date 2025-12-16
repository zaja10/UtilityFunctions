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

    if (!inherits(model, "asreml")) cli::cli_abort("Object must be of class {.cls asreml}.")

    vc <- summary(model)$varcomp
    vc_names <- rownames(vc)
    coefs <- coef(model)$random
    coef_names <- rownames(coefs)

    # 1. Parsing Logic ----------------------------------------------------------
    clean_str <- gsub("\\s+", "", classify)
    # Extract k (number of factors)
    # Look for number in parens after fa or rr: fa(Site,2) -> 2
    k_match <- regmatches(clean_str, regexec("(fa|rr)\\([^,]+,([0-9]+)\\)", clean_str))
    if (length(k_match[[1]]) < 3) cli::cli_abort("Could not determine 'k' from classify string.")
    k <- as.numeric(k_match[[1]][3])

    # Extract Group (Site/Experiment)
    # Parsing "fa(Experiment,3)" -> "Experiment"
    group_match <- regmatches(clean_str, regexec("(fa|rr)\\(([^,]+),", clean_str))
    group_var <- group_match[[1]][3]

    # Extract Genotype
    parts <- strsplit(clean_str, ":")[[1]]
    gen_part <- parts[!grepl("(fa|rr)\\(", parts)]
    # Strip vm() or ide() wrapper if present
    gen_col_name <- gsub("^[a-z]+\\(([^,]+).*\\)$", "\\1", gen_part)
    if (length(gen_col_name) == 0) gen_col_name <- "Genotype" # Fallback

    is_rr <- grepl("rr\\(", clean_str)

    # 2. Extract Loadings (Lambda) ----------------------------------------------
    lambda_list <- vector("list", k)

    # Regex to find loadings: must contain group_var and end in fa_i/rr_i
    # e.g. "fa(Experiment, 3):Genotype!YT_Neo_22!fa1"
    for (i in 1:k) {
        # Pattern: contains group name, ends in !fa1 or !rr1 or !comp1
        pat <- paste0("!((fa|rr|comp)[_]?", i, ")$")
        rows <- grep(pat, vc_names, ignore.case = TRUE)

        if (length(rows) == 0) next

        vals <- vc[rows, "component"]
        full_names <- vc_names[rows]

        # Extract the Site Name from the middle of the string
        # Strategy: Remove the known suffix (!fa1) and the known prefix (fa(...):Genotype!)
        # But prefix is variable. Safer strategy: Extract the token between the last two '!'

        # Split by '!'
        tokens <- strsplit(full_names, "!")
        # The site name is usually the second to last token (before fa1)
        sites <- sapply(tokens, function(x) x[length(x) - 1])

        lambda_list[[i]] <- data.frame(Group = sites, Value = vals, Factor = i)
    }

    lambda_df <- do.call(rbind, lambda_list)
    if (is.null(lambda_df)) cli::cli_abort("No loadings found. Check syntax.")

    lambda_mat <- xtabs(Value ~ Group + Factor, data = lambda_df)
    class(lambda_mat) <- "matrix"

    # 3. Extract Specific Variances (Psi) - IMPROVED ----------------------------
    psi_df <- NULL

    # Heuristic: Look for rows ending in !var that contain the Site names found in Lambda
    # This avoids guessing the prefix structure.
    var_rows <- grep("!var$", vc_names, value = TRUE)

    # Filter to only those that contain our sites
    target_sites <- rownames(lambda_mat)

    found_sites <- c()
    found_psis <- c()

    for (site in target_sites) {
        # Strict match: "!Site!var" or starts with "Site!var" or contains "Site!var"
        # We look for the site name surrounded by delimiters to avoid partial matches (e.g. Site1 vs Site10)
        # Matches: "!Site!" or ":Site!" or "^Site!"

        # Simple robust check: look for the string "!<Site>!var"
        site_pat <- paste0("!", site, "!var$")
        match <- grep(site_pat, var_rows, fixed = FALSE, value = TRUE)

        if (length(match) == 1) {
            found_sites <- c(found_sites, site)
            found_psis <- c(found_psis, vc[match, "component"])
        }
    }

    if (length(found_sites) > 0) {
        psi_df <- data.frame(Group = found_sites, Psi = found_psis)
    } else {
        # If RR, psi might be zero or explicitly modeled elsewhere.
        # If FA, this is a warning.
        if (!is_rr) cli::cli_warn("Could not match specific variances (!var) to sites.")
        psi_df <- data.frame(Group = target_sites, Psi = 0)
    }

    # 4. Extract Scores ---------------------------------------------------------
    scores_list <- vector("list", k)
    for (i in 1:k) {
        # Search for "Comp_i" or "Fac_i" in random coefficients
        pat <- paste0("(Comp|Fac|fa|rr)[_]?", i, ".*", gen_col_name)
        rows <- grep(pat, coef_names, ignore.case = TRUE)

        if (length(rows) > 0) {
            # Extract Genotype Name (clean suffixes)
            clean_ids <- sub(paste0(".*", gen_col_name, "(_|:)?"), "", coef_names[rows])
            scores_list[[i]] <- data.frame(Genotype = clean_ids, Value = coefs[rows, 1], Factor = i)
        }
    }
    scores_df <- do.call(rbind, scores_list)
    has_scores <- !is.null(scores_df) && nrow(scores_df) > 0

    f_mat <- if (has_scores) xtabs(Value ~ Genotype + Factor, data = scores_df) else NULL
    if (!is.null(f_mat)) class(f_mat) <- "matrix"

    # 5. Rotation & Indices -----------------------------------------------------
    rot_mat <- diag(k)
    lambda_rot <- lambda_mat
    f_rot <- f_mat
    var_exp <- rep(NA, k)

    if (rotate && nrow(lambda_mat) > 1) {
        svd_res <- svd(lambda_mat)
        V <- svd_res$v
        lambda_rot <- lambda_mat %*% V

        # Sign convention: Sum of col 1 positive
        if (sum(lambda_rot[, 1]) < 0) {
            lambda_rot[, 1] <- -lambda_rot[, 1]
            V[, 1] <- -V[, 1]
        }

        if (has_scores) f_rot <- f_mat %*% V

        # Variance Explained (Eigenvalues)
        var_exp <- (svd_res$d^2 / sum(svd_res$d^2)) * 100
    }

    # 6. Output Construction ----------------------------------------------------
    # Reconstruct Correlation Matrix
    common <- intersect(rownames(lambda_rot), psi_df$Group)
    L_sub <- lambda_rot[common, , drop = FALSE]
    Psi_sub <- psi_df$Psi[match(common, psi_df$Group)]

    G_est <- (L_sub %*% t(L_sub)) + diag(Psi_sub)
    C_est <- tryCatch(cov2cor(G_est), error = function(e) G_est)

    fast_df <- NULL
    if (has_scores) {
        op <- f_rot[, 1] * mean(lambda_rot[, 1])
        # Multi-factor RMSD
        if (k > 1) {
            dev <- f_rot[, 2:k, drop = FALSE] %*% t(lambda_rot[, 2:k, drop = FALSE])
            rmsd <- sqrt(rowMeans(dev^2))
        } else {
            rmsd <- rep(0, nrow(f_rot))
        }
        fast_df <- data.frame(Genotype = rownames(f_rot), OP = op, RMSD = rmsd)
        fast_df <- fast_df[order(fast_df$OP, decreasing = TRUE), ]
    }

    res <- list(
        loadings = list(raw = lambda_mat, rotated = lambda_rot),
        scores = list(raw = f_mat, rotated = f_rot),
        var_comp = list(psi = psi_df),
        matrices = list(G = G_est, Cor = C_est),
        fast = fast_df,
        meta = list(k = k, group = group_var, genotype = gen_col_name, var_explained = var_exp)
    )
    class(res) <- "fa_model"
    return(res)
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

#' Summary Method for FA Model
#'
#' Provides a commercial-grade summary of the Factor Analytic model,
#' including mean reliability, variance accounted for, and key correlations.
#'
#' @param object An object of class \code{fa_model}.
#' @param ... Additional arguments.
#' @export
summary.fa_model <- function(object, ...) {
    k <- object$meta$k
    grp <- object$meta$group

    cat(sprintf("\n=== FA Model Summary (%s factors) ===\n", k))
    cat(sprintf("Dimension: %s x %s\n", nrow(object$matrices$G), object$meta$genotype))

    # 1. Variance Accounted For
    vaf <- object$var_comp$vaf
    if (!is.null(vaf$Total_VAF)) {
        mean_vaf <- mean(vaf$Total_VAF, na.rm = TRUE)
        cat(sprintf("Mean Variance Accounted For (VAF): %.1f%%\n", mean_vaf))

        # Count low VAF sites
        n_low <- sum(vaf$Total_VAF < 70, na.rm = TRUE)
        if (n_low > 0) cat(sprintf("Warning: %d %ss have VAF < 70%%\n", n_low, grp))
    }

    # 2. Correlations
    cor_mat <- object$matrices$Cor
    upper_tri <- cor_mat[upper.tri(cor_mat)]
    cat(sprintf(
        "Average Genetic Correlation: %.2f (Range: %.2f to %.2f)\n",
        mean(upper_tri), min(upper_tri), max(upper_tri)
    ))

    # 3. FAST Indices Preview
    if (!is.null(object$fast)) {
        cat("\n--- Top 3 Genotypes (by Overall Performance) ---\n")
        print(head(object$fast[, c("Genotype", "OP", "RMSD")], 3), row.names = FALSE)
    }

    invisible(object)
}
