#' Extract and Rotate Factor Analytic Model Parameters
#'
#' A universal parser for ASReml-R Factor Analytic (FA) and Reduced Rank (RR) models.
#' Calculates FAST indices (OP, RMSD), Global Factor Variance, and Site VAF.
#'
#' @param model A fitted object of class \code{asreml}.
#' @param classify Character. The term to extract. E.g., \code{"fa(Site, 2):Genotype"}.
#' @param psi_term Character (Optional). Specific variance term for RR models.
#' @param rotate Logical. Perform SVD rotation to PC solution? Default TRUE.
#' @param annotation Optional dataframe or named list to annotate sites immediately.
#'        Passed to \code{\link{annotate_model}}.
#'
#' @return An object of class \code{fa_model} with 'loadings', 'scores', 'fast', and 'variance'.
#' @importFrom stats coef cov2cor sd
#' @import cli
#' @export
fit_fa_model <- function(model, classify, psi_term = NULL, rotate = TRUE, annotation = NULL) {
    cli::cli_h1("Extracting FA/RR Model Parameters")

    if (!inherits(model, "asreml")) cli::cli_abort("Object must be of class {.cls asreml}.")

    vc <- summary(model)$varcomp
    vc_names <- rownames(vc)
    coefs <- coef(model)$random
    coef_names <- rownames(coefs)

    # 0. Input Validation -------------------------------------------------------
    # Basic check if the classify term is plausible in the model
    # We strip whitespace and parens to see if the main variables exist
    clean_str <- gsub("\\s+", "", classify)

    # 1. Parsing Logic (Improved Regex) -----------------------------------------
    # Extract k (number of factors) - Handles spaces after comma now
    k_match <- regmatches(clean_str, regexec("(fa|rr)\\([^,]+,\\s*([0-9]+)\\)", clean_str))
    if (length(k_match[[1]]) < 3) cli::cli_abort("Could not determine 'k' from classify string: {.val {classify}}")
    k <- as.numeric(k_match[[1]][3])

    # Extract Group (Site/Experiment)
    group_match <- regmatches(clean_str, regexec("(fa|rr)\\(([^,]+),", clean_str))
    group_var <- group_match[[1]][3]

    # Extract Genotype
    parts <- strsplit(clean_str, ":")[[1]]
    gen_part <- parts[!grepl("(fa|rr)\\(", parts)]
    gen_col_name <- gsub("^[a-z]+\\(([^,]+).*\\)$", "\\1", gen_part)
    if (length(gen_col_name) == 0) gen_col_name <- "Genotype"

    is_rr <- grepl("rr\\(", clean_str)

    # 2. Extract Loadings (Lambda) - Refactored to lapply -----------------------
    lambda_list <- lapply(1:k, function(i) {
        # Pattern: contains group name, ends in !fa1 or !rr1
        pat <- paste0("!((fa|rr|comp)[_]?", i, ")$")
        rows <- grep(pat, vc_names, ignore.case = TRUE)

        if (length(rows) == 0) {
            return(NULL)
        }

        vals <- vc[rows, "component"]
        full_names <- vc_names[rows]

        # Robust Site Name Extraction
        tokens <- strsplit(full_names, "!")
        # The site name is usually the second to last token (before fa1)
        sites <- sapply(tokens, function(x) x[length(x) - 1])

        data.frame(Group = sites, Value = vals, Factor = i)
    })

    # Filter NULLs and combine
    lambda_list <- lambda_list[!sapply(lambda_list, is.null)]
    lambda_df <- do.call(rbind, lambda_list)

    if (is.null(lambda_df)) cli::cli_abort("No loadings found. Check syntax or model convergence.")

    lambda_mat <- xtabs(Value ~ Group + Factor, data = lambda_df)
    class(lambda_mat) <- "matrix"

    # Ensure all factors present
    if (ncol(lambda_mat) < k) {
        pad <- matrix(0, nrow(lambda_mat), k - ncol(lambda_mat))
        lambda_mat <- cbind(lambda_mat, pad)
    }

    # 3. Extract Specific Variances (Psi) ---------------------------------------
    psi_df <- NULL
    var_rows <- grep("!var$", vc_names, value = TRUE)
    target_sites <- rownames(lambda_mat)

    found_sites <- c()
    found_psis <- c()

    for (site in target_sites) {
        # Look for "!Site!var" to avoid partial matches
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
        # FIX-03: Abort if FA model but no variances found
        if (!is_rr) {
            cli::cli_abort(c(
                "No specific variances (!var) found for sites in FA model.",
                "i" = "If this is a Reduced Rank (RR) model, ensure 'rr(' is used in classify.",
                "x" = "If this is an FA model, check model specification."
            ))
        }
        psi_df <- data.frame(Group = target_sites, Psi = 0)
    }

    # 4. Extract Scores ---------------------------------------------------------
    scores_list <- vector("list", k)
    for (i in 1:k) {
        pat <- paste0("(Comp|Fac|fa|rr)[_]?", i, ".*", gen_col_name)
        rows <- grep(pat, coef_names, ignore.case = TRUE)

        if (length(rows) > 0) {
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

        # Variance Explained (Eigenvalues) for Plots
        var_exp <- (svd_res$d^2 / sum(svd_res$d^2)) * 100

        lambda_rot <- lambda_mat %*% V

        if (sum(lambda_rot[, 1]) < 0) {
            lambda_rot[, 1] <- -lambda_rot[, 1]
            V[, 1] <- -V[, 1]
        }

        if (has_scores) f_rot <- f_mat %*% V
    }

    colnames(lambda_rot) <- paste0("Fac", 1:k)
    if (has_scores) colnames(f_rot) <- paste0("Fac", 1:k)

    # 6. Reconstruct & VAF ------------------------------------------------------
    common <- intersect(rownames(lambda_rot), psi_df$Group)
    if (length(common) == 0) cli::cli_abort("No matching sites between Loadings and Variances.")

    L_sub <- lambda_rot[common, , drop = FALSE]
    Psi_sub <- psi_df$Psi[match(common, psi_df$Group)]

    # Genetic Covariance Matrix
    G_est <- (L_sub %*% t(L_sub)) + diag(Psi_sub)
    C_est <- tryCatch(cov2cor(G_est), error = function(e) G_est)

    # Calculate VAF per site
    G_diag <- diag(G_est)
    vaf_df <- data.frame(Group = common, stringsAsFactors = FALSE)
    total_vaf <- numeric(length(common))

    for (i in 1:k) {
        # Variance explained by factor i
        v_fac <- L_sub[, i]^2
        # % of Total Genetic Variance
        pct <- ifelse(G_diag > 1e-9, (v_fac / G_diag) * 100, 0)
        vaf_df[[paste0("VAF_Fac", i)]] <- round(pct, 2)
        total_vaf <- total_vaf + pct
    }
    vaf_df$Total_VAF <- round(total_vaf, 2)

    # 7. FAST Indices -----------------------------------------------------------
    fast_df <- NULL
    if (has_scores) {
        op <- f_rot[, 1] * mean(lambda_rot[, 1])
        # Note: RMSD uses rowMeans (simple average).
        # Weighted mean might be better but equals simple mean if factors are balanced.
        rmsd <- if (k > 1) sqrt(rowMeans((f_rot[, 2:k, drop = FALSE] %*% t(lambda_rot[, 2:k, drop = FALSE]))^2)) else rep(0, nrow(f_rot))
        fast_df <- data.frame(Genotype = rownames(f_rot), OP = op, RMSD = rmsd)
        fast_df <- fast_df[order(fast_df$OP, decreasing = TRUE), ]
    }

    res <- list(
        loadings = list(raw = lambda_mat, rotated = lambda_rot),
        scores = list(raw = f_mat, rotated = f_rot),
        var_comp = list(psi = psi_df, vaf = vaf_df),
        matrices = list(G = G_est, Cor = C_est),
        fast = fast_df,
        meta = list(k = k, group = group_var, genotype = gen_col_name, var_explained = var_exp)
    )
    class(res) <- "fa_model"

    # 8. Apply Annotation if Provided -------------------------------------------
    if (!is.null(annotation)) {
        # If it's a dataframe, pass as df, else pass as ...
        if (is.data.frame(annotation)) {
            res <- annotate_model(res, df = annotation)
        } else if (is.list(annotation)) {
            # Convert list to args
            # We can't easily pass a list to ... in a wrapper without do.call
            # So we assume user passes a named list or we use checks
            res <- do.call(annotate_model, c(list(object = res), annotation))
        }
    }

    cli::cli_alert_success("Extraction Complete (k={k}).")
    return(res)
}
