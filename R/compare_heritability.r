#' Compare Site-Specific Heritability Metrics for FA Models
#'
#' @description
#' Calculates and compares three standard heritability metrics (Cullis, Oakey, Standard)
#' for each environment within a Factor Analytic (FA) or Reduced Rank (RR) model.
#' It leverages the heterogeneous genetic variances ($\sigma^2_{g_j}$) estimated by the FA model.
#'
#' @param model A fitted \code{asreml} object.
#' @param fa_object An object of class \code{fa_asreml} produced by \code{fa.asreml()}.
#' @param grm Optional. The genomic relationship matrix (inverse not required).
#'        Required only for the "Oakey" method.
#' @param classify Character string. The prediction term. If NULL, attempts to detect
#'        automatically from \code{fa_object$meta$classify}.
#' @param methods Character vector. Metrics to calculate: c("Cullis", "Oakey", "Standard").
#'
#' @return A data.frame summarizing H2, Genetic Variance (Vg), and Error Metrics per Site.
#'
#' @details
#' \strong{Metrics:}
#' \itemize{
#'   \item \strong{Standard:} \eqn{1 - (PEV / V_g)}. Simple but often biased in unbalanced trials.
#'   \item \strong{Cullis:} \eqn{1 - (\bar{v}_{\Delta}^{BLUP} / 2V_g)}. The standard for phenotypic trials.
#'   \item \strong{Oakey:} \eqn{1 - tr(C^{22}G^{-1})/n}. The gold standard for genomic trials.
#' }
#'
#' \strong{Optimizations:}
#' This function implements two key optimizations for large-scale genomic trials:
#' \enumerate{
#'   \item \strong{Split Predictions:} Predictions are run per-site to avoid calculating massive
#'   PEV matrices (e.g., 20k x 20k) that would exhaust memory.
#'   \item \strong{Trace Trick (Oakey):} Uses \code{sum(A * B)} instead of \code{trace(A \%*\% B)}
#'   to calculate the trace of the product of two symmetric matrices, drastically reducing
#'   computation time.
#' }
#'
#' @importFrom asreml predict.asreml
#' @importFrom MASS ginv
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
compare_h2 <- function(model, fa_object, grm = NULL,
                       classify = NULL,
                       methods = c("Cullis", "Oakey", "Standard")) {
    if (!inherits(fa_object, "fa_asreml")) stop("fa_object must be of class 'fa_asreml'.")

    # Auto-detect classify if NULL
    if (is.null(classify)) {
        if (!is.null(fa_object$meta$classify)) {
            classify <- fa_object$meta$classify
            cat(sprintf("-> Detected prediction term: '%s'\n", classify))
        } else {
            stop("Classify term not provided and not found in fa_object metadata.")
        }
    }

    # 1. SETUP: Extract Site-Specific Genetic Variances
    # -------------------------------------------------
    # We use the diagonal of the reconstructed G matrix (Lambda*Lambda' + Psi)
    G_mat <- fa_object$matrices$G
    site_stats <- data.frame(Site = rownames(G_mat), Vg = diag(G_mat))
    sites <- as.character(site_stats$Site)

    # Pre-calculate GRM Inverse for Oakey if needed
    grm_inv <- NULL
    if ("Oakey" %in% methods) {
        if (is.null(grm)) {
            warning("GRM missing. Skipping Oakey method.")
            methods <- setdiff(methods, "Oakey")
        } else {
            cat("-> Inverting GRM for Oakey calculation...\n")
            # Use tryCatch for robust inversion (singular matrices common in markers)
            grm_inv <- tryCatch(solve(grm), error = function(e) MASS::ginv(grm))
        }
    }

    results_list <- list()

    # Parse site column name from classify string (e.g. "Site:Genotype" -> "Site")
    # Assuming the first term is the site term for splitting
    term_parts <- unlist(strsplit(classify, ":"))

    # CLEAN CLASSIFY FOR PREDICT:
    # Transform "fa(Site,2):Gen" -> "Site:Gen"
    # Transform "rr(Site,2):Gen" -> "Site:Gen"
    # This ensures predict() sums BOTH the rr() term AND the diag() term for RR models.
    clean_factors <- sapply(term_parts, function(x) {
        # Remove function calls parentheses e.g. fa(Site,2) -> Site
        # Regex: remove everything from ( to end, and descriptors?
        # ASReml syntax is tricky. "fa(Site, k)" -> we want "Site".
        # "gn(Gen)" -> "Gen"
        # Simple heuristic: extract the variable name inside the first parens if present.
        if (grepl("\\(", x)) {
            sub("^[a-z]+\\(([^,]+).*", "\\1", x) # Extract first arg of function
        } else {
            x
        }
    })
    classify_for_predict <- paste(clean_factors, collapse = ":")

    site_term <- clean_factors[1] # Heuristic: First term is usually Site/Env

    cat(sprintf("-> Running split predictions for '%s' (Term: %s) across %d sites...\n", classify_for_predict, classify, length(sites)))

    # Progress Bar
    pb <- utils::txtProgressBar(min = 0, max = length(sites), style = 3)

    # 2. LOOP BY SITE (Optimization: Split Prediction)
    # ------------------------------------------------

    for (i in seq_along(sites)) {
        site <- sites[i]
        utils::setTxtProgressBar(pb, i)

        # Run predict for this specific site ONLY.
        # We restrict the level of the site factor to just the current site.
        # This forces ASReml to compute PEV only for the relevant subset, saving RAM.

        # Construct levels list dynamically
        levels_list <- list()
        levels_list[[site_term]] <- site

        tryCatch(
            {
                preds <- asreml::predict.asreml(model,
                    classify = classify_for_predict,
                    levels = levels_list,
                    only = classify_for_predict, # Minimize output size
                    sed = "Cullis" %in% methods, # Needed for Cullis
                    vcov = any(c("Oakey", "Standard") %in% methods), # Needed for Oakey/Standard
                    trace = FALSE
                )

                pvals <- preds$pvals
                total_sed <- preds$sed
                total_vcov <- preds$vcov

                # Filter indices (should be all of them since we filtered by levels, but double check)
                # pvals usually contains the columns specified in classify
                # We find the column that matches the site name
                # Note: pvals column names might match term_parts

                # Since we predicted for ONE site, all rows in pvals should belong to that site
                # (or ASReml might return all with NAs, depending on version, but levels usually restricts it)

                # Let's find the valid rows for this site
                site_col_idx <- which(sapply(pvals, function(c) any(c == site)))
                if (length(site_col_idx) > 0) {
                    idx <- which(pvals[[site_col_idx[1]]] == site)
                } else {
                    # Fallback: assume all rows if prediction worked for single level
                    idx <- 1:nrow(pvals)
                }

                n_gen <- length(idx)

                if (n_gen == 0) {
                    warning(sprintf("No predictions found for site %s", site))
                    next
                }

                # Get local Genetic Variance
                Vg <- site_stats$Vg[site_stats$Site == site]

                # --- A. CULLIS (Phenotypic Standard) ---
                if ("Cullis" %in% methods && !is.null(total_sed)) {
                    # Extract SED sub-matrix
                    sub_sed <- total_sed[idx, idx]

                    # Calculation: Mean Variance of Difference / 2*Vg
                    # Square SED to get Variances of Difference
                    vd_mat <- sub_sed^2
                    avg_vd <- mean(vd_mat[upper.tri(vd_mat, diag = FALSE)], na.rm = TRUE)

                    h_cullis <- 1 - (avg_vd / (2 * Vg))
                    results_list[[paste0(site, "_Cullis")]] <- data.frame(
                        Site = site, Method = "Cullis", Value = h_cullis, Reliability = h_cullis, Vg = Vg
                    )
                }

                # --- B. OAKEY (Genomic Standard) ---
                if ("Oakey" %in% methods && !is.null(grm_inv) && !is.null(total_vcov)) {
                    # Extract PEV (VCOV) sub-matrix
                    sub_vcov <- total_vcov[idx, idx]

                    # We must match the GRM dimensions to the sub_vcov dimensions.
                    # Genotype column: identify the column that is NOT the site column
                    gen_col_name <- setdiff(names(pvals), site_term)[1] # Heuristic
                    if (is.na(gen_col_name)) gen_col_name <- names(pvals)[2] # Fallback

                    current_gens <- as.character(pvals[[gen_col_name]][idx])

                    if (all(current_gens %in% rownames(grm_inv))) {
                        sub_grm_inv <- grm_inv[current_gens, current_gens]

                        # Optimization: Trace(A %*% B) = Sum(Elementwise(A * t(B)))
                        # Since both are symmetric: Sum(A * B)
                        trace_val <- sum(sub_vcov * sub_grm_inv)

                        # Oakey Formula: 1 - (Trace / (n * Vg_local))
                        # Note: We scale G_inv by (1/Vg) effectively
                        h_oakey <- 1 - (trace_val / (n_gen * Vg))

                        results_list[[paste0(site, "_Oakey")]] <- data.frame(
                            Site = site, Method = "Oakey", Value = h_oakey, Reliability = h_oakey, Vg = Vg
                        )
                    } else {
                        # warning(sprintf("Genotype mismatch in Oakey calc for site %s", site))
                    }
                }

                # --- C. STANDARD (Line Mean) ---
                if ("Standard" %in% methods && !is.null(total_vcov)) {
                    # Mean PEV
                    sub_pev <- diag(total_vcov)[idx]
                    mean_pev <- mean(sub_pev, na.rm = TRUE)

                    h_std <- 1 - (mean_pev / Vg)
                    results_list[[paste0(site, "_Std")]] <- data.frame(
                        Site = site, Method = "Standard", Value = h_std, Reliability = h_std, Vg = Vg
                    )
                }
            },
            error = function(e) {
                warning(sprintf("Prediction failed for site %s: %s", site, e$message))
            }
        )
    }

    close(pb)

    # Compile Results
    final_df <- do.call(rbind, results_list)
    rownames(final_df) <- NULL
    return(final_df)
}
