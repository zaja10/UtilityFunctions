#' Compare Site-Specific Repeatability Metrics
#'
#' @description
#' Calculates the Model-Based Reliability ($R$) for each environment in the trial.
#' Unlike simple repeatability based on variance components, this function uses
#' the **Prediction Error Variance (PEV)** from the BLUPs to provide a rigorous
#' quality metric that accounts for experimental design imbalance and explicit
#' specific variances ($\Psi$) in Factor Analytic or Reduced Rank models.
#'
#' Metric: \deqn{R = 1 - \frac{\text{PEV}_{avg}}{V_g}}
#'
#' @param model A fitted \code{asreml} object.
#' @param fa_object An object of class \code{fa_asreml} produced by \code{fa.asreml()}.
#'
#' @return A data.frame summarizing Repeatability, Genetic Variance (Vg), and Residual Variance (Ve = NA) per Site.
#'
#' @details
#' This function performs a split-prediction strategy to compute the average PEV
#' for each site efficiently. It automatically cleans the \code{classify} term
#' (e.g., transforming \code{rr(Site,2):Gen} to \code{Site:Gen}) to ensure that
#' \code{predict()} captures the Total Genetic Effect (Factor + Specific/Diag).
#'
#' @export
compare_repeatability <- function(model, fa_object) {
    if (!inherits(fa_object, "fa_asreml")) stop("fa_object must be of class 'fa_asreml'.")

    # 1. SETUP: Extract Site-Specific Genetic Variances
    G_mat <- fa_object$matrices$G
    site_stats <- data.frame(Site = rownames(G_mat), Vg = diag(G_mat))
    sites <- as.character(site_stats$Site)

    # 2. CALCULATE REPEATABILITY (Rigorous: 1 - PEV/Vg)
    # -------------------------------------------------------------------------


    # 3. CALCULATE REPEATABILITY
    results_list <- list()

    # 3. CALCULATE REPEATABILITY (Rigorous: 1 - PEV/Vg)
    # The user requested using predict() to get PEV, which accounts for the specific variances
    # and experimental design (unbalance) more accurately than raw Vc.

    cat("-> Calculating Repeatability using PEV (Model-Based Reliability)...\n")

    # Clean Classify Logic (Same as compare_h2)
    # -------------------------------------------------------------------------
    # Auto-detect classify if NULL (Need to add classify arg to function signature first?)
    # We will try to find it in fa_object$meta$classify
    classify <- fa_object$meta$classify
    if (is.null(classify)) stop("Could not determine 'classify' string from fa_object metadata.")

    term_parts <- unlist(strsplit(classify, ":"))
    clean_factors <- sapply(term_parts, function(x) {
        if (grepl("\\(", x)) sub("^[a-z]+\\(([^,]+).*", "\\1", x) else x
    })
    classify_for_predict <- paste(clean_factors, collapse = ":")
    site_term <- clean_factors[1]

    results_list <- list()
    pb <- utils::txtProgressBar(min = 0, max = length(sites), style = 3)

    for (i in seq_along(sites)) {
        site <- sites[i]
        utils::setTxtProgressBar(pb, i)

        Vg <- site_stats$Vg[site_stats$Site == site]
        levels_list <- list()
        levels_list[[site_term]] <- site

        R_val <- NA

        tryCatch(
            {
                # Run prediction for this site
                preds <- asreml::predict.asreml(model,
                    classify = classify_for_predict,
                    levels = levels_list,
                    only = classify_for_predict,
                    vcov = TRUE,
                    trace = FALSE
                )

                # Extract PEV
                # Find valid indices for this site
                pvals <- preds$pvals
                total_vcov <- preds$vcov

                # Identify rows for this site
                site_col_idx <- which(sapply(pvals, function(c) any(c == site)))
                if (length(site_col_idx) > 0) {
                    idx <- which(pvals[[site_col_idx[1]]] == site)
                } else {
                    idx <- 1:nrow(pvals)
                }

                sub_pev <- diag(total_vcov)[idx]
                mean_pev <- mean(sub_pev, na.rm = TRUE)

                # Calculate Reliability (standard formula)
                # This serves as our "Rigorous Repeatability"
                if (Vg > 1e-8) {
                    R_val <- 1 - (mean_pev / Vg)
                } else {
                    R_val <- 0
                }
            },
            error = function(e) {
                # Fallback to NA if predict fails
                warning(sprintf("Prediction failed for site %s: %s", site, e$message))
            }
        )

        results_list[[site]] <- data.frame(
            Site = site,
            Vg = Vg,
            Ve = NA, # PEV based, no single 'Ve'
            Repeatability = R_val
        )
    }

    close(pb)

    final_df <- do.call(rbind, results_list)
    rownames(final_df) <- NULL

    # Sort by Repeatability descending
    final_df <- final_df[order(-final_df$Repeatability), ]

    return(final_df)
}
