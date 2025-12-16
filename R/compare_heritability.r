#' Compare Site-Specific Heritability Metrics
#'
#' @description
#' Calculates and compares three standard heritability metrics (Cullis, Oakey, Standard)
#' for each environment. Supports Factor Analytic (FA), Unstructured (US), Diagonal,
#' and Single-Site models.
#'
#' @param model A fitted \code{asreml} object.
#' @param fa_object Optional. An object of class \code{fa_asreml} produced by \code{fa.asreml()}.
#'        If provided, Vg is extracted from the FA reconstruction. If NULL, Vg is extracted
#'        directly from the model summary (supports US, Diag, and scalar components).
#' @param grm Optional. The genomic relationship matrix (inverse not required).
#'        Required only for the "Oakey" method.
#' @param classify Character string. The prediction term (e.g. "Site:Genotype" or "Genotype").
#'        If NULL, attempts to detect automatically from \code{fa_object}.
#' @param id_var Character string. The name of the genotype factor (trace term).
#'        Defaults to "Genotype". Used when parsing variance components if fa_object is NULL.
#' @param methods Character vector. Metrics to calculate: c("Cullis", "Oakey", "Standard").
#'
#' @return A data.frame summarizing H2, Genetic Variance (Vg), and Error Metrics per Site.
#'
#' @details
#' \strong{Metrics:}
#' \itemize{
#'   \item \strong{Standard (Line Mean):} \eqn{1 - (PEV / V_g)}. Often referred to as Line Mean Reliability or Repeatability.
#'   \item \strong{Cullis (Generalized Repeatability):} \eqn{1 - (\bar{v}_{\Delta}^{BLUP} / 2V_g)}. Generalized Heritability/Repeatability for phenotypic trials.
#'   \item \strong{Oakey (Genomic Reliability):} \eqn{1 - tr(C^{22}G^{-1})/n}. Gold standard for genomic prediction accuracy.
#' }
#'
#' \strong{Note on Terminology:}
#' While often loosely called "Heritability", these PEV-based metrics are technically forms of
#' \strong{Generalized Repeatability} or \strong{Reliability} (Squared Accuracy) of the
#' predictions. True heritability is defined as the slope of the regression of True Genetic Values
#' on Phenotypes (or BLUEs on BLUPs), which requires knowing the true effects or fitting dual models.
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
compare_h2 <- function(model, fa_object = NULL, grm = NULL,
                       classify = NULL, id_var = "Genotype",
                       methods = c("Cullis", "Oakey", "Standard")) {
    if (inherits(model, "fa_asreml")) {
        stop("First argument must be an asreml model object. Did you mean to pass 'fa_object'?")
    }

    # --- STEP 1: VARIANCE EXTRACTION (SAME AS ABOVE) ---
    site_stats <- NULL
    sites <- NULL
    site_term <- NULL
    is_multisite <- FALSE
    site_term_name <- NULL # The actual column name

    if (!is.null(fa_object)) {
        if (!inherits(fa_object, "fa_asreml")) stop("fa_object must be of class 'fa_asreml'.")
        if (is.null(classify) && !is.null(fa_object$meta$classify)) classify <- fa_object$meta$classify

        cat("-> Extracting Vg from FA object...\n")
        G_mat <- fa_object$matrices$G
        site_stats <- data.frame(Site = rownames(G_mat), Vg = diag(G_mat))
        sites <- as.character(site_stats$Site)
        is_multisite <- TRUE

        # Heuristic
        term_parts <- unlist(strsplit(classify, ":"))
        site_term_name <- term_parts[1]
    } else {
        if (is.null(classify)) stop("The 'classify' argument is required when fa_object is NULL.")
        cat("-> Extracting Vg from Model Summary...\n")

        vc <- summary(model)$varcomp
        vc_names <- rownames(vc)
        term_parts <- unlist(strsplit(classify, ":"))

        if (length(term_parts) > 1) {
            site_term_name <- term_parts[1]
            id_var <- term_parts[2]
            is_multisite <- TRUE

            # Try to get sites from data
            data_name <- as.character(model$call$data)
            if (exists(data_name)) {
                df <- get(data_name)
                if (site_term_name %in% names(df)) {
                    sites <- as.character(unique(df[[site_term_name]]))
                    sites <- sites[!is.na(sites)]
                }
            }
            if (is.null(sites)) stop("Could not determine site levels (data not found).")

            vg_list <- list()
            # Basic grep for Vg
            for (s in sites) {
                # Strategy: Convert site name to regex safe pattern?
                # Assuming ASReml component naming: "Term!Level"
                # Or "us(Term):Gen!Level:Level"

                # 1. Find components containing the site name
                matches <- grep(s, vc_names, fixed = TRUE, value = TRUE)
                # 2. Must contain genotype term
                matches <- grep(id_var, matches, fixed = TRUE, value = TRUE)

                val <- NA
                if (length(matches) > 0) {
                    # Take the first one that looks like a variance (diagonal)
                    # If US, expect "s:s"
                    # If Diag/At, expect just "s" at end or "at(..., s)"

                    # Preference logic:
                    us_diag <- grep(paste0(s, ":", s), matches, fixed = TRUE)
                    if (length(us_diag) > 0) {
                        val <- vc[matches[us_diag[1]], "component"]
                    } else {
                        val <- vc[matches[1], "component"]
                    } # Fallback
                }

                if (!is.na(val)) vg_list[[s]] <- val
            }
            if (length(vg_list) == 0) stop("Could not extract Vg. Check model terms.")
            site_stats <- data.frame(Site = names(vg_list), Vg = as.numeric(vg_list))
            sites <- as.character(site_stats$Site)
        } else {
            cat("-> Detected Single-Site/Global model.\n")
            is_multisite <- FALSE
            sites <- "Global"
            # Scalar Vg
            idx <- grep(id_var, vc_names)
            if (length(idx) == 0) stop("No Vg found.")
            Vg <- vc[idx[1], "component"]
            site_stats <- data.frame(Site = "Global", Vg = Vg)
        }
    }

    # GRM setup
    grm_inv <- NULL
    if ("Oakey" %in% methods) {
        if (is.null(grm)) {
            warning("GRM missing. Skipping Oakey.")
            methods <- setdiff(methods, "Oakey")
        } else {
            grm_inv <- tryCatch(solve(grm), error = function(e) MASS::ginv(grm))
        }
    }

    results_list <- list()

    if (is_multisite) {
        cat(sprintf("-> Running split predictions across %d sites...\n", length(sites)))
        pb <- utils::txtProgressBar(min = 0, max = length(sites), style = 3)
    }

    for (i in seq_along(sites)) {
        site_label <- sites[i]
        if (is_multisite) utils::setTxtProgressBar(pb, i)

        # Setup prediction args
        pred_args <- list(
            object = model, classify = classify, only = classify,
            sed = "Cullis" %in% methods,
            vcov = any(c("Oakey", "Standard") %in% methods),
            trace = FALSE
        )

        # If multisite, restrict levels
        if (is_multisite) {
            lvl <- list()
            lvl[[site_term_name]] <- site_label
            pred_args$levels <- lvl
        }

        # Execute
        tryCatch(
            {
                preds <- do.call(asreml::predict.asreml, pred_args)

                pvals <- preds$pvals
                total_sed <- preds$sed
                total_vcov <- preds$vcov

                # Identify valid indices
                if (is_multisite) {
                    # Find column matching site_term_name
                    if (!site_term_name %in% names(pvals)) {
                        # Warning: Classify string might not match column names in pvals?
                        # Case: Classify "Site:Genotype", pvals cols "Site", "Genotype"
                        # Just use logic: find column that has value 'site_label'
                        # But safer to filter by rows
                    }
                    # Because we used 'levels=', hopefully only relevant rows returned?
                    # ASReml sometimes returns all rows with NAs.
                    # Filter where Site == site_label
                    idx <- which(pvals[[site_term_name]] == site_label)
                } else {
                    idx <- 1:nrow(pvals)
                }

                n_gen <- length(idx)
                if (n_gen == 0) next

                Vg <- site_stats$Vg[site_stats$Site == site_label]

                # --- METRICS ---
                if ("Cullis" %in% methods && !is.null(total_sed)) {
                    sub_sed <- total_sed[idx, idx]
                    vd_mat <- sub_sed^2
                    avg_vd <- mean(vd_mat[upper.tri(vd_mat, diag = FALSE)], na.rm = TRUE)
                    if (is.nan(avg_vd)) avg_vd <- 0 # Single genotype case?

                    results_list[[paste0(site_label, "_Cullis")]] <- data.frame(
                        Site = site_label, Method = "Cullis", Value = 1 - (avg_vd / (2 * Vg)), Vg = Vg
                    )
                }

                if ("Standard" %in% methods && !is.null(total_vcov)) {
                    sub_pev <- diag(total_vcov)[idx]
                    mean_pev <- mean(sub_pev, na.rm = TRUE)
                    results_list[[paste0(site_label, "_Std")]] <- data.frame(
                        Site = site_label, Method = "Standard", Value = 1 - (mean_pev / Vg), Vg = Vg
                    )
                }

                if ("Oakey" %in% methods && !is.null(grm_inv)) {
                    # Match genotypes
                    # Find genotype col
                    gen_cols <- setdiff(names(pvals), c(site_term_name, "predicted.value", "standard.error", "status"))
                    # Heuristic: usually last term in classify
                    gen_col_name <- tail(unlist(strsplit(classify, ":")), 1)

                    current_gens <- as.character(pvals[[gen_col_name]][idx])
                    if (all(current_gens %in% rownames(grm_inv))) {
                        sub_vcov <- total_vcov[idx, idx]
                        sub_grm <- grm_inv[current_gens, current_gens]
                        trace_val <- sum(sub_vcov * sub_grm)
                        results_list[[paste0(site_label, "_Oakey")]] <- data.frame(
                            Site = site_label, Method = "Oakey", Value = 1 - (trace_val / (n_gen * Vg)), Vg = Vg
                        )
                    }
                }
            },
            error = function(e) warning(paste("Error site", site_label, ":", e$message))
        )
    }

    if (is_multisite) close(pb)

    final_df <- do.call(rbind, results_list)
    rownames(final_df) <- NULL
    class(final_df) <- c("h2_comparison", "data.frame")
    return(final_df)
}
