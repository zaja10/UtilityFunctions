#' Force Model Convergence
#'
#' Iteratively updates an ASReml model until the variance components stabilize.
#'
#' @param model An ASReml model object.
#' @param max_tries Integer. Maximum number of iterations to attempt convergence. Default is 20.
#' @param tolerance Numeric. Percentage change threshold for variance components to consider converged. Default is 1.0.
#'
#' @return A converged ASReml model object, or the last iteration if convergence failed.
#' @export
force_convergence <- function(model, max_tries = 20, tolerance = 1.0) {
    if (!requireNamespace("asreml", quietly = TRUE)) {
        stop("The 'asreml' package is required for this function.")
    }

    try_count <- 1
    converged <- FALSE

    # Check initial status
    if (!is.null(summary(model)$varcomp)) {
        pct_chg <- summary(model)$varcomp[, "%ch"]
        if (!any(pct_chg > tolerance, na.rm = TRUE)) {
            converged <- TRUE
        }
    }

    while (!converged && try_count < max_tries) {
        try_count <- try_count + 1
        # suppress warnings during update loop to avoid clutter
        model <- suppressWarnings(try(update(model), silent = TRUE))

        if (inherits(model, "try-error")) {
            warning("Model update failed during convergence loop.")
            break
        }

        pct_chg <- summary(model)$varcomp[, "%ch"]

        if (!any(pct_chg > tolerance, na.rm = TRUE)) {
            converged <- TRUE
        }
    }

    return(model)
}

#' Analyze Single Variate Single Trial
#'
#' Performs single trial analysis for multiple traits across multiple studies/environments.
#' Automates data cleaning, spatial grid completion, spatial model selection (BIC),
#' and extraction of BLUEs and heritability metrics.
#'
#' @param data A data frame containing the trial data.
#' @param trait Character. Column name of the trait to analyze.
#' @param study_col Character. Column name for the study/environment identifier.
#' @param genotype_col Character. Column name for the genotype identifier.
#' @param row_col Character. Column name for the row coordinate.
#' @param col_col Character. Column name for the column coordinate.
#' @param rep_col Character. Column name for the replicate identifier.
#' @param block_col Character. Column name for the block identifier.
#' @param spatial_models A named list of residual formulas to test. Default tests independent and AR1xAR1 errors.
#' @param extra_random Character. Additional random terms to add to the model formulation (specifically for AR1xAR1 model). Default is NULL.
#' @param plot_residuals Logical. Whether to generate residual plots (currently returned as data/objects, not saved to disk).
#'
#' @return A list containing:
#' \item{blues}{Data frame of Best Linear Unbiased Estimators with columns: genotype, trait, study, predicted.value, std.error, reliability.}
#' \item{model_stats}{Data frame of model statistics including BIC, selected spatial model, and convergence status.}
#' \item{correlations}{Phenotypic correlation matrix of BLUEs (if applicable).}
#'
#' @importFrom dplyr filter mutate select arrange group_by summarize bind_rows everything
#' @importFrom tidyr expand replace_na pivot_wider
#' @importFrom stats as.formula update predict var cor
#' @export
analyze_single_trial <- function(data,
                                 trait,
                                 study_col = "studyName",
                                 genotype_col = "germplasmName",
                                 row_col = "rowNumber",
                                 col_col = "colNumber",
                                 rep_col = "replicate",
                                 block_col = "blockNumber",
                                 spatial_models = list(
                                     "Independent" = "~id(col_f):id(row_f)",
                                     "AR1xAR1" = "~ar1(col_f):ar1(row_f)"
                                 ),
                                 extra_random = NULL,
                                 plot_residuals = FALSE) {
    if (!requireNamespace("asreml", quietly = TRUE)) {
        stop("The 'asreml' package is required for this function.")
    }

    # Ensure required columns exist
    req_cols <- c(study_col, genotype_col, row_col, col_col, rep_col, block_col, trait)
    missing_cols <- req_cols[!req_cols %in% colnames(data)]
    if (length(missing_cols) > 0) {
        stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))
    }

    # Standardize -9 to NA if present in trait
    for (trt in trait) {
        if (any(data[[trt]] == -9, na.rm = TRUE)) {
            data[[trt]][data[[trt]] == -9] <- NA
        }
    }

    results_blues <- list()
    results_stats <- list()

    studies <- unique(as.character(data[[study_col]]))

    for (std in studies) {
        # message(paste("Analyzing Study:", std))

        # 1. Prepare Data for Study
        sub_data <- data[data[[study_col]] == std, ]

        # Create complete spatial grid
        # ASReml spatial analysis requires a complete grid (row x col)
        # We use tidyr::expand to generate all combinations and join back
        # Note: We need to handle potential duplicates if multiple records exist per coord (shouldn't for single trial yield, but good to be safe)

        # Ensure coordinates are integers for range calculation
        r_min <- min(sub_data[[row_col]], na.rm = TRUE)
        r_max <- max(sub_data[[row_col]], na.rm = TRUE)
        c_min <- min(sub_data[[col_col]], na.rm = TRUE)
        c_max <- max(sub_data[[col_col]], na.rm = TRUE)

        # Create grid
        grid_df <- expand.grid(
            rowNumber_temp = seq(r_min, r_max),
            colNumber_temp = seq(c_min, c_max)
        )
        colnames(grid_df) <- c(row_col, col_col)

        # Merge grid with data
        spatial_data <- merge(grid_df, sub_data, by = c(row_col, col_col), all.x = TRUE)

        # Convert design variables to factors
        # We create internal names to avoid conflicts in formula
        spatial_data$row_f <- factor(spatial_data[[row_col]], levels = seq(r_min, r_max))
        spatial_data$col_f <- factor(spatial_data[[col_col]], levels = seq(c_min, c_max))
        spatial_data$geno_f <- factor(spatial_data[[genotype_col]])
        spatial_data$rep_f <- factor(spatial_data[[rep_col]])
        spatial_data$blk_f <- factor(spatial_data[[block_col]])


        spatial_data <- spatial_data[order(spatial_data[[row_col]], spatial_data[[col_col]]), ]

        for (trt in trait) {
            # message(paste("  Trait:", trt))

            # Skip if trait is entirely missing or has 0 variance
            if (all(is.na(spatial_data[[trt]])) || var(spatial_data[[trt]], na.rm = TRUE) == 0) {
                warning(paste("Skipping trait", trt, "in study", std, "due to no variance or all NA."))
                next
            }

            # Prepare storage for model comparison
            model_comp <- data.frame(
                resid_model = names(spatial_models),
                BIC = NA,
                converged = FALSE,
                mean_reliability = 0,
                stringsAsFactors = FALSE
            )

            best_model_obj <- NULL
            best_res_name <- NULL
            min_bic <- Inf

            # 2. Iterate through spatial models (Random Genotype Model for Selection)
            for (mod_name in names(spatial_models)) {
                res_formula_str <- spatial_models[[mod_name]]

                # Base random model: Genotype + Block + Residual
                # Note: Fixed intercept only ("~1")
                # Genotype is random here for variance component estimation & spatial selection

                fixed_form <- as.formula(paste(trt, "~ 1"))
                random_form <- ~ geno_f + blk_f
                resid_form <- as.formula(res_formula_str)

                # Add extra random terms if AR1 model and specified
                if (mod_name == "AR1xAR1" && !is.null(extra_random)) {
                    random_form <- update(random_form, paste("~ . +", extra_random))
                }

                # Fit model
                # Using try() to catch ASReml errors
                # na.action included to match user script logic
                mod_fit <- try(
                    asreml::asreml(
                        fixed = fixed_form,
                        random = random_form,
                        residual = resid_form,
                        data = spatial_data,
                        na.action = asreml::na.method(y = "include", x = "include"),
                        trace = FALSE
                    ),
                    silent = TRUE
                )

                if (!inherits(mod_fit, "try-error")) {
                    # Attempt convergence
                    mod_fit <- try(force_convergence(mod_fit), silent = TRUE)

                    if (!inherits(mod_fit, "try-error")) {
                        # Extract stats
                        summ <- summary(mod_fit)
                        bic_val <- tryCatch(infoCriteria(mod_fit)$BIC, error = function(e) NA)

                        # Calculate reliability
                        # rel = 1 - PEV/Vg
                        # We need predictions for genotypes
                        pvals <- try(predict(mod_fit, classify = "geno_f", ignore = c("(Intercept)"))$pvals, silent = TRUE)

                        mean_rel <- 0
                        if (!inherits(pvals, "try-error")) {
                            pev <- pvals[["std.error"]]^2
                            # Extract Genetic Variance (Vg)
                            # Dependable on component name, usually "geno_f"
                            vg <- summ$varcomp["geno_f", "component"]

                            if (!is.null(vg) && vg > 0) {
                                rels <- 1 - (pev / vg)
                                mean_rel <- mean(rels, na.rm = TRUE)
                            }
                        }

                        # Store results
                        idx <- which(model_comp$resid_model == mod_name)
                        model_comp$BIC[idx] <- bic_val
                        model_comp$converged[idx] <- mod_fit$converge
                        model_comp$mean_reliability[idx] <- mean_rel

                        # Update best model
                        if (!is.na(bic_val) && bic_val < min_bic) {
                            min_bic <- bic_val
                            best_res_name <- mod_name
                            # We don't save the random model object for BLUEs, we refit as Fixed below
                        }
                    }
                }
            } # End spatial model loop

            # 3. Fit Final Model (Fixed Genotype) with Best Spatial Structure
            if (!is.null(best_res_name)) {
                final_res_form <- as.formula(spatial_models[[best_res_name]])

                # Fixed Genotype Model
                fixed_form_final <- as.formula(paste(trt, "~ 1 + geno_f"))
                random_form_final <- ~blk_f

                if (best_res_name == "AR1xAR1" && !is.null(extra_random)) {
                    random_form_final <- update(random_form_final, paste("~ . +", extra_random))
                }

                final_mod <- try(
                    asreml::asreml(
                        fixed = fixed_form_final,
                        random = random_form_final,
                        residual = final_res_form,
                        data = spatial_data,
                        na.action = asreml::na.method(y = "include", x = "include"),
                        trace = FALSE
                    ),
                    silent = TRUE
                )

                if (!inherits(final_mod, "try-error")) {
                    final_mod <- force_convergence(final_mod)

                    # Predict BLUEs
                    # classify="geno_f"
                    blues_pred <- predict(final_mod, classify = "geno_f", pworkspace = 64e6)$pvals

                    if (!inherits(blues_pred, "try-error")) {
                        # Construct result row
                        # Extract reliability from the *Random* model phase (common practice is to report H2/Rel from random model)
                        best_rel_idx <- which(model_comp$resid_model == best_res_name)
                        final_rel <- model_comp$mean_reliability[best_rel_idx]

                        # Format BLUEs dataframe
                        df_blues <- data.frame(
                            genotype = blues_pred$geno_f,
                            trait = trt,
                            study = std,
                            predicted.value = blues_pred$predicted.value,
                            std.error = blues_pred$std.error,
                            status = blues_pred$status,
                            reliability = final_rel,
                            resid_model = best_res_name,
                            stringsAsFactors = FALSE
                        )

                        # Add original genotype names if needed (geno_f is factor)
                        # Actually geno_f levels should match original if data wasn't subset weirdly.
                        # But let's ensure we return valid strings.
                        df_blues$genotype <- as.character(df_blues$genotype)

                        results_blues[[paste(std, trt, sep = "_")]] <- df_blues

                        # Save stats
                        stat_row <- model_comp[model_comp$resid_model == best_res_name, ]
                        stat_row$study <- std
                        stat_row$trait <- trt
                        results_stats[[paste(std, trt, sep = "_")]] <- stat_row
                    }
                }
            } else {
                warning(paste("No successful models for trait", trt, "in study", std))
            }
        } # End trait loop
    } # End study loop

    # Combine results
    all_blues <- do.call(rbind, results_blues)
    all_stats <- do.call(rbind, results_stats)

    # Calculate Correlations if requested and data permits
    pheno_cor <- NULL
    if (!is.null(all_blues) && nrow(all_blues) > 0) {
        # Pivot wider for correlation
        # We want correlation of predicted values between traits (or consistency across studies?)
        # Usually phenotypic correlation matrix is trait x trait (averaged across studies or within study)
        # The user script pivoted by trait+study.

        wide_df <- all_blues %>%
            select(genotype, trait, study, predicted.value) %>%
            pivot_wider(names_from = c(trait, study), values_from = predicted.value)

        # Calculate cor
        if (ncol(wide_df) > 2) { # genotype + at least 2 cols
            mat <- wide_df %>%
                select(-genotype) %>%
                as.matrix()
            pheno_cor <- round(cor(mat, use = "pairwise.complete.obs"), 2)
        }
    }

    return(list(
        blues = all_blues,
        model_stats = all_stats,
        correlations = pheno_cor
    ))
}
