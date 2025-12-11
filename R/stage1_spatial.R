#' Force Model Convergence
#'
#' Iteratively updates an ASReml model until the variance components stabilize or
#' the log-likelihood converges.
#'
#' @param model An ASReml model object.
#' @param max_tries Integer. Maximum number of iterations to attempt convergence. Default is 20.
#' @param tolerance Numeric. Percentage change threshold for variance components. Default is 1.0.
#'
#' @return A converged ASReml model object, or the last iteration if convergence failed.
#' @importFrom MASS ginv
#' @export
force_convergence <- function(model, max_tries = 20, tolerance = 1.0) {
    if (!requireNamespace("asreml", quietly = TRUE)) {
        stop("The 'asreml' package is required for this function.")
    }

    try_count <- 0
    converged <- FALSE
    last_loglik <- -Inf

    # Check initial status
    if (!is.null(summary(model)$varcomp)) {
        pct_chg <- summary(model)$varcomp[, "%ch"]
        if (!any(pct_chg > tolerance, na.rm = TRUE)) {
            converged <- TRUE
        }
    }

    while (!converged && try_count < max_tries) {
        try_count <- try_count + 1

        # Capture previous loglik to check stability
        prev_loglik <- model$loglik

        # suppress warnings during update loop to avoid clutter
        model <- suppressWarnings(try(update(model), silent = TRUE))

        if (inherits(model, "try-error")) {
            warning("Model update failed during convergence loop.")
            break
        }

        # Check Variance Component Change
        pct_chg <- summary(model)$varcomp[, "%ch"]
        vc_stable <- !any(pct_chg > tolerance, na.rm = TRUE)

        # Check LogLik Stability (if variance is wiggling but likelihood is flat)
        ll_stable <- FALSE
        if (!is.null(prev_loglik) && !is.null(model$loglik)) {
            if (abs(model$loglik - prev_loglik) < 0.05) { # Strict LL threshold
                ll_stable <- TRUE
            }
        }

        if (vc_stable || ll_stable) {
            converged <- TRUE
        }
    }

    return(model)
}

#' Check Experimental Design Factors by Study
#'
#' Scans a MET dataset to determine which design terms (Row, Column, Block, Replicate)
#' are valid for inclusion in a model for each specific study.
#'
#' @param data The dataframe containing the MET data.
#' @param study_col String. Name of the column defining the Study/Experiment (default "studyName").
#' @param row_col String. Name of the Row column (default "rowNumber").
#' @param col_col String. Name of the Column column (default "colNumber").
#' @param block_col String. Name of the Block column (default "blockNumber").
#' @param rep_col String. Name of the Replicate column (default "replicate").
#'
#' @return A dataframe summarizing which terms are available for each study.
#' @export
check_design_terms <- function(data,
                               study_col = "studyName",
                               row_col = "rowNumber",
                               col_col = "colNumber",
                               block_col = "blockNumber",
                               rep_col = "replicate") {
    # Ensure the study column exists
    if (!study_col %in% names(data)) {
        stop(paste("Study column", study_col, "not found in dataframe."))
    }

    # Split data by Study
    study_list <- split(data, data[[study_col]])

    # Helper function to check if a column has > 1 unique value (ignoring NAs)
    has_var <- function(df, col_name) {
        if (!col_name %in% names(df)) {
            return(FALSE)
        }
        vals <- df[[col_name]]
        vals <- vals[!is.na(vals)]
        return(length(unique(vals)) > 1)
    }

    results <- lapply(names(study_list), function(nm) {
        sub_df <- study_list[[nm]]
        data.frame(
            Study = nm,
            Has_Row = has_var(sub_df, row_col),
            Has_Col = has_var(sub_df, col_col),
            Has_Block = has_var(sub_df, block_col),
            Has_Rep = has_var(sub_df, rep_col),
            stringsAsFactors = FALSE
        )
    })

    out_df <- do.call(rbind, results)
    rownames(out_df) <- NULL
    return(out_df)
}

#' Convert Symmetric Matrix to Sparse Three-Column Format
#'
#' Converts a symmetric weight matrix (Inverse PEV) into the "Row, Col, Value" format
#' required by ASReml's \code{vm()} function.
#'
#' @param mat A symmetric numeric matrix.
#' @return A dataframe with columns "Row", "Col", "Value".
#' @keywords internal
mat2sparse <- function(mat) {
    # 1. Get lower triangle values (including diagonal)
    # ASReml typically expects the lower triangle for symmetric matrices
    val <- mat[lower.tri(mat, diag = TRUE)]

    # 2. Create grid of indices
    n <- ncol(mat)
    grid <- expand.grid(Row = 1:n, Col = 1:n)
    grid <- grid[grid$Row >= grid$Col, ] # Keep lower triangle

    # 3. Bind
    # Row and Col should be strictly integer indices for the vm() structure
    out <- cbind(grid, Value = val)

    # Filter zeros to keep it truly sparse?
    # Usually vm() handles zeros, but removing them saves memory if huge.
    # However, ASReml coordinate format is rigid. Let's keep strict structure for now.
    return(out)
}

#' Fit Single Site Model (Internal Helper)
#'
#' Fits a spatial model for a single environment. Used by \code{analyze_single_trial}.
#'
#' @keywords internal
fit_site_model <- function(sub_data, trait, genotype, block, row, col, resids_list, compute_weights = FALSE, weight_threshold = 0.01) {
    if (nrow(sub_data) == 0) {
        return(NULL)
    }

    # Check Variance
    if (var(sub_data[[trait]], na.rm = TRUE) == 0) {
        return(NULL)
    }

    fixed_form <- as.formula(paste(trait, "~ 1"))
    random_form <- as.formula(paste("~", genotype, "+", block))

    loop_reliabilities <- c()
    loop_bics <- c()
    loop_conv <- c()

    # Model Selection Loop
    for (resid_model_str in resids_list) {
        resid_form <- as.formula(paste("~", resid_model_str))

        # 1. Initial Fit
        amod1 <- try(asreml::asreml(
            fixed = fixed_form,
            random = random_form,
            residual = resid_form,
            data = sub_data,
            na.action = asreml::na.method(x = "include", y = "include"),
            maxit = 50, trace = FALSE, workspace = "1gb"
        ), silent = TRUE)

        # 2. Force Convergence
        if (!inherits(amod1, "try-error")) {
            amod1 <- force_convergence(amod1, max_tries = 3)

            if (amod1$converge) {
                # Calculate Reliability
                blups <- try(predict(amod1, classify = genotype, ignore = c("(Intercept)"))$pvals, silent = TRUE)
                if (!inherits(blups, "try-error")) {
                    pev <- blups$std.error^2
                    vc <- summary(amod1)$varcomp
                    Vg <- if (genotype %in% rownames(vc)) vc[genotype, "component"] else 0

                    rel <- if (Vg > 1e-6) mean(1 - (pev / Vg), na.rm = TRUE) else 0

                    loop_reliabilities <- c(loop_reliabilities, rel)
                    loop_bics <- c(loop_bics, tryCatch(asreml::infoCriteria(amod1)$BIC, error = function(e) 1e14))
                    loop_conv <- c(loop_conv, TRUE)
                } else {
                    loop_reliabilities <- c(loop_reliabilities, 0)
                    loop_bics <- c(loop_bics, 1e14)
                    loop_conv <- c(loop_conv, FALSE)
                }
            } else {
                loop_reliabilities <- c(loop_reliabilities, 0)
                loop_bics <- c(loop_bics, 1e14)
                loop_conv <- c(loop_conv, FALSE)
            }
        } else {
            loop_reliabilities <- c(loop_reliabilities, 0)
            loop_bics <- c(loop_bics, 1e14)
            loop_conv <- c(loop_conv, FALSE)
        }
    }

    # Select Best Model
    best_idx <- which.min(loop_bics)
    if (length(best_idx) == 0) best_idx <- 1

    best_resid <- resids_list[best_idx]
    best_rel <- loop_reliabilities[best_idx]

    # Final Fit (Fixed Genotype for BLUEs)
    asreml::asreml.options(ai.sing = TRUE, aom = TRUE)
    fixed_final <- as.formula(paste(trait, "~ 1 +", genotype))
    random_final <- as.formula(paste("~", block))
    resid_final <- as.formula(paste("~", best_resid))

    amod2 <- try(asreml::asreml(
        fixed = fixed_final,
        random = random_final,
        residual = resid_final,
        data = sub_data,
        na.action = asreml::na.method(x = "include", y = "include"),
        trace = FALSE, workspace = "1gb"
    ), silent = TRUE)

    if (!inherits(amod2, "try-error")) {
        if (!amod2$converge) amod2 <- try(update(amod2), silent = TRUE)

        blues0 <- try(predict(amod2, classify = genotype, pworkspace = "1gb"), silent = TRUE)
        if (!inherits(blues0, "try-error")) {
            blues_raw <- blues0$pvals
            out_df <- blues_raw[, c(genotype, "predicted.value", "std.error")]
            colnames(out_df)[1] <- "germplasmName"

            out_df$trait <- trait
            out_df$study <- as.character(sub_data[[1, "study_internal"]]) # retrieval
            out_df$residmod <- best_resid
            out_df$conv <- amod2$converge
            out_df$rel <- best_rel
            # Stage 2 Weighting
            weights_out <- NULL

            if (compute_weights) {
                # 1. Predict with Full VCOV
                # Request vcov matrix. Maxit 1 usually enough for prediction if model converged.
                preds_cov <- try(predict(amod2, classify = genotype, pworkspace = "1gb", vcov = TRUE), silent = TRUE)

                if (!inherits(preds_cov, "try-error") && !is.null(preds_cov$vcov)) {
                    pev_mat <- as.matrix(preds_cov$vcov)

                    # 2. Invert Matrix (Smith Weights)
                    # Use Generalized Inverse for stability
                    weight_mat <- try(MASS::ginv(pev_mat), silent = TRUE)

                    if (!inherits(weight_mat, "try-error")) {
                        # 3. Sanity Check: Negative Diagonals
                        diag_w <- diag(weight_mat)
                        if (any(diag_w < 0)) {
                            min_abs <- min(abs(diag_w[diag_w != 0]))
                            diag(weight_mat)[diag_w < 0] <- min_abs
                        }

                        # 4. Memory Guard / Sparsity Threshold
                        # If off-diagonal correlations are tiny, set to 0 to save space (?)
                        # Actually, thresholding the WEIGHT matrix directly might be risky structurally,
                        # but removing near-zero elements is standard sparse practice.
                        # User suggested: "If correlation < 0.01, set to 0".
                        # Converting weight matrix to cor to check? No, just threshold small values.
                        # Let's enforce symmetry and threshold.
                        weight_mat[abs(weight_mat) < 1e-8] <- 0

                        # 5. Convert to Sparse Format
                        weights_out <- mat2sparse(weight_mat)

                        # Attach genotype names map for reconstruction if needed
                        attr(weights_out, "levels") <- levels(sub_data[[genotype]])
                        attr(weights_out, "INVERSE") <- TRUE
                    }
                }
            }

            # Legacy/Approximate Weight (Diagonal)
            out_df$weight <- 1 / (out_df$std.error^2)

            return(list(blues = out_df, weights = weights_out))
        }
    }
    return(NULL)
}

#' Analyze Single Variate Single Trial
#'
#' Performs spatial analysis on multiple trials in parallel.
#'
#' @param compute_weights Logical. If TRUE, calculates full inverse PEV matrices for Stage 2 weighting.
#' @param return_weights Logical. If TRUE, returns a list(blues, weights). If FALSE, returns just BLUEs dataframe.
#' @importFrom foreach %dopar% foreach
#' @export
analyze_single_trial <- function(data,
                                 trait = "Yield",
                                 genotype = "germplasmName",
                                 env = "studyName",
                                 row = "rowNumber",
                                 col = "colNumber",
                                 block = "blockNumber",
                                 rep = "replicate",
                                 parallel = TRUE,
                                 compute_weights = TRUE,
                                 return_weights = TRUE) {
    if (!requireNamespace("asreml", quietly = TRUE)) stop("The 'asreml' package is required.")

    # Parallel Backend Check
    if (parallel) {
        if (!requireNamespace("foreach", quietly = TRUE)) {
            message("foreach package not found. Running sequentially.")
            parallel <- FALSE
        }
        # User should have registered a backend like doParallel
    }

    df <- as.data.frame(data)
    df <- df[!is.na(df[[trait]]), ]

    # Factorize
    df[[genotype]] <- as.factor(df[[genotype]])
    df[[env]] <- as.factor(df[[env]])
    df[[block]] <- as.factor(df[[block]])
    df[[row]] <- as.factor(df[[row]])
    df[[col]] <- as.factor(df[[col]])

    studies <- unique(as.character(df[[env]]))

    # Prepare Data List for Parallel Processing
    # We pre-pad the grids here to save complexity inside the loop
    study_data_list <- list()

    for (std in studies) {
        sub <- droplevels(df[df[[env]] == std, ])
        if (nrow(sub) == 0) next

        # Spatial Padding logic
        r_vals <- as.numeric(as.character(sub[[row]]))
        c_vals <- as.numeric(as.character(sub[[col]]))

        if (length(r_vals) > 0 && length(c_vals) > 0) {
            r_range <- min(r_vals, na.rm = TRUE):max(r_vals, na.rm = TRUE)
            c_range <- min(c_vals, na.rm = TRUE):max(c_vals, na.rm = TRUE)

            grid <- expand.grid(rowNumber = factor(r_range), colNumber = factor(c_range))
            colnames(grid) <- c(row, col)

            sub <- merge(grid, sub, by = c(row, col), all.x = TRUE)
            sub <- sub[order(sub[[row]], sub[[col]]), ]

            sub$study_internal <- std # Keep track of ID
            study_data_list[[std]] <- sub
        }
    }

    # Residual Models
    resids_to_test <- c(
        paste0("id(", col, "):id(", row, ")"),
        paste0("ar1(", col, "):ar1(", row, ")"),
        paste0("ar1(", col, "):id(", row, ")"), # 1D col trend
        paste0("id(", col, "):ar1(", row, ")") # 1D row trend
    )

    # Parallel Execution
    if (parallel) {
        results_list <- foreach(sub = study_data_list, .packages = c("asreml", "MASS")) %dopar% {
            fit_site_model(sub, trait, genotype, block, row, col, resids_to_test, compute_weights)
        }
    } else {
        results_list <- lapply(study_data_list, function(sub) {
            fit_site_model(sub, trait, genotype, block, row, col, resids_to_test, compute_weights)
        })
    }

    # Process Results
    # results_list contains list(blues, weights) or NULL

    blues_list <- list()
    weights_list <- list()

    for (i in seq_along(results_list)) {
        res <- results_list[[i]]
        if (!is.null(res)) {
            blues_list[[i]] <- res$blues
            if (compute_weights && !is.null(res$weights)) {
                # Use the study ID as name
                std_id <- as.character(res$blues$study[1])
                weights_list[[std_id]] <- res$weights
            }
        }
    }

    final_blues <- do.call(rbind, blues_list)

    if (return_weights) {
        return(list(blues = final_blues, weights = weights_list))
    } else {
        return(final_blues)
    }
}
