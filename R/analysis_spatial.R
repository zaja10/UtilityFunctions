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
#' @importFrom cli cli_alert_warning cli_alert_info
#'
force_convergence <- function(model, max_tries = 20, tolerance = 1.0) {
    if (!check_asreml_availability(action = "warning")) {
        return(model)
    }

    try_count <- 0
    converged <- FALSE

    while (!converged && try_count < max_tries) {
        try_count <- try_count + 1
        prev_loglik <- model$loglik

        # Safe update
        model_upd <- suppressWarnings(try(update(model), silent = TRUE))

        if (inherits(model_upd, "try-error") || is.null(model_upd)) {
            cli::cli_alert_warning("Model update failed. Returning previous iteration.")
            break
        }

        model <- model_upd

        # Check Convergence Flags first (fastest)
        if (isTRUE(model$converge)) {
            converged <- TRUE
        } else {
            # Fallback: Check LogLik stability if 'converge' flag is loose
            if (!is.null(prev_loglik) && !is.null(model$loglik)) {
                if (abs(model$loglik - prev_loglik) < 0.05) converged <- TRUE
            }
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
#' @param study_col String. Name of the column defining the Study/Experiment.
#' @param row_col String. Name of the Row column.
#' @param col_col String. Name of the Column column.
#' @param block_col String. Name of the Block column.
#' @param rep_col String. Name of the Replicate column.
#'
#' @return A dataframe summarizing which terms are available for each study.
check_design_terms <- function(data,
                               study_col = "studyName",
                               row_col = "rowNumber",
                               col_col = "colNumber",
                               block_col = "blockNumber",
                               rep_col = "replicate") {
    if (!study_col %in% names(data)) stop(paste("Study column", study_col, "not found in dataframe."))

    study_list <- split(data, data[[study_col]])

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
#' @importFrom Matrix Matrix
#' @keywords internal
mat2sparse <- function(mat) {
    # Force lower triangle (zero out upper) to ensure strictly lower triangular sparse matrix
    # (assuming symmetric input)
    mat[upper.tri(mat)] <- 0

    # Use Matrix package for efficient storage
    mat_sparse <- Matrix::Matrix(mat, sparse = TRUE)

    # Convert to Triplet format (i, j, x)
    # TsparseMatrix is safe for extraction
    # Force to 'general' (dgCMatrix) first to ensure explicit storage (handles ddiMatrix issue)
    dsT <- as(as(mat_sparse, "dgCMatrix"), "TsparseMatrix")

    # Convert to 1-based indices
    out <- data.frame(
        Row = dsT@i + 1,
        Col = dsT@j + 1,
        Value = dsT@x
    )

    # Ensure Row >= Col (though our zeroing of upper tri should handle this)
    out <- out[out$Row >= out$Col, ]

    return(out)
}

# --- Helper Functions for fit_site_model ---

#' Calculate Reliability from Model
#' @keywords internal
calculate_reliability <- function(model, genotype) {
    blups <- try(predict(model, classify = genotype, ignore = c("(Intercept)"))$pvals, silent = TRUE)
    if (!inherits(blups, "try-error")) {
        pev <- blups$std.error^2
        vc <- summary(model)$varcomp
        Vg <- if (genotype %in% rownames(vc)) vc[genotype, "component"] else 0
        if (Vg > 1e-6) {
            return(mean(1 - (pev / Vg), na.rm = TRUE))
        }
    }
    return(0)
}

#' Select Best Spatial Model
#' @keywords internal
select_spatial_model <- function(sub_data, trait, genotype, block, resids_list) {
    fixed_form <- as.formula(paste(trait, "~ 1"))
    random_form <- as.formula(paste("~", genotype, "+", block))

    best_bic <- 1e14
    best_resid <- resids_list[1]
    best_rel <- 0

    for (resid_model_str in resids_list) {
        resid_form <- as.formula(paste("~", resid_model_str))

        # Initial Fit
        amod1 <- try(asreml::asreml(
            fixed = fixed_form, random = random_form, residual = resid_form,
            data = sub_data, na.action = asreml::na.method(x = "include", y = "include"),
            maxit = 50, trace = FALSE, workspace = "1gb" # Parametrized workspace could be better
        ), silent = TRUE)

        if (!inherits(amod1, "try-error")) {
            amod1 <- force_convergence(amod1, max_tries = 3)

            if (amod1$converge) {
                bic <- tryCatch(asreml::infoCriteria(amod1)$BIC, error = function(e) 1e14)
                if (bic < best_bic) {
                    best_bic <- bic
                    best_resid <- resid_model_str
                    best_rel <- calculate_reliability(amod1, genotype)
                }
            }
        }
    }

    return(list(resid = best_resid, rel = best_rel, bic = best_bic))
}

#' Compute Stage 2 Weights form Model
#' @keywords internal
compute_weights_from_model <- function(model, genotype, threshold = 1e-8) {
    preds_cov <- try(predict(model, classify = genotype, pworkspace = "1gb", vcov = TRUE), silent = TRUE)

    if (!inherits(preds_cov, "try-error") && !is.null(preds_cov$vcov)) {
        pev_mat <- as.matrix(preds_cov$vcov)
        weight_mat <- try(MASS::ginv(pev_mat), silent = TRUE)

        if (!inherits(weight_mat, "try-error")) {
            # Sanity Check
            diag_w <- diag(weight_mat)
            if (any(diag_w < 0)) {
                min_abs <- min(abs(diag_w[diag_w != 0]))
                diag(weight_mat)[diag_w < 0] <- min_abs
            }

            # Thresholding for sparsity
            weight_mat[abs(weight_mat) < threshold] <- 0

            weights_out <- mat2sparse(weight_mat)
            return(weights_out)
        }
    }
    return(NULL)
}

#' Fit Single Site Model (Internal Helper)
#'
#' Fits a spatial model for a single environment. Used by \code{analyze_single_trial}.
#'
#' @keywords internal
fit_site_model <- function(sub_data, trait, genotype, block, row, col, resids_list, compute_weights = FALSE, weight_threshold = 0.01) {
    check_asreml_availability()

    if (nrow(sub_data) == 0) {
        return(NULL)
    }
    if (var(sub_data[[trait]], na.rm = TRUE) == 0) {
        return(NULL)
    }

    # 1. Model Selection
    selection <- select_spatial_model(sub_data, trait, genotype, block, resids_list)
    best_resid <- selection$resid
    best_rel <- selection$rel

    # 2. Final Fit (Fixed Genotype for BLUEs)
    asreml::asreml.options(ai.sing = TRUE, aom = TRUE)
    fixed_final <- as.formula(paste(trait, "~ 1 +", genotype))
    random_final <- as.formula(paste("~", block))
    resid_final <- as.formula(paste("~", best_resid))

    amod2 <- try(asreml::asreml(
        fixed = fixed_final, random = random_final, residual = resid_final,
        data = sub_data, na.action = asreml::na.method(x = "include", y = "include"),
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
            out_df$study <- as.character(sub_data[[1, "study_internal"]])
            out_df$residmod <- best_resid
            out_df$conv <- amod2$converge
            out_df$rel <- best_rel

            # Stage 2 Weighting
            weights_out <- NULL
            if (compute_weights) {
                weights_out <- compute_weights_from_model(amod2, genotype)
                if (!is.null(weights_out)) {
                    attr(weights_out, "levels") <- levels(sub_data[[genotype]])
                    attr(weights_out, "INVERSE") <- TRUE
                }
            }

            # Legacy/Approximate Weight
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
    check_asreml_availability()

    # Parallel Backend Check
    if (parallel) {
        if (!requireNamespace("foreach", quietly = TRUE)) {
            message("foreach package not found. Running sequentially.")
            parallel <- FALSE
        }
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
    study_data_list <- list()

    # Pre-processing (Padding) loop
    for (std in studies) {
        sub <- droplevels(df[df[[env]] == std, ])
        if (nrow(sub) == 0) next

        r_vals <- as.numeric(as.character(sub[[row]]))
        c_vals <- as.numeric(as.character(sub[[col]]))

        if (length(r_vals) > 0 && length(c_vals) > 0) {
            r_range <- min(r_vals, na.rm = TRUE):max(r_vals, na.rm = TRUE)
            c_range <- min(c_vals, na.rm = TRUE):max(c_vals, na.rm = TRUE)

            grid <- expand.grid(rowNumber = factor(r_range), colNumber = factor(c_range))
            colnames(grid) <- c(row, col)

            sub <- merge(grid, sub, by = c(row, col), all.x = TRUE)
            sub <- sub[order(sub[[row]], sub[[col]]), ]
            sub$study_internal <- std
            study_data_list[[std]] <- sub
        }
    }

    resids_to_test <- c(
        paste0("id(", col, "):id(", row, ")"),
        paste0("ar1(", col, "):ar1(", row, ")"),
        paste0("ar1(", col, "):id(", row, ")"),
        paste0("id(", col, "):ar1(", row, ")")
    )

    if (parallel) {
        results_list <- foreach(sub = study_data_list, .packages = c("asreml", "MASS", "Matrix", "UtilityFunctions")) %dopar% {
            fit_site_model(sub, trait, genotype, block, row, col, resids_to_test, compute_weights)
        }
    } else {
        results_list <- lapply(study_data_list, function(sub) {
            fit_site_model(sub, trait, genotype, block, row, col, resids_to_test, compute_weights)
        })
    }

    blues_list <- list()
    weights_list <- list()

    for (i in seq_along(results_list)) {
        res <- results_list[[i]]
        if (!is.null(res)) {
            blues_list[[i]] <- res$blues
            if (compute_weights && !is.null(res$weights)) {
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
