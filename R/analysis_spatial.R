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
#' Pad Trial Layout for Spatial Analysis (MET-Ready)
#'
#' @description Prepares a dataset for spatial analysis by ensuring a complete
#' Row x Column grid. Missing plots are filled with NAs.
#'
#' For Multi-Environment Trials (METs), specify `group_cols` to pad each
#' environment independently.
#'
#' @param data A data.frame containing the trial data.
#' @param row_col Character string. Name of the Row column.
#' @param col_col Character string. Name of the Column column.
#' @param group_cols Character vector. Columns identifying unique trials (e.g. c("Site", "Year")).
#' The function will split data by these columns, pad each subset, and combine them.
#'
#' @return An object of class `padded_trial` (inherits from data.frame).
#' @export
pad_trial_layout <- function(data, row_col = "Row", col_col = "Column", group_cols = NULL) {
    # --- 1. Input Validation ---
    req_cols <- c(row_col, col_col, group_cols)
    if (!all(req_cols %in% names(data))) {
        missing <- req_cols[!req_cols %in% names(data)]
        stop(paste("Columns not found in data:", paste(missing, collapse = ", ")))
    }

    # Ensure Row/Col are numeric
    if (is.factor(data[[row_col]])) data[[row_col]] <- as.numeric(as.character(data[[row_col]]))
    if (is.factor(data[[col_col]])) data[[col_col]] <- as.numeric(as.character(data[[col_col]]))

    # --- 2. Helper Function: Pad Single Trial ---
    pad_single <- function(sub_data) {
        # Determine grid extents for THIS specific trial
        r_range <- range(sub_data[[row_col]], na.rm = TRUE)
        c_range <- range(sub_data[[col_col]], na.rm = TRUE)

        # Create perfect grid
        full_grid <- expand.grid(
            Row = seq(r_range[1], r_range[2]),
            Col = seq(c_range[1], c_range[2])
        )
        names(full_grid) <- c(row_col, col_col)

        # Preserve grouping info (e.g. if this is Site A, put "Site A" in the new rows)
        if (!is.null(group_cols)) {
            for (grp in group_cols) {
                # We take the first value because they are all identical in this subset
                full_grid[[grp]] <- sub_data[[grp]][1]
            }
        }

        # Merge to keep data
        merged <- merge(full_grid, sub_data, by = c(row_col, col_col, group_cols), all.x = TRUE)
        return(merged)
    }

    # --- 3. Execution Logic (Single vs MET) ---

    if (is.null(group_cols)) {
        # Scenario A: Single Trial (No grouping)
        padded_data <- pad_single(data)
        n_added <- nrow(padded_data) - nrow(data)
    } else {
        # Scenario B: MET (Split - Apply - Combine)

        # Create a splitting factor based on all group columns
        split_fac <- interaction(data[group_cols], drop = TRUE)
        data_split <- split(data, split_fac)

        # Apply padding to each list element
        padded_list <- lapply(data_split, pad_single)

        # Combine back into one dataframe
        padded_data <- do.call(rbind, padded_list)
        rownames(padded_data) <- NULL # Clean row names

        n_added <- nrow(padded_data) - nrow(data)
    }

    # --- 4. Diagnostics & Attributes ---
    sparsity <- n_added / nrow(padded_data)

    attr(padded_data, "sparsity") <- sparsity
    attr(padded_data, "n_added") <- n_added
    attr(padded_data, "is_met") <- !is.null(group_cols)
    attr(padded_data, "coords") <- c(row = row_col, col = col_col)

    class(padded_data) <- c("padded_trial", "data.frame")

    return(padded_data)
}

#' Plot Trial Layout
#' @export
plot.padded_trial <- function(x, fill_col = NULL, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plotting.")

    row_name <- attr(x, "coords")["row"]
    col_name <- attr(x, "coords")["col"]

    # Helper to find first non-structural column to detect missingness
    # We skip row, col, and any attributes
    candidate_cols <- names(x)[!names(x) %in% c(row_name, col_name)]
    # Also skip attributes if they act like columns? No, attributes are separate.

    # Usually any trait column will be NA in padded rows
    # Heuristic: Find a column that has NAs where we padded.
    # But we don't know the traits.
    # Let's use the first non-coordinate column as a proxy for "Data Present"
    check_col <- candidate_cols[1]

    x$Status <- ifelse(is.na(x[[check_col]]),
        "Missing (Padded)", "Data Present"
    )

    p <- ggplot2::ggplot(x, ggplot2::aes_string(x = col_name, y = row_name)) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Trial Layout Integrity")

    # If it's a MET, we MUST facet by the grouping columns
    if (attr(x, "is_met")) {
        # Logic to identify group cols (heuristically or passed attributes)
        # Ideally, store group_cols in attributes in the main function.
        # For now, we assume the user handles faceting externally or we create a generic facet
        # Here is a generic approach:
        p <- p + ggplot2::facet_wrap(~., scales = "free")
    }

    if (!is.null(fill_col)) {
        p <- p + ggplot2::geom_tile(ggplot2::aes_string(fill = fill_col), color = "grey90")
    } else {
        p <- p + ggplot2::geom_tile(ggplot2::aes(fill = Status), color = "white") +
            ggplot2::scale_fill_manual(values = c(
                "Data Present" = "steelblue",
                "Missing (Padded)" = "firebrick"
            ))
    }

    return(p)
}
#' Find Missing Row/Column Combinations in MET Data
#'
#' Identifies missing plots required to complete a rectangular grid for spatial analysis.
#' Useful for diagnosing ASREml "residual model implies X" errors.
#'
#' @param data Dataframe containing the trial data.
#' @param experiment String. Column name for the Experiment/Trial ID.
#' @param row String. Column name for Row.
#' @param col String. Column name for Column.
#'
#' @return A dataframe listing the Experiment, Row, and Column of missing plots.
#' @export
find_missing_plots <- function(data, experiment = "Experiment", row = "Row", col = "Column") {
    if (!all(c(experiment, row, col) %in% names(data))) {
        stop("Specified columns not found in dataframe.")
    }

    # Ensure strict types for merging
    data[[experiment]] <- as.factor(data[[experiment]])

    # Iterate through each experiment
    missing_list <- lapply(levels(data[[experiment]]), function(exp_id) {
        # Subset data for this experiment
        sub_df <- data[data[[experiment]] == exp_id, ]

        if (nrow(sub_df) == 0) {
            return(NULL)
        }

        r_vals <- as.numeric(as.character(sub_df[[row]]))
        c_vals <- as.numeric(as.character(sub_df[[col]]))

        min_r <- min(r_vals, na.rm = TRUE)
        max_r <- max(r_vals, na.rm = TRUE)
        min_c <- min(c_vals, na.rm = TRUE)
        max_c <- max(c_vals, na.rm = TRUE)

        expected_grid <- expand.grid(
            Row_Num = seq(min_r, max_r),
            Col_Num = seq(min_c, max_c)
        )

        sub_df$key <- paste(r_vals, c_vals, sep = "_")
        expected_grid$key <- paste(expected_grid$Row_Num, expected_grid$Col_Num, sep = "_")

        missing_keys <- setdiff(expected_grid$key, sub_df$key)

        if (length(missing_keys) > 0) {
            coords <- do.call(rbind, strsplit(missing_keys, "_"))
            return(data.frame(
                Experiment = exp_id,
                Row = as.numeric(coords[, 1]),
                Column = as.numeric(coords[, 2]),
                Status = "MISSING",
                stringsAsFactors = FALSE
            ))
        } else {
            return(NULL)
        }
    })

    out_df <- do.call(rbind, missing_list)

    if (is.null(out_df)) {
        message("Success: All experiments have complete rectangular grids.")
        return(invisible(NULL))
    } else {
        message(paste("Found", nrow(out_df), "missing plots across", length(unique(out_df$Experiment)), "experiments."))
        return(out_df)
    }
}

#' Convert Yield Units
#'
#' Converts yield from Bushels per Acre (bu/ac) to Tonnes per Hectare (t/ha).
#'
#' @param yield Numeric vector of yield values.
#' @param crop String. Crop type to determine bushel weight. Options: "wheat", "soybeans", "corn", "barley", "oats".
#' @param lbs_per_bu Numeric. Custom weight per bushel (overrides crop).
#'
#' @return Numeric vector of yield in t/ha.
#' @export
convert_buac_to_tha <- function(yield, crop = "wheat", lbs_per_bu = NULL) {
    weights <- c(
        "wheat" = 60, "soybeans" = 60, "soy" = 60,
        "corn" = 56, "maize" = 56, "barley" = 48, "oats" = 32
    )

    if (!is.null(lbs_per_bu)) {
        weight <- lbs_per_bu
    } else {
        crop <- tolower(crop)
        if (crop %in% names(weights)) {
            weight <- weights[[crop]]
        } else {
            warning("Crop not recognized. Defaulting to 60 lbs/bu (Wheat). Use 'lbs_per_bu' to specify.")
            weight <- 60
        }
    }

    # Factor = (weight_lbs * 0.453592 kg/lb) / (0.404686 ha/ac * 1000 kg/t)
    factor <- weight * 0.00112085116
    return(yield * factor)
}

#' Flexible Trial Map
#'
#' Plots the spatial layout of a specific trial. Adapts to missing spatial coordinates
#' by falling back to a 1D visualization if necessary.
#'
#' @param data A data.frame.
#' @param trial_val Value. The specific trial ID to filter.
#' @param trial_col String. Column name for Trial ID.
#' @param row String. Column name for Row.
#' @param col String. Column name for Column.
#' @param val String. Variable to fill heatmap (e.g. "Yield").
#'
#' @importFrom grDevices hcl.colors
#' @importFrom graphics image layout par box axis points grid plot
#' @export
plot_trial_map <- function(data, trial_val = NULL, trial_col = "Trial",
                           row = "Row", col = "Column", val = "Yield") {
    # Filter Data
    if (!is.null(trial_val) && trial_col %in% names(data)) {
        plot_data <- data[data[[trial_col]] == trial_val, ]
        title_sub <- paste(trial_val)
    } else {
        plot_data <- data
        title_sub <- "Single Site"
    }

    # Check if Row/Col exist
    has_spatial <- all(c(row, col) %in% names(plot_data))

    # Setup Palette
    vals <- plot_data[[val]]
    if (is.null(vals)) stop(paste("Variable", val, "not found."))

    cols <- hcl.colors(20, "Spectral", rev = TRUE)

    if (has_spatial) {
        # Padded Map
        padded <- pad_trial_layout(plot_data, row, col)

        r_vals <- as.numeric(padded[[row]])
        c_vals <- as.numeric(padded[[col]])

        mat <- matrix(NA, nrow = max(r_vals, na.rm = T), ncol = max(c_vals, na.rm = T))
        for (i in 1:nrow(padded)) {
            if (!is.na(r_vals[i]) & !is.na(c_vals[i])) {
                mat[r_vals[i], c_vals[i]] <- padded[[val]][i]
            }
        }

        # Image Plot
        layout(matrix(1:2, ncol = 2), widths = c(4, 1))
        par(mar = c(4, 4, 3, 1))
        image(1:ncol(mat), 1:nrow(mat), t(mat),
            col = cols,
            main = paste("Trial Map:", title_sub),
            xlab = col, ylab = row
        )
        box()

        # Legend
        par(mar = c(4, 1, 3, 2))
        image(1, seq(min(vals, na.rm = T), max(vals, na.rm = T), length.out = 20),
            t(as.matrix(1:20)),
            axes = F, xlab = "", ylab = "", col = cols
        )
        axis(4, las = 1)
        layout(1)
    } else {
        # 1D Fallback Plan
        warning("Spatial columns missing. Plotting linear sequence.")
        plot(vals, type = "n", main = paste("Linear Plot:", title_sub), ylab = val, xlab = "Index")
        points(vals, pch = 21, bg = cols[cut(vals, 20)], cex = 1.5)
        grid()
    }
}

#' Plot MET Trends
#'
#' Visualizes performance trends over time using a combined Boxplot (distribution)
#' and Linear Regression (rate of change) approach.
#'
#' @param data A data.frame.
#' @param x String. Categorical variable for X-axis (e.g. "Year"). Must be convertible to numeric for trend analysis.
#' @param y String. Response variable (e.g. "Yield").
#' @param main String. Plot title prefix.
#' @param ... Additional arguments passed to `plot` or `boxplot`.
#'
#' @importFrom stats lm coef
#' @importFrom graphics boxplot abline legend
#' @export
plot_met_trend <- function(data, x = "Year", y = "Yield", main = "Yield Trend", ...) {
    if (!all(c(x, y) %in% names(data))) stop("Variables not found in data.")

    raw_x <- data[[x]]
    if (is.numeric(raw_x)) {
        x_num <- raw_x
    } else {
        x_num <- suppressWarnings(as.numeric(as.character(raw_x)))
        if (all(is.na(x_num))) stop(paste("Column", x, "must be numeric (e.g. Year) for trend analysis."))
    }
    x_fac <- as.factor(raw_x)

    # Fits
    tmp_df <- data.frame(Y = data[[y]], X = x_num)
    model <- lm(Y ~ X, data = tmp_df)
    slope <- coef(model)["X"]
    r_sq <- summary(model)$r.squared

    # Visuals
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    layout(matrix(1:2, nrow = 2))

    # Boxplot
    par(mar = c(3, 5, 3, 2))
    boxplot(data[[y]] ~ x_fac,
        col = "lightgreen", border = "darkgreen",
        main = paste(main, "(Distribution)"),
        ylab = y, las = 1, ...
    )

    # Line Plot
    par(mar = c(5, 5, 2, 2))
    means <- tapply(tmp_df$Y, tmp_df$X, mean, na.rm = TRUE)
    years_plot <- as.numeric(names(means))

    plot(years_plot, means,
        type = "b", pch = 19, lwd = 2, col = "blue",
        main = paste(main, "(Rate of Change)"),
        xlab = x, ylab = paste("Mean", y),
        xaxt = "n", las = 1
    )
    axis(1, at = years_plot, labels = years_plot)
    grid()
    abline(model, col = "red", lwd = 2)

    legend("topleft",
        legend = c(paste("Rate:", round(slope, 4), "unit/yr"), paste("R2:", round(r_sq, 3))),
        bty = "n", cex = 0.9, text.col = "black"
    )

    layout(1)
}
