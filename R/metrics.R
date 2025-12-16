#' Calculate Heritability from a List of ASReml Models
#'
#' This function iterates through a list of fitted asreml model objects and
#' calculates broad-sense heritability (H2) based on the Cullis method
#' (1 - mean(PEV) / Vg), using the variance-covariance matrix of predictions.
#' Note: This function uses `predict(..., vcov = TRUE)` and may fail with a
#' "long vector not supported" error if the genetic factor has many levels,
#' leading to a very large prediction variance-covariance matrix. A safer
#' alternative uses standard errors from predict.
#'
#' @param model_list A named list containing fitted `asreml` objects.
#' @param id A character string specifying the name of the factor representing
#'   the genetic random effect (e.g., "Genotype", "Variety"). Defaults to "Genotype".
#'
#' @return A named list or vector containing the calculated heritability (H2)
#'   for each valid and converged model in the input list. Returns NA for
#'   models that are invalid, did not converge, or encountered errors during
#'   calculation.
#' Genetic Metrics
#'
#' Functions to calculate heritability and reliability.
#'
#' @name genetic_metrics
NULL

#' @describeIn genetic_metrics Calculate Cullis H2
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' fm <- list() # Placeholder for a list of asreml models
#' # h2_values <- calc_h2_cullis(fm, id = "Variety")
#' }
calc_h2_cullis <- function(model_list, id = "Genotype") {
    # Check if input is a list
    if (!is.list(model_list)) {
        stop("Input 'model_list' must be a list.")
    }

    # Initialize an empty list to store the results
    h2_results <- list()

    # Loop through each model in the list 'model_list'
    for (model_name in names(model_list)) {
        message("Processing model: ", model_name)
        asr <- model_list[[model_name]]

        # --- Basic Checks ---
        # Check if the element is actually an asreml object and converged
        if (is.null(asr) || !inherits(asr, "asreml")) {
            message("Skipping ", model_name, " - Not a valid asreml object.")
            h2_results[[model_name]] <- NA # Store NA for invalid objects
            next # Go to the next iteration
        }
        if (!asr$converge) {
            message("Skipping ", model_name, " - Model did not converge.")
            h2_results[[model_name]] <- NA # Store NA for non-converged models
            next # Go to the next iteration
        }

        # --- Calculation within tryCatch ---
        h2_calc <- tryCatch(
            {
                # Extract response variable name (safer way)
                response <- NA
                # Try different ways to access the response variable name based on asreml versions/call structure
                if (!is.null(asr$call$fixed) && length(as.list(asr$call$fixed)) >= 2) {
                    response <- as.character(as.list(asr$call$fixed)[[2]])
                } else if (!is.null(asr$call[[2]]) && length(asr$call[[2]]) >= 2) {
                    # Fallback for older style or direct formula input potentially
                    response <- as.character(asr$call[[2]][[2]])
                }
                if (is.na(response)) {
                    warning(paste("Could not determine response variable name for", model_name))
                    response <- "h2" # Default name if extraction fails
                }

                # Run predict with vcov=TRUE (this is where the "long vector" error might occur)
                # Note: 'only=id' is removed as it wasn't in the final user script and can sometimes cause issues
                mypred <- predict(asr, classify = id, maxit = 1, vcov = TRUE) # Added average, sed

                # Get the prediction vcov matrix
                if (!"vcov" %in% names(mypred)) {
                    stop("vcov component missing from predict output.")
                }
                my.vcov <- mypred$vcov

                # Find the genetic variance component index
                summary_asr <- summary(asr, coef = FALSE) # Get summary efficiently
                var_comp_table <- summary_asr$varcomp

                # Improved pattern to match id!id or vm(id, ...) more reliably
                pattern <- paste0("^(", id, ")!|^(vm\\(", id, "[ ,])") # Matches id! or vm(id, or vm(id<space>
                which.vc <- grep(pattern, rownames(var_comp_table))

                # Check if exactly one component was found
                if (length(which.vc) != 1) {
                    # Try a simpler grep if the first failed, maybe simpler naming used
                    which.vc <- grep(paste0("^", id, "$"), rownames(var_comp_table)) # Match exact term name
                    if (length(which.vc) != 1) {
                        stop(paste0(
                            "Could not uniquely identify variance component for '", id,
                            "'. Found: ", length(which.vc), " matches. Check component names: ",
                            paste(rownames(var_comp_table), collapse = ", ")
                        ))
                    }
                }
                Vg <- var_comp_table[which.vc, "component"]

                # Check if Vg is valid
                if (is.na(Vg) || Vg <= 1e-8) {
                    stop(paste("Genetic variance component for '", id, "' is zero, negative, or NA (", Vg, ")."))
                }

                # Calculate h2 using diag() - this might fail for large matrices
                h2_value <- 1 - sum(diag(my.vcov)) / (Vg * nrow(my.vcov))

                # Ensure h2 is within bounds [0, 1]
                h2_value <- max(0, min(1, h2_value))

                # Name the result
                names(h2_value) <- response

                # Return the calculated value
                h2_value
            },
            error = function(e) {
                # If any error occurred in the try block
                message("Error processing ", model_name, ": ", e$message)
                return(NA) # Return NA if calculation fails
            }
        )

        # Store the result (either the h2 value or NA)
        h2_results[[model_name]] <- h2_calc
    }

    # --- Combine results ---
    # Unlist to return a named vector
    h2_vector <- unlist(h2_results)

    # Print summary of NAs at the end
    message("------------------------------------")
    message("Heritability calculation summary:")
    message("Number of models processed: ", length(model_list))
    message("Number of successful calculations: ", sum(!is.na(h2_vector)))
    message("Number of failed calculations (NA): ", sum(is.na(h2_vector)))
    message("------------------------------------")

    return(h2_vector)
}

#' @export
calculate_h2_from_list <- function(...) {
    warning("This function is deprecated. Please use 'calc_h2_cullis()' instead.")
    calc_h2_cullis(...)
}

# --- How to use the function ---
# Assuming 'fm' is your list of models:
# h2_values <- calc_h2_cullis(fm, id = "Genotype")
# print(h2_values)
#' Calculate Realized Genetic Gain
#'
#' Estimates the rate of genetic gain by regressing breeding values (OP) against year of release.
#' Excludes long-term checks from the regression to avoid bias.
#'
#' @param data A dataframe containing genotype performance.
#' @param year_col String. Column for Year of Release (numeric).
#' @param value_col String. Column for Breeding Value / OP.
#' @param check_list Character vector. Names of checks to exclude from slope calculation.
#'
#' @return A list containing the regression model, gain percentage, and summary stats.
#' @export
calculate_realized_gain <- function(data, year_col = "Year", value_col = "OP", check_list = NULL) {
    df <- data[!is.na(data[[year_col]]) & !is.na(data[[value_col]]), ]

    # Filter breeding population (exclude checks)
    breed_pop <- df
    if (!is.null(check_list)) {
        breed_pop <- df[!df$Genotype %in% check_list, ]
    }

    if (nrow(breed_pop) < 5) stop("Not enough data points for regression.")

    form <- as.formula(paste(value_col, "~", year_col))
    mod <- lm(form, data = breed_pop)

    slope <- coef(mod)[2]
    intercept <- coef(mod)[1]

    # Calculate gain relative to 1990 baseline (or min year)
    base_year <- min(breed_pop[[year_col]])
    base_val <- predict(mod, newdata = data.frame(setNames(list(base_year), year_col)))

    pct_gain <- (slope / base_val) * 100

    return(list(
        model = mod,
        slope = slope,
        r_squared = summary(mod)$r.squared,
        pct_gain_per_year = pct_gain,
        base_year = base_year
    ))
}

#' Plot Genetic Trend (Era Plot)
#'
#' Visualizes the genetic trend over time, highlighting checks and the regression line.
#'
#' @param data A dataframe containing genotype performance.
#' @param year_col String. Column for Year of Release.
#' @param value_col String. Column for Breeding Value / OP.
#' @param check_list Character vector. Names of checks to highlight.
#' @param ... Additional arguments to plot.
#' @export
plot_genetic_trend <- function(data, year_col = "Year", value_col = "OP", check_list = NULL, ...) {
    df <- data[!is.na(data[[year_col]]) & !is.na(data[[value_col]]), ]

    # Plot Base
    plot(df[[year_col]], df[[value_col]],
        type = "n",
        xlab = "Year of Release", ylab = "Genetic Value (BLUP)",
        main = "Genetic Trend (Era Plot)", ...
    )
    grid()

    # Add Regression (Breeding Pop only)
    breed_pop <- df
    if (!is.null(check_list)) {
        breed_pop <- df[!df$Genotype %in% check_list, ]
    }

    if (nrow(breed_pop) > 5) {
        mod <- lm(as.formula(paste(value_col, "~", year_col)), data = breed_pop)
        abline(mod, col = "blue", lwd = 2)

        # Add Gain Text
        slope <- coef(mod)[2]
        base_val <- predict(mod, newdata = data.frame(setNames(list(min(breed_pop[[year_col]])), year_col)))
        pct <- round((slope / base_val) * 100, 2)
        legend("topleft", legend = paste("Gain:", pct, "% / yr"), bty = "n", text.col = "blue")
    }

    # Points
    # Breeding Pop
    points(breed_pop[[year_col]], breed_pop[[value_col]], pch = 19, col = alpha("grey", 0.6))

    # Checks
    if (!is.null(check_list)) {
        checks <- df[df$Genotype %in% check_list, ]
        if (nrow(checks) > 0) {
            points(checks[[year_col]], checks[[value_col]], pch = 17, col = "red", cex = 1.2)
            text(checks[[year_col]], checks[[value_col]], labels = checks$Genotype, pos = 3, cex = 0.7, col = "red")
        }
    }
}

# Helper for manual alpha if scales not loaded
alpha <- function(col, alpha_val) {
    rgb_vals <- col2rgb(col)
    rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], max = 255, alpha = alpha_val * 255)
}
#' Compare Site-Specific Heritability Metrics
#'
#' @description
#' Calculates and compares three standard heritability metrics (Cullis, Oakey, Standard)
#' for each environment. Supports Factor Analytic (FA), Unstructured (US), Diagonal,
#' and Single-Site models.
#'
#' @param model A fitted \code{asreml} object.
#' @param fa_object Optional. An object of class \code{fa_model} produced by \code{fit_fa_model()}.
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

#' @importFrom MASS ginv
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @describeIn genetic_metrics Compare H2 across models
#' @export
compare_h2 <- function(model, fa_object = NULL, grm = NULL,
                       classify = NULL, id_var = "Genotype",
                       methods = c("Cullis", "Oakey", "Standard")) {
    if (inherits(model, "fa_model")) {
        stop("First argument must be an asreml model object. Did you mean to pass 'fa_object'?")
    }

    # --- STEP 1: VARIANCE EXTRACTION (SAME AS ABOVE) ---
    site_stats <- NULL
    sites <- NULL
    site_term <- NULL
    is_multisite <- FALSE
    site_term_name <- NULL # The actual column name

    if (!is.null(fa_object)) {
        if (!inherits(fa_object, "fa_model")) stop("fa_object must be of class 'fa_model'.")
        if (is.null(classify) && !is.null(fa_object$meta$classify)) classify <- fa_object$meta$classify

        message("-> Extracting Vg from FA object...")
        G_mat <- fa_object$matrices$G
        site_stats <- data.frame(Site = rownames(G_mat), Vg = diag(G_mat))
        sites <- as.character(site_stats$Site)
        is_multisite <- TRUE

        # Heuristic
        term_parts <- unlist(strsplit(classify, ":"))
        site_term_name <- term_parts[1]
    } else {
        if (is.null(classify)) stop("The 'classify' argument is required when fa_object is NULL.")
        message("-> Extracting Vg from Model Summary...")

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
            message("-> Detected Single-Site/Global model.")
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
        message(sprintf("-> Running split predictions across %d sites...", length(sites)))
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
                # FIX: do.call fails with asreml::predict because it passes the object by value.
                # We need to construct the call using the symbol 'model' and evaluate it.
                call_args <- list(quote(asreml::predict.asreml), object = quote(model))

                # Append other arguments
                call_args <- c(call_args, pred_args[names(pred_args) != "object"])

                # Evaluate
                preds <- eval(as.call(call_args))

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
#' Calculate Multi-Trait Stability Index
#'
#' Computes a composite selection index across multiple traits, incorporating
#' both performance (OP) and stability (RMSD) for each trait.
#'
#' @param fa_objects_list A named list of \code{fa_model} objects (one per trait).
#' @param weights Named numeric vector of economic weights (w) for each trait.
#' @param penalties Named numeric vector of stability penalties (lambda) for each trait.
#'
#' @return A dataframe with standardized scores per trait and the final Index.
#' @export
calculate_mt_index <- function(fa_objects_list, weights, penalties) {
    if (is.null(names(fa_objects_list))) stop("fa_objects_list must be named (Trait names).")

    traits <- names(fa_objects_list)
    common_genos <- NULL

    # 1. Extract and Standardize per trait
    res_list <- list()

    for (tr in traits) {
        obj <- fa_objects_list[[tr]]
        if (is.null(obj$fast)) {
            # Try to calculate if missing
            obj <- tryCatch(calculate_op(obj), error = function(e) NULL)
        }

        if (is.null(obj) || is.null(obj$selection)) {
            warning(paste("Skipping trait", tr, "- no FAST indices available."))
            next
        }

        df <- obj$selection
        if (is.null(common_genos)) {
            common_genos <- df$Genotype
        } else {
            common_genos <- intersect(common_genos, df$Genotype)
        }

        # Standardize (Z-score)
        df$OP_Z <- (df$OP - mean(df$OP, na.rm = TRUE)) / sd(df$OP, na.rm = TRUE)
        df$RMSD_Z <- (df$RMSD - mean(df$RMSD, na.rm = TRUE)) / sd(df$RMSD, na.rm = TRUE)

        res_list[[tr]] <- df
    }

    if (length(common_genos) == 0) stop("No common genotypes found across traits.")

    # 2. Combine
    final_df <- data.frame(Genotype = common_genos)
    final_df$Index <- 0

    for (tr in traits) {
        w <- if (tr %in% names(weights)) weights[[tr]] else 1.0
        lam <- if (tr %in% names(penalties)) penalties[[tr]] else 0.0

        sub <- res_list[[tr]]
        match_idx <- match(final_df$Genotype, sub$Genotype)

        op_z <- sub$OP_Z[match_idx]
        rmsd_z <- sub$RMSD_Z[match_idx]

        # Component Index: w * (OP - lambda * RMSD)
        # Note: RMSD is "bad", so -RMSD is "good" stability
        # Formula: w * (OP_Z - lambda * RMSD_Z)

        comp_val <- w * (op_z - (lam * rmsd_z))

        final_df[[paste0(tr, "_OP_Z")]] <- round(op_z, 2)
        final_df[[paste0(tr, "_RMSD_Z")]] <- round(rmsd_z, 2)
        final_df[[paste0(tr, "_Comp")]] <- comp_val

        final_df$Index <- final_df$Index + comp_val
    }

    final_df <- final_df[order(final_df$Index, decreasing = TRUE), ]
    final_df$Rank <- 1:nrow(final_df)

    return(final_df)
}

#' Predict Genetic Gain
#'
#' Estimates expected genetic gain (R) using the Breeder's Equation: \eqn{R = i r \sigma_g}.
#' Useful for evaluating the ROI of the breeding cycle.
#'
#' @param object An object of class \code{fa_model}.
#' @param selection_percent Numeric. Top percentage to select (e.g., 5 for 5\%).
#' @return A dataframe of predicted gains per environment/trait.
#' @export
predict_gain <- function(object, selection_percent = 5) {
    if (selection_percent <= 0 || selection_percent >= 100) stop("Selection % must be between 0 and 100")

    # 1. Selection Intensity (i)
    # Standardized selection differential for normal distribution
    p <- selection_percent / 100
    i_val <- dnorm(qnorm(1 - p)) / p

    # 2. Genetic Variance (Sigma_g)
    # Diagonal of the G matrix (Genetic Variance per site)
    G <- object$matrices$G
    sigma_g <- sqrt(diag(G))

    # 3. Prediction
    # R = i * sigma_g (assuming r=1 for potential gain)
    pred_gain <- i_val * sigma_g

    res <- data.frame(
        Group = names(pred_gain),
        Sigma_G = round(sigma_g, 3),
        Intensity_i = round(i_val, 2),
        Potential_Gain = round(pred_gain, 3),
        row.names = NULL
    )

    return(res)
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

    cli::cli_h1("FA Model Report")
    cli::cli_alert_info("Dimension: {nrow(object$matrices$G)} Sites x {nrow(object$scores$rotated)} Genotypes")
    cli::cli_alert_info("Factors: {k} | Rotation: SVD")

    # NEW: Annotation Report
    if (!is.null(object$meta$annotation)) {
        cli::cli_text("\n")
        cli::cli_h2("Site Annotations")

        tags <- object$meta$annotation$Tag
        counts <- table(tags)

        # Print distinct count
        cli::cli_ul(paste0(names(counts), ": ", counts))

        # If "Excluded" or similar negative tags exist, highlight them
        bad_tags <- c("Frost", "Drought", "Exclude", "Error")
        found_bad <- intersect(names(counts), bad_tags)

        if (length(found_bad) > 0) {
            cli::cli_alert_warning("Note: {sum(counts[found_bad])} sites flagged with issues ({paste(found_bad, collapse=', ')}).")
        }
    }

    # 1. Variance Accounted For
    vaf <- object$var_comp$vaf
    # Check if Total_VAF exists (it should from fit_fa_model)
    if (!is.null(vaf$Total_VAF)) {
        mean_vaf <- mean(vaf$Total_VAF, na.rm = TRUE)
        cli::cli_alert_success("Mean Variance Accounted For (VAF): {round(mean_vaf, 1)}%")

        # Identify poor sites
        poor_sites <- vaf[vaf$Total_VAF < 70, 1]
        if (length(poor_sites) > 0) {
            cli::cli_alert_warning("Low Quality Sites (VAF < 70%): {.val {head(poor_sites, 3)}} {if(length(poor_sites)>3) '...' else ''}")
        }
    }

    # 2. Correlations
    cor_mat <- object$matrices$Cor
    upper_tri <- cor_mat[upper.tri(cor_mat)]
    avg_cor <- mean(upper_tri, na.rm = TRUE)

    cat("\n")
    cli::cli_text("Genetic Correlation Summary:")
    cli::cli_ul(c(
        "Average: {round(avg_cor, 2)}",
        "Range: {round(min(upper_tri), 2)} to {round(max(upper_tri), 2)}"
    ))

    # 3. FAST Indices Preview
    if (!is.null(object$fast)) {
        cat("\n")
        cli::cli_text("Top 3 Genotypes (Overall Performance):")
        print(head(object$fast[, c("Genotype", "OP", "RMSD")], 3), row.names = FALSE)
    }

    invisible(object)
}

#' Calculate Interaction Classes (Generalized)
#'
#' Classifies environments/genotypes based on the sign of their Factor Loadings/Scores.
#' Supports any number of factors (k).
#'
#' @param object An object of class \code{fa_model}.
#' @return A list with site classes and genotype classes.
#' @export
calculate_i_classes <- function(object) {
    lam <- object$loadings$rotated
    sco <- object$scores$rotated
    k <- object$meta$k

    # Helper to convert sign vector to string (e.g. "+-+")
    get_sign_class <- function(mat) {
        apply(mat, 1, function(x) paste0(ifelse(x > 0, "+", "-"), collapse = ""))
    }

    site_classes <- get_sign_class(lam)
    geno_classes <- if (!is.null(sco)) get_sign_class(sco) else NULL

    # Calculate Mean Performance of each class
    # (Optional: You could calc mean yield per class if needed)

    return(list(
        site_classes = data.frame(Group = names(site_classes), Class = site_classes),
        geno_classes = if (!is.null(geno_classes)) data.frame(Genotype = names(geno_classes), Class = geno_classes) else NULL
    ))
}

#' Calculate Reliability (Deprecated)
#'
#' @export
calculate_reliability <- function(...) {
    cli::cli_alert_warning("This function is deprecated/unreliable and will be removed.")
    return(NULL)
}
