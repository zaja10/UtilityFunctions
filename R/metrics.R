#' Genetic Metrics Suite
#'
#' A unified collection of functions for calculating heritability, reliabilities,
#' and genetic gain from ASReml models and phenotypic data.
#'
#' @name genetic_metrics
NULL

#' Calculate Heritability and Reliability
#'
#' Computes broad-sense heritability metrics from ASReml models.
#' Supports both "Cullis" (Generalized Heritability) and "Standard" (Line Mean Reliability) methods.
#'
#' @param model A fitted \code{asreml} object or a list of such objects.
#' @param method Character string. "Cullis" (default) or "Standard".
#' @param id_var Character string. Name of the genetic factor (default "Genotype").
#' @param vcov Logical. If TRUE, forces usage of the variance-covariance matrix (required for Standard).
#'   For large number of genotypes (>2000), defaults to FALSE to avoid memory exhaustion unless explicitly set.
#'
#' @details
#' \strong{Methods:}
#' \itemize{
#'   \item \strong{Cullis:} \eqn{1 - \frac{\bar{v}_{\Delta}^{BLUP}}{2 \sigma^2_g}}. Based on the average standard error of difference (SED). Recommended for unbalanced trials.
#'   \item \strong{Standard:} \eqn{1 - \frac{PEV}{\sigma^2_g}}. Based on the prediction error variance.
#' }
#'
#' @return A named vector of heritability values.
#' @importFrom stats predict density
#' @export
calculate_heritability <- function(model, method = "Cullis", id_var = "Genotype", vcov = NULL) {
    # Recursion for lists
    if (is.list(model) && !inherits(model, "asreml")) {
        return(sapply(model, function(m) calculate_heritability(m, method, id_var, vcov)))
    }

    if (!inherits(model, "asreml")) {
        return(NA)
    }
    if (!model$converge) {
        warning("Model did not converge.")
        return(NA)
    }

    # 1. Extract Genetic Variance (Vg)
    summ <- summary(model)$varcomp
    # Regex to find the variance component for the ID
    # Matches "Genotype", "Genotype!Genotype", "vm(Genotype", "ide(Genotype"
    pat <- paste0("^(", id_var, ")(!|$)|vm\\(", id_var, "|ide\\(", id_var)

    idx <- grep(pat, rownames(summ))
    if (length(idx) == 0) {
        warning(paste("No variance component found for", id_var))
        return(NA)
    }

    # Preference: Take the first match that isn't bound at boundary (optional refinement)
    Vg <- summ[idx[1], "component"]

    if (Vg <= 1e-6) {
        warning("Genetic variance is effectively zero.")
        return(0)
    }

    # 2. Logic Control
    n_gen <- length(levels(model$model[[id_var]]))

    # Auto-switch to lighter method if N is large
    if ((is.null(vcov) && n_gen > 2000) || method == "Cullis") {
        use_vcov <- FALSE
    } else {
        use_vcov <- TRUE
    }

    if (method == "Standard") use_vcov <- TRUE # Mandatory for Standard

    # 3. Prediction
    # We suppress warnings about "aliasing" which are common in prediction
    # Predict only the genetic term.

    preds <- tryCatch(
        {
            suppressWarnings(
                predict(model, classify = id_var, vcov = use_vcov, sed = !use_vcov)
            )
        },
        error = function(e) {
            warning("Prediction failed: ", e$message)
            return(NULL)
        }
    )

    if (is.null(preds)) {
        return(NA)
    }

    # 4. Calculation
    H2 <- NA

    if (method == "Cullis") {
        if (!is.null(preds$sed)) {
            # Use SED square
            av_sed_sq <- mean(preds$sed^2, na.rm = TRUE) # SED matrix is actually provided as dense or we calculate from pvals?
            # ASReml predict returns $sed as a matrix usually
            # But wait, predict(sed=TRUE) returns the matrix in $sed

            # Note: For large N, we might only have Average SED if not full matrix?
            # asreml returns full matrix for sed=TRUE usually.

            if (is.matrix(preds$sed)) {
                vd <- preds$sed^2
                av_sed_diff <- mean(vd[upper.tri(vd)], na.rm = TRUE)
                H2 <- 1 - (av_sed_diff / (2 * Vg))
            } else {
                # Fallback if SED is just a value? Unlikely in asreml-r
                warning("SED matrix not returned.")
            }
        }
    } else if (method == "Standard") {
        if (!is.null(preds$vcov)) {
            pev <- diag(preds$vcov)
            mean_pev <- mean(pev, na.rm = TRUE)
            H2 <- 1 - (mean_pev / Vg)
        }
    }

    return(max(0, min(1, H2))) # Clamp to 0-1
}


#' Calculate Genetic Gain
#'
#' Estimates genetic gain using either "Realized" (Historical Regression) or
#' "Predicted" (Breeder's Equation) methods.
#'
#' @param input Either a dataframe (for "Realized") or an `fa_model`/`asreml` object (for "Predicted").
#' @param ... arguments handling specifics.
#' \describe{
#'   \item{Realized}{Requires `year_col` (Year) and `value_col` (Yield/BV).}
#'   \item{Predicted}{Requires `selection_percent` (default 5).}
#' }
#'
#' @return
#' For Realized: A list with slope, R2, and pct_gain.
#' For Predicted: A dataframe of potential gain per trait/site.
#' @importFrom stats lm coef dnorm qnorm as.formula
#' @export
calculate_genetic_gain <- function(input, ...) {
    args <- list(...)

    # MODE 1: PREDICTED (Model Object)
    if (inherits(input, c("asreml", "fa_model"))) {
        sel_pct <- if ("selection_percent" %in% names(args)) args$selection_percent else 5
        if (sel_pct <= 0 || sel_pct >= 100) stop("Selection % must be 0-100")

        # Intensity
        p <- sel_pct / 100
        i_val <- dnorm(qnorm(1 - p)) / p

        # Sigma G
        if (inherits(input, "fa_model")) {
            Vg <- diag(input$matrices$G)
            lbl <- names(Vg)
        } else {
            # Simple scalar model assumption
            summ <- summary(input)$varcomp
            Vg <- summ[grep("Genotype", rownames(summ))[1], "component"]
            lbl <- "Global"
        }

        param_gain <- i_val * sqrt(Vg) # Assuming r=1

        return(data.frame(
            Group = lbl,
            Sigma_G = sqrt(Vg),
            Intensity = i_val,
            Predicted_Gain = param_gain
        ))
    }

    # MODE 2: REALIZED (Dataframe)
    else if (is.data.frame(input)) {
        year_col <- if ("year_col" %in% names(args)) args$year_col else "Year"
        val_col <- if ("value_col" %in% names(args)) args$value_col else "Yield"
        checks <- if ("check_list" %in% names(args)) args$check_list else NULL

        if (!all(c(year_col, val_col) %in% names(input))) stop("Columns missing.")

        df <- input
        if (!is.null(checks)) df <- df[!df$Genotype %in% checks, ]

        form <- as.formula(paste(val_col, "~", year_col))
        mod <- lm(form, data = df)

        slope <- coef(mod)[2]
        base_yr <- min(df[[year_col]], na.rm = TRUE)
        base_val <- predict(mod, newdata = data.frame(setNames(list(base_yr), year_col)))

        return(list(
            Method = "Realized (Regression)",
            Slope = slope,
            Pct_Gain_Per_Year = (slope / base_val) * 100,
            R2 = summary(mod)$r.squared
        ))
    } else {
        stop("Input must be a dataframe or asreml/fa_model object.")
    }
}

#' Calculate Heritability from List
#'
#' Robustly calculates heritability for a list of ASReml models.
#'
#' @param model_list List of ASReml objects.
#' @param id Character. The genetic term name (default "Genotype").
#'
#' @return A named vector of H2 values.
#' @export
calculate_h2_from_list <- function(model_list, id = "Genotype") {
    # Check if input is a list
    if (!is.list(model_list)) {
        stop("Input 'model_list' must be a list.")
    }

    # Initialize an empty list to store the results
    h2_results <- list()

    # Loop through each model in the list 'model_list'
    for (model_name in names(model_list)) {
        message(paste("Processing model:", model_name))
        asr <- model_list[[model_name]]

        # --- Basic Checks ---
        # Check if the element is actually an asreml object and converged
        if (is.null(asr) || !inherits(asr, "asreml")) {
            message(paste("Skipping", model_name, "- Not a valid asreml object."))
            h2_results[[model_name]] <- NA # Store NA for invalid objects
            next # Go to the next iteration
        }
        if (!asr$converge) {
            message(paste("Skipping", model_name, "- Model did not converge."))
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
                message(paste("Error processing", model_name, ":", e$message))
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
    message(paste("Number of models processed:", length(model_list)))
    message(paste("Number of successful calculations:", sum(!is.na(h2_vector))))
    message(paste("Number of failed calculations (NA):", sum(is.na(h2_vector))))
    message("------------------------------------")

    return(h2_vector)
}
