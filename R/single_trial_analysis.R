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
#' @export
analyze_single_trial <- function(data,
                                 trait = "Yield",
                                 genotype = "germplasmName",
                                 env = "studyName",
                                 row = "rowNumber",
                                 col = "colNumber",
                                 block = "blockNumber",
                                 rep = "replicate") {
    if (!requireNamespace("asreml", quietly = TRUE)) {
        stop("The 'asreml' package is required for this function.")
    }

    cat(paste0("--- Processing Trait: ", trait, " ---\n"))

    # 1. Setup Data
    df <- as.data.frame(data)

    # Check if trait/env exist
    if (!trait %in% colnames(df)) stop(paste("Trait", trait, "not found."))
    if (!env %in% colnames(df)) stop(paste("Env col", env, "not found."))

    df <- df[!is.na(df[[trait]]), ] # Remove missing trait data

    # Ensure factors
    # We use get/assign or just [[ formatting
    df[[genotype]] <- as.factor(df[[genotype]])
    df[[env]] <- as.factor(df[[env]])
    df[[block]] <- as.factor(df[[block]])
    df[[row]] <- as.factor(df[[row]])
    df[[col]] <- as.factor(df[[col]])
    df[[rep]] <- as.factor(df[[rep]])

    # Store for binding later
    all_studies_list <- list()

    # Unique Environments
    studies <- unique(as.character(df[[env]]))

    for (i in 1:length(studies)) {
        current_study <- studies[i]
        cat(paste0("  > Analyzing Study: ", current_study, "\n"))

        # Subset and Drop Levels
        sub <- droplevels(df[df[[env]] == current_study, ])

        # Check if data exists for study
        if (nrow(sub) == 0) next

        # --- 2. Spatial Grid Padding (Matches your script logic) ---
        r_vals <- as.numeric(as.character(sub[[row]]))
        c_vals <- as.numeric(as.character(sub[[col]]))

        if (length(r_vals) == 0 || length(c_vals) == 0) {
            warning(paste("No coordinates for study", current_study))
            next
        }

        r_range <- min(r_vals, na.rm = TRUE):max(r_vals, na.rm = TRUE)
        c_range <- min(c_vals, na.rm = TRUE):max(c_vals, na.rm = TRUE)

        grid <- expand.grid(rowNumber = factor(r_range), colNumber = factor(c_range))
        colnames(grid) <- c(row, col)

        # Merge and Sort (Order is critical for AR1)
        sub <- merge(grid, sub, by = c(row, col), all.x = TRUE)
        sub <- sub[order(sub[[row]]), ]
        sub <- sub[order(sub[[col]]), ] # Match script order: Row then Col

        # Check Variance
        if (var(sub[[trait]], na.rm = TRUE) > 0) {
            # Define Residual Models to Test (Matches your script)
            # We construct the formula strings dynamically
            resids_to_test <- c(
                paste0("id(", col, "):id(", row, ")"),
                paste0("ar1(", col, "):ar1(", row, ")")
            )

            # Storage for Model Selection Loop
            loop_reliabilities <- c()
            loop_bics <- c()
            loop_conv <- c()

            # Formulas (Matches your script)
            fixed_form <- as.formula(paste(trait, "~ 1"))
            random_form <- as.formula(paste("~", genotype, "+", block))

            # --- 3. MODEL SELECTION LOOP (ROBUST) ---
            for (k in 1:length(resids_to_test)) {
                resid_form <- as.formula(paste("~", resids_to_test[k]))

                # 1. Initial Fit with higher maxit
                # Note: asreml namespace call
                amod1 <- try(asreml::asreml(
                    fixed = fixed_form,
                    random = random_form,
                    residual = resid_form,
                    data = sub,
                    na.action = asreml::na.method(x = "include", y = "include"),
                    # INCREASED MAXIT HERE
                    maxit = 50,
                    trace = FALSE, workspace = "4gb"
                ), silent = TRUE)

                # 2. Force Convergence Loop
                # If it didn't crash, but didn't converge, update it up to 3 times
                if (!inherits(amod1, "try-error")) {
                    attempt <- 1
                    while (!amod1$converge && attempt <= 3) {
                        cat(paste0("    ... Updating model (Attempt ", attempt, ")\n"))
                        amod1 <- try(update(amod1, maxit = 50), silent = TRUE)
                        attempt <- attempt + 1
                    }

                    # 3. Calculate Metrics (Only if it exists)
                    if (!inherits(amod1, "try-error") && amod1$converge) {
                        # Predict Random to get PEV (ignore Intercept)
                        blups <- try(predict(amod1, classify = genotype, ignore = c("(Intercept)"))$pvals, silent = TRUE)

                        if (!inherits(blups, "try-error")) {
                            pev <- blups$std.error^2

                            # Extract Vg safely
                            vc <- summary(amod1)$varcomp
                            if (genotype %in% rownames(vc)) {
                                Vg <- vc[genotype, "component"]
                            } else {
                                Vg <- 0
                            }

                            # Reliability Calculation
                            if (Vg > 1e-6) {
                                rel_calc <- 1 - (pev / Vg)
                                loop_reliabilities <- append(loop_reliabilities, mean(rel_calc, na.rm = TRUE))
                            } else {
                                loop_reliabilities <- append(loop_reliabilities, 0)
                            }

                            # Use asreml::infoCriteria
                            bic_val <- tryCatch(asreml::infoCriteria(amod1)$BIC, error = function(e) 1e14)
                            loop_bics <- append(loop_bics, bic_val)
                            loop_conv <- append(loop_conv, TRUE)
                        } else {
                            # Prediction failed
                            loop_reliabilities <- append(loop_reliabilities, 0)
                            loop_bics <- append(loop_bics, 1e14)
                            loop_conv <- append(loop_conv, FALSE)
                        }
                    } else {
                        # Still failed after updates
                        loop_reliabilities <- append(loop_reliabilities, 0)
                        loop_bics <- append(loop_bics, 1e14)
                        loop_conv <- append(loop_conv, FALSE) # Mark as NOT converged
                    }
                } else {
                    # Crashed immediately
                    loop_reliabilities <- append(loop_reliabilities, 0)
                    loop_bics <- append(loop_bics, 1e14)
                    loop_conv <- append(loop_conv, FALSE)
                }
            }
            # --- 4. SELECT BEST MODEL ---
            best_idx <- which.min(loop_bics)

            # If all models failed (BICs are all huge), default to index 1
            if (length(best_idx) == 0) best_idx <- 1

            best_resid_model <- resids_to_test[best_idx]
            best_rel <- loop_reliabilities[best_idx]
            best_conv <- loop_conv[best_idx]

            # --- 5. FIT FINAL MODEL (Fixed Genotype for BLUEs) ---
            # Matches script: fixed includes genotype, random is just block
            asreml::asreml.options(ai.sing = TRUE, aom = TRUE)

            fixed_form_final <- as.formula(paste(trait, "~ 1 +", genotype))
            random_form_final <- as.formula(paste("~", block))
            resid_form_final <- as.formula(paste("~", best_resid_model))

            amod2 <- try(asreml::asreml(
                fixed = fixed_form_final,
                random = random_form_final,
                residual = resid_form_final,
                data = sub,
                na.action = asreml::na.method(x = "include", y = "include"),
                trace = FALSE, workspace = "4gb"
            ), silent = TRUE)

            if (!inherits(amod2, "try-error")) {
                # Ensure convergence
                if (!amod2$converge) amod2 <- try(update(amod2), silent = TRUE)

                # Predict BLUEs
                blues0 <- try(predict(amod2, classify = genotype, pworkspace = "4gb"), silent = TRUE)

                if (!inherits(blues0, "try-error")) {
                    blues_raw <- blues0$pvals

                    # Construct Output Dataframe
                    out_df <- blues_raw[, c(genotype, "predicted.value", "std.error")]
                    colnames(out_df)[1] <- "germplasmName" # Ensure name match

                    # Add Metrics
                    out_df$trait <- trait
                    out_df$study <- current_study
                    out_df$residmod <- best_resid_model
                    out_df$conv <- amod2$converge
                    out_df$rel <- best_rel
                    out_df$weight <- 1 / (out_df$std.error^2)

                    all_studies_list[[i]] <- out_df
                }
            } else {
                cat(paste0("    ! Final model failed for: ", current_study, "\n"))
            }
        }
    }

    # Bind Results
    final_blues <- do.call(rbind, all_studies_list)
    return(final_blues)
}
