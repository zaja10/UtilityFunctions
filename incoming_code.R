# ---------------------------------------------------------
# Incoming Code / Scratchpad
# ---------------------------------------------------------
# Paste your functions, snippets, or scripts below.
# I will review this file, understand the code, and integrate
# it into the appropriate locations within the package.
#
# Once a function is successfully integrated, we can clear
# it from this file.
# ---------------------------------------------------------


# Function to safely ensure ASReml model convergence
mkConv <- function(mod) {
    pctchg <- summary(mod)$varcomp[, c("%ch")]
    tries <- 1
    while (any(pctchg > 1, na.rm = TRUE) & tries < 20) {
        tries <- tries + 1
        mod <- suppressWarnings(update(mod))
        pctchg <- summary(mod)$varcomp[, c("%ch")]
    }
    return(mod)
}

# Function to fit and compare spatial models
fit_asreml_comparison <- function(fixed_f, random_f, residual_f, model_name, df) {
    tryCatch(
        {
            m <- asreml(
                fixed = fixed_f,
                random = random_f,
                residual = residual_f,
                data = df,
                trace = FALSE,
                maxiter = 50,
                na.action = na.method(y = "include", x = "include")
            )
            if (!m$converge) m <- mkConv(m) # Using your custom convergence function

            aic_val <- summary(m)$aic
            data.frame(
                Model = model_name,
                AIC = aic_val,
                LogLik = m$loglik,
                Converged = m$converge,
                stringsAsFactors = FALSE
            )
        },
        error = function(e) {
            data.frame(Model = model_name, AIC = Inf, LogLik = NA, Converged = FALSE)
        }
    )
}

# Function to run the final best model and extract metrics
run_final_analysis <- function(params, model_name, df, trial_name) {
    m_final <- asreml(
        fixed = params$f,
        random = params$r,
        residual = params$res,
        data = df,
        trace = FALSE,
        maxiter = 50,
        na.action = na.method(y = "include", x = "include")
    )
    if (!m_final$converge) m_final <- mkConv(m_final)

    vars <- m_final$sigma2 * m_final$vparameters

    # Extract Genetic Variance from the gkeep term
    Vg <- vars[grep("vm\\(gkeep", names(vars))]
    if (length(Vg) == 0) Vg <- NA

    # Extract Error and Spatial Variances
    Ve <- vars[grep("colNumber:rowNumber!R", names(vars))]
    Vr <- if (grepl("AxAe", model_name)) vars[grep("units", names(vars))] else 0

    Vtotal <- Ve + Vr
    h2 <- if (!is.na(Vg)) Vg / (Vg + Vtotal / 2) else NA

    tryCatch(
        {
            # Predict using the gkeep column to map to the G.inv
            pred <- predict(m_final, classify = "gkeep", only = "vm(gkeep, G.inv)", vcov = TRUE, sed = TRUE)
            mean_pev <- mean(diag(pred$vcov))
            avsed <- pred$avsed[2]^2 / 2
            GH_c <- if (!is.na(Vg)) (1 - avsed / Vg) else NA
            Acc <- if (!is.na(GH_c) && GH_c > 0) sqrt(GH_c) else NA
        },
        error = function(e) {
            mean_pev <<- NA
            avsed <<- NA
            GH_c <<- NA
            Acc <<- NA
        }
    )

    return(data.frame(
        Trial = trial_name, Model = model_name, rG = NA, rE = NA, rR = NA,
        Vg = Vg, Ve = Ve, Vr = Vr, h2 = h2, PEV = mean_pev, SED = avsed,
        GH_c = GH_c, Acc = Acc, LogLik = m_final$loglik, stringsAsFactors = FALSE
    ))
}

#' Extract Standardized Residuals and Identify Outliers
#'
#' This function updates a fitted asreml model to calculate the Average Outlier
#' Measure (AOM), extracts the standardized conditional residuals, and filters
#' for potential outliers based on a defined threshold.
#'
#' @param best_model A converged asreml model object
#' @param trial_data The original dataframe used for the trial
#' @param threshold Numeric value defining the outlier cutoff (default is 4)
#' @return A dataframe containing only the flagged outliers
extract_asreml_outliers <- function(best_model, trial_data, threshold = 4) {
    # 1. Update the model to compute outlier statistics
    # suppressWarnings is used to keep the console clean from minor update messages
    out.fm <- suppressWarnings(update.asreml(best_model, aom = TRUE))

    # 2. Initialize a new column for the residuals
    trial_data$scres <- NA

    # 3. Safely map residuals back to the dataframe
    # ASReml's working dataframe (out.fm$mf) only contains rows without missing target values.
    # We extract the row indices of the working data to place the residuals correctly.
    valid_rows <- as.numeric(rownames(out.fm$mf))

    # Extract the standardized conditional residuals (column 2 of the AOM R matrix)
    trial_data$scres[valid_rows] <- out.fm$aom$R[, 2]

    # 4. Filter the dataframe for true outliers using absolute value
    outliers <- trial_data[!is.na(trial_data$scres) & abs(trial_data$scres) > threshold, ]

    # 5. Keep only the most relevant columns to make the review easy
    outliers_clean <- subset(outliers, select = c(
        studyName, colNumber, rowNumber,
        replicate, germplasmName,
        grain_yield, scres
    ))

    return(outliers_clean)
}


# 2. Define a helper function to evaluate matches for a single row
match_germplasm <- function(name, syn_string, valid_ids) {
    # Convert inputs to plain characters (since they are currently factors)
    name <- as.character(name)
    syn_string <- as.character(syn_string)

    # Check 1: Does the primary germplasmName exist in the genomic matrix?
    if (name %in% valid_ids) {
        return(name)
    }

    # Check 2: If no primary match, check the synonyms
    if (!is.na(syn_string) && syn_string != "") {
        # Split the string by comma and remove any accidental blank spaces
        syn_list <- str_trim(unlist(str_split(syn_string, ",")))

        # Identify which synonyms exist in the genomic data
        matched_syns <- syn_list[syn_list %in% valid_ids]

        # If at least one synonym matches, return the first valid one
        if (length(matched_syns) > 0) {
            return(matched_syns[1])
        }
    }

    # Return NA if absolutely nothing matches
    return(NA_character_)
}

# 3. Apply the matching logic to your phenotypic trial dataframe
YT.df <- YT.df %>%
    rowwise() %>%
    mutate(
        # Create the gkeep column using our custom matching function
        gkeep = match_germplasm(germplasmName, germplasmSynonyms, genotyped_ids),

        # Create the gdrop column based on the result of gkeep
        gdrop = if_else(is.na(gkeep), as.character(germplasmName), NA_character_)
    ) %>%
    ungroup() # Always ungroup after rowwise operations to maintain performance later


# Get a list of unique trials
trial_list <- unique(as.character(YT.df$studyName))

for (trial_name_str in trial_list) {
    message("\n==================================================")
    message(paste("Processing Trial:", trial_name_str))
    message("==================================================")

    # Define expected output files
    selection_file <- file.path(tables_dir, paste0(trial_name_str, "_Spatial_Model_Selection.csv"))
    metrics_file <- file.path(tables_dir, paste0(trial_name_str, "_Final_Metrics.csv"))
    plot_file_path <- file.path(diagnostics_dir, paste0(trial_name_str, "_Diagnostics.pdf"))

    # 1. FULL SKIP: If the diagnostics plot already exists, the trial is 100% complete
    if (file.exists(plot_file_path)) {
        message("  [✓] Trial fully processed (Models and Diagnostics complete). Skipping.")
        next
    }

    tryCatch(
        {
            # --- Data Preparation ---
            trial_data <- YT.df %>% filter(studyName == trial_name_str)
            if (nrow(trial_data) == 0) stop("No data found for this trial.")

            analysis_data <- trial_data %>%
                mutate(Yield = as.numeric(scale(grain_yield))) %>%
                mutate(across(c(rowNumber, colNumber, any_of("blockNumber")), as.factor)) %>%
                arrange(colNumber, rowNumber)

            analysis_data$gkeep <- factor(analysis_data$gkeep, levels = attr(G.inv, "rowNames"))
            analysis_data$gdrop <- as.factor(analysis_data$gdrop)

            # --- Build Random Terms ---
            has_blocks <- "blockNumber" %in% names(analysis_data) && n_distinct(analysis_data$blockNumber, na.rm = TRUE) > 1
            has_ungenotyped <- sum(!is.na(analysis_data$gdrop)) > 0

            rand_terms <- c("vm(gkeep, G.inv)", "colNumber", "rowNumber")
            if (has_ungenotyped) rand_terms <- c(rand_terms, "gdrop")
            if (has_blocks) rand_terms <- c(rand_terms, "blockNumber")

            base_rand_str <- paste("~", paste(rand_terms, collapse = " + "))

            # Define available models
            model_list <- list(
                "RCD"  = list(f = Yield ~ 1, r = as.formula(base_rand_str), res = ~ colNumber:rowNumber),
                "AxA"  = list(f = Yield ~ 1, r = as.formula(base_rand_str), res = ~ ar1(colNumber):ar1(rowNumber)),
                "AxAe" = list(f = Yield ~ 1, r = as.formula(paste(base_rand_str, "+ units")), res = ~ ar1(colNumber):ar1(rowNumber))
            )

            # 2. SMART BYPASS OR FULL RUN
            if (file.exists(selection_file) && file.exists(metrics_file)) {
                # BYPASS: Read the existing file to find the best model
                message("  [✓] Model selection files found. Loading best model...")
                saved_results <- read.csv(selection_file)
                best_model_name <- saved_results$Model[1] # Assumes CSV is sorted by AIC
                best_params <- model_list[[best_model_name]]
            } else {
                # FULL RUN: Fit all models, compare, and save
                message("  Fitting Univariate Spatial Models for comparison...")
                asreml_results <- purrr::imap_dfr(model_list, ~ fit_asreml_comparison(.x$f, .x$r, .x$res, .y, analysis_data)) %>% arrange(AIC)

                best_model_name <- asreml_results$Model[1]
                best_params <- model_list[[best_model_name]]

                univ_results <- run_final_analysis(best_params, paste0(best_model_name, "_Univ_Yield"), analysis_data, trial_name_str)
                all_model_results <- bind_rows(univ_results %>% mutate(across(everything(), as.character)))

                write.csv(asreml_results, file = selection_file, row.names = FALSE)
                write.csv(all_model_results, file = metrics_file, row.names = FALSE)
                message(paste("  [✓] Model selection saved. Best Model:", best_model_name))
            }

            # 3. FIT THE BEST MODEL FOR DIAGNOSTICS
            message(paste("  Refitting best model (", best_model_name, ") for diagnostics..."))
            best_m_fit <- asreml(
                fixed = best_params$f, random = best_params$r,
                residual = best_params$res, data = analysis_data,
                trace = FALSE, maxiter = 50,
                na.action = na.method(y = "include", x = "include")
            )

            # Ensure it converged before extracting residuals
            if (!best_m_fit$converge) best_m_fit <- mkConv(best_m_fit)

            # --- Outlier Detection & Plotting ---
            message("  Extracting outliers and saving diagnostic plots...")

            trial_outliers <- extract_asreml_outliers(
                best_model = best_m_fit,
                trial_data = analysis_data,
                threshold = 4
            )

            if (nrow(trial_outliers) > 0) {
                write.csv(trial_outliers,
                    file = file.path(tables_dir, paste0(trial_name_str, "_Flagged_Outliers.csv")),
                    row.names = FALSE
                )
                message(paste("  [!] Found", nrow(trial_outliers), "outliers. Saved to review list."))
            } else {
                message("  [✓] No outliers > 4 detected.")
            }

            # Update model to generate Average Outlier Measure (AOM) for plotting
            out.fm <- suppressWarnings(update.asreml(best_m_fit, aom = TRUE))

            pdf(file = plot_file_path, width = 8, height = 6)

            tryCatch(
                {
                    plot(out.fm, res = "stdCond", main = paste("Standardised Residuals:", trial_name_str))
                },
                error = function(e) message("    [!] Could not plot residuals: ", e$message)
            )

            tryCatch(
                {
                    met.asreml(best_m_fit)
                },
                error = function(e) message("    [!] Could not plot variogram: ", e$message)
            )

            dev.off()
            message("  [✓] Diagnostic plots saved to folder.")
        },
        error = function(e) {
            message(paste("  [!] Error processing trial:", e$message))
        }
    )
}


run_vi_mlr_pipeline <- function(data, target_yield, vi_names, method_name = "MLR_UAS") {
    # 1. CLEAN DATA
    # Ensure we only use rows with yield and at least some VI data
    valid_data <- data %>% filter(!is.na(!!sym(target_yield)))

    # 2. CONSTRUCT FORMULA
    # Format: Yield ~ VI_1 + VI_2 + ...
    mlr_formula <- as.formula(paste(target_yield, "~", paste(vi_names, collapse = " + ")))

    # 3. EXTRACT GLOBAL IMPORTANCE
    # Using the absolute t-statistic as an importance weight
    global_fit <- lm(mlr_formula, data = valid_data)
    importance_df <- tidy(global_fit) %>%
        filter(!term %in% c("(Intercept)")) %>%
        mutate(
            Method = method_name,
            Importance_Weight = abs(statistic),
            P_Value = p.value
        ) %>%
        select(Method, Predictor = term, Importance_Weight, P_Value) %>%
        arrange(desc(Importance_Weight))

    # 4. RUN LEAVE-ONE-COLUMN-OUT (LOCO) CV
    # This tests if the VIs predict yield in "new" parts of the field
    data$Imputed_Yield <- NA_real_
    unique_cols <- unique(data$colNumber[!is.na(data$colNumber)])

    for (k in unique_cols) {
        train_data <- data %>% filter(colNumber != k & !is.na(!!sym(target_yield)))
        test_data <- data %>% filter(colNumber == k)

        if (nrow(train_data) < length(vi_names)) next

        tryCatch(
            {
                # Fit model on everything EXCEPT column k
                fit <- lm(mlr_formula, data = train_data)
                # Predict for column k
                data$Imputed_Yield[data$colNumber == k] <- predict(fit, newdata = test_data)
            },
            error = function(e) {
                message(paste("Column", k, "skipped due to rank deficiency."))
            }
        )
    }

    # 5. CALCULATE ACCURACIES
    predictions_df <- data %>%
        select(plotNumber, colNumber, rowNumber, germplasmName, Imputed_Yield)

    # Validation stats
    valid_res <- data %>% filter(!is.na(!!sym(target_yield)) & !is.na(Imputed_Yield))
    acc_r <- if (nrow(valid_res) > 2) cor(valid_res[[target_yield]], valid_res$Imputed_Yield) else NA

    col_acc_df <- valid_res %>%
        group_by(colNumber) %>%
        summarize(
            N_Plots = n(),
            Column_Acc_r = if (n() > 2 && sd(Imputed_Yield) > 0) cor(!!sym(target_yield), Imputed_Yield) else NA_real_,
            .groups = "drop"
        ) %>%
        mutate(Method = method_name)

    return(list(
        Acc = acc_r,
        Importance = importance_df,
        Col_Acc = col_acc_df,
        Predictions = predictions_df
    ))
}


prepare_trial_grm <- function(master_grm, trial_data) {
    trial_lines <- as.character(unique(trial_data$germplasmName))
    genotyped_lines <- intersect(trial_lines, rownames(master_grm))
    ungenotyped_lines <- setdiff(trial_lines, genotyped_lines)

    if (length(genotyped_lines) == 0) stop("No genotyped lines found in the GRM for this trial.")

    grm_sub <- master_grm[genotyped_lines, genotyped_lines, drop = FALSE]

    # Scale the GRM
    scale_factor <- mean(diag(grm_sub))
    grm_sub <- grm_sub / scale_factor

    # Pad ungenotyped lines
    if (length(ungenotyped_lines) > 0) {
        pad_mat <- diag(1, nrow = length(ungenotyped_lines))
        rownames(pad_mat) <- colnames(pad_mat) <- ungenotyped_lines
        grm_final <- as.matrix(bdiag(grm_sub, pad_mat))
        rownames(grm_final) <- colnames(grm_final) <- c(genotyped_lines, ungenotyped_lines)
    } else {
        grm_final <- grm_sub
    }

    diag(grm_final) <- diag(grm_final) + 1e-4
    ordered_levels <- levels(factor(trial_data$germplasmName))

    return(grm_final[ordered_levels, ordered_levels])
}


fit_asreml_comparison <- function(fixed_f, random_f, residual_f, model_name, df) {
    df <- df %>% mutate(Yield = as.numeric(scale(Yield)))
    tryCatch(
        {
            m <- asreml(fixed = fixed_f, random = random_f, residual = residual_f, data = df, trace = FALSE, maxiter = 50, na.action = na.method(y = "include", x = "include"))
            if (!m$converge) m <- update.asreml(m)
            aic_val <- summary(m)$aic
            data.frame(Model = model_name, AIC = aic_val, LogLik = m$loglik, Converged = m$converge, stringsAsFactors = FALSE)
        },
        error = function(e) data.frame(Model = model_name, AIC = Inf, LogLik = NA, Converged = FALSE)
    )
}

run_final_analysis <- function(params, model_name, df, trial_name) {
    m_final <- asreml(fixed = params$f, random = params$r, residual = params$res, data = df, trace = FALSE, maxiter = 50, na.action = na.method(y = "include", x = "include"))
    if (!m_final$converge) m_final <- update.asreml(m_final, maxit = 100)

    vars <- m_final$sigma2 * m_final$vparameters
    Vg <- vars[grep("vm\\(germplasmName", names(vars))]
    if (length(Vg) == 0) Vg <- NA
    Ve <- vars[grep("colNumber:rowNumber!R", names(vars))]
    Vr <- if (grepl("AxAe", model_name)) vars[grep("units", names(vars))] else 0
    Vtotal <- Ve + Vr
    h2 <- if (!is.na(Vg)) Vg / (Vg + Vtotal / 2) else NA

    tryCatch(
        {
            pred <- predict(m_final, classify = "germplasmName", only = "vm(germplasmName, trial_ginv)", vcov = TRUE, sed = TRUE)
            mean_pev <- mean(diag(pred$vcov))
            avsed <- pred$avsed[2]^2 / 2
            GH_c <- if (!is.na(Vg)) (1 - avsed / Vg) else NA
            Acc <- if (!is.na(GH_c) && GH_c > 0) sqrt(GH_c) else NA
        },
        error = function(e) {
            mean_pev <<- NA
            avsed <<- NA
            GH_c <<- NA
            Acc <<- NA
        }
    )

    return(data.frame(
        Trial = trial_name, Model = model_name, rG = NA, rE = NA, rR = NA,
        Vg = Vg, Ve = Ve, Vr = Vr, h2 = h2, PEV = mean_pev, SED = avsed,
        GH_c = GH_c, Acc = Acc, LogLik = m_final$loglik, stringsAsFactors = FALSE
    ))
}


run_bivariate_analysis <- function(secondary_trait, df, rand_f, res_f, model_name, trial_name) {
    pair_data <- df %>% filter(Trait %in% c("Yield", secondary_trait))
    sec_trait_str <- paste(secondary_trait, collapse = "_")
    full_model_name <- paste0(model_name, "_Biv_", sec_trait_str)
    res_f_str <- paste(res_f, collapse = ", ")

    tryCatch(
        {
            m_biv <- asreml(fixed = Value_Scaled ~ Trait, random = rand_f, residual = res_f, data = pair_data, trace = FALSE, maxiter = 50, na.action = na.method(y = "include", x = "include"))
            if (!m_biv$converge) m_biv <- update.asreml(m_biv, maxit = 100)

            vars <- m_biv$sigma2 * m_biv$vparameters

            Vg_matches <- vars[grep("Trait:vm\\(germplasmName.*!Trait_Yield:Yield", names(vars))]
            Vg <- if (length(Vg_matches) > 0) Vg_matches[1] else NA

            rG_mat <- summary(m_biv, vparameters = TRUE)$vparameters$`Trait:vm(germplasmName, trial_ginv)`
            rG <- if (!is.null(rG_mat) && is.matrix(rG_mat)) cov2cor(rG_mat)[2, 1] else NA

            if (grepl("us\\(", deparse(res_f))) {
                vparams <- summary(m_biv, vparameters = TRUE)$vparameters
                rE_mat <- vparams$`Trait:colNumber:rowNumber`$Trait
                rE <- if (!is.null(rE_mat) && is.matrix(rE_mat)) cov2cor(rE_mat)[2, 1] else NA

                Ve_matches <- vars[grep("Trait:colNumber:rowNumber!Trait_Yield:Yield", names(vars))]
                Ve <- if (length(Ve_matches) > 0) Ve_matches[1] else NA
            } else {
                rE <- NA
                Ve_matches <- vars[grep("Trait_Yield!R", names(vars))]
                Ve <- if (length(Ve_matches) > 0) Ve_matches[1] else NA
            }

            if (grepl("AxAe", model_name)) {
                Vr_matches <- vars[grep("Trait:units!Trait_Yield(:Yield)?$", names(vars))]
                Vr <- if (length(Vr_matches) > 0) Vr_matches[1] else 0

                # Using string match on model_name to avoid formula parsing bug
                if (grepl("usUnits", model_name)) {
                    rR_mat <- summary(m_biv, vparameters = TRUE)$vparameters$`Trait:units`
                    rR <- if (!is.null(rR_mat) && is.matrix(rR_mat)) cov2cor(rR_mat)[2, 1] else NA
                } else {
                    rR <- NA
                }
            } else {
                Vr <- 0
                rR <- NA
            }

            Vtotal <- Ve + Vr
            h2 <- if (!is.na(Vg) && !is.na(Vtotal)) Vg / (Vg + Vtotal / 2) else NA

            pred <- predict(m_biv, classify = "Trait:germplasmName", only = "us(Trait):vm(germplasmName, trial_ginv)", levels = list(Trait = "Yield"), pworkspace = 1e8, vcov = TRUE, sed = TRUE)
            yield_idx <- which(pred$pvals$Trait == "Yield")
            mean_pev <- mean(diag(pred$vcov[yield_idx, yield_idx]))
            mean_sed <- mean(pred$sed[yield_idx, yield_idx][upper.tri(pred$sed[yield_idx, yield_idx])])
            avsed <- mean_sed^2 / 2
            GH_c <- if (!is.na(Vg)) (1 - avsed / Vg) else NA
            Acc <- if (!is.na(GH_c) && GH_c > 0) sqrt(GH_c) else NA

            return(data.frame(Trial = trial_name, Model = full_model_name, rG = paste(rG, collapse = ", "), rE = paste(rE, collapse = ", "), rR = paste(rR, collapse = ", "), Vg = paste(Vg, collapse = ", "), Ve = paste(Ve, collapse = ", "), Vr = paste(Vr, collapse = ", "), h2 = paste(h2, collapse = ", "), PEV = paste(mean_pev, collapse = ", "), SED = paste(avsed, collapse = ", "), GH_c = paste(GH_c, collapse = ", "), Acc = paste(Acc, collapse = ", "), residual = res_f_str, LogLik = m_biv$loglik, stringsAsFactors = FALSE))
        },
        error = function(e) {
            message(paste("  [!] Skipping", full_model_name, "- Variance structure issue."))
            return(data.frame(Trial = trial_name, Model = full_model_name, rG = NA_character_, rE = NA_character_, rR = NA_character_, Vg = NA_character_, Ve = NA_character_, Vr = NA_character_, h2 = NA_character_, PEV = NA_character_, SED = NA_character_, GH_c = NA_character_, Acc = NA_character_, residual = res_f_str, LogLik = NA_real_, stringsAsFactors = FALSE))
        }
    )
}


for (current_file in processed_files) {
    # Extract trial name from the filename
    trial_name_str <- str_remove(basename(current_file), "_AnalysisReady\\.parquet$")

    message("\n==================================================")
    message(paste("Processing Trial:", trial_name_str))
    message("==================================================")

    # Check if already processed
    expected_output <- file.path(tables_dir, paste0(trial_name_str, "_Final_Metrics.csv"))
    if (file.exists(expected_output)) {
        message("  [✓] Already analyzed. Skipping.")
        next
    }

    # Wrap the execution in a tryCatch so one bad trial doesn't stop the loop
    tryCatch(
        {
            # Load Data
            final_trial_clean <- read_parquet(current_file)

            # Map Protocol & Load GRM
            current_protocol <- protocol_map %>%
                filter(study_name == trial_name_str) %>%
                pull(genoProtocolName) %>%
                unique()
            if (length(current_protocol) == 0 || is.na(current_protocol)) stop("Protocol not found in CSV.")

            safe_protocol_name <- gsub(" ", "_", current_protocol)
            master_grm_path <- here("data", "processed", "grms", paste0(safe_protocol_name, "_grm.rds"))
            if (!file.exists(master_grm_path)) stop(paste("GRM missing:", master_grm_path))

            master_grm <- readRDS(master_grm_path)

            # Assign the GRM globally so ASReml can find it
            trial_ginv <<- prepare_trial_grm(master_grm, final_trial_clean)

            analysis_data <- final_trial_clean %>%
                mutate(Yield = as.numeric(scale(Yield))) %>%
                mutate(across(c(rowNumber, colNumber, any_of("blockNumber"), germplasmName), as.factor)) %>%
                arrange(colNumber, rowNumber)

            # --- Univariate Models ---
            message("  Fitting Univariate Spatial Models...")
            has_blocks <- "blockNumber" %in% names(analysis_data) && n_distinct(analysis_data$blockNumber, na.rm = TRUE) > 1
            rand_terms <- c("vm(germplasmName, trial_ginv)", "colNumber", "rowNumber")
            if (has_blocks) rand_terms <- c(rand_terms, "blockNumber")
            base_rand_str <- paste("~", paste(rand_terms, collapse = " + "))

            model_list <- list(
                "RCD"  = list(f = Yield ~ 1, r = as.formula(base_rand_str), res = ~ colNumber:rowNumber),
                "AxA"  = list(f = Yield ~ 1, r = as.formula(base_rand_str), res = ~ ar1(colNumber):ar1(rowNumber)),
                "AxAe" = list(f = Yield ~ 1, r = as.formula(paste(base_rand_str, "+ units")), res = ~ ar1(colNumber):ar1(rowNumber))
            )

            asreml_results <- imap_dfr(model_list, ~ fit_asreml_comparison(.x$f, .x$r, .x$res, .y, analysis_data)) %>% arrange(AIC)
            best_model_name <- asreml_results$Model[1]
            best_params <- model_list[[best_model_name]]

            univ_results <- run_final_analysis(best_params, paste0(best_model_name, "_Univ_Yield"), analysis_data, trial_name_str)

            # --- Bivariate Models ---
            message("  Fitting Bivariate Models...")
            bivariate_data_long <- final_trial_clean %>%
                pivot_longer(cols = any_of(c("Yield", "PC1", "PC2", "Imputed_Yield", "TestWeight")), names_to = "Trait", values_to = "Value") %>%
                group_by(Trait) %>%
                mutate(Value_Scaled = as.numeric(scale(Value))) %>%
                ungroup() %>%
                mutate(Trait = as.factor(Trait)) %>%
                arrange(Trait, colNumber, rowNumber)

            rand_terms_base <- c("us(Trait):vm(germplasmName, trial_ginv)", "diag(Trait):colNumber", "diag(Trait):rowNumber")
            if (has_blocks) rand_terms_base <- c(rand_terms_base, "diag(Trait):blockNumber")

            winning_res <- if (best_model_name == "RCD") "id(colNumber):id(rowNumber)" else "ar1(colNumber):ar1(rowNumber)"
            biv_res_f_dsum <- as.formula(paste("~ dsum(~", winning_res, "| Trait)"))
            biv_res_f_us <- as.formula(paste("~ us(Trait):", winning_res))

            # Dynamically inject the 4 configurations for AxAe
            if (best_model_name == "AxAe") {
                biv_configs <- list(
                    "dsum_diagUnits" = list(res_f = biv_res_f_dsum, rand_f = as.formula(paste("~", paste(c(rand_terms_base, "diag(Trait):units"), collapse = " + ")))),
                    "dsum_usUnits"   = list(res_f = biv_res_f_dsum, rand_f = as.formula(paste("~", paste(c(rand_terms_base, "us(Trait):units"), collapse = " + ")))),
                    "us_diagUnits"   = list(res_f = biv_res_f_us, rand_f = as.formula(paste("~", paste(c(rand_terms_base, "diag(Trait):units"), collapse = " + ")))),
                    "us_usUnits"     = list(res_f = biv_res_f_us, rand_f = as.formula(paste("~", paste(c(rand_terms_base, "us(Trait):units"), collapse = " + "))))
                )
            } else {
                base_rand_f <- as.formula(paste("~", paste(rand_terms_base, collapse = " + ")))
                biv_configs <- list(
                    "dsum" = list(res_f = biv_res_f_dsum, rand_f = base_rand_f),
                    "us"   = list(res_f = biv_res_f_us, rand_f = base_rand_f)
                )
            }

            available_traits <- intersect(names(final_trial_clean), c("PC1", "PC2", "Imputed_Yield", "TestWeight"))
            secondary_traits <- list()
            if ("PC1" %in% available_traits) secondary_traits <- append(secondary_traits, list("PC1"))
            if ("PC1" %in% available_traits && "PC2" %in% available_traits) secondary_traits <- append(secondary_traits, list(c("PC1", "PC2")))
            remaining_traits <- setdiff(available_traits, c("PC1", "PC2"))
            if (length(remaining_traits) > 0) secondary_traits <- append(secondary_traits, as.list(remaining_traits))

            execution_grid <- expand_grid(trait_var = secondary_traits, config_name = names(biv_configs))
            biv_results <- pmap_dfr(execution_grid, function(trait_var, config_name) {
                config <- biv_configs[[config_name]]
                run_bivariate_analysis(trait_var, bivariate_data_long, config$rand_f, config$res_f, paste0(best_model_name, "_", config_name), trial_name_str)
            })

            # --- Save Outputs ---
            all_model_results <- bind_rows(
                univ_results %>% mutate(across(everything(), as.character)),
                biv_results %>% mutate(across(everything(), as.character))
            )

            write.csv(asreml_results, file = file.path(tables_dir, paste0(trial_name_str, "_Spatial_Model_Selection.csv")), row.names = FALSE)
            write.csv(all_model_results, file = file.path(tables_dir, paste0(trial_name_str, "_Final_Metrics.csv")), row.names = FALSE)

            message("  [✓] Successfully saved all outputs.")
        },
        error = function(e) {
            message(paste("  [!] Error processing trial:", e$message))
        }
    )
}
vccmake.mfxlm <- function(obj = asr.sv, data = set1.df, Rownam = "Row!cor",
                          Colnam = "Col!cor", Resnam = "!R$", Envnam = "Env", Trialnam = "Trial") {
    temp <- obj$vparameters.table
    sigmas <- temp[grep(Resnam, temp$Component), ]
    colphis <- temp[grep(Colnam, temp$Component), ]
    rowphis <- temp[grep(Rownam, temp$Component), ]
    ns <- dim(sigmas)[1]
    ncp <- dim(colphis)[1]
    nrp <- dim(rowphis)[1]
    temp.ss <- rbind(sigmas, colphis, rowphis)
    order <- temp$Component[temp$Component %in% temp.ss$Component]
    temp.ss$Trial <- sapply(strsplit(as.character(temp.ss$Component), split = c("_")), function(x) x[2])
    temp.ss$Trial <- sapply(strsplit(temp.ss$Trial, split = c("!")), function(x) x[1])
    temp.ss$Env <- factor(tapply(as.character(data[[Envnam]]), data[[Trialnam]], unique)[temp.ss$Trial])
    temp.ss$Partype <- factor(rep(c("sig", "colphi", "rowphi"), times = c(ns, ncp, nrp)))
    names(temp.ss)[1] <- "Vparameter"
    if ((ncp + nrp) == 0) {
        if (length(levels(temp.ss$Env)) == 1) {
            M <- vcm.lm(~1, data = temp.ss) # note: vcm.lm is an asreml function
        } else {
            M <- vcm.lm(~Env, data = temp.ss)
        }
    } else {
        if (length(levels(temp.ss$Env)) == 1) {
            M <- vcm.lm(~Partype, data = temp.ss)
        } else {
            M <- vcm.lm(~ Env:Partype, data = temp.ss)
        }
    }
    attr(M, "assign") <- NULL
    attr(M, "contrasts") <- NULL
    if (ns > 0) {
        envmeans <- with(subset(temp.ss, Partype == "sig"), tapply(Value, Env, mean))
        temp.ss$Value[temp.ss$Partype == "sig"] <- envmeans[as.character(temp.ss$Env[temp.ss$Partype == "sig"])]
    }
    if (ncp > 0) {
        envmeans <- with(subset(temp.ss, Partype == "colphi"), tapply(Value, Env, mean))
        temp.ss$Value[temp.ss$Partype == "colphi"] <- envmeans[as.character(temp.ss$Env[temp.ss$Partype == "colphi"])]
    }
    if (nrp > 0) {
        envmeans <- with(subset(temp.ss, Partype == "rowphi"), tapply(Value, Env, mean))
        temp.ss$Value[temp.ss$Partype == "rowphi"] <- envmeans[as.character(temp.ss$Env[temp.ss$Partype == "rowphi"])]
    }
    names(temp.ss)[1] <- "Component"
    rownames(temp.ss) <- temp.ss$Component
    Mcc <- cbind((col(M) * M) %*% rep(1, dim(M)[2]), rep(1, dim(M)[1]))
    # x <- as.vector(t(matrix(1:231, nrow=77, ncol=3)))
    # M2 <- M[x,]
    temp.ss <- temp.ss[order, ]
    Mvcm <- M[order, ]
    return(list(gammas = temp.ss, Mcc = Mcc))
}
faidesum.mfxlm <- function(obj.asr, data.df, ide.k = 1, Gfac = "name", Efac = "location", trim = FALSE) {
    ide.sum <- list()
    snam <- as.character(levels(data.df[[Efac]]))
    1
    gnam <- as.character(levels(data.df[[Gfac]]))
    ss <- summary(obj.asr)$varcomp
    ss$param <- row.names(ss)
    # get all variance parameters associated with Gfac (rr+diag) in model
    pp <- ss[grep(Gfac, ss$param), c("component", "bound", "param")]
    ide.gam.temp <- loadsum.mfxlm(param = pp, snam, k = ide.k, trim)
    ide.sum$gammas <- ide.gam.temp$out.list
    # now for the blups
    cc <- obj.asr$coef$random
    vc <- obj.asr$vcoeff$random
    # fix problems with funny characters '\xa0' so can use grep
    dimnames(cc)[[1]] <- enc2native(dimnames(cc)[[1]])
    names(vc) <- enc2native(dimnames(cc)[[1]])
    cc.ide <- cc[grep(Gfac, dimnames(cc)[[1]]), 1]
    vc.ide <- vc[grep(Gfac, names(vc))]
    ide.blup.temp <- blupsum.mfxlm(cc.ide, vc.ide, ide.k, gnam, snam, data.df, Gfac, Efac,
        ide.gam.temp,
        Excludevm = NULL, trim
    )
    ide.sum$blups <- ide.blup.temp$blup.df
    ide.sum$scores <- ide.blup.temp$score.df
    # add accuracies to CVE blups and re-name regblup as CVE
    names(ide.sum$blups)[4:5] <- c("CVE", "seCVE")
    Lam <- ide.sum$gammas$`rotated loads`
    Dmat <- ide.sum$gammas$Dmat
    gvar <- diag(Lam %*% Dmat %*% t(Lam))
    names(gvar) <- as.character(dimnames(Lam)[[1]])
    blups <- ide.sum$blups
    blups$gvar <- gvar[blups[[Efac]]]
    ide.sum$blups$accCVE <- sqrt(1 - ide.sum$blups$seCVE^2 / blups$gvar)
    ide.sum
}
loadsum.mfxlm <- function(param, snam, k, trim) {
    nsite <- length(snam)
    # get any GxE parameters with .var in name. this will include dummy specific variances
    # (all fixed to 0) from rr term and actual specific variances from diag term
    # psi <- param[grep('.var',param$param),]
    # get rid of specific variances associated with rr term
    psi <- param[grep("rr\\(", param$param, invert = TRUE), "component"]
    # 24/10/20 extra to allow trimming of LOF
    if (trim) {
        ss.psi <- param[grep("rr\\(", param$param, invert = TRUE), ]
        # re-order so in alphabetical order of environments (assumes snam is alpha order)
        # ss.psi will have names of at(Efac, ??)... where ?? is an element of snam
        # so alpha order of names will put in alpha order of snam
        ss.psi <- ss.psi[order(ss.psi$param), ]
        psi <- ss.psi$component
    }
    names(psi) <- snam
    Lam <- matrix(NA, nrow = nsite, ncol = k)
    Lam[, 1] <- param[grep(".fa1", param$param), ]$component
    if (k > 1) {
        for (i in 2:k) {
            Lam[, i] <- param[grep(paste(".fa", i, sep = ""), param$param), ]$component
        }
    }
    Lam.orig <- Lam
    svdL <- svd.mfxlm(Lam.orig)
    Lam <- svdL$u
    Dmat <- svdL$d2
    2
    dimnames(Lam) <- list(snam, paste("fac", 1:k, sep = "_"))
    dimnames(Dmat) <- list(paste("fac", 1:k, sep = "_"), paste("fac", 1:k, sep = "_"))
    Gmat <- Lam %*% Dmat %*% t(Lam) + diag(psi)
    Cmat <- cov2cor(Gmat)
    paf.site <- matrix(0, nrow = nsite, ncol = k)
    dimnames(paf.site) <- list(snam, paste("fac", 1:k, sep = "_"))
    for (i in 1:k) {
        paf.site[, i] <- 100 * (Dmat[i, i] * diag(Lam[, i] %*% t(Lam[, i]))) / diag(Gmat)
    }
    if (k > 1) {
        all <- 100 * diag(Lam %*% Dmat %*% t(Lam)) / diag(Gmat)
        paf.site <- cbind(paf.site, all)
    }
    paf.mod <- 100 * sum(diag(Lam %*% Dmat %*% t(Lam))) / sum(diag(Gmat))
    paf.fac <- 100 * diag(Dmat) / sum(diag(Gmat))
    dd <- 1 / sqrt(diag(Gmat))
    Lamc <- diag(dd) %*% Lam
    dimnames(Lamc) <- dimnames(Lam.orig) <- dimnames(Lam)
    out.list <- list(Gmat, Cmat, paf.site, paf.fac, paf.mod, Lam, psi, Lamc, Lam.orig, Dmat)
    names(out.list) <- c(
        "Gmat", "Cmat", "site %vaf", "factor %vaf", "total %vaf",
        "rotated loads", "specific var", "rotated loads- c", "raw loads", "Dmat"
    )
    list(out.list = out.list, svdL = svdL)
}
svd.mfxlm <- function(Lam) {
    # note can use this for nfac=1. make sure Lam is a matrix first
    svdL <- svd(Lam)
    u1neg <- 100 * length(svdL$u[, 1][svdL$u[, 1] < 0]) / dim(Lam)[[1]]
    if (u1neg > 50) {
        svdL$u <- -1 * svdL$u
        svdL$v <- -1 * svdL$v
    }
    # make svdL$d a diagonal matrix here- especially important
    # for nfac=1 case- also save squared singular values
    svdL$d2 <- diag(svdL$d^2, nrow = length(svdL$d))
    svdL$d <- diag(svdL$d, nrow = length(svdL$d))
    svdL
}
blupsum.mfxlm <- function(coef, vcoef, nfac, gnam.ped, snam, data.df, Gfac, Efac, gam.temp, Excludevm, trim) {
    nsite <- length(snam)
    ngeno.ped <- length(gnam.ped) # genos in pedigree
    # get regblups and pev directly from rr term
    coef.all <- coef[grep("rr\\(", names(coef))]
    # check length
    if (length(coef.all) != (nsite + nfac) * ngeno.ped) {
        cat("\nWarning: rr BLUPs are incorrect length")
    }
    coef.rr <- coef.all[1:(nsite * ngeno.ped)] # take first nsite x ngeno.ped blups. rest are scores
    vcoef.rr <- vcoef[grep("rr\\(", names(vcoef))][1:(nsite * ngeno.ped)]
    # get delta blups from diag term (ie. not rr term)
    # 24/10/20 only do this if not using trim.
    if (trim == FALSE) {
        coef.dd <- coef[grep("rr\\(", names(coef), invert = TRUE)]
        if (!is.null(Excludevm)) coef.dd <- coef.dd[grep(Excludevm, names(coef.dd), invert = TRUE)] # NEW LINE!!!
        # check length
        if (length(coef.dd) != nsite * ngeno.ped) {
            cat("\nWarning: delta BLUPs are incorrect length")
        }
        3
    }
    if (trim == TRUE) {
        coef.dd <- rep(NA, length(coef.rr))
    }
    blup.df <- data.frame(blup = coef.rr + coef.dd, regblup = coef.rr, seregblup = sqrt(vcoef.rr))
    blup.df[[Gfac]] <- rep(gnam.ped, nsite)
    blup.df[[Efac]] <- rep(snam, each = ngeno.ped)
    blup.df$GE <- paste(blup.df[[Gfac]], blup.df[[Efac]], sep = "x")
    pp <- table(data.df[[Gfac]], data.df[[Efac]])
    gnam <- dimnames(pp)[[1]]
    ngeno <- length(gnam)
    reps <- as.numeric(pp)
    names(reps) <- paste(rep(gnam, nsite), rep(snam, each = ngeno), sep = "x")
    blup.df$reps <- reps[blup.df$GE]
    blup.df$reps[is.na(blup.df$reps)] <- 0
    blup.df <- blup.df[, c(5, 4, 1, 2, 3, 7)]
    # now for scores ...
    scores <- coef.all[(nsite * ngeno.ped + 1):length(coef.all)]
    score.mat <- matrix(scores, ncol = nfac)
    score.mat <- score.mat %*% gam.temp$svdL$v %*% gam.temp$svdL$d
    score.df <- data.frame(score = as.vector(score.mat))
    score.df$comp <- paste("Comp", rep(1:nfac, each = ngeno.ped), sep = "")
    score.df[[Gfac]] <- rep(gnam.ped, nfac)
    score.df <- score.df[, c(2, 3, 1)]
    # check on regblups.
    regblup <- as.vector(score.mat %*% t(gam.temp$out.list$"rotated loads"))
    pcdiff <- 100 * abs((blup.df$regblup - regblup) / regblup)
    pcdiff[is.na(pcdiff)] <- 0
    mm <- mean(pcdiff)
    # if (mm > 0.5) {
    cat("\nMean % difference between regblup coeffs and (rotated) loads * scores", mm)
    # }
    list(blup.df = blup.df, score.df = score.df)
}
update.mfxlm <- function(obj, sv, diff, boundfix = TRUE) {
    ##########
    ## Function to update from one model to another
    ## should work for any general problem with vm and ide
    ## companion to pikk.update.coloc
    ## start with R gammas
    ##
    ## we have a pikk.coloc
    ## eg move from fak no spatial to fak spatial
    ## input obj- target asreml object for spatial terms can be fa or diag,
    ## input sv- template for new model
    ## input- optional for colocated trials as the output from vcc
    ##
    vparams.tab <- sv$vparameters.table
    oldgams <- summary(obj)$varcomp
    oldgams$Component <- dimnames(oldgams)[[1]]
    oldvpars <- summary(obj, vparameters = TRUE)$vparameters
    rpar.rows <- grep(
        "\\!R",
        substring(
            vparams.tab$Component, nchar(vparams.tab$Component) - 1,
            nchar(vparams.tab$Component)
        )
    )
    rpar.rows <- seq(from = min(rpar.rows), to = nrow(vparams.tab))
    vv.vcc <- vparams.tab[rpar.rows, ]
    temp <- merge(vv.vcc, oldgams, by = "Component")
    temp <- temp[, c("Component", "component", "Constraint")]
    names(temp)[2] <- "Value"
    R.sv <- temp
    4
    #################
    # now do the genetic and non-genetic G
    # drop R terms which are in co-located
    # I will need another hook if not colocated
    # best in another function
    # one less argument
    # the way forward is to have a list which does a match on Component &
    # names(oldpars) eg here we have 'rr.*vm'
    # as the string match and then we want to replace
    # fa1 and fa2 in component with what we provide
    # e.g.
    # diff=list('rr.*vm'=list('fa1'=NULL,'fa2'=c(0,rep(.02,nsite-1))))
    # as the argument to the function
    # this takes the diff and uses the elements to replace NAs in the component in the merge
    # which are those variance parameters which are in terms which have changed
    # replace them with oldgams if is.null, else replace them with values in the
    # named vector where the name of the vector is used to identify the variance parameters in the
    # term to be replaced!!
    temp <- subset(vparams.tab, !is.element(Component, vv.vcc$Component))
    temp$index <- 1:nrow(temp)
    temp <- merge(temp, oldgams, by = "Component", all.x = TRUE)
    temp <- temp[order(temp$index), ]
    G.sv <- subset(temp, !is.na(component)) # this stores those with a match
    ####################
    # now process those which are not a match and have NA in component
    #
    temp <- subset(temp, is.na(component))
    cat("Number of non matched variance parameters", nrow(temp))
    # diff <-list('rr.*vm'=list('fa1'=NULL,'fa2'=c(0,rep(.02,nsite-1))))
    for (tt in names(diff)) {
        cat(
            "...and those which match the diff name", tt, " is",
            nrow(temp[grep(pattern = tt, x = temp$Component), ]), "\n"
        )
        if (any(sapply(diff[[tt]], is.null))) thisgam <- oldvpars[[grep(pattern = tt, x = names(oldvpars))]]
        for (pp in names(diff[[tt]])) {
            xx <- diff[[tt]][[pp]]
            print(xx)
            cat(paste(tt, pp, sep = ".*"))
            if (is.null(xx)) {
                temp[grep(pattern = paste(tt, pp, sep = ".*"), x = temp$Component), "component"] <
                    thisgam[grep(pattern = paste(tt, pp, sep = ".*"), x = names(thisgam))]
            } else {
                temp[grep(pattern = paste(tt, pp, sep = ".*"), x = temp$Component), "component"] <- xx
            }
        }
    }
    temp[is.na(temp$component), "component"] <- temp[is.na(temp$component), "Value"]
    G.sv <- rbind(G.sv, temp)
    G.sv <- G.sv[, c("Component", "component", "Constraint")]
    names(G.sv)[2] <- "Value"
    ################
    ## now tidy up the starting values (Values) which are too small for Constraint=='P'
    ## in G level terms
    if (boundfix) G.sv$Value[G.sv$Constraint == "P" & G.sv$Value < 1e-5] <- 1e-4
    return(list(G.sv = G.sv, R.sv = R.sv))
}
