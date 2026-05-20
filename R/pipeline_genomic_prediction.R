#' @title Run Univariate Spatial Pipeline
#' @description Iterates over all trials in the dataset to select the best spatial model, calculate metrics, and output diagnostic plots.
#' @param data Processed dataframe containing phenotype data. Must include `studyName`, `colNumber`, `rowNumber`, `grain_yield`, `gkeep`, and `gdrop`
#' @param grm Genomic relationship matrix inverse
#' @param output_dir Path to the output directory (default "outputs")
#' @return NULL (Writes output to tables/ and diagnostics/ subdirectories)
#' @importFrom dplyr filter mutate arrange across any_of bind_rows everything %>%
#' @importFrom rlang .data
#' @importFrom purrr imap_dfr
#' @importFrom asreml asreml na.method update.asreml
#' @importFrom ASExtras4 met.asreml
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot
#' @importFrom utils write.csv read.csv
#' @importFrom stats as.formula
#' @export
run_univariate_pipeline <- function(data, grm, output_dir = "outputs") {
  tables_dir <- file.path(output_dir, "tables")
  diagnostics_dir <- file.path(output_dir, "diagnostics")
  
  if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
  if (!dir.exists(diagnostics_dir)) dir.create(diagnostics_dir, recursive = TRUE)
  
  trial_list <- unique(as.character(data$studyName))
  
  for (trial_name_str in trial_list) {
    message(sprintf("\n==================================================\nProcessing Trial: %s\n==================================================", trial_name_str))
    
    selection_file <- file.path(tables_dir, paste0(trial_name_str, "_Spatial_Model_Selection.csv"))
    metrics_file   <- file.path(tables_dir, paste0(trial_name_str, "_Final_Metrics.csv"))
    plot_file_path <- file.path(diagnostics_dir, paste0(trial_name_str, "_Diagnostics.pdf"))
    
    if (file.exists(plot_file_path)) {
      message("  [✓] Trial fully processed (Models and Diagnostics complete). Skipping.")
      next
    }
    
    tryCatch({
      trial_data <- data %>% dplyr::filter(.data$studyName == trial_name_str)
      if(nrow(trial_data) == 0) stop("No data found for this trial.")
      
      analysis_data <- trial_data %>%
        dplyr::mutate(Yield = as.numeric(scale(.data$grain_yield))) %>% 
        dplyr::mutate(dplyr::across(c("rowNumber", "colNumber", dplyr::any_of("blockNumber")), as.factor)) %>%
        dplyr::arrange(.data$colNumber, .data$rowNumber)
      
      analysis_data$gkeep <- factor(analysis_data$gkeep, levels = attr(grm, "rowNames"))
      analysis_data$gdrop <- as.factor(analysis_data$gdrop)
      
      has_blocks <- "blockNumber" %in% names(analysis_data) && length(unique(analysis_data$blockNumber[!is.na(analysis_data$blockNumber)])) > 1
      has_ungenotyped <- sum(!is.na(analysis_data$gdrop)) > 0
      
      rand_terms <- c("vm(gkeep, grm)", "colNumber", "rowNumber")
      if (has_ungenotyped) rand_terms <- c(rand_terms, "gdrop")
      if (has_blocks) rand_terms <- c(rand_terms, "blockNumber")
      
      base_rand_str <- paste("~", paste(rand_terms, collapse = " + "))
      
      model_list <- list(
        "RCD"  = list(f = Yield ~ 1, r = stats::as.formula(base_rand_str), res = ~ colNumber:rowNumber),
        "AxA"  = list(f = Yield ~ 1, r = stats::as.formula(base_rand_str), res = ~ ar1(colNumber):ar1(rowNumber)),
        "AxAe" = list(f = Yield ~ 1, r = stats::as.formula(paste(base_rand_str, "+ units")), res = ~ ar1(colNumber):ar1(rowNumber))
      )
      
      if (file.exists(selection_file) && file.exists(metrics_file)) {
        message("  [✓] Model selection files found. Loading best model...")
        saved_results <- utils::read.csv(selection_file)
        best_model_name <- saved_results$Model[1]
        best_params <- model_list[[best_model_name]]
      } else {
        message("  Fitting Univariate Spatial Models for comparison...")
        asreml_results <- purrr::imap_dfr(model_list, ~ fit_asreml_comparison(.x$f, .x$r, .x$res, .y, analysis_data)) %>% dplyr::arrange(.data$AIC)
        
        best_model_name <- asreml_results$Model[1]
        best_params <- model_list[[best_model_name]]
        
        univ_results <- run_final_analysis(best_params, paste0(best_model_name, "_Univ_Yield"), analysis_data, trial_name_str, grm)
        all_model_results <- dplyr::bind_rows(univ_results %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.character)))
        
        utils::write.csv(asreml_results, file = selection_file, row.names = FALSE)
        utils::write.csv(all_model_results, file = metrics_file, row.names = FALSE)
        message(paste("  [✓] Model selection saved. Best Model:", best_model_name))
      }
      
      message(paste("  Refitting best model (", best_model_name, ") for diagnostics..."))
      best_m_fit <- asreml::asreml(fixed = best_params$f, random = best_params$r, 
                           residual = best_params$res, data = analysis_data, 
                           trace = FALSE, maxiter = 50, 
                           na.action = asreml::na.method(y = "include", x = "include"))
      
      if (!best_m_fit$converge) best_m_fit <- mkConv(best_m_fit)
      
      message("  Extracting outliers and saving diagnostic plots...")
      trial_outliers <- extract_asreml_outliers(best_model = best_m_fit, 
                                                trial_data = analysis_data, 
                                                threshold = 4)
      
      if (nrow(trial_outliers) > 0) {
        utils::write.csv(trial_outliers, file = file.path(tables_dir, paste0(trial_name_str, "_Flagged_Outliers.csv")), row.names = FALSE)
        message(paste("  [!] Found", nrow(trial_outliers), "outliers. Saved to review list."))
      } else {
        message("  [✓] No outliers > 4 detected.")
      }
      
      out.fm <- suppressWarnings(asreml::update.asreml(best_m_fit, aom = TRUE))
      
      grDevices::pdf(file = plot_file_path, width = 8, height = 6)
      tryCatch({
        graphics::plot(out.fm, res = 'stdCond', main = paste("Standardised Residuals:", trial_name_str))
      }, error = function(e) message("    [!] Could not plot residuals: ", e$message))
      
      tryCatch({
        ASExtras4::met.asreml(best_m_fit) 
      }, error = function(e) message("    [!] Could not plot variogram: ", e$message))
      grDevices::dev.off()
      message("  [✓] Diagnostic plots saved to folder.")
      
    }, error = function(e) {
      message(paste("  [!] Error processing trial:", e$message))
    })
  }
}

#' @title Generate Diagnostic Review Template
#' @description Reads the selection files and creates an Excel template for human review of models and outliers.
#' @param selection_files_dir Path to the directory containing spatial model selection CSVs
#' @param full_data Dataframe of the full trial data
#' @param output_path Path to save the Excel template
#' @return NULL (Writes an Excel file)
#' @importFrom dplyr filter %>%
#' @importFrom rlang .data
#' @importFrom purrr map_dfr
#' @importFrom writexl write_xlsx
#' @importFrom utils read.csv
#' @export
generate_diagnostic_template <- function(selection_files_dir, full_data, output_path) {
  message("\n==================================================\nGenerating Trial Terms Template for Diagnostic Review\n==================================================")
  
  selection_files <- list.files(selection_files_dir, pattern = "_Spatial_Model_Selection\\.csv$", full.names = TRUE)
  
  build_template_row <- function(file_path, full_data) {
    file_name <- basename(file_path)
    trial_name <- gsub("_Spatial_Model_Selection\\.csv$", "", file_name)
    
    model_data <- utils::read.csv(file_path)
    best_model <- as.character(model_data$Model[1])
    
    trial_data <- full_data %>% dplyr::filter(.data$studyName == trial_name)
    
    has_blocks <- "blockNumber" %in% names(trial_data) && length(unique(trial_data$blockNumber[!is.na(trial_data$blockNumber)])) > 1
    has_ungenotyped <- sum(!is.na(trial_data$gdrop)) > 0
    
    rand_terms <- c("vm(gkeep, grm)", "colNumber", "rowNumber")
    if (has_ungenotyped) rand_terms <- c(rand_terms, "gdrop")
    if (has_blocks) rand_terms <- c(rand_terms, "blockNumber")
    if (!is.na(best_model) && best_model == "AxAe") rand_terms <- c(rand_terms, "units")
    
    base_fixed <- "Yield ~ 1"
    base_random <- paste("~", paste(rand_terms, collapse = " + "))
    
    data.frame(
      studyName = trial_name,
      best_model = best_model,
      fixed_terms = base_fixed,
      random_terms = base_random,
      remove_outliers = FALSE,  
      stringsAsFactors = FALSE
    )
  }
  
  if (length(selection_files) > 0) {
    template_df <- purrr::map_dfr(selection_files, ~build_template_row(.x, full_data))
    writexl::write_xlsx(template_df, output_path)
    message(paste("  [✓] Excel template generated at:", output_path))
  } else {
    message("  [!] No spatial model selection files found. Run the main loop first.")
  }
}

#' @title Run Modified Trials
#' @description Re-evaluates trials based on the human-modified template, optionally masking outliers and refitting custom formulas.
#' @param data Processed dataframe containing phenotype data
#' @param template_path Path to the human-reviewed Diagnostic Review Template Excel file
#' @param grm Genomic relationship matrix inverse
#' @param output_dir Path to the output directory (default "outputs")
#' @return NULL (Updates existing output files and creates final summary)
#' @importFrom readxl read_xlsx
#' @importFrom dplyr filter mutate left_join select arrange bind_rows across everything %>%
#' @importFrom rlang .data
#' @importFrom asreml asreml na.method
#' @importFrom utils read.csv write.csv
#' @importFrom writexl write_xlsx
#' @importFrom purrr map_dfr
#' @importFrom stats as.formula
#' @export
run_modified_trials <- function(data, template_path, grm, output_dir = "outputs") {
  message("\n==================================================\nEvaluating Human-Modified (MOD) Spatial Models\n==================================================")
  
  tables_dir <- file.path(output_dir, "tables")
  
  if (!file.exists(template_path)) {
    stop("Excel template not found. Please run the template generation step first.")
  }
  
  mod_data <- readxl::read_xlsx(template_path)
  
  for (i in 1:nrow(mod_data)) {
    trial_name_str <- mod_data$studyName[i]
    base_model     <- mod_data$best_model[i]
    rem_outliers   <- as.logical(mod_data$remove_outliers[i])
    
    message(paste("\nTesting MOD for Trial:", trial_name_str))
    
    trial_data <- data %>% dplyr::filter(.data$studyName == trial_name_str)
    if(nrow(trial_data) == 0) next
    
    if (isTRUE(rem_outliers)) {
      outlier_file <- file.path(tables_dir, paste0(trial_name_str, "_Flagged_Outliers.csv"))
      if (file.exists(outlier_file)) {
        outliers_df <- utils::read.csv(outlier_file)
        message(paste("  [!] Removing", nrow(outliers_df), "flagged outliers..."))
        
        trial_data <- trial_data %>%
          dplyr::mutate(col_match = as.character(.data$colNumber), row_match = as.character(.data$rowNumber)) %>%
          dplyr::left_join(
            outliers_df %>% 
              dplyr::mutate(col_match = as.character(.data$colNumber), row_match = as.character(.data$rowNumber)) %>%
              dplyr::select("col_match", "row_match", is_outlier = "scres"),
            by = c("col_match", "row_match")
          ) %>%
          dplyr::mutate(grain_yield = ifelse(!is.na(.data$is_outlier), NA, .data$grain_yield)) %>%
          dplyr::select(-"col_match", -"row_match", -"is_outlier")
      } else {
        message("  [!] remove_outliers is TRUE, but no Outlier CSV was found. Proceeding normally.")
      }
    }
    
    analysis_data <- trial_data %>%
      dplyr::mutate(Yield = as.numeric(scale(.data$grain_yield))) %>% 
      dplyr::mutate(dplyr::across(c("rowNumber", "colNumber", dplyr::any_of("blockNumber")), as.factor)) %>%
      dplyr::arrange(.data$colNumber, .data$rowNumber)
    
    analysis_data$gkeep <- factor(analysis_data$gkeep, levels = attr(grm, "rowNames"))
    analysis_data$gdrop <- as.factor(analysis_data$gdrop)
    
    tryCatch({
      f_form <- stats::as.formula(mod_data$fixed_terms[i])
      r_form <- stats::as.formula(mod_data$random_terms[i])
      
      if (grepl("AxA", base_model)) {
        res_form <- ~ ar1(colNumber):ar1(rowNumber)
      } else {
        res_form <- ~ colNumber:rowNumber
      }
      
      m_mod <- asreml::asreml(
        fixed = f_form, 
        random = r_form, 
        residual = res_form, 
        data = analysis_data, 
        trace = FALSE, 
        maxiter = 50, 
        na.action = asreml::na.method(y = "include", x = "include")
      )
      if (!m_mod$converge) m_mod <- mkConv(m_mod)
      
      sel_file <- file.path(tables_dir, paste0(trial_name_str, "_Spatial_Model_Selection.csv"))
      if (file.exists(sel_file)) {
        sel_df <- utils::read.csv(sel_file) %>% dplyr::filter(.data$Model != "MOD")
        mod_row <- data.frame(
          Model = "MOD", 
          AIC = summary(m_mod)$aic, 
          LogLik = m_mod$loglik, 
          Converged = m_mod$converge,
          stringsAsFactors = FALSE
        )
        sel_df <- dplyr::bind_rows(sel_df, mod_row) %>% dplyr::arrange(.data$AIC)
        utils::write.csv(sel_df, sel_file, row.names = FALSE)
        message(paste("  [✓] Updated AIC selection. MOD AIC:", round(mod_row$AIC, 2)))
      }
      
      met_file <- file.path(tables_dir, paste0(trial_name_str, "_Final_Metrics.csv"))
      if (file.exists(met_file)) {
        met_df <- utils::read.csv(met_file) %>% 
          dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>% 
          dplyr::filter(!grepl("MOD", .data$Model))
        
        mod_params <- list(f = f_form, r = r_form, res = res_form)
        mod_metrics <- run_final_analysis(mod_params, "MOD_Univ_Yield", analysis_data, trial_name_str, grm)
        mod_metrics <- mod_metrics %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
        met_df <- dplyr::bind_rows(met_df, mod_metrics)
        
        utils::write.csv(met_df, met_file, row.names = FALSE)
        message("  [✓] Updated variance components and metrics.")
      }
    }, error = function(e) {
      message(paste("  [!] Error evaluating MOD model:", e$message))
    })
  }
  
  message("\n==================================================\nGenerating Final Model Summary Report\n==================================================")
  selection_files <- list.files(tables_dir, pattern = "_Spatial_Model_Selection\\.csv$", full.names = TRUE)
  
  get_best_model_info <- function(file_path, full_data, template_df) {
    trial_name <- gsub("_Spatial_Model_Selection\\.csv$", "", basename(file_path))
    model_data <- utils::read.csv(file_path)
    best_row <- model_data[1, ]
    final_mod_name <- best_row$Model
    
    trial_template <- template_df %>% dplyr::filter(.data$studyName == trial_name)
    base_model_from_template <- trial_template$best_model[1]
    
    outliers_rem <- if(nrow(trial_template) > 0 && !is.na(trial_template$remove_outliers[1])) {
      as.logical(trial_template$remove_outliers[1])
    } else {
      FALSE
    }
    
    if (final_mod_name == "MOD") {
      f_term <- trial_template$fixed_terms[1]
      r_term <- trial_template$random_terms[1]
      res_term <- ifelse(grepl("AxA", base_model_from_template), 
                         "~ ar1(colNumber):ar1(rowNumber)", 
                         "~ colNumber:rowNumber")
    } else {
      trial_data_sub <- full_data %>% dplyr::filter(.data$studyName == trial_name)
      has_blocks <- "blockNumber" %in% names(trial_data_sub) && length(unique(trial_data_sub$blockNumber[!is.na(trial_data_sub$blockNumber)])) > 1
      has_ungenotyped <- sum(!is.na(trial_data_sub$gdrop)) > 0
      
      f_term <- "Yield ~ 1"
      rand_list <- c("vm(gkeep, grm)", "colNumber", "rowNumber")
      if (has_ungenotyped) rand_list <- c(rand_list, "gdrop")
      if (has_blocks) rand_list <- c(rand_list, "blockNumber")
      if (final_mod_name == "AxAe") rand_list <- c(rand_list, "units")
      
      r_term <- paste("~", paste(rand_list, collapse = " + "))
      res_term <- ifelse(grepl("AxA", final_mod_name), 
                         "~ ar1(colNumber):ar1(rowNumber)", 
                         "~ colNumber:rowNumber")
    }
    
    data.frame(
      studyName = trial_name,
      final_model = final_mod_name,
      outliers_removed = outliers_rem,
      fixed_terms = f_term,
      random_terms = r_term,
      residual_terms = res_term,
      AIC = best_row$AIC,
      Converged = best_row$Converged,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(selection_files) > 0) {
    final_summary <- purrr::map_dfr(selection_files, ~get_best_model_info(.x, data, mod_data)) %>%
      dplyr::arrange(.data$studyName)
    output_path <- file.path(output_dir, "Final_Model_Summary.xlsx")
    writexl::write_xlsx(final_summary, output_path)
    message(paste("  [✓] Final summary report generated at:", output_path))
  } else {
    message("  [!] No spatial model selection files found.")
  }
}

#' @title Construct and Run 1-Stage MET Model (Phase 1)
#' @description Pre-processes the dataset masking defined outliers and dynamically constructs fixed, random, and residual equations to fit the step-up stage 1 MET model.
#' @param data Dataframe containing phenotype data
#' @param grm Genomic relationship matrix inverse
#' @param final_summary_path Path to the Final Model Summary Excel file
#' @param tables_dir Path to the output tables directory containing outlier CSVs (default "outputs/tables")
#' @param save_dir Path to directory to save the output RDS (default "outputs")
#' @return A fitted asreml MET model object
#' @importFrom readxl read_xlsx
#' @importFrom dplyr filter mutate case_when across any_of arrange %>% row_number
#' @importFrom rlang .data
#' @importFrom asreml asreml na.method
#' @importFrom stats as.formula
#' @importFrom utils read.csv
#' @export
run_met_pipeline <- function(data, grm, final_summary_path, tables_dir = "outputs/tables", save_dir = "outputs") {
  message("\n==================================================\nConstructing 1-Stage MET Model (Step-Up Phase 1)\n==================================================")
  
  if (!file.exists(final_summary_path)) {
    stop("Final_Model_Summary.xlsx not found. Please ensure the univariate pipeline is complete.")
  }
  
  final_summary <- readxl::read_xlsx(final_summary_path)
  
  met_data <- data %>%
    dplyr::filter(.data$studyName %in% final_summary$studyName) %>%
    dplyr::mutate(
      Yield = .data$grain_yield,
      rowNum_num = as.numeric(as.character(.data$rowNumber)),
      colNum_num = as.numeric(as.character(.data$colNumber))
    )
  
  for (i in 1:nrow(final_summary)) {
    t_name <- final_summary$studyName[i]
    rem_outliers <- as.logical(final_summary$outliers_removed[i])
    
    if (isTRUE(rem_outliers)) {
      outlier_file <- file.path(tables_dir, paste0(t_name, "_Flagged_Outliers.csv"))
      if (file.exists(outlier_file)) {
        outliers <- utils::read.csv(outlier_file)
        message(paste("  [-] Masking", nrow(outliers), "outliers for trial:", t_name))
        
        met_data <- met_data %>%
          dplyr::mutate(
            Yield = dplyr::case_when(
              .data$studyName == t_name & paste(.data$colNumber, .data$rowNumber, sep="_") %in% 
                paste(outliers$colNumber, outliers$rowNumber, sep="_") ~ NA_real_,
              TRUE ~ .data$Yield
            )
          )
      }
    }
  }
  
  met_data <- met_data %>%
    dplyr::mutate(
      units = as.factor(dplyr::row_number()), 
      dplyr::across(dplyr::any_of(c("studyName", "blockNumber", "rowNumber", "colNumber")), as.factor)
    )
  
  met_data$gkeep <- factor(met_data$gkeep, levels = attr(grm, "rowNames"))
  met_data$gdrop <- as.factor(met_data$gdrop)
  
  message("  Building dynamic MET equations...")
  
  fixed_str <- "Yield ~ studyName"
  trials_lin_row <- final_summary$studyName[grepl("lin\\(rowNumber\\)", final_summary$fixed_terms)]
  trials_lin_col <- final_summary$studyName[grepl("lin\\(colNumber\\)", final_summary$fixed_terms)]
  
  for (t in trials_lin_row) fixed_str <- paste0(fixed_str, " + at(studyName, '", t, "'):rowNum_num")
  for (t in trials_lin_col) fixed_str <- paste0(fixed_str, " + at(studyName, '", t, "'):colNum_num")
  
  fixed_eq <- stats::as.formula(fixed_str)
  
  random_str <- "~ diag(studyName):vm(gkeep, grm) + diag(studyName):gdrop"
  
  for (i in 1:nrow(final_summary)) {
    t <- final_summary$studyName[i]
    r_term <- final_summary$random_terms[i]
    
    if (grepl("blockNumber", r_term)) random_str <- paste0(random_str, " + at(studyName, '", t, "'):blockNumber")
    if (grepl("colNumber", r_term))   random_str <- paste0(random_str, " + at(studyName, '", t, "'):colNumber")
    if (grepl("rowNumber", r_term))   random_str <- paste0(random_str, " + at(studyName, '", t, "'):rowNumber")
    if (grepl("units", r_term))       random_str <- paste0(random_str, " + at(studyName, '", t, "'):units")
  }
  
  random_eq <- stats::as.formula(random_str)
  residual_eq <- ~ dsum(~ ar1(colNumber):ar1(rowNumber) | studyName)
  
  message("\n  Executing Base MET Model (diag GxE)...")
  message(paste("  Fixed Eq:", fixed_str))
  
  met_data <- met_data %>% dplyr::arrange(.data$studyName, .data$colNumber, .data$rowNumber)
  
  tryCatch({
    met_step1 <- asreml::asreml(
      fixed = fixed_eq, 
      random = random_eq, 
      residual = residual_eq, 
      data = met_data,
      workspace = "46gb", 
      pworkspace = "32gb",
      trace = TRUE, 
      maxiter = 50, 
      na.action = asreml::na.method(y = "include", x = "include")
    )
    
    if (!met_step1$converge) {
      message("  [!] Initial MET model did not converge. Attempting mkConv update...")
      met_step1 <- mkConv(met_step1)
    }
    
    message("  [✓] MET Step 1 Model (DIAG) converged successfully.")
    saveRDS(met_step1, file.path(save_dir, "MET_Step1_Diag_Model.rds"))
    
    return(met_step1)
  }, error = function(e) {
    message("  [!] Error executing MET model: ", e$message)
    return(NULL)
  })
}
