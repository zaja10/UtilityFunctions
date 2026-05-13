#' Fit and Compare Spatial ASReml Models
#'
#' Fits a list of ASReml models and returns a comparison dataframe sorted by AIC.
#'
#' @param model_list A named list. Each element must be a list containing `fixed`, `random`, and `residual` formula objects.
#' @param data A dataframe containing the trial data.
#' @return A dataframe comparing the models (Model, AIC, LogLik, Converged), sorted by AIC.
#' @importFrom stats as.formula
#' @export
compare_spatial_models <- function(model_list, data) {
  check_asreml_availability()
  
  # Ensure mkConv is available
  if (!exists("mkConv", mode = "function")) {
    stop("mkConv helper function is required but not found.")
  }

  results <- purrr::imap_dfr(model_list, function(params, model_name) {
    tryCatch({
      # Dynamically evaluate the formulas within the package environment
      m <- asreml::asreml(
        fixed = params$fixed,
        random = params$random,
        residual = params$residual,
        data = data,
        trace = FALSE,
        maxiter = 50,
        na.action = asreml::na.method(y = "include", x = "include")
      )
      
      if (!m$converge) {
        m <- mkConv(m)
      }
      
      aic_val <- summary(m)$aic
      
      data.frame(
        Model = model_name,
        AIC = aic_val,
        LogLik = m$loglik,
        Converged = m$converge,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        Model = model_name,
        AIC = Inf,
        LogLik = NA,
        Converged = FALSE,
        stringsAsFactors = FALSE
      )
    })
  })
  
  results <- results[order(results$AIC), ]
  return(results)
}

#' Extract Variance Components and Key Metrics
#'
#' Extracts genetic, environmental, and spatial variances from an ASReml model,
#' and computes heritability, mean PEV, and accuracy.
#'
#' @param model A converged ASReml model.
#' @param model_name Character. Name of the model (for output tracking).
#' @param trial_name Character. Name of the trial (for output tracking).
#' @param genotype_term Character. The exact name of the genotype term in the vparameters (e.g., "vm(Genotype, G.inv)").
#' @param error_term Character. The term name for the residual error (e.g., "colNumber:rowNumber!R").
#' @param spatial_term Character. The term name for the spatial error (e.g., "units"). Defaults to NULL.
#' @return A one-row dataframe with the extracted metrics.
#' @export
extract_variance_components <- function(model, model_name = "Model", trial_name = "Trial", 
                                        genotype_term = NULL, error_term = NULL, spatial_term = NULL) {
  
  check_asreml_availability()
  
  if (!inherits(model, "asreml")) stop("Input must be an asreml object.")
  
  vars <- model$sigma2 * model$vparameters
  
  # Extract Genetic Variance
  Vg <- NA
  if (!is.null(genotype_term)) {
    match_idx <- grep(genotype_term, names(vars), fixed = TRUE)
    if (length(match_idx) > 0) Vg <- vars[match_idx[1]]
  }
  
  # Extract Error Variance
  Ve <- NA
  if (!is.null(error_term)) {
    match_idx <- grep(error_term, names(vars), fixed = TRUE)
    if (length(match_idx) > 0) Ve <- vars[match_idx[1]]
  }
  
  # Extract Spatial Variance (e.g., units)
  Vr <- 0
  if (!is.null(spatial_term)) {
    match_idx <- grep(spatial_term, names(vars), fixed = TRUE)
    if (length(match_idx) > 0) Vr <- vars[match_idx[1]]
  }
  
  Vtotal <- Ve + Vr
  h2 <- if (!is.na(Vg) && !is.na(Vtotal)) Vg / (Vg + Vtotal / 2) else NA
  
  # Extract PEV and SED if classification is possible
  mean_pev <- NA; avsed <- NA; GH_c <- NA; Acc <- NA
  
  if (!is.na(Vg) && !is.null(genotype_term)) {
    tryCatch({
      # Parse classify term from the genotype term string 
      # e.g. "vm(Genotype, G.inv)" -> "Genotype"
      classify_str <- gsub("vm\\(([^,]+),.*\\)", "\\1", genotype_term)
      
      pred <- predict(model, classify = classify_str, only = genotype_term, vcov = TRUE, sed = TRUE)
      mean_pev <- mean(diag(pred$vcov))
      avsed <- pred$avsed[2]^2 / 2
      GH_c <- (1 - avsed / Vg)
      if (GH_c > 0) Acc <- sqrt(GH_c)
    }, error = function(e) {
      # Ignore prediction errors silently
    })
  }
  
  return(data.frame(
    Trial = trial_name, 
    Model = model_name, 
    Vg = Vg, 
    Ve = Ve, 
    Vr = Vr, 
    h2 = h2, 
    PEV = mean_pev, 
    SED = avsed, 
    GH_c = GH_c, 
    Acc = Acc, 
    LogLik = model$loglik, 
    stringsAsFactors = FALSE
  ))
}

#' Extract Standardized Residuals and Identify Outliers
#'
#' Updates a fitted ASReml model to compute the Average Outlier Measure (AOM),
#' extracts standardized conditional residuals, and filters for potential outliers.
#'
#' @param best_model A converged ASReml model object.
#' @param trial_data The original dataframe used for the trial.
#' @param threshold Numeric. Absolute value cutoff for identifying outliers (default 4).
#' @param cols_to_keep Character vector. Column names to retain in the output summary.
#' @return A dataframe containing only the flagged outliers and their residuals (`scres`).
#' @export
extract_asreml_outliers <- function(best_model, trial_data, threshold = 4, cols_to_keep = NULL) {
  check_asreml_availability()
  
  # Update model to compute outlier statistics
  out.fm <- suppressWarnings(update(best_model, aom = TRUE))
  
  trial_data$scres <- NA
  
  # Map residuals back using valid row indices
  valid_rows <- as.numeric(rownames(out.fm$mf))
  trial_data$scres[valid_rows] <- out.fm$aom$R[, 2]
  
  # Filter outliers
  outliers <- trial_data[!is.na(trial_data$scres) & abs(trial_data$scres) > threshold, ]
  
  if (nrow(outliers) > 0 && !is.null(cols_to_keep)) {
    cols_to_keep <- intersect(cols_to_keep, names(outliers))
    outliers <- subset(outliers, select = c(cols_to_keep, "scres"))
  }
  
  return(outliers)
}
