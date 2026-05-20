#' @title Safely ensure ASReml model convergence
#' @description This function checks the percentage change in variance components and updates the model if it hasn't converged, up to a maximum number of tries.
#' @param mod A fitted asreml model object
#' @return A converged asreml model object
#' @importFrom asreml asreml summary update
#' @export
mkConv <- function(mod) {
  pctchg <- summary(mod)$varcomp[, c('%ch')]
  tries <- 1
  while(any(pctchg > 1, na.rm = TRUE) & tries < 20) {
    tries <- tries + 1
    mod <- suppressWarnings(update(mod))
    pctchg <- summary(mod)$varcomp[, c('%ch')]
  }
  return(mod)
}

#' @title Extract Standardized Residuals and Identify Outliers
#' @description This function updates a fitted asreml model to calculate the Average Outlier Measure (AOM), extracts the standardized conditional residuals, and filters for potential outliers based on a defined threshold.
#' @param best_model A converged asreml model object
#' @param trial_data The original dataframe used for the trial
#' @param threshold Numeric value defining the outlier cutoff (default is 4)
#' @return A dataframe containing only the flagged outliers
#' @importFrom asreml update.asreml
#' @export
extract_asreml_outliers <- function(best_model, trial_data, threshold = 4) {
  out.fm <- suppressWarnings(asreml::update.asreml(best_model, aom = TRUE))
  
  trial_data$scres <- NA
  valid_rows <- as.numeric(rownames(out.fm$mf))
  trial_data$scres[valid_rows] <- out.fm$aom$R[, 2]
  
  outliers <- trial_data[!is.na(trial_data$scres) & abs(trial_data$scres) > threshold, ]
  
  # Use base subsetting for package-safe column extraction
  outliers_clean <- outliers[, c("studyName", "colNumber", "rowNumber", 
                                 "replicate", "germplasmName", 
                                 "grain_yield", "scres")]
  
  return(outliers_clean)
}
