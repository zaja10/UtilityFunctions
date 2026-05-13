#' Fit Bivariate ASReml Model and Extract Metrics
#'
#' Fits a bivariate ASReml model and extracts genetic and environmental 
#' correlations, along with heritability for the primary trait.
#'
#' @param data Dataframe in long format with a `Trait` column and a `Value_Scaled` column.
#' @param primary_trait Character. Name of the primary trait (e.g., "Yield").
#' @param secondary_trait Character. Name of the secondary trait.
#' @param fixed_f Formula. Fixed effects (usually `Value_Scaled ~ Trait`).
#' @param rand_f Formula. Random effects (e.g., `~ us(Trait):vm(Genotype, G.inv)`).
#' @param res_f Formula. Residual structure (e.g., `~ dsum(~ ar1(col):ar1(row) | Trait)`).
#' @param model_name Character. Name of the model.
#' @param trial_name Character. Name of the trial.
#' @param genotype_term Character. Exact name of genotype term in `vparameters` for covariance extraction.
#' @param spatial_term Character. Optional term name for spatial units.
#' @return A one-row dataframe with extracted variances and correlations.
#' @importFrom stats as.formula cov2cor
#' @export
fit_bivariate_models <- function(data, primary_trait, secondary_trait, 
                                 fixed_f = Value_Scaled ~ Trait, rand_f, res_f, 
                                 model_name = "Biv_Model", trial_name = "Trial", 
                                 genotype_term = "Trait:vm(Genotype, G.inv)", 
                                 spatial_term = "Trait:units") {
  
  check_asreml_availability()
  
  pair_data <- data[data$Trait %in% c(primary_trait, secondary_trait), ]
  sec_trait_str <- paste(secondary_trait, collapse = "_")
  full_model_name <- paste0(model_name, "_Biv_", sec_trait_str)
  res_f_str <- paste(deparse(res_f), collapse = ", ")
  
  tryCatch({
    m_biv <- asreml::asreml(
      fixed = fixed_f, 
      random = rand_f, 
      residual = res_f, 
      data = pair_data, 
      trace = FALSE, 
      maxiter = 50, 
      na.action = asreml::na.method(y = "include", x = "include")
    )
    
    if (!m_biv$converge && exists("mkConv", mode = "function")) {
      m_biv <- mkConv(m_biv)
    }
    
    vars <- m_biv$sigma2 * m_biv$vparameters
    vparams_sum <- summary(m_biv, vparameters = TRUE)$vparameters
    
    # 1. Genetic Variance (Primary Trait)
    vg_regex <- paste0(gsub("\\(", "\\\\(", genotype_term), ".*!Trait_", primary_trait, ":", primary_trait)
    Vg_matches <- vars[grep(vg_regex, names(vars))]
    Vg <- if (length(Vg_matches) > 0) Vg_matches[1] else NA
    
    # 2. Genetic Correlation (rG)
    rG_mat <- vparams_sum[[genotype_term]]
    rG <- if (!is.null(rG_mat) && is.matrix(rG_mat)) cov2cor(rG_mat)[2, 1] else NA
    
    # 3. Environmental Variance and Correlation (Ve, rE)
    Ve <- NA; rE <- NA
    if (grepl("us\\(", deparse(res_f))) {
      # Identify the residual term (heuristic)
      res_names <- names(vparams_sum)
      res_match <- res_names[grep("Trait:.*!R", res_names)]
      if (length(res_match) > 0) {
        rE_mat <- vparams_sum[[res_match[1]]]$Trait
        rE <- if (!is.null(rE_mat) && is.matrix(rE_mat)) cov2cor(rE_mat)[2, 1] else NA
      }
      
      ve_regex <- paste0(".*!R.*!Trait_", primary_trait, ":", primary_trait)
      Ve_matches <- vars[grep(ve_regex, names(vars))]
      if (length(Ve_matches) > 0) Ve <- Ve_matches[1]
    } else {
      ve_regex <- paste0("Trait_", primary_trait, "!R")
      Ve_matches <- vars[grep(ve_regex, names(vars))]
      if (length(Ve_matches) > 0) Ve <- Ve_matches[1]
    }
    
    # 4. Spatial Variance and Correlation (Vr, rR)
    Vr <- 0; rR <- NA
    if (!is.null(spatial_term) && any(grepl(spatial_term, names(vars), fixed = TRUE))) {
      vr_regex <- paste0(gsub("\\(", "\\\\(", spatial_term), ".*!Trait_", primary_trait, "(:", primary_trait, ")?$")
      Vr_matches <- vars[grep(vr_regex, names(vars))]
      if (length(Vr_matches) > 0) Vr <- Vr_matches[1]
      
      rR_mat <- vparams_sum[[spatial_term]]
      if (!is.null(rR_mat) && is.matrix(rR_mat)) rR <- cov2cor(rR_mat)[2, 1]
    }
    
    Vtotal <- Ve + Vr
    h2 <- if (!is.na(Vg) && !is.na(Vtotal)) Vg / (Vg + Vtotal / 2) else NA
    
    # Predict PEV and Acc for Primary Trait
    mean_pev <- NA; avsed <- NA; GH_c <- NA; Acc <- NA
    tryCatch({
      # Dynamically extract classify term from genotype_term (e.g., "Trait:vm(Genotype, G.inv)" -> "Trait:Genotype")
      classify_str <- gsub("vm\\(([^,]+),.*\\)", "\\1", genotype_term)
      
      # Prepare levels list for prediction (e.g., list(Trait = "Yield"))
      pred_levels <- list()
      pred_levels[["Trait"]] <- primary_trait
      
      pred <- predict(m_biv, classify = classify_str, only = paste0("us(Trait):", gsub("Trait:", "", genotype_term)), 
                      levels = pred_levels, pworkspace = 1e8, vcov = TRUE, sed = TRUE)
      
      # Find index where Trait == primary_trait in pvals
      yield_idx <- which(pred$pvals$Trait == primary_trait)
      if (length(yield_idx) > 0) {
        mean_pev <- mean(diag(pred$vcov[yield_idx, yield_idx, drop=FALSE]))
        sed_mat <- pred$sed[yield_idx, yield_idx, drop=FALSE]
        mean_sed <- mean(sed_mat[upper.tri(sed_mat)]) 
        avsed <- mean_sed^2 / 2
        GH_c <- if (!is.na(Vg)) (1 - avsed / Vg) else NA 
        Acc <- if (!is.na(GH_c) && GH_c > 0) sqrt(GH_c) else NA
      }
    }, error = function(e) {
      # Prediction failed, values remain NA
    })
    
    return(data.frame(
      Trial = trial_name, 
      Model = full_model_name, 
      rG = paste(rG, collapse = ", "), 
      rE = paste(rE, collapse = ", "), 
      rR = paste(rR, collapse = ", "), 
      Vg = paste(Vg, collapse = ", "), 
      Ve = paste(Ve, collapse = ", "), 
      Vr = paste(Vr, collapse = ", "), 
      h2 = paste(h2, collapse = ", "), 
      PEV = paste(mean_pev, collapse = ", "), 
      SED = paste(avsed, collapse = ", "), 
      GH_c = paste(GH_c, collapse = ", "), 
      Acc = paste(Acc, collapse = ", "), 
      residual = res_f_str, 
      LogLik = m_biv$loglik, 
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    message(paste("  [!] Skipping", full_model_name, "- Variance structure issue."))
    return(data.frame(
      Trial = trial_name, Model = full_model_name, rG = NA_character_, 
      rE = NA_character_, rR = NA_character_, Vg = NA_character_, 
      Ve = NA_character_, Vr = NA_character_, h2 = NA_character_, 
      PEV = NA_character_, SED = NA_character_, GH_c = NA_character_, 
      Acc = NA_character_, residual = res_f_str, LogLik = NA_real_, 
      stringsAsFactors = FALSE
    ))
  })
}
