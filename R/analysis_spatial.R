#' @title Fit and Compare Spatial Models
#' @description Fits a specified asreml spatial model and returns its AIC, LogLikelihood, and convergence status.
#' @param fixed_f Formula for fixed effects
#' @param random_f Formula for random effects
#' @param residual_f Formula for residual effects
#' @param model_name Character string for the model name
#' @param df Dataframe containing the trial data
#' @return A dataframe containing Model name, AIC, LogLik, and Converged status
#' @importFrom asreml asreml na.method
#' @importFrom stats as.formula update terms
#' @export
fit_asreml_comparison <- function(fixed_f, random_f, residual_f, model_name, df) {
  tryCatch({
    m <- asreml::asreml(
      fixed = fixed_f, 
      random = random_f, 
      residual = residual_f, 
      data = df, 
      trace = FALSE, 
      maxiter = 50, 
      na.action = asreml::na.method(y = "include", x = "include")
    )
    if (!m$converge) m <- mkConv(m)
    
    aic_val  <- summary(m)$aic
    data.frame(
      Model = model_name, 
      AIC = aic_val, 
      LogLik = m$loglik, 
      Converged = m$converge, 
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(Model = model_name, AIC = Inf, LogLik = NA, Converged = FALSE)
  })
}

#' @title Run Final Best Model and Extract Metrics
#' @description Runs the best selected spatial model and extracts key quantitative genetics metrics like Vg, Ve, Vr, h2, and Accuracy.
#' @param params List containing formulas: f (fixed), r (random), res (residual)
#' @param model_name Name of the final model
#' @param df Dataframe containing the trial data
#' @param trial_name Name of the trial being analyzed
#' @param grm Genomic relationship matrix inverse
#' @return A dataframe of extracted variance components and metrics
#' @importFrom asreml asreml na.method
#' @importFrom stats predict terms
#' @export
run_final_analysis <- function(params, model_name, df, trial_name, grm) {
  m_final <- asreml::asreml(
    fixed = params$f, 
    random = params$r, 
    residual = params$res, 
    data = df, 
    trace = FALSE, 
    maxiter = 50, 
    na.action = asreml::na.method(y = "include", x = "include")
  )
  if (!m_final$converge) m_final <- mkConv(m_final)
  
  vars <- m_final$sigma2 * m_final$vparameters
  
  Vg <- vars[grep("vm\\(gkeep", names(vars))]
  if (length(Vg) == 0) Vg <- NA 
  
  Ve <- vars[grep("colNumber:rowNumber!R", names(vars))]
  Vr <- if (grepl("AxAe", model_name)) vars[grep("units", names(vars))] else 0
  
  Vtotal <- Ve + Vr
  h2 <- if(!is.na(Vg)) Vg / (Vg + Vtotal/2) else NA
  
  mean_pev <- NA; avsed <- NA; GH_c <- NA; Acc <- NA
  tryCatch({
    fixed_terms <- attr(stats::terms(params$f), "term.labels")
    terms_to_ignore <- fixed_terms[grepl("lin\\(|spl\\(", fixed_terms)]
    
    if (length(terms_to_ignore) > 0) {
      pred <- stats::predict(m_final, classify = "gkeep", only = "vm(gkeep, grm)", 
                      ignore = terms_to_ignore, vcov = TRUE, sed = TRUE)
    } else {
      pred <- stats::predict(m_final, classify = "gkeep", only = "vm(gkeep, grm)", 
                      vcov = TRUE, sed = TRUE)
    }
    
    mean_pev <- mean(diag(pred$vcov))
    avsed <- pred$avsed[2]^2 / 2
    GH_c <- if(!is.na(Vg)) (1 - avsed / Vg) else NA
    Acc <- if(!is.na(GH_c) && GH_c > 0) sqrt(GH_c) else NA
  }, error = function(e) {
  })
  
  return(data.frame(
    Trial = trial_name, Model = model_name, rG = NA, rE = NA, rR = NA, 
    Vg = Vg, Ve = Ve, Vr = Vr, h2 = h2, PEV = mean_pev, SED = avsed, 
    GH_c = GH_c, Acc = Acc, LogLik = m_final$loglik, stringsAsFactors = FALSE
  ))
}
