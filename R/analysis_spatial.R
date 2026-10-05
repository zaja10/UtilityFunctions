#' Fit and Compare a Spatial Model
#'
#' Fits one \code{asreml} model and returns its fit statistics. Errors are
#' caught and reported in the \code{Error} column rather than aborting.
#'
#' @param fixed_f Formula for fixed effects.
#' @param random_f Formula for random effects.
#' @param residual_f Formula for residual effects.
#' @param model_name Character string for the model name.
#' @param df Data frame containing the trial data.
#' @param grm Optional relationship matrix bound to the name \code{grm} in
#'   the formulae.
#' @param env Environment in which other objects named in the formulae
#'   (e.g. a user's \code{G.inv}) are found (default: the caller's).
#' @return A one-row data frame with \code{Model}, \code{AIC}, \code{LogLik},
#'   \code{Converged} and \code{Error} (\code{NA} when the fit succeeded).
#' @export
fit_asreml_comparison <- function(fixed_f, random_f, residual_f, model_name, df, grm = NULL, env = parent.frame()) {
  tryCatch(
    {
      m <- fit_asreml(fixed_f, random_f, residual_f, data = df, env = model_env(grm, env))
      data.frame(
        Model = model_name, AIC = summary(m)$aic, LogLik = m$loglik,
        Converged = isTRUE(m$converge), Error = NA_character_,
        stringsAsFactors = FALSE
      )
    },
    error = function(e) {
      data.frame(
        Model = model_name, AIC = Inf, LogLik = NA_real_, Converged = FALSE,
        Error = conditionMessage(e), stringsAsFactors = FALSE
      )
    }
  )
}

#' Compare Several Spatial Models
#'
#' All models must share the same response, data and fixed effects for their
#' REML likelihoods (and hence AIC) to be comparable; this is checked.
#'
#' @param models Named list. Each element is a list with \code{fixed},
#'   \code{random} and \code{residual} formulae (or the short names \code{f},
#'   \code{r}, \code{res}).
#' @param data Data frame.
#' @param grm,env See \code{\link{fit_asreml_comparison}}.
#' @return Data frame of fit statistics ordered by \code{\link{select_best_model}}
#'   (converged models first, then by AIC).
#' @export
compare_spatial_models <- function(models, data, grm = NULL, env = parent.frame()) {
  if (is.null(names(models)) || any(!nzchar(names(models)))) {
    cli::cli_abort("{.arg models} must be a named list.")
  }
  models <- lapply(models, normalise_model_spec)
  fixed_txt <- unique(vapply(models, function(m) paste(deparse(m$fixed), collapse = ""), ""))
  if (length(fixed_txt) > 1) {
    cli::cli_abort(c(
      "x" = "Models have different fixed effects; REML AIC is not comparable across them.",
      "i" = "Found: {.val {fixed_txt}}"
    ))
  }
  res <- do.call(rbind, Map(function(m, nm) {
    fit_asreml_comparison(m$fixed, m$random, m$residual, nm, data, grm = grm, env = env)
  }, models, names(models)))
  rownames(res) <- NULL
  rank_model_table(res)
}

#' Choose the Best Model from a Comparison Table
#'
#' Picks the lowest-AIC model among converged fits. If nothing converged, the
#' lowest-AIC model is returned with a warning.
#'
#' @param results Data frame with \code{Model}, \code{AIC} and \code{Converged}.
#' @return The name of the selected model.
#' @export
select_best_model <- function(results) {
  ranked <- rank_model_table(results)
  if (!isTRUE(as.logical(ranked$Converged[1]))) {
    cli::cli_warn("No model converged; selecting {.val {ranked$Model[1]}} by AIC anyway.")
  }
  ranked$Model[1]
}

# Converged first, then ascending AIC.
rank_model_table <- function(results) {
  conv <- as.logical(results$Converged)
  conv[is.na(conv)] <- FALSE
  results[order(!conv, results$AIC), , drop = FALSE]
}

normalise_model_spec <- function(m) {
  list(
    fixed = m$fixed %||% m$f,
    random = m$random %||% m$r,
    residual = m$residual %||% m$res
  )
}

#' Extract Variance Components, Heritability and Reliability
#'
#' @param model A fitted \code{asreml} object.
#' @param genotype_term Character. Exact name of the genetic variance
#'   component, e.g. \code{"vm(gkeep, grm)"} or \code{"Genotype"}.
#' @param error_term Character. Exact name of the residual variance component.
#'   If \code{NULL}, it is detected: the single residual section's
#'   \code{"!R"} component, multiplied by its variance parameter when the
#'   structure carries one (e.g. \code{units!R} = 1 and \code{units!units}
#'   for \code{residual = ~ idv(units)}); correlations are ignored.
#' @param nugget_term Character. Name of the nugget component (default
#'   \code{"units"}); contributes 0 when absent from the model.
#' @param n_reps Numeric. Effective number of replicates per genotype for the
#'   entry-mean heritability (default 1 = plot basis).
#' @param vg_scale Numeric. Multiplier converting the genetic variance
#'   parameter to the genetic variance of the population. For genomic models
#'   use \code{mean(diag(G)) - mean(G)} (Legarra 2016); 1 otherwise.
#' @param ignore Character vector of fixed terms to ignore in prediction
#'   (e.g. \code{lin()}/\code{spl()} trend terms).
#' @param env Environment in which the model's data and \code{vm()} matrices
#'   are found (default: the caller's).
#' @return One-row data frame: \code{Vg}, \code{Ve}, \code{Vr}, \code{h2},
#'   \code{PEV} (mean prediction error variance), \code{vBLUP} (mean variance
#'   of a difference between two BLUPs), \code{Gen_Reliability} (generalised
#'   reliability \eqn{\bar r^2_g}), \code{Acc}
#'   (\eqn{\bar r_g = \sqrt{\bar r^2_g}}), \code{LogLik}, \code{Converged}.
#' @details
#' \deqn{h^2 = \frac{s V_g}{s V_g + (V_e + V_r)/n}, \qquad
#'   \bar r^2_g = 1 - \frac{\bar v_{BLUP}}{2 s V_g}}
#' where \eqn{s} is \code{vg_scale}. \eqn{\bar r^2_g} is the generalised
#' heritability of Cullis et al. (2006), referred to here as a generalised
#' reliability: the average reliability \eqn{r^2_g} of the genotype
#' predictions (exactly so for independent genotypes and balanced data).
#' \eqn{\bar v_{BLUP}} is computed exactly as the mean of all pairwise \eqn{\mathrm{SED}^2} from the prediction error
#' variance matrix, not as the square of the mean SED (which is biased low).
#' @references Cullis, B. R., Smith, A. B. & Coombes, N. E. (2006). On the
#'   design of early generation variety trials with correlated data.
#'   \emph{JABES} 11, 381-393.
#' @export
extract_variance_components <- function(model, genotype_term, error_term = NULL,
                                        nugget_term = "units", n_reps = 1,
                                        vg_scale = 1, ignore = NULL, env = parent.frame()) {
  if (!inherits(model, "asreml")) cli::cli_abort("{.arg model} must be an asreml object.")
  check_asreml_availability()

  vc <- asreml_components(model)
  if (!genotype_term %in% names(vc)) {
    cli::cli_abort(c(
      "x" = "Genetic term {.val {genotype_term}} not found among variance components.",
      "i" = "Available: {.val {names(vc)}}"
    ))
  }
  Vg <- unname(vc[genotype_term]) * vg_scale

  if (is.null(error_term)) {
    Ve <- residual_variance(vc)
  } else {
    if (!error_term %in% names(vc)) cli::cli_abort("Residual term {.val {error_term}} not found.")
    Ve <- unname(vc[error_term])
  }
  Vr <- if (!is.null(nugget_term) && nugget_term %in% names(vc)) unname(vc[nugget_term]) else 0

  h2 <- Vg / (Vg + (Ve + Vr) / n_reps)

  classify <- sub("^vm\\(\\s*([^,]+),.*$", "\\1", genotype_term)
  pev <- vblup <- rel <- Acc <- NA_real_
  pred <- tryCatch(
    {
      args <- list(classify = classify, only = genotype_term, vcov = TRUE)
      if (length(ignore) > 0) args$ignore <- ignore
      suppressWarnings(asreml_eval(quote(stats::predict), model, args, env))
    },
    error = function(e) {
      cli::cli_warn("Prediction for {.val {classify}} failed: {conditionMessage(e)}")
      NULL
    }
  )
  if (!is.null(pred)) {
    C <- as.matrix(pred$vcov)
    pev <- mean(diag(C))
    vblup <- mean_pairwise_vd(C)
    rel <- 1 - vblup / (2 * Vg)
    Acc <- if (is.finite(rel) && rel > 0) sqrt(rel) else NA_real_
  }

  data.frame(
    Vg = Vg, Ve = Ve, Vr = Vr, h2 = h2, PEV = pev, vBLUP = vblup,
    Gen_Reliability = rel, Acc = Acc, LogLik = model$loglik,
    Converged = isTRUE(model$converge), stringsAsFactors = FALSE
  )
}

# Residual variance of a single-section residual model. ASReml reports a
# scale "<section>!R" plus any structure parameters "<section>!<name>"; for
# idv()-type structures the scale is fixed at 1 and the variance is a
# separate parameter. Correlation parameters ("!cor", "!corr") are ignored.
residual_variance <- function(vc) {
  r_terms <- grep("!R$", names(vc), value = TRUE)
  if (length(r_terms) != 1) {
    cli::cli_abort("Found {length(r_terms)} residual sections ({.val {r_terms}}); specify {.arg error_term}.")
  }
  section <- sub("!R$", "", r_terms)
  params <- names(vc)[startsWith(names(vc), paste0(section, "!")) & names(vc) != r_terms]
  var_params <- params[!grepl("!(cor|corr|corgh|pow)[^!]*$", params)]
  if (length(var_params) > 1) {
    cli::cli_abort("Residual section {.val {section}} has several variance parameters ({.val {var_params}}); specify {.arg error_term}.")
  }
  unname(vc[r_terms] * if (length(var_params) == 1) vc[var_params] else 1)
}

# Mean over i < j of Var(u_i - u_j) = C_ii + C_jj - 2 C_ij.
mean_pairwise_vd <- function(C) {
  n <- nrow(C)
  if (n < 2) {
    return(NA_real_)
  }
  2 * (n * sum(diag(C)) - sum(C)) / (n * (n - 1))
}

#' Run Final Best Model and Extract Metrics
#'
#' Fits the selected spatial model and extracts variance components,
#' heritability and accuracy via \code{\link{extract_variance_components}}.
#'
#' @param params List of formulas: \code{f} (fixed), \code{r} (random),
#'   \code{res} (residual).
#' @param model_name Name of the model (used only as a label).
#' @param df Data frame containing the trial data.
#' @param trial_name Name of the trial being analysed.
#' @param grm Genomic relationship matrix (inverse) referenced as \code{grm}
#'   inside the formulas.
#' @param genotype_term Character. Genetic variance component name
#'   (default \code{"vm(gkeep, grm)"}).
#' @param vg_scale See \code{\link{extract_variance_components}}.
#' @return A one-row data frame of metrics, including \code{n_reps}, the
#'   harmonic mean number of observed plots per genotype used for \code{h2}.
#' @export
run_final_analysis <- function(params, model_name, df, trial_name, grm,
                               genotype_term = "vm(gkeep, grm)", vg_scale = 1) {
  params <- normalise_model_spec(params)
  env <- model_env(grm, parent.frame())
  m_final <- fit_asreml(params$fixed, params$random, params$residual, data = df, env = env)

  fixed_terms <- attr(stats::terms(params$fixed), "term.labels")
  ignore <- fixed_terms[grepl("^(lin|spl)\\(", fixed_terms)]

  response <- all.vars(params$fixed)[1]
  geno_col <- sub("^vm\\(\\s*([^,]+),.*$", "\\1", genotype_term)
  n_reps <- harmonic_reps(df[[geno_col]][!is.na(df[[response]])])

  vcs <- extract_variance_components(
    m_final, genotype_term,
    n_reps = n_reps, vg_scale = vg_scale, ignore = ignore, env = env
  )
  cbind(
    data.frame(Trial = trial_name, Model = model_name, rG = NA_real_, rE = NA_real_, rR = NA_real_, stringsAsFactors = FALSE),
    vcs,
    n_reps = n_reps
  )
}

# Harmonic mean of replicate counts over genotypes that were observed.
harmonic_reps <- function(geno) {
  counts <- table(geno[!is.na(geno)])
  counts <- counts[counts > 0]
  if (length(counts) == 0) {
    return(NA_real_)
  }
  length(counts) / sum(1 / counts)
}
