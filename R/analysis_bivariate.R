#' Fit Bivariate ASReml Model and Extract Metrics
#'
#' Fits a bivariate model for a primary and a secondary trait and extracts
#' the genetic, residual and nugget correlations together with the primary
#' trait's heritability, prediction error variance, generalised reliability
#' (\eqn{\bar r^2_g}; see \code{\link{extract_variance_components}}) and accuracy.
#'
#' @param data Long-format data frame with a \code{Trait} column and a
#'   response column (default \code{Value_Scaled}).
#' @param primary_trait Character. Primary trait (e.g. \code{"Yield"}).
#' @param secondary_trait Character. Secondary trait.
#' @param fixed_f Formula. Fixed effects (default \code{Value_Scaled ~ Trait}).
#' @param rand_f Formula. Random effects, e.g.
#'   \code{~ us(Trait):vm(Genotype, G.inv)}.
#' @param res_f Formula. Residual, e.g.
#'   \code{~ dsum(~ ar1(col):ar1(row) | Trait)} or \code{~ units:us(Trait)}.
#' @param model_name Character. Model label.
#' @param trial_name Character. Trial label.
#' @param genotype_term Character. Variance-component prefix of the genetic
#'   term as ASReml names it, e.g. \code{"Trait:vm(Genotype, G.inv)"} for
#'   \code{us(Trait):vm(Genotype, G.inv)}.
#' @param spatial_term Character. Prefix of a trait-specific nugget term
#'   (e.g. \code{"Trait:units"}), or \code{NULL}.
#' @param grm Optional relationship matrix bound to the name \code{grm} in the
#'   formulas.
#' @return A one-row data frame with numeric columns \code{rG}, \code{rE},
#'   \code{rR}, \code{Vg}, \code{Ve}, \code{Vr}, \code{h2}, \code{PEV},
#'   \code{vBLUP}, \code{Gen_Reliability}, \code{Acc}, \code{LogLik}, \code{Converged},
#'   plus \code{residual} (formula text) and \code{Error} (\code{NA} on success).
#' @export
fit_bivariate_models <- function(data, primary_trait, secondary_trait,
                                 fixed_f = Value_Scaled ~ Trait, rand_f, res_f,
                                 model_name = "Biv_Model", trial_name = "Trial",
                                 genotype_term = "Trait:vm(Genotype, G.inv)",
                                 spatial_term = "Trait:units", grm = NULL) {
  check_asreml_availability()
  caller <- parent.frame()
  if (!"Trait" %in% names(data)) cli::cli_abort("{.arg data} must have a {.field Trait} column.")

  traits <- c(primary_trait, secondary_trait)
  full_name <- paste0(model_name, "_Biv_", paste(secondary_trait, collapse = "_"))
  res_txt <- paste(deparse(res_f), collapse = " ")

  out <- data.frame(
    Trial = trial_name, Model = full_name, rG = NA_real_, rE = NA_real_, rR = NA_real_,
    Vg = NA_real_, Ve = NA_real_, Vr = 0, h2 = NA_real_, PEV = NA_real_, vBLUP = NA_real_,
    Gen_Reliability = NA_real_, Acc = NA_real_, residual = res_txt, LogLik = NA_real_, Converged = FALSE,
    Error = NA_character_, stringsAsFactors = FALSE
  )

  tryCatch(
    {
      d <- as.data.frame(data[data$Trait %in% traits, , drop = FALSE])
      d$Trait <- factor(as.character(d$Trait), levels = traits)
      d <- d[order(d$Trait), , drop = FALSE]
      rownames(d) <- NULL

      env <- model_env(grm, caller)
      m <- fit_asreml(fixed_f, rand_f, res_f, data = d, env = env)
      vc <- asreml_components(m)
      p <- primary_trait
      s <- secondary_trait

      g <- us_block(vc, genotype_term, p, s)
      if (is.null(g)) cli::cli_abort("Genetic us() block {.val {genotype_term}} not found.")
      out$Vg <- g[["var_p"]]
      out$rG <- g[["cor"]]

      prefixes <- setdiff(unique(sub("!.*$", "", names(vc))), c(genotype_term, spatial_term))
      r_blocks <- Filter(Negate(is.null), lapply(prefixes, us_block, vc = vc, p = p, s = s))
      if (length(r_blocks) > 0) {
        out$Ve <- r_blocks[[1]][["var_p"]]
        out$rE <- r_blocks[[1]][["cor"]]
      } else {
        ve <- vc[paste0("Trait_", p, "!R")]
        if (is.na(ve)) cli::cli_abort("Residual variance for {.val {p}} not found.")
        out$Ve <- unname(ve)
      }

      if (!is.null(spatial_term)) {
        r <- us_block(vc, spatial_term, p, s)
        if (!is.null(r)) {
          out$Vr <- r[["var_p"]]
          out$rR <- r[["cor"]]
        }
      }

      geno_col <- sub("^.*vm\\(\\s*([^,]+),.*$", "\\1", genotype_term)
      n_reps <- harmonic_reps(d[[geno_col]][d$Trait == p & !is.na(d[[all.vars(fixed_f)[1]]])])
      out$h2 <- out$Vg / (out$Vg + (out$Ve + out$Vr) / n_reps)

      only_term <- grep(geno_col, attr(stats::terms(rand_f), "term.labels"), value = TRUE, fixed = TRUE)[1]
      pred <- tryCatch(
        suppressWarnings(asreml_eval(quote(stats::predict), m, list(
          classify = paste0("Trait:", geno_col), only = only_term,
          levels = list(Trait = p), vcov = TRUE
        ), env)),
        error = function(e) {
          cli::cli_warn("Prediction failed: {conditionMessage(e)}")
          NULL
        }
      )
      if (!is.null(pred)) {
        keep <- which(pred$pvals$Trait == p)
        C <- as.matrix(pred$vcov)[keep, keep, drop = FALSE]
        out$PEV <- mean(diag(C))
        out$vBLUP <- mean_pairwise_vd(C)
        out$Gen_Reliability <- 1 - out$vBLUP / (2 * out$Vg)
        out$Acc <- if (is.finite(out$Gen_Reliability) && out$Gen_Reliability > 0) sqrt(out$Gen_Reliability) else NA_real_
      }
      out$LogLik <- m$loglik
      out$Converged <- isTRUE(m$converge)
      out
    },
    error = function(e) {
      cli::cli_alert_danger("Skipping {full_name}: {conditionMessage(e)}")
      out$Error <- conditionMessage(e)
      out
    }
  )
}

# Variances and correlation from an ASReml us() block named
# "<prefix>!Trait_<a>:<b>".
us_block <- function(vc, prefix, p, s) {
  key <- function(a, b) paste0(prefix, "!Trait_", a, ":", b)
  vp <- vc[key(p, p)]
  vs <- vc[key(s, s)]
  cv <- vc[key(s, p)]
  if (is.na(cv)) cv <- vc[key(p, s)]
  if (anyNA(c(vp, vs, cv))) {
    return(NULL)
  }
  c(var_p = unname(vp), var_s = unname(vs), cor = unname(cv / sqrt(vp * vs)))
}
