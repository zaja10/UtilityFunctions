#' Abort unless ASReml-R is installed
#'
#' @return Invisibly \code{TRUE}.
#' @keywords internal
#' @noRd
check_asreml_availability <- function() {
  if (!requireNamespace("asreml", quietly = TRUE)) {
    cli::cli_abort(c(
      "x" = "Package {.pkg asreml} is required for this function.",
      "i" = "ASReml-R is licensed software; see {.url https://vsni.co.uk/software/asreml-r}."
    ))
  }
  invisible(TRUE)
}

#' Named vector of estimated variance components
#'
#' Reads \code{summary(model)$varcomp[, "component"]}, which is already on the
#' variance scale (unlike \code{model$vparameters}, which are relative to
#' \code{sigma2} for some parameterisations).
#'
#' @param model A fitted \code{asreml} object.
#' @return Named numeric vector.
#' @keywords internal
#' @noRd
asreml_components <- function(model) {
  vc <- summary(model)$varcomp
  stats::setNames(vc[, "component"], rownames(vc))
}

#' Environment for evaluating ASReml calls
#'
#' ASReml-R resolves \code{vm()} matrices and the model's \code{data} by
#' name in the environment that \emph{called} \code{asreml()},
#' \code{update()} or \code{predict()}. All package calls are therefore
#' evaluated in a fresh environment that inherits from the caller's (so user
#' objects such as \code{G.inv} resolve) and optionally binds \code{grm}.
#'
#' @keywords internal
#' @noRd
model_env <- function(grm = NULL, parent = parent.frame()) {
  env <- new.env(parent = parent)
  if (!is.null(grm)) assign("grm", grm, envir = env)
  env
}

#' Call update()/predict() on an asreml model inside `env`
#' @keywords internal
#' @noRd
asreml_eval <- function(fun, model, args = list(), env = parent.frame()) {
  e <- new.env(parent = env)
  assign(".uf_model", model, envir = e)
  eval(as.call(c(list(fun, quote(.uf_model)), args)), e)
}

#' Fit an ASReml model with the package's standard options
#'
#' @param fixed,random,residual Model formulae (\code{random}/\code{residual}
#'   may be \code{NULL}).
#' @param data Data frame.
#' @param maxiter Maximum REML iterations for the initial fit.
#' @param converge Logical. Run \code{\link{mkConv}} if not converged.
#' @param env Environment in which matrices named in the formulae are found
#'   (see \code{model_env()}).
#' @param ... Further arguments passed to \code{asreml::asreml()}
#'   (e.g. \code{workspace}).
#' @return A fitted \code{asreml} object.
#' @keywords internal
#' @noRd
fit_asreml <- function(fixed, random = NULL, residual = NULL, data,
                       maxiter = 50, converge = TRUE, trace = FALSE, env = parent.frame(), ...) {
  check_asreml_availability()
  args <- list(
    fixed = fixed, data = data, trace = trace, maxiter = maxiter,
    na.action = asreml::na.method(y = "include", x = "include"), ...
  )
  if (!is.null(random)) args$random <- random
  if (!is.null(residual)) args$residual <- residual
  call_env <- new.env(parent = env)
  m <- eval(as.call(c(list(quote(asreml::asreml)), args)), call_env)
  if (converge && !isTRUE(m$converge)) m <- mkConv(m, env = env)
  m
}

#' Safely ensure ASReml model convergence
#'
#' Re-runs \code{update()} on a fitted model until the log-likelihood has
#' converged and no variance parameter changed by more than \code{tol} percent
#' on the last iteration, or until \code{max_tries} updates have been made.
#'
#' @param mod A fitted \code{asreml} model object.
#' @param max_tries Integer. Maximum number of updates (default 20).
#' @param tol Numeric. Maximum acceptable absolute percentage change in any
#'   variance parameter (default 1).
#' @param env Environment in which the model's data and any \code{vm()}
#'   matrices can be found (default: the caller's).
#' @return The updated \code{asreml} object. A warning is raised if it still
#'   has not converged after \code{max_tries} updates.
#' @export
mkConv <- function(mod, max_tries = 20, tol = 1, env = parent.frame()) {
  check_asreml_availability()
  is_stable <- function(m) {
    pct <- summary(m)$varcomp[, "%ch"]
    isTRUE(m$converge) && !any(abs(pct) > tol, na.rm = TRUE)
  }
  tries <- 0L
  while (!is_stable(mod) && tries < max_tries) {
    tries <- tries + 1L
    mod <- suppressWarnings(asreml_eval(quote(stats::update), mod, env = env))
  }
  if (!is_stable(mod)) {
    cli::cli_warn(c(
      "!" = "Model not converged after {max_tries} update{?s}.",
      "i" = "Treat its estimates and likelihood with caution."
    ))
  }
  mod
}

#' Extract Standardized Residuals and Identify Outliers
#'
#' Refits a fitted \code{asreml} model with \code{aom = TRUE} (alternative
#' outlier model; Gogel et al. 2018), attaches the standardised conditional
#' residuals to the data and returns the rows whose absolute value exceeds
#' \code{threshold}.
#'
#' @param best_model A converged \code{asreml} model object.
#' @param trial_data The data frame the model was fitted to, in the same row
#'   order.
#' @param threshold Numeric. Outlier cut-off for |standardised residual|
#'   (default 4).
#' @param aom_fit Optional. An already-computed \code{update(best_model, aom = TRUE)}
#'   fit, to avoid refitting.
#' @param env Environment in which the model's data and any \code{vm()}
#'   matrices can be found (default: the caller's).
#' @return The flagged rows of \code{trial_data} with an added \code{scres}
#'   column (all original columns are kept).
#' @export
extract_asreml_outliers <- function(best_model, trial_data, threshold = 4, aom_fit = NULL, env = parent.frame()) {
  check_asreml_availability()
  if (is.null(aom_fit)) {
    aom_fit <- suppressWarnings(asreml_eval(quote(stats::update), best_model, list(aom = TRUE), env))
  }

  R <- aom_fit$aom$R
  res_col <- if ("stdCondRes" %in% colnames(R)) "stdCondRes" else 2L
  scres <- R[, res_col]

  if (length(scres) == nrow(trial_data)) {
    trial_data$scres <- scres
  } else {
    # Rows were dropped by the fit: align on row names, never on position.
    idx <- match(rownames(aom_fit$mf), rownames(trial_data))
    if (anyNA(idx)) {
      cli::cli_abort("Cannot align AOM residuals to {.arg trial_data}: row names differ from the fitted data.")
    }
    trial_data$scres <- NA_real_
    trial_data$scres[idx] <- scres
  }

  trial_data[!is.na(trial_data$scres) & abs(trial_data$scres) > threshold, , drop = FALSE]
}
