#' Calculate Genetic Gain
#'
#' \strong{Predicted} (model input): the breeder's equation
#' \eqn{R = i \, r \, \sigma_G / L}, with selection intensity \eqn{i} for the
#' top \code{selection_percent}, selection accuracy \eqn{r} and cycle length
#' \eqn{L}.
#'
#' \strong{Realized} (data frame input): the slope of \code{value_col} on
#' \code{year_col}. Raw phenotypes confound genetic and environmental trend;
#' supply genotype BLUPs against year of first testing, or supply
#' \code{check_list}: checks are genetically constant, so their slope
#' estimates the environmental trend, which is subtracted.
#'
#' @param input A data frame (realized) or an \code{fa_model}/\code{asreml}
#'   object (predicted).
#' @param selection_percent Numeric in (0, 100). Percentage selected (default 5).
#' @param accuracy Numeric in (0, 1]. Selection accuracy \eqn{r} (scalar or
#'   one per group). If \code{NULL}, \eqn{r = 1} is used and a warning is
#'   raised because the result is then an upper bound.
#' @param cycle_length Numeric. Years per breeding cycle (default 1).
#' @param varcomp_term Character. For \code{asreml} input, the exact name of
#'   the genetic variance component.
#' @param year_col,value_col,genotype_col Column names for realized gain.
#' @param check_list Character vector of check genotypes (realized gain).
#'
#' @return Predicted: data frame with \code{Group}, \code{Sigma_G},
#'   \code{Intensity}, \code{Accuracy}, \code{Predicted_Gain} (per year).
#'   Realized: list with \code{Method}, \code{Slope}, \code{Env_Slope},
#'   \code{Pct_Gain_Per_Year} and \code{R2}.
#' @export
calculate_genetic_gain <- function(input, selection_percent = 5, accuracy = NULL, cycle_length = 1,
                                   varcomp_term = NULL, year_col = "Year", value_col = "Yield",
                                   genotype_col = "Genotype", check_list = NULL) {
  if (inherits(input, c("asreml", "fa_model"))) {
    if (!is.numeric(selection_percent) || selection_percent <= 0 || selection_percent >= 100) {
      cli::cli_abort("{.arg selection_percent} must be in (0, 100).")
    }
    p <- selection_percent / 100
    i_val <- stats::dnorm(stats::qnorm(1 - p)) / p

    if (inherits(input, "fa_model")) {
      Vg <- diag(input$matrices$G)
      lbl <- names(Vg) %||% rownames(input$matrices$G)
    } else {
      check_asreml_availability()
      vc <- asreml_components(input)
      if (is.null(varcomp_term) || !varcomp_term %in% names(vc)) {
        cli::cli_abort(c(
          "x" = "Supply {.arg varcomp_term}, the exact genetic variance component name.",
          "i" = "Available: {.val {names(vc)}}"
        ))
      }
      Vg <- unname(vc[varcomp_term])
      lbl <- varcomp_term
    }

    if (is.null(accuracy)) {
      cli::cli_warn("{.arg accuracy} not supplied: assuming r = 1, so {.field Predicted_Gain} is an upper bound.")
      accuracy <- 1
    }
    if (any(accuracy <= 0 | accuracy > 1)) cli::cli_abort("{.arg accuracy} must be in (0, 1].")

    return(data.frame(
      Group = lbl, Sigma_G = sqrt(Vg), Intensity = i_val, Accuracy = accuracy,
      Predicted_Gain = i_val * accuracy * sqrt(Vg) / cycle_length,
      stringsAsFactors = FALSE, row.names = NULL
    ))
  }

  if (!is.data.frame(input)) cli::cli_abort("{.arg input} must be a data frame or an asreml/fa_model object.")
  if (!all(c(year_col, value_col) %in% names(input))) {
    cli::cli_abort("Columns {.val {c(year_col, value_col)}} are required.")
  }
  df <- input[!is.na(input[[value_col]]) & !is.na(input[[year_col]]), , drop = FALSE]
  df[[year_col]] <- suppressWarnings(as.numeric(as.character(df[[year_col]])))
  if (anyNA(df[[year_col]])) cli::cli_abort("{.val {year_col}} must be numeric (or a factor of numbers).")

  form <- stats::as.formula(paste0("`", value_col, "` ~ `", year_col, "`"))
  env_slope <- 0
  if (!is.null(check_list)) {
    if (!genotype_col %in% names(df)) cli::cli_abort("Column {.val {genotype_col}} is required with {.arg check_list}.")
    is_check <- df[[genotype_col]] %in% check_list
    if (length(unique(df[[year_col]][is_check])) < 2) {
      cli::cli_abort("Checks must be observed in at least two years to estimate the environmental trend.")
    }
    env_slope <- unname(stats::coef(stats::lm(form, data = df[is_check, , drop = FALSE]))[2])
    df <- df[!is_check, , drop = FALSE]
  }
  mod <- stats::lm(form, data = df)
  slope <- unname(stats::coef(mod)[2]) - env_slope
  base <- unname(stats::predict(mod, newdata = stats::setNames(data.frame(min(df[[year_col]])), year_col)))

  list(
    Method = if (is.null(check_list)) "Realized (regression)" else "Realized (check-adjusted regression)",
    Slope = slope, Env_Slope = env_slope,
    Pct_Gain_Per_Year = 100 * slope / base,
    R2 = summary(mod)$r.squared
  )
}

#' Calculate Multi-Environment Heritability
#'
#' Plot-basis heritability per trial from a MET model with \code{diag}
#' genetic terms.
#'
#' Two designs are supported:
#' \describe{
#'   \item{\code{disjoint = TRUE} (default; the package pipeline)}{\code{gkeep}
#'     holds genotyped lines (\code{vm()} term) and \code{gdrop} holds the
#'     \emph{other}, ungenotyped lines. Each plot carries only one of the two
#'     variances, so they are never summed: \eqn{h^2 = V_a / (V_a + V_e)} for
#'     genotyped lines and \eqn{H^2_{ungenotyped} = V_{na} / (V_{na} + V_e)}.}
#'   \item{\code{disjoint = FALSE}}{Both terms apply to the same lines
#'     (additive + non-additive): \eqn{H^2 = (V_a + V_{na}) / (V_a + V_{na} + V_e)}
#'     and \eqn{h^2 = V_a / (V_a + V_{na} + V_e)}.}
#' }
#' \eqn{V_e} is the trial's spatial residual plus nugget variance.
#'
#' @param model A fitted \code{asreml} object.
#' @param gdrop_term Character. Non-additive/ungenotyped genetic term in
#'   \code{vparameters} (default \code{"studyName:gdrop"}).
#' @param gkeep_term Character. Additive/genotyped genetic term (default
#'   \code{"studyName:vm(gkeep, G.inv)"}).
#' @param env_var Character. Environment factor name (default \code{"studyName"}).
#' @param vg_scale Numeric. Multiplier for the additive component, e.g.
#'   \code{mean(diag(G)) - mean(G)} (default 1).
#' @param disjoint Logical. Do the two genetic terms cover different lines
#'   (default \code{TRUE})?
#'
#' @return A data frame with \code{Trial}, \code{Vg_add}, \code{Vg_na},
#'   \code{V_spatial}, \code{V_units} and the heritabilities above
#'   (\code{h2}, \code{H2_ungenotyped} when \code{disjoint}; otherwise
#'   \code{Vg_total}, \code{H2}, \code{h2}). Components are matched to trials
#'   by name; a trial without an ungenotyped-line component gets
#'   \code{Vg_na = 0}.
#' @export
calculate_met_h2 <- function(model, gdrop_term = "studyName:gdrop", gkeep_term = "studyName:vm(gkeep, G.inv)",
                             env_var = "studyName", vg_scale = 1, disjoint = TRUE) {
  if (!inherits(model, "asreml")) cli::cli_abort("Input must be a valid asreml model.")
  check_asreml_availability()
  varp <- summary(model, vparameters = TRUE)$vparameters
  if (is.null(varp)) cli::cli_abort("Could not extract vparameters from the model summary.")
  met_h2_from_vparameters(varp, gdrop_term, gkeep_term, env_var, vg_scale, disjoint)
}

met_h2_from_vparameters <- function(varp, gdrop_term, gkeep_term, env_var, vg_scale = 1, disjoint = TRUE) {
  if (!gkeep_term %in% names(varp)) {
    cli::cli_abort(c("x" = "Term {.val {gkeep_term}} not found.", "i" = "Available: {.val {names(varp)}}"))
  }
  clean <- function(n) sub(paste0("^", env_var, "_"), "", sub(".*!", "", n))
  va <- varp[[gkeep_term]]
  va <- stats::setNames(as.numeric(va), clean(names(va)))
  vna <- if (gdrop_term %in% names(varp)) {
    v <- varp[[gdrop_term]]
    stats::setNames(as.numeric(v), clean(names(v)))
  } else {
    cli::cli_alert_info("Term {.val {gdrop_term}} not found; Vg_na set to 0.")
    numeric(0)
  }

  trials <- names(va)
  extra <- setdiff(names(vna), trials)
  if (length(extra) > 0) cli::cli_warn("Trials with only an ungenotyped component are ignored: {.val {extra}}.")

  all_names <- names(varp)
  first_match <- function(pattern) {
    hit <- grep(pattern, all_names, fixed = TRUE)
    if (length(hit) == 0) NA_real_ else as.numeric(varp[[hit[1]]][1])
  }
  v_sp <- vapply(trials, function(tr) {
    if (tr %in% all_names) as.numeric(varp[[tr]][1]) else first_match(paste0(env_var, "_", tr, "!R"))
  }, numeric(1))
  v_un <- vapply(trials, function(tr) {
    u <- first_match(paste0("'", tr, "'):units"))
    if (is.na(u)) u <- first_match(paste0("\"", tr, "\"):units"))
    if (is.na(u)) 0 else u
  }, numeric(1))
  if (anyNA(v_sp)) cli::cli_warn("No residual variance found for {.val {trials[is.na(v_sp)]}}; heritability is NA there.")

  res <- data.frame(
    Trial = trials,
    Vg_add = unname(va) * vg_scale,
    Vg_na = unname(ifelse(trials %in% names(vna), vna[trials], 0)),
    V_spatial = unname(v_sp), V_units = unname(v_un),
    stringsAsFactors = FALSE
  )
  v_res <- res$V_spatial + res$V_units
  if (disjoint) {
    # Genotyped and ungenotyped lines are different sets: never add their variances.
    res$h2 <- res$Vg_add / (res$Vg_add + v_res)
    res$H2_ungenotyped <- ifelse(res$Vg_na > 0, res$Vg_na / (res$Vg_na + v_res), NA_real_)
    return(res[, c("Trial", "Vg_add", "Vg_na", "V_spatial", "V_units", "h2", "H2_ungenotyped")])
  }
  res$Vg_total <- res$Vg_add + res$Vg_na
  res$H2 <- res$Vg_total / (res$Vg_total + v_res)
  res$h2 <- res$Vg_add / (res$Vg_total + v_res)
  res[, c("Trial", "Vg_add", "Vg_na", "Vg_total", "V_spatial", "V_units", "H2", "h2")]
}
