# -----------------------------------------------------------------------------
# Shared helpers for the univariate -> review -> MET pipeline.
#
# Expected data columns: studyName, colNumber, rowNumber, gkeep (genotyped ID,
# NA when ungenotyped), gdrop (ungenotyped ID, NA when genotyped), optional
# blockNumber, plus the trait column (default "grain_yield").
# -----------------------------------------------------------------------------

PIPELINE_COLS <- c("studyName", "colNumber", "rowNumber", "gkeep", "gdrop")
BASE_MODELS <- c("RCD", "AxA", "AxAe")
BASE_FIXED <- "Yield ~ 1"

validate_pipeline_data <- function(data, trait) {
  missing <- setdiff(c(PIPELINE_COLS, trait), names(data))
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "x" = "Pipeline data is missing required column{?s}: {.val {missing}}.",
      "i" = "Required: {.val {c(PIPELINE_COLS, trait)}}."
    ))
  }
  if (!is.numeric(data[[trait]])) cli::cli_abort("Trait column {.val {trait}} must be numeric.")
  invisible(TRUE)
}

# Factor whose levels follow numeric order when the values are numeric.
numeric_factor <- function(x) {
  x_chr <- as.character(x)
  num <- suppressWarnings(as.numeric(x_chr))
  if (!anyNA(num[!is.na(x_chr)])) {
    factor(x_chr, levels = as.character(sort(unique(num[!is.na(num)]))))
  } else {
    factor(x_chr)
  }
}

has_blocks <- function(d) {
  "blockNumber" %in% names(d) && length(unique(stats::na.omit(d$blockNumber))) > 1
}

has_ungenotyped <- function(d) any(!is.na(d$gdrop))

# Prepare one trial: standardise the trait, set factors, align gkeep to the GRM.
prepare_trial_data <- function(trial_data, trait, grm, scale_trait = TRUE) {
  d <- as.data.frame(trial_data)
  d$Yield <- if (scale_trait) as.numeric(scale(d[[trait]])) else d[[trait]]
  d$colNumber <- numeric_factor(d$colNumber)
  d$rowNumber <- numeric_factor(d$rowNumber)
  if ("blockNumber" %in% names(d)) d$blockNumber <- factor(d$blockNumber)

  ids <- grm_ids(grm)
  raw <- as.character(d$gkeep)
  d$gkeep <- factor(raw, levels = ids)
  lost <- sum(!is.na(raw) & is.na(d$gkeep))
  n_geno <- sum(!is.na(raw))
  if (n_geno > 0 && lost == n_geno) {
    cli::cli_abort("None of the {n_geno} genotyped plot{?s} match the GRM IDs; check {.field gkeep} against the GRM.")
  }
  if (lost > 0) cli::cli_warn("{lost} of {n_geno} genotyped plot{?s} have IDs absent from the GRM and will be treated as missing.")
  d$gdrop <- factor(d$gdrop)

  d <- pad_grid(d)
  d <- d[order(d$colNumber, d$rowNumber), , drop = FALSE]
  rownames(d) <- NULL
  d
}

# Separable AR1 x AR1 residuals need every row x column position present:
# add empty plots (missing response and genotype) for any gaps.
pad_grid <- function(d) {
  d$.padded <- FALSE
  dup <- duplicated(d[, c("colNumber", "rowNumber")])
  if (any(dup)) cli::cli_abort("{sum(dup)} duplicated row/column position{?s} in trial {.val {d$studyName[1]}}.")
  grid <- expand.grid(colNumber = levels(d$colNumber), rowNumber = levels(d$rowNumber), stringsAsFactors = FALSE)
  key <- paste(d$colNumber, d$rowNumber)
  gaps <- grid[!paste(grid$colNumber, grid$rowNumber) %in% key, , drop = FALSE]
  if (nrow(gaps) == 0) {
    return(d)
  }
  pad <- d[rep(NA_integer_, nrow(gaps)), , drop = FALSE]
  pad$studyName <- d$studyName[1]
  pad$colNumber <- factor(gaps$colNumber, levels = levels(d$colNumber))
  pad$rowNumber <- factor(gaps$rowNumber, levels = levels(d$rowNumber))
  pad$.padded <- TRUE
  cli::cli_alert_info("Padded {nrow(gaps)} empty plot position{?s} to complete the grid.")
  rbind(d, pad)
}

spatial_random_terms <- function(d, model_name) {
  terms <- c("vm(gkeep, grm)", "colNumber", "rowNumber")
  if (has_ungenotyped(d)) terms <- c(terms, "gdrop")
  if (has_blocks(d)) terms <- c(terms, "blockNumber")
  if (identical(model_name, "AxAe")) terms <- c(terms, "units")
  paste("~", paste(terms, collapse = " + "))
}

spatial_residual <- function(model_name) {
  if (model_name %in% c("AxA", "AxAe")) "~ ar1(colNumber):ar1(rowNumber)" else "~ colNumber:rowNumber"
}

build_model_list <- function(d) {
  stats::setNames(lapply(BASE_MODELS, function(nm) {
    list(
      fixed = stats::as.formula(BASE_FIXED, env = globalenv()),
      random = stats::as.formula(spatial_random_terms(d, nm), env = globalenv()),
      residual = stats::as.formula(spatial_residual(nm), env = globalenv())
    )
  }), BASE_MODELS)
}

# Canonical form of a formula string for change detection.
formula_key <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) {
    return("")
  }
  f <- stats::as.formula(x)
  tl <- sort(gsub("\\s+", "", attr(stats::terms(f), "term.labels")))
  lhs <- if (length(f) == 3) deparse(f[[2]]) else ""
  paste(lhs, paste(tl, collapse = "+"), sep = "~")
}

trial_file <- function(dir, trial, suffix) file.path(dir, paste0(trial, suffix))

read_selection <- function(path) {
  sel <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!"Comparable" %in% names(sel)) sel$Comparable <- TRUE
  sel
}

#' Choose the final model for a trial
#'
#' A reviewer-requested MOD model that changes the data (outlier removal) or
#' the fixed effects cannot be compared to the base models by AIC; if it
#' converged it is adopted as an explicit reviewer decision. Otherwise the
#' best converged model by AIC among comparable models is chosen.
#' @return list(model, basis)
#' @keywords internal
#' @noRd
choose_final_model <- function(sel) {
  comparable <- if (is.null(sel$Comparable)) rep(TRUE, nrow(sel)) else as.logical(sel$Comparable)
  comparable[is.na(comparable)] <- TRUE
  mod <- sel[sel$Model == "MOD" & !comparable, , drop = FALSE]
  if (nrow(mod) > 0 && isTRUE(as.logical(mod$Converged[1]))) {
    return(list(model = "MOD", basis = "Reviewer override (not AIC-comparable)"))
  }
  list(model = select_best_model(sel[comparable, , drop = FALSE]), basis = "AIC")
}

write_diagnostics_pdf <- function(path, aom_fit, model, trial) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = tmp, width = 8, height = 6)
  dev <- grDevices::dev.cur()
  on.exit(if (dev %in% grDevices::dev.list()) grDevices::dev.off(dev), add = TRUE)

  graphics::plot(aom_fit, res = "stdCond", main = paste("Standardised Residuals:", trial))
  if (requireNamespace("ASExtras4", quietly = TRUE)) {
    tryCatch(
      getExportedValue("ASExtras4", "met.asreml")(model),
      error = function(e) cli::cli_warn("Variogram failed for {.val {trial}}: {conditionMessage(e)}")
    )
  }
  grDevices::dev.off(dev)
  if (!file.copy(tmp, path, overwrite = TRUE)) cli::cli_abort("Could not write {.file {path}}.")
  unlink(tmp)
  invisible(path)
}

report_status <- function(status, what) {
  failed <- status[status$Status == "error", , drop = FALSE]
  if (nrow(failed) > 0) {
    cli::cli_warn(c(
      "!" = "{what}: {nrow(failed)} trial{?s} failed.",
      stats::setNames(paste0(failed$Trial, ": ", failed$Message), rep("x", nrow(failed)))
    ))
  }
  invisible(status)
}

status_row <- function(trial, status, model = NA_character_, message = NA_character_) {
  data.frame(Trial = trial, Status = status, Model = model, Message = message, stringsAsFactors = FALSE)
}

#' Run Univariate Spatial Pipeline
#'
#' For every trial: fits RCD, AR1xAR1 and AR1xAR1 + nugget spatial models,
#' selects the best converged model by AIC, extracts metrics, flags AOM
#' outliers and writes residual diagnostics. A trial is skipped when its
#' diagnostics PDF already exists; the PDF is only written once every other
#' step has succeeded, so interrupted trials are re-run.
#'
#' @param data Phenotype data frame with columns \code{studyName},
#'   \code{colNumber}, \code{rowNumber}, \code{gkeep}, \code{gdrop}, optional
#'   \code{blockNumber}, and the trait.
#' @param grm Genomic relationship matrix inverse (ASReml sparse format, see
#'   \code{\link{prepare_asreml_grm}}, or a dense matrix with row names).
#' @param output_dir Output directory (default \code{"outputs"}).
#' @param trait Character. Trait column (default \code{"grain_yield"}). It is
#'   standardised within trial before fitting.
#' @param outlier_threshold Numeric. |standardised residual| cut-off (default 4).
#' @return Invisibly, a data frame with one row per trial: \code{Trial},
#'   \code{Status} (\code{"done"}, \code{"skipped"}, \code{"error"}),
#'   \code{Model} and \code{Message}. Failures are also summarised in a warning.
#' @export
run_univariate_pipeline <- function(data, grm, output_dir = "outputs",
                                    trait = "grain_yield", outlier_threshold = 4) {
  check_asreml_availability()
  validate_pipeline_data(data, trait)
  tables_dir <- file.path(output_dir, "tables")
  diagnostics_dir <- file.path(output_dir, "diagnostics")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ASExtras4", quietly = TRUE)) {
    cli::cli_inform(c("i" = "{.pkg ASExtras4} not installed; variograms will be omitted from diagnostics."))
  }

  trials <- unique(as.character(data$studyName))
  status <- lapply(trials, function(trial) {
    cli::cli_h2("Trial: {trial}")
    selection_file <- trial_file(tables_dir, trial, "_Spatial_Model_Selection.csv")
    metrics_file <- trial_file(tables_dir, trial, "_Final_Metrics.csv")
    outlier_file <- trial_file(tables_dir, trial, "_Flagged_Outliers.csv")
    plot_file <- trial_file(diagnostics_dir, trial, "_Diagnostics.pdf")

    if (file.exists(plot_file)) {
      cli::cli_alert_success("Already processed; skipping.")
      return(status_row(trial, "skipped"))
    }

    tryCatch(
      {
        d <- prepare_trial_data(data[data$studyName == trial, , drop = FALSE], trait, grm)
        if (nrow(d) == 0) cli::cli_abort("No data for this trial.")
        models <- build_model_list(d)
        env <- model_env(grm)

        if (file.exists(selection_file) && file.exists(metrics_file)) {
          sel <- read_selection(selection_file)
          best <- select_best_model(sel[sel$Model %in% BASE_MODELS, , drop = FALSE])
          cli::cli_alert_info("Using saved model selection: {.val {best}}.")
        } else {
          sel <- compare_spatial_models(models, d, grm = grm)
          sel$Comparable <- TRUE
          best <- select_best_model(sel)
          metrics <- run_final_analysis(models[[best]], paste0(best, "_Univ_Yield"), d, trial, grm)
          utils::write.csv(sel, selection_file, row.names = FALSE)
          utils::write.csv(metrics, metrics_file, row.names = FALSE)
          cli::cli_alert_success("Best model: {.val {best}}.")
        }

        p <- models[[best]]
        fit <- fit_asreml(p$fixed, p$random, p$residual, data = d, env = env)
        aom_fit <- suppressWarnings(asreml_eval(quote(stats::update), fit, list(aom = TRUE), env))
        outliers <- extract_asreml_outliers(fit, d, threshold = outlier_threshold, aom_fit = aom_fit)
        if (nrow(outliers) > 0) {
          utils::write.csv(outliers, outlier_file, row.names = FALSE)
          cli::cli_alert_warning("{nrow(outliers)} outlier{?s} flagged.")
        } else if (file.exists(outlier_file)) {
          unlink(outlier_file)
        }

        write_diagnostics_pdf(plot_file, aom_fit, fit, trial)
        status_row(trial, "done", best)
      },
      error = function(e) {
        cli::cli_alert_danger("Failed: {conditionMessage(e)}")
        status_row(trial, "error", message = conditionMessage(e))
      }
    )
  })
  report_status(do.call(rbind, status), "Univariate pipeline")
}

#' Generate Diagnostic Review Template
#'
#' Writes an Excel template, one row per trial, listing the selected model's
#' fixed, random and residual terms for a reviewer to edit. Set
#' \code{remove_outliers} to \code{TRUE} to mask flagged outliers.
#'
#' @param selection_files_dir Directory containing the
#'   \code{*_Spatial_Model_Selection.csv} files.
#' @param full_data Data frame of the full trial data.
#' @param output_path Path of the Excel file to write.
#' @return Invisibly, the template data frame (\code{NULL} if no selection
#'   files were found).
#' @export
generate_diagnostic_template <- function(selection_files_dir, full_data, output_path) {
  files <- list.files(selection_files_dir, pattern = "_Spatial_Model_Selection\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    cli::cli_warn("No spatial model selection files found in {.path {selection_files_dir}}. Run {.fn run_univariate_pipeline} first.")
    return(invisible(NULL))
  }

  template <- do.call(rbind, lapply(files, function(path) {
    trial <- sub("_Spatial_Model_Selection\\.csv$", "", basename(path))
    sel <- read_selection(path)
    best <- select_best_model(sel[sel$Model %in% BASE_MODELS, , drop = FALSE])
    d <- full_data[full_data$studyName == trial, , drop = FALSE]
    data.frame(
      studyName = trial,
      best_model = best,
      fixed_terms = BASE_FIXED,
      random_terms = spatial_random_terms(d, best),
      residual_terms = spatial_residual(best),
      remove_outliers = FALSE,
      stringsAsFactors = FALSE
    )
  }))

  writexl::write_xlsx(template, output_path)
  cli::cli_alert_success("Template written to {.file {output_path}}.")
  invisible(template)
}

#' Run Modified Trials
#'
#' Refits each trial whose template row differs from its selected model
#' ("MOD"). A MOD that only changes random/residual structure is compared to
#' the base models by AIC. A MOD that removes outliers or changes the fixed
#' effects is fitted to different data or a different fixed model, so its
#' REML likelihood is not comparable; it is recorded with
#' \code{Comparable = FALSE} and, if it converged, adopted as an explicit
#' reviewer decision. A final summary workbook is then written.
#'
#' @param data Phenotype data frame (see \code{\link{run_univariate_pipeline}}).
#' @param template_path Path to the reviewed template (Excel).
#' @param grm Genomic relationship matrix inverse.
#' @param output_dir Output directory (default \code{"outputs"}).
#' @param trait Character. Trait column (default \code{"grain_yield"}).
#' @return Invisibly, a per-trial status data frame.
#' @export
run_modified_trials <- function(data, template_path, grm, output_dir = "outputs", trait = "grain_yield") {
  check_asreml_availability()
  validate_pipeline_data(data, trait)
  if (!file.exists(template_path)) {
    cli::cli_abort("Template {.file {template_path}} not found. Run {.fn generate_diagnostic_template} first.")
  }
  tables_dir <- file.path(output_dir, "tables")
  template <- as.data.frame(readxl::read_xlsx(template_path))
  if (!"residual_terms" %in% names(template)) template$residual_terms <- NA_character_

  status <- lapply(seq_len(nrow(template)), function(i) {
    row <- template[i, , drop = FALSE]
    trial <- as.character(row$studyName)
    cli::cli_h2("MOD for trial: {trial}")
    sel_file <- trial_file(tables_dir, trial, "_Spatial_Model_Selection.csv")
    met_file <- trial_file(tables_dir, trial, "_Final_Metrics.csv")

    tryCatch(
      {
        raw <- data[data$studyName == trial, , drop = FALSE]
        if (nrow(raw) == 0) cli::cli_abort("No data for this trial.")
        base_model <- as.character(row$best_model)
        rem <- isTRUE(as.logical(row$remove_outliers))
        fixed_str <- as.character(row$fixed_terms)
        random_str <- as.character(row$random_terms)
        residual_str <- if (is.na(row$residual_terms)) spatial_residual(base_model) else as.character(row$residual_terms)

        same_fixed <- formula_key(fixed_str) == formula_key(BASE_FIXED)
        unchanged <- !rem && same_fixed &&
          formula_key(random_str) == formula_key(spatial_random_terms(raw, base_model)) &&
          formula_key(residual_str) == formula_key(spatial_residual(base_model))

        if (unchanged) {
          drop_mod_rows(sel_file, met_file)
          cli::cli_alert_info("Template unchanged; no MOD model fitted.")
          return(status_row(trial, "unchanged", base_model))
        }

        if (rem) raw <- mask_outliers(raw, trait, trial_file(tables_dir, trial, "_Flagged_Outliers.csv"))
        d <- prepare_trial_data(raw, trait, grm)
        params <- list(
          fixed = stats::as.formula(fixed_str, env = globalenv()),
          random = stats::as.formula(random_str, env = globalenv()),
          residual = stats::as.formula(residual_str, env = globalenv())
        )

        comparable <- !rem && same_fixed
        mod_row <- fit_asreml_comparison(params$fixed, params$random, params$residual, "MOD", d, grm = grm)
        mod_row$Comparable <- comparable
        if (!is.na(mod_row$Error)) cli::cli_abort("MOD fit failed: {mod_row$Error}")

        if (file.exists(sel_file)) {
          sel <- read_selection(sel_file)
          sel <- sel[sel$Model != "MOD", , drop = FALSE]
          sel <- rank_model_table(dplyr::bind_rows(sel, mod_row))
          utils::write.csv(sel, sel_file, row.names = FALSE)
        }
        if (file.exists(met_file)) {
          met <- utils::read.csv(met_file, stringsAsFactors = FALSE)
          met <- met[!grepl("^MOD", met$Model), , drop = FALSE]
          mod_metrics <- run_final_analysis(params, "MOD_Univ_Yield", d, trial, grm)
          utils::write.csv(dplyr::bind_rows(met, mod_metrics), met_file, row.names = FALSE)
        }
        if (!comparable) {
          cli::cli_alert_info("MOD changes the data or fixed effects: AIC not comparable, recorded as reviewer choice.")
        }
        cli::cli_alert_success("MOD AIC: {round(mod_row$AIC, 2)} (converged: {mod_row$Converged}).")
        status_row(trial, "done", "MOD")
      },
      error = function(e) {
        cli::cli_alert_danger("Failed: {conditionMessage(e)}")
        status_row(trial, "error", message = conditionMessage(e))
      }
    )
  })
  status <- report_status(do.call(rbind, status), "Modified trials")

  write_final_summary(tables_dir, data, template, output_dir)
  invisible(status)
}

drop_mod_rows <- function(sel_file, met_file) {
  if (file.exists(sel_file)) {
    sel <- read_selection(sel_file)
    if (any(sel$Model == "MOD")) utils::write.csv(sel[sel$Model != "MOD", , drop = FALSE], sel_file, row.names = FALSE)
  }
  if (file.exists(met_file)) {
    met <- utils::read.csv(met_file, stringsAsFactors = FALSE)
    if (any(grepl("^MOD", met$Model))) utils::write.csv(met[!grepl("^MOD", met$Model), , drop = FALSE], met_file, row.names = FALSE)
  }
}

# Set the trait to NA for plots listed in a trial's outlier file.
mask_outliers <- function(trial_data, trait, outlier_file) {
  if (!file.exists(outlier_file)) {
    cli::cli_warn("{.arg remove_outliers} is TRUE but {.file {basename(outlier_file)}} does not exist; nothing masked.")
    return(trial_data)
  }
  out <- utils::read.csv(outlier_file, stringsAsFactors = FALSE)
  key <- paste(trial_data$colNumber, trial_data$rowNumber, sep = "_")
  hit <- key %in% unique(paste(out$colNumber, out$rowNumber, sep = "_"))
  trial_data[[trait]][hit] <- NA
  cli::cli_alert_info("Masked {sum(hit)} outlier plot{?s}.")
  trial_data
}

write_final_summary <- function(tables_dir, data, template, output_dir) {
  files <- list.files(tables_dir, pattern = "_Spatial_Model_Selection\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    cli::cli_warn("No spatial model selection files found; final summary not written.")
    return(invisible(NULL))
  }
  summary_df <- do.call(rbind, lapply(files, function(path) {
    trial <- sub("_Spatial_Model_Selection\\.csv$", "", basename(path))
    sel <- read_selection(path)
    choice <- choose_final_model(sel)
    best_row <- sel[sel$Model == choice$model, , drop = FALSE][1, ]
    tmpl <- template[template$studyName == trial, , drop = FALSE]

    if (choice$model == "MOD") {
      f_term <- tmpl$fixed_terms[1]
      r_term <- tmpl$random_terms[1]
      res_term <- if (is.na(tmpl$residual_terms[1])) spatial_residual(tmpl$best_model[1]) else tmpl$residual_terms[1]
      outliers <- isTRUE(as.logical(tmpl$remove_outliers[1]))
    } else {
      f_term <- BASE_FIXED
      r_term <- spatial_random_terms(data[data$studyName == trial, , drop = FALSE], choice$model)
      res_term <- spatial_residual(choice$model)
      outliers <- FALSE
    }
    data.frame(
      studyName = trial, final_model = choice$model, selection_basis = choice$basis,
      outliers_removed = outliers, fixed_terms = f_term, random_terms = r_term,
      residual_terms = res_term, AIC = best_row$AIC, Converged = best_row$Converged,
      stringsAsFactors = FALSE
    )
  }))
  summary_df <- summary_df[order(summary_df$studyName), , drop = FALSE]
  out <- file.path(output_dir, "Final_Model_Summary.xlsx")
  writexl::write_xlsx(summary_df, out)
  cli::cli_alert_success("Final summary written to {.file {out}}.")
  invisible(summary_df)
}

#' Construct and Run the 1-Stage MET Model (Step-Up Phase 1)
#'
#' Builds a single-stage MET model from the per-trial decisions in the final
#' summary: trial-specific linear row/column trends, blocks, row/column and
#' nugget effects, and a residual structure that matches each trial's
#' selected model (AR1xAR1 or independent), with \code{diag} genetic
#' variance across trials.
#'
#' @param data Phenotype data frame (see \code{\link{run_univariate_pipeline}}).
#' @param grm Genomic relationship matrix inverse.
#' @param final_summary_path Path to \code{Final_Model_Summary.xlsx}.
#' @param tables_dir Directory with outlier CSVs (default \code{"outputs/tables"}).
#' @param save_dir Directory to save the fitted model RDS (default \code{"outputs"}).
#' @param trait Character. Trait column (default \code{"grain_yield"}).
#' @param trace Logical. Print ASReml iteration trace (default \code{TRUE}).
#' @param ... Further arguments for \code{asreml::asreml()}, e.g.
#'   \code{workspace = "8gb"}.
#' @return The fitted \code{asreml} model. Errors are raised, not swallowed.
#' @export
run_met_pipeline <- function(data, grm, final_summary_path, tables_dir = "outputs/tables",
                             save_dir = "outputs", trait = "grain_yield", trace = TRUE, ...) {
  check_asreml_availability()
  validate_pipeline_data(data, trait)
  if (!file.exists(final_summary_path)) {
    cli::cli_abort("{.file {final_summary_path}} not found. Run {.fn run_modified_trials} first.")
  }
  fs <- as.data.frame(readxl::read_xlsx(final_summary_path))
  trials <- as.character(fs$studyName)

  met <- as.data.frame(data[data$studyName %in% trials, , drop = FALSE])
  for (i in which(as.logical(fs$outliers_removed) %in% TRUE)) {
    t <- trials[i]
    rows <- met$studyName == t
    met[rows, ] <- mask_outliers(met[rows, , drop = FALSE], trait, trial_file(tables_dir, t, "_Flagged_Outliers.csv"))
  }

  met <- do.call(rbind, lapply(split(met, factor(met$studyName, levels = trials)), function(td) {
    td$colNumber <- numeric_factor(td$colNumber)
    td$rowNumber <- numeric_factor(td$rowNumber)
    td <- pad_grid(td)
    td$colNumber <- as.character(td$colNumber)
    td$rowNumber <- as.character(td$rowNumber)
    td
  }))

  met$Yield <- met[[trait]]
  met$rowNum_num <- as.numeric(as.character(met$rowNumber))
  met$colNum_num <- as.numeric(as.character(met$colNumber))
  met$studyName <- factor(met$studyName, levels = trials)
  met$colNumber <- numeric_factor(met$colNumber)
  met$rowNumber <- numeric_factor(met$rowNumber)
  if ("blockNumber" %in% names(met)) met$blockNumber <- factor(met$blockNumber)
  met$gkeep <- factor(as.character(met$gkeep), levels = grm_ids(grm))
  met$gdrop <- factor(met$gdrop)
  met <- met[order(met$studyName, met$colNumber, met$rowNumber), , drop = FALSE]
  met$units <- factor(seq_len(nrow(met)))
  rownames(met) <- NULL

  eqs <- build_met_equations(fs, has_gdrop = has_ungenotyped(met))
  cli::cli_inform(c("i" = "Fixed: {eqs$fixed}", "i" = "Residual: {eqs$residual}"))

  model <- tryCatch(
    fit_asreml(
      stats::as.formula(eqs$fixed, env = globalenv()),
      stats::as.formula(eqs$random, env = globalenv()),
      stats::as.formula(eqs$residual, env = globalenv()),
      data = met, trace = trace, env = model_env(grm), ...
    ),
    error = function(e) cli::cli_abort("MET model failed.", parent = e)
  )

  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(model, file.path(save_dir, "MET_Step1_Diag_Model.rds"))
  cli::cli_alert_success("MET step 1 (diag) model saved.")
  model
}

# Build MET formula strings from the final summary table.
build_met_equations <- function(fs, has_gdrop = TRUE) {
  at_term <- function(t, term) sprintf("at(studyName, %s):%s", encodeString(t, quote = "\""), term)
  labels_of <- function(x) {
    if (is.null(x) || is.na(x) || !nzchar(x)) {
      return(character(0))
    }
    gsub("\\s+", "", attr(stats::terms(stats::as.formula(x)), "term.labels"))
  }

  fixed <- "Yield ~ studyName"
  random <- c("diag(studyName):vm(gkeep, grm)", if (has_gdrop) "diag(studyName):gdrop")
  supported_fixed <- c("lin(rowNumber)", "lin(colNumber)")
  supported_random <- c("vm(gkeep,grm)", "gdrop", "blockNumber", "colNumber", "rowNumber", "units")

  for (i in seq_len(nrow(fs))) {
    t <- as.character(fs$studyName[i])
    fl <- labels_of(fs$fixed_terms[i])
    rl <- labels_of(fs$random_terms[i])
    bad <- c(setdiff(fl, supported_fixed), setdiff(rl, supported_random))
    if (length(bad) > 0) cli::cli_warn("Trial {.val {t}}: term{?s} {.val {bad}} not carried into the MET model.")

    if ("lin(rowNumber)" %in% fl) fixed <- paste(fixed, "+", at_term(t, "rowNum_num"))
    if ("lin(colNumber)" %in% fl) fixed <- paste(fixed, "+", at_term(t, "colNum_num"))
    for (term in intersect(c("blockNumber", "colNumber", "rowNumber", "units"), rl)) {
      random <- c(random, at_term(t, term))
    }
  }

  is_ar1 <- grepl("ar1", fs$residual_terms %||% rep("ar1", nrow(fs)))
  is_ar1[is.na(is_ar1)] <- TRUE
  lv <- function(x) paste0("c(", paste(encodeString(as.character(x), quote = "\""), collapse = ", "), ")")
  residual <- if (all(is_ar1)) {
    "~ dsum(~ ar1(colNumber):ar1(rowNumber) | studyName)"
  } else if (!any(is_ar1)) {
    "~ dsum(~ id(colNumber):id(rowNumber) | studyName)"
  } else {
    paste0(
      "~ dsum(~ ar1(colNumber):ar1(rowNumber) | studyName, levels = ", lv(fs$studyName[is_ar1]), ")",
      " + dsum(~ id(colNumber):id(rowNumber) | studyName, levels = ", lv(fs$studyName[!is_ar1]), ")"
    )
  }

  list(fixed = fixed, random = paste("~", paste(random, collapse = " + ")), residual = residual)
}
