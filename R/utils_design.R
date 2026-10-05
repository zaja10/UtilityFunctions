#' Diagnose Trial Design Structure
#'
#' Summarises a MET: genotypes and plots per trial, replication, and the
#' genotype overlap between trials.
#'
#' @param data Data frame containing the trial data.
#' @param genotype Character. Genotype column.
#' @param trial Character. Trial column.
#' @param rep Character. Optional replicate/block column.
#' @return An object of class \code{design_diagnosis}: a list with
#'   \code{counts} (overall genotypes \code{G}, trials \code{S}, plots),
#'   \code{per_trial} (data frame), \code{connectivity} (trial x trial shared
#'   genotype counts; the diagonal is the number of genotypes per trial) and
#'   \code{n_clusters}.
#' @export
diagnose_design <- function(data, genotype, trial, rep = NULL) {
  cols <- c(genotype, trial, rep)
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) cli::cli_abort("Column{?s} {.val {missing}} not found in {.arg data}.")

  g <- as.character(data[[genotype]])
  t <- as.character(data[[trial]])
  ok <- !is.na(g) & !is.na(t)
  d <- data[ok, , drop = FALSE]
  g <- g[ok]
  t <- t[ok]

  per_trial <- do.call(rbind, lapply(split(seq_along(t), factor(t, levels = unique(t))), function(i) {
    reps_per_g <- table(g[i])
    data.frame(
      Trial = t[i[1]], Plots = length(i), Genotypes = length(reps_per_g),
      Mean_Reps = mean(reps_per_g), Prop_Replicated = mean(reps_per_g > 1),
      N_Reps = if (!is.null(rep)) length(unique(stats::na.omit(d[[rep]][i]))) else NA_integer_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(per_trial) <- NULL

  conn <- suppressMessages(check_met_connectivity(
    data.frame(T = t, G = g),
    trial_col = "T", gen_col = "G", warn_below = 0
  ))

  structure(
    list(
      counts = list(G = length(unique(g)), S = length(unique(t)), plots = length(g)),
      per_trial = per_trial,
      connectivity = conn$shared,
      n_clusters = conn$n_clusters
    ),
    class = "design_diagnosis"
  )
}

#' @export
print.design_diagnosis <- function(x, ...) {
  cli::cli_h1("Design diagnosis")
  cli::cli_text("{x$counts$G} genotypes, {x$counts$S} trials, {x$counts$plots} plots; {x$n_clusters} connected cluster{?s}.")
  print(x$per_trial, row.names = FALSE, digits = 3)
  off <- x$connectivity[upper.tri(x$connectivity)]
  if (length(off) > 0) {
    cli::cli_text("Shared genotypes between trial pairs: min {min(off)}, median {stats::median(off)}, max {max(off)}.")
  }
  invisible(x)
}

#' Convert Factor to Numeric safely
#'
#' Converts via \code{as.character()}, so factor labels (not level codes)
#' are used.
#'
#' @param df Vector or column to convert.
#' @return Numeric vector.
#' @export
makeNm <- function(df) {
  as.numeric(as.character(df))
}

#' Find Missing Plots
#'
#' Identifies gaps in each trial's rectangular row x column grid.
#'
#' @param data Data frame containing the trial layout.
#' @param experiment Character. Experiment/trial column.
#' @param row Character. Row column (numeric or numeric-valued factor).
#' @param col Character. Column column.
#' @return A data frame of missing positions named after \code{experiment},
#'   \code{row} and \code{col}, or \code{NULL} if every grid is complete.
#' @export
find_missing_plots <- function(data, experiment = "Experiment", row = "Row", col = "Column") {
  missing <- setdiff(c(experiment, row, col), names(data))
  if (length(missing) > 0) cli::cli_abort("Column{?s} {.val {missing}} not found in {.arg data}.")

  out <- do.call(rbind, lapply(unique(data[[experiment]]), function(id) {
    sub <- data[data[[experiment]] == id & !is.na(data[[experiment]]), , drop = FALSE]
    r <- makeNm(sub[[row]])
    cc <- makeNm(sub[[col]])
    ok <- !is.na(r) & !is.na(cc)
    if (!any(ok)) {
      return(NULL)
    }
    grid <- expand.grid(R = seq(min(r[ok]), max(r[ok])), C = seq(min(cc[ok]), max(cc[ok])))
    gap <- !paste(grid$R, grid$C) %in% paste(r[ok], cc[ok])
    if (!any(gap)) {
      return(NULL)
    }
    stats::setNames(data.frame(id, grid$R[gap], grid$C[gap], stringsAsFactors = FALSE), c(experiment, row, col))
  }))

  if (is.null(out)) cli::cli_alert_success("All grids are complete.")
  out
}
