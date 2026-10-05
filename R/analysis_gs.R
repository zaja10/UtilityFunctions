#' Evaluate Genomic Selection Predictions from ASReml FA Models
#'
#' Extracts the FA solution with \code{\link{fa.asreml}} (SVD rotation, FAST
#' indices of Smith & Cullis 2018), labels genotypes as training or target,
#' computes a selection index and, optionally, target-to-training genomic
#' relatedness, and assembles diagnostic plots.
#'
#' @param model An \code{asreml} object.
#' @param target_ids Character vector. IDs of the target (unphenotyped) lines.
#' @param training_ids Character vector. IDs of the training (phenotyped) lines.
#' @param target_term Character. The FA term, e.g.
#'   \code{"fa(studyName, 2):vm(gkeep, G.inv)"} (a trailing \code{_CompN} is
#'   ignored).
#' @param g_inv Optional. The inverse GRM used in the model; only used to
#'   diagnose target IDs missing from the predictions.
#' @param grm Optional. The GRM \strong{G}. Used for reliabilities
#'   (\eqn{r^2_g = 1 - PEV/g_{ii}}) and target-to-training relatedness.
#' @param output_prefix Character. Prefix for output files.
#' @param top_pct Numeric. Quantile of target OP marked on the selection plot
#'   (default 0.95).
#' @param save_output Logical. Write \code{<prefix>_BLUPS.csv} and
#'   \code{<prefix>_Figure.pdf} (default \code{FALSE}).
#' @param weight_stability Numeric. Penalty on RMSD in the selection index
#'   \eqn{OP - b \cdot RMSD} (default 1; see \code{\link{rank_genotypes}}).
#'
#' @return An object of class \code{gs_predictions}: a list with \code{blups}
#'   (training and target lines), \code{missing_targets}, \code{plots}
#'   (a patchwork figure) and \code{fa} (the \code{fa_model}).
#' @export
evaluate_gs_predictions <- function(model, target_ids, training_ids, target_term, g_inv = NULL,
                                    grm = NULL, output_prefix = "GS_Predictions", top_pct = 0.95,
                                    save_output = FALSE, weight_stability = 1) {
  if (!inherits(model, "asreml")) cli::cli_abort("Model must be an object of class {.cls asreml}.")
  base_term <- sub("_Comp[0-9]+", "", target_term)
  fa <- fa.asreml(model, classify = base_term, rotate = "svd", grm = grm)
  if (is.null(fa$fast)) cli::cli_abort("No genotype scores could be extracted for {.val {base_term}}.")
  k <- fa$meta$k

  fast <- fa$fast
  blups <- data.frame(
    id = fast$Genotype, GEBV_OP = fast$OP, PEV_OP = fast$PEV_OP %||% NA_real_,
    Reliability = fast$Reliability %||% NA_real_, Stability_RMSD = fast$RMSD,
    stringsAsFactors = FALSE
  )

  target_ids <- as.character(target_ids)
  training_ids <- as.character(training_ids)
  overlap <- intersect(target_ids, training_ids)
  if (length(overlap) > 0) {
    cli::cli_warn("{length(overlap)} ID{?s} listed as both training and target; treated as training.")
  }
  blups$status <- NA_character_
  blups$status[blups$id %in% target_ids] <- "Target"
  blups$status[blups$id %in% training_ids] <- "Training"

  missing_targets <- setdiff(target_ids, blups$id)
  if (length(missing_targets) > 0) diagnose_missing_targets(missing_targets, g_inv)

  blups <- blups[!is.na(blups$status), , drop = FALSE]
  if (!any(blups$status == "Target")) cli::cli_abort("No target IDs were matched in the model predictions.")
  blups$Selection_Index <- blups$GEBV_OP - weight_stability * blups$Stability_RMSD
  blups$status <- factor(blups$status, levels = c("Training", "Target"))

  if (!is.null(grm)) blups <- add_relatedness(blups, grm)
  blups <- blups[order(blups$status, -blups$Selection_Index), , drop = FALSE]
  rownames(blups) <- NULL

  plots <- gs_plots(blups, k, top_pct)

  if (save_output) {
    csv_file <- paste0(output_prefix, "_BLUPS.csv")
    utils::write.csv(blups, csv_file, row.names = FALSE)
    pdf_file <- paste0(output_prefix, "_Figure.pdf")
    ggplot2::ggsave(pdf_file, plots, device = "pdf", width = 11, height = 9)
    cli::cli_alert_success("Saved {.file {csv_file}} and {.file {pdf_file}}.")
  }

  structure(
    list(blups = blups, missing_targets = missing_targets, plots = plots, fa = fa, k = k),
    class = c("gs_predictions", "list")
  )
}

diagnose_missing_targets <- function(missing, g_inv) {
  cli::cli_alert_warning("{length(missing)} target ID{?s} missing from the model predictions.")
  if (is.null(g_inv)) {
    return(invisible())
  }
  ids <- attr(g_inv, "rowNames") %||% rownames(g_inv)
  if (is.null(ids)) {
    return(invisible())
  }
  in_g <- intersect(missing, ids)
  if (length(in_g) == 0) {
    cli::cli_alert_danger("All missing targets are absent from the G inverse.")
  } else {
    cli::cli_alert_danger("{length(in_g)} missing target{?s} {?is/are} in the G inverse but not in the predictions; check factor levels.")
  }
  invisible()
}

# Mean and maximum genomic relationship of each line to the training set
# (excluding self-relationships for training lines).
add_relatedness <- function(blups, grm) {
  grm <- as.matrix(grm)
  ids <- rownames(grm)
  train <- intersect(as.character(blups$id[blups$status == "Training"]), ids)
  target <- intersect(as.character(blups$id[blups$status == "Target"]), ids)
  if (length(train) == 0 || length(target) == 0) {
    cli::cli_alert_warning("Insufficient GRM overlap with training/target IDs; relatedness not computed.")
    return(blups)
  }
  Gtt <- grm[target, train, drop = FALSE]
  Grr <- grm[train, train, drop = FALSE]
  diag(Grr) <- NA
  safe_max <- function(m) apply(m, 1, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
  rel <- rbind(
    data.frame(id = target, Mean_Rel = rowMeans(Gtt), Max_Rel = safe_max(Gtt), stringsAsFactors = FALSE),
    data.frame(id = train, Mean_Rel = rowMeans(Grr, na.rm = TRUE), Max_Rel = safe_max(Grr), stringsAsFactors = FALSE)
  )
  merge(blups, rel, by = "id", all.x = TRUE, sort = FALSE)
}

gs_plots <- function(blups, k, top_pct) {
  cols <- c(Training = "#56B4E9", Target = "#E69F00")
  means <- tapply(blups$GEBV_OP, blups$status, mean, na.rm = TRUE)
  target <- blups[blups$status == "Target", , drop = FALSE]
  top <- utils::head(target[order(-target$Selection_Index), , drop = FALSE], 5)
  cutoff <- stats::quantile(target$GEBV_OP, top_pct, na.rm = TRUE)

  pA <- ggplot2::ggplot(blups, ggplot2::aes(x = .data$GEBV_OP, fill = .data$status)) +
    ggplot2::geom_density(alpha = 0.6, color = "black", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = means[["Training"]], linetype = "dashed", color = cols[["Training"]], linewidth = 1) +
    ggplot2::geom_vline(xintercept = means[["Target"]], linetype = "dashed", color = cols[["Target"]], linewidth = 1) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = "A. Overall Performance (OP)", x = "OP (rotated factor 1)", y = "Density", fill = "Cohort") +
    theme_genetics()

  x_var <- if (k > 1) "Stability_RMSD" else "Reliability"
  x_lab <- if (k > 1) "Crossover GxE (RMSD)" else "Reliability"
  rel_lab <- quote(bold(Reliability ~ (r[g]^2)))
  x_axis <- if (k > 1) x_lab else rel_lab
  pB <- ggplot2::ggplot(blups, ggplot2::aes(x = .data$status, y = .data[[x_var]], fill = .data$status)) +
    ggplot2::geom_violin(alpha = 0.6, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.2, fill = "white", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = paste("B.", x_lab), x = "Cohort", y = x_axis) +
    theme_genetics() +
    ggplot2::theme(legend.position = "none")

  pC <- ggplot2::ggplot(blups, ggplot2::aes(x = .data[[x_var]], y = .data$GEBV_OP)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$status), alpha = 0.7, size = 2.5, stroke = 0) +
    ggplot2::geom_hline(yintercept = cutoff, linetype = "dotted", linewidth = 1) +
    ggrepel::geom_label_repel(
      data = top, ggplot2::aes(label = .data$id), size = 3, fontface = "bold",
      box.padding = 0.5, segment.color = "grey50", max.overlaps = Inf
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(
      title = "C. Selection Space", x = x_axis, y = "OP",
      caption = sprintf("Dotted line: target OP %g%% quantile; labels: top 5 targets by selection index", 100 * top_pct)
    ) +
    theme_genetics() +
    ggplot2::theme(legend.position = "none")

  title <- patchwork::plot_annotation(
    title = "Genomic Prediction Validation and Selection",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
  )

  if ("Max_Rel" %in% names(blups) && any(!is.na(blups$Max_Rel) & !is.na(blups$Reliability))) {
    pD <- ggplot2::ggplot(blups, ggplot2::aes(x = .data$Max_Rel, y = .data$Reliability, color = .data$status)) +
      ggplot2::geom_point(alpha = 0.7, size = 2, stroke = 0) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1) +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::labs(title = "D. Relatedness vs Reliability", x = "Max genomic relationship to training", y = rel_lab) +
      theme_genetics() +
      ggplot2::theme(legend.position = "none")
    (pA | pD) / (pB | pC) + title
  } else {
    pA / (pB | pC) + title
  }
}

#' Print Method for Genomic Selection Predictions
#'
#' @param x A \code{gs_predictions} object.
#' @param n Integer. Number of top targets to show (default 5).
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.gs_predictions <- function(x, n = 5, ...) {
  b <- x$blups
  cli::cli_h1("Genomic Selection Predictions")
  cli::cli_text("Cohorts: {sum(b$status == 'Training')} training | {sum(b$status == 'Target')} target | k = {x$k}")
  if (length(x$missing_targets) > 0) cli::cli_alert_warning("{length(x$missing_targets)} target line{?s} missing from the model.")
  cli::cli_h2("Top {n} targets by selection index")
  tgt <- b[b$status == "Target", , drop = FALSE]
  cols <- intersect(c("id", "GEBV_OP", "Stability_RMSD", "Selection_Index", "Reliability"), names(tgt))
  print(utils::head(tgt[order(-tgt$Selection_Index), cols, drop = FALSE], n), row.names = FALSE)
  invisible(x)
}

#' Summary Method for Genomic Selection Predictions
#'
#' @param object A \code{gs_predictions} object.
#' @param ... Ignored.
#' @return An object of class \code{summary.gs_predictions}: a data frame of
#'   cohort means.
#' @export
summary.gs_predictions <- function(object, ...) {
  b <- object$blups
  cols <- intersect(c("GEBV_OP", "Reliability", "Stability_RMSD", "Max_Rel"), names(b))
  out <- do.call(rbind, lapply(split(b, b$status), function(s) {
    data.frame(Cohort = as.character(s$status[1]), N = nrow(s), as.list(colMeans(s[, cols, drop = FALSE], na.rm = TRUE)))
  }))
  rownames(out) <- NULL
  structure(out, class = c("summary.gs_predictions", "data.frame"))
}

#' @export
print.summary.gs_predictions <- function(x, ...) {
  cli::cli_h1("Genomic Selection Validation Summary (cohort means)")
  print(structure(x, class = "data.frame"), row.names = FALSE, digits = 3)
  invisible(x)
}

#' Plot Method for Genomic Selection Predictions
#'
#' @param x A \code{gs_predictions} object.
#' @param ... Ignored.
#' @return The \code{patchwork} figure.
#' @export
plot.gs_predictions <- function(x, ...) x$plots
