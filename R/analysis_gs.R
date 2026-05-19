#' Evaluate Genomic Selection Predictions from ASReml FA Models
#'
#' Extracts latent scores and prediction error variances (PEVs) from a fitted
#' ASReml factor analytic (FA) model. Following the FAST framework (Smith & Cullis, 2018),
#' this function applies a Principal Component (SVD) rotation to the raw FA parameters.
#' Overall Performance (OP) is defined by the rotated 1st Factor, while Stability (RMSD)
#' is calculated from the higher-order rotated factors (k >= 2) representing crossover GxE.
#'
#' @param model An \code{asreml} object.
#' @param target_ids Character vector. IDs of the target (unphenotyped) lines.
#' @param training_ids Character vector. IDs of the training (phenotyped) lines.
#' @param target_term Character. The base random term to extract from the ASReml model
#'   (e.g., \code{"fa(studyName, 2):vm(gkeep, G.inv)"}).
#' @param g_inv Optional. The inverse Genomic Relationship Matrix (GRM) used in the model,
#'   used for diagnosing missing targets. Expected to have dimnames corresponding to IDs.
#' @param grm Optional. The Genomic Relationship Matrix (GRM). If provided, it is used to calculate
#'   Genetic Relatedness (row means) for plotting.
#' @param output_prefix Character. Prefix for the output CSV and PDF files. Default \code{"GS_Predictions"}.
#' @param top_pct Numeric. The quantile threshold for highlighting top candidates (e.g., 0.95 for top 5 percent).
#' @param save_output Logical. Whether to save the CSV and PDF outputs. Default \code{TRUE}.
#'
#' @return A list containing the following elements:
#' * \code{blups}: A dataframe of scaled GEBVs, PEVs, Reliability, and Stability for all lines.
#' * \code{missing_targets}: A character vector of target IDs that were not found in the model.
#' * \code{plots}: A combined \code{patchwork} plot object.
#'
#' @import ggplot2
#' @import patchwork
#' @import ggrepel
#' @import dplyr
#' @import cli
#' @importFrom stats quantile
#' @importFrom utils head write.csv
#' @export
evaluate_gs_predictions <- function(
  model, target_ids, training_ids, target_term, g_inv = NULL,
  grm = NULL, output_prefix = "GS_Predictions", top_pct = 0.95,
  save_output = TRUE
) {
  cli::cli_h1("Evaluating Genomic Selection Predictions (FAST Framework)")
  if (!inherits(model, "asreml")) {
    cli::cli_abort("Model must be an object of class {.cls asreml}")
  }

  base_term <- sub("_Comp[0-9]+", "", target_term)
  cli::cli_alert_info("Analyzing FA base term: {.val {base_term}}")

  # --- BUG FIX: Strict extraction of FA parameters for THIS term only ---
  vparams <- model$vparameters
  all_fa_keys <- grep("!fa[0-9]+$", names(vparams), value = TRUE)

  # Remove spaces for robust matching (protects against "fa(studyName, 2)" vs "fa(studyName,2)")
  clean_base <- gsub(" ", "", base_term)
  clean_keys <- gsub(" ", "", all_fa_keys)

  # Only keep keys that belong to the specific target_term
  fa_keys <- all_fa_keys[grepl(clean_base, clean_keys, fixed = TRUE)]

  if (length(fa_keys) == 0) {
    cli::cli_abort("Could not find any factor loadings (!fa) specifically for the term {.val {base_term}}.")
  }
  # ----------------------------------------------------------------------

  factors_detected <- unique(sub(".*!(fa[0-9]+)$", "\\1", fa_keys))
  k <- length(factors_detected)
  cli::cli_alert_info("Detected FA model with k = {k} factor(s). Applying SVD PC Rotation.")

  # The fixed regex that handles exclamation marks
  env_names <- unique(sub(".*[:!]([^:!]+)!fa[0-9]+$", "\\1", fa_keys))
  if (length(env_names) == 0) {
    env_names <- paste0("Env", seq_along(grep("!fa1$", fa_keys)))
  }

  raw_loadings_list <- list()
  for (i in seq_len(k)) {
    fac_label <- paste0("fa", i)
    keys_k <- grep(paste0("!", fac_label, "$"), fa_keys, value = TRUE)
    raw_loadings_list[[i]] <- vparams[keys_k]
  }

  Lambda_raw <- do.call(cbind, raw_loadings_list)
  rownames(Lambda_raw) <- env_names

  summ_coef <- summary(model, coef = TRUE)$coef.random
  rand_names <- rownames(summ_coef)
  term_parts <- unlist(strsplit(base_term, ":"))
  target_idx <- seq_along(rand_names)
  for (part in term_parts) {
    target_idx <- intersect(target_idx, grep(part, rand_names, fixed = TRUE))
  }

  if (length(target_idx) == 0) {
    cli::cli_abort("Could not find any random coefficients matching the base term {.val {base_term}}.")
  }

  coef_sub <- summ_coef[target_idx, , drop = FALSE]
  names_sub <- rownames(coef_sub)
  raw_scores_list <- list()
  raw_pev_list <- list()
  clean_ids_master <- NULL

  for (i in seq_len(k)) {
    comp_pat <- paste0("_Comp", i)
    idx_k <- grep(comp_pat, names_sub, fixed = TRUE)
    if (length(idx_k) == 0) {
      next
    }
    sub_n <- names_sub[idx_k]
    sub_c <- coef_sub[idx_k, , drop = FALSE]
    clean_ids <- sub(paste0(".*", comp_pat, "(_|:)?"), "", sub_n)
    clean_ids <- sub("^:", "", clean_ids)
    clean_ids <- sub("vm\\([^)]+\\)_?", "", clean_ids)
    if (is.null(clean_ids_master)) {
      clean_ids_master <- clean_ids
    }
    sol_col <- grep("solution|value", colnames(sub_c), ignore.case = TRUE, value = TRUE)[1]
    se_col <- grep("std", colnames(sub_c), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(sol_col) || is.na(se_col)) {
      cli::cli_abort("Could not identify solution or std error columns in model coefficients.")
    }
    raw_scores_list[[i]] <- as.numeric(sub_c[, sol_col])
    raw_pev_list[[i]] <- as.numeric(sub_c[, se_col])^2
  }

  if (length(raw_scores_list) == 0) {
    cli::cli_abort("Failed to process any _Comp elements. Please check the target_term format.")
  }

  F_raw <- do.call(cbind, raw_scores_list)
  PEV_raw <- do.call(cbind, raw_pev_list)
  rownames(F_raw) <- clean_ids_master
  Lambda_rot <- Lambda_raw
  F_rot <- F_raw
  V <- diag(k)

  if (k > 1) {
    svd_res <- svd(Lambda_raw)
    V <- svd_res$v
    Lambda_rot <- Lambda_raw %*% V
    F_rot <- F_raw %*% V
  }

  if (mean(Lambda_rot[, 1]) < 0) {
    Lambda_rot[, 1] <- -Lambda_rot[, 1]
    F_rot[, 1] <- -F_rot[, 1]
    V[, 1] <- -V[, 1]
  }

  mean_load_rot <- colMeans(Lambda_rot)
  OP <- F_rot[, 1] * mean_load_rot[1]
  Var_A_OP <- mean_load_rot[1]^2
  PEV_F1_rot <- rowSums(PEV_raw * (matrix(V[, 1]^2, nrow = nrow(PEV_raw), ncol = k, byrow = TRUE)))
  PEV_OP <- PEV_F1_rot * Var_A_OP

  blups_gs <- data.frame(
    id = clean_ids_master, GEBV_OP = OP,
    PEV_OP = PEV_OP, Reliability = 1 - (PEV_OP / Var_A_OP),
    Stability_RMSD = 0, stringsAsFactors = FALSE
  )

  blups_gs$Reliability <- ifelse(blups_gs$Reliability < 0, 0, blups_gs$Reliability)

  if (k > 1) {
    cli::cli_alert_info("Calculating crossover Stability (RMSD) from rotated higher-order factors...")
    L_int <- Lambda_rot[, 2:k, drop = FALSE]
    F_int <- F_rot[, 2:k, drop = FALSE]
    I_mat <- F_int %*% t(L_int)
    blups_gs$Stability_RMSD <- sqrt(rowMeans(I_mat^2))
  }

  target_ids <- as.character(target_ids)
  training_ids <- as.character(training_ids)
  blups_gs$status <- NA_character_
  blups_gs$status[blups_gs$id %in% training_ids] <- "Training"
  blups_gs$status[blups_gs$id %in% target_ids] <- "Target"
  missing_targets <- setdiff(target_ids, blups_gs$id)

  if (length(missing_targets) > 0) {
    cli::cli_alert_warning("Found {length(missing_targets)} target IDs missing from the model predictions.")
    if (!is.null(g_inv)) {
      grm_names <- rownames(g_inv)
      if (is.null(grm_names) && !is.null(attr(g_inv, "dimnames"))) {
        grm_names <- attr(g_inv, "dimnames")[[1]]
      }
      if (!is.null(grm_names)) {
        missing_from_g <- setdiff(missing_targets, grm_names)
        if (length(missing_from_g) == length(missing_targets)) {
          cli::cli_alert_danger("All missing lines are completely missing from the G.inv matrix.")
        } else {
          cli::cli_alert_danger("Some missing lines ARE in G.inv but failed to match ASReml output names.")
        }
      }
    }
  }

  blups_assigned <- blups_gs[!is.na(blups_gs$status), ]

  if (nrow(blups_assigned[blups_assigned$status == "Target", ]) == 0) {
    cli::cli_abort("No target IDs were matched in the model predictions.")
  }

  if (k > 1) {
    blups_assigned$Selection_Index <- blups_assigned$GEBV_OP - blups_assigned$Stability_RMSD
  } else {
    blups_assigned$Selection_Index <- blups_assigned$GEBV_OP
  }

  blups_target <- blups_assigned[blups_assigned$status == "Target", ]
  blups_target <- blups_target[order(blups_target$Selection_Index, decreasing = TRUE), ]

  mean_rel_tr <- mean(blups_assigned$Reliability[blups_assigned$status == "Training"], na.rm = TRUE)
  mean_rel_ta <- mean(blups_target$Reliability, na.rm = TRUE)

  cli::cli_h2("Genomic Selection Diagnostics")
  cli::cli_text("Mean OP Reliability (Training): {.val {round(mean_rel_tr, 3)}}")
  cli::cli_text("Mean OP Reliability (Target)  : {.val {round(mean_rel_ta, 3)}}")
  cli::cli_text("Top 5 Candidates for Advancement:")

  if (k > 1) {
    print(head(blups_target[, c("id", "GEBV_OP", "Stability_RMSD", "Selection_Index", "Reliability")], 5))
  } else {
    print(head(blups_target[, c("id", "GEBV_OP", "PEV_OP", "Reliability")], 5))
  }

  if (save_output) {
    csv_file <- paste0(output_prefix, "_BLUPS.csv")
    write.csv(blups_target, csv_file, row.names = FALSE)
    cli::cli_alert_success("Saved target BLUPs to {.file {csv_file}}")
  }

  cli::cli_alert_info("Generating enhanced diagnostic plots...")
  blups_assigned$status <- factor(blups_assigned$status, levels = c("Training", "Target"))
  mean_train <- mean(blups_assigned$GEBV_OP[blups_assigned$status == "Training"], na.rm = TRUE)
  mean_target <- mean(blups_assigned$GEBV_OP[blups_assigned$status == "Target"], na.rm = TRUE)

  if (!is.null(grm)) {
    cli::cli_alert_info("Partitioning GRM to calculate targeted target-to-training relatedness metrics...")

    grm_ids <- rownames(grm)
    valid_train <- intersect(training_ids, grm_ids)
    valid_target <- intersect(blups_assigned$id[blups_assigned$status == "Target"], grm_ids)

    if (length(valid_train) == 0 || length(valid_target) == 0) {
      cli::cli_alert_warning("Insufficient GRM overlap with training/target IDs. Skipping relatedness plots.")
    } else {
      # 1. Isolate the target-to-training submatrix
      G_target_train <- grm[valid_target, valid_train, drop = FALSE]

      # 2. Calculate the metrics per target line
      rel_target <- data.frame(
        id = valid_target,
        Mean_Rel = rowMeans(G_target_train, na.rm = TRUE),
        Max_Rel = apply(G_target_train, 1, max, na.rm = TRUE),
        stringsAsFactors = FALSE
      )

      # For plotting context, calculate training lines' relation to the REST of the training set
      G_train_train <- grm[valid_train, valid_train, drop = FALSE]
      diag(G_train_train) <- NA # Remove self-relatedness (1.0) to find true max relative

      rel_train <- data.frame(
        id = valid_train,
        Mean_Rel = rowMeans(G_train_train, na.rm = TRUE),
        Max_Rel = apply(G_train_train, 1, max, na.rm = TRUE),
        stringsAsFactors = FALSE
      )

      rel_df <- rbind(rel_target, rel_train)
      blups_assigned <- merge(blups_assigned, rel_df, by = "id", all.x = TRUE)
    }
  }

  plot_A <- ggplot(blups_assigned, aes(x = .data$GEBV_OP, fill = .data$status)) +
    geom_density(alpha = 0.6, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c(Training = "#56B4E9", Target = "#E69F00")) +
    geom_vline(xintercept = mean_train, linetype = "dashed", color = "#56B4E9", linewidth = 1) +
    geom_vline(xintercept = mean_target, linetype = "dashed", color = "#E69F00", linewidth = 1) +
    labs(title = "A. Distribution of Overall Performance (OP)", x = "Genomic Estimated Breeding Value (Rotated Factor 1)", y = "Density", fill = "Cohort") +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")

  top_cutoff <- quantile(blups_assigned$GEBV_OP[blups_assigned$status == "Target"], top_pct, na.rm = TRUE)
  top_candidates <- head(blups_target, 5)

  if (k > 1) {
    plot_B <- ggplot(blups_assigned, aes(x = .data$status, y = .data$Stability_RMSD, fill = .data$status)) +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.8) +
      scale_fill_manual(values = c(Training = "#56B4E9", Target = "#E69F00")) +
      labs(title = paste("B. Stability (RMSD of Factors 2 to", k, ")"), x = "Cohort", y = "RMSD (Lower = More Stable)") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")

    plot_C <- ggplot(blups_assigned, aes(x = .data$Stability_RMSD, y = .data$GEBV_OP)) +
      geom_point(aes(color = .data$status), alpha = 0.7, size = 2.5, stroke = 0) +
      scale_color_manual(values = c(Training = "#56B4E9", Target = "#E69F00")) +
      geom_hline(yintercept = top_cutoff, linetype = "dotted", color = "black", linewidth = 1) +
      ggrepel::geom_label_repel(data = top_candidates, aes(label = .data$id), size = 3, fontface = "bold", box.padding = 0.5, point.padding = 0.3, segment.color = "grey50", max.overlaps = Inf) +
      labs(title = "C. Breeder's Selection Space", x = "Crossover GxE (Stability RMSD)", y = "Genomic Merit (GEBV_OP)", color = "Cohort") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")
  } else {
    plot_B <- ggplot(blups_assigned, aes(x = .data$status, y = .data$Reliability, fill = .data$status)) +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.8) +
      scale_fill_manual(values = c(Training = "#56B4E9", Target = "#E69F00")) +
      labs(title = "B. Prediction Reliability", x = "Cohort", y = "Reliability (Accuracy²)") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")

    plot_C <- ggplot(blups_assigned, aes(x = .data$Reliability, y = .data$GEBV_OP)) +
      geom_point(aes(color = .data$status), alpha = 0.7, size = 2.5, stroke = 0) +
      scale_color_manual(values = c(Training = "#56B4E9", Target = "#E69F00")) +
      geom_hline(yintercept = top_cutoff, linetype = "dotted", color = "black", linewidth = 1) +
      ggrepel::geom_label_repel(data = top_candidates, aes(label = .data$id), size = 3, fontface = "bold", box.padding = 0.5, point.padding = 0.3, segment.color = "grey50", max.overlaps = Inf) +
      labs(title = "C. Breeder's Selection Space", x = "Prediction Reliability", y = "Genomic Merit (GEBV_OP)", color = "Cohort") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")
  }

  has_plot_D <- FALSE
  if ("Max_Rel" %in% names(blups_assigned)) {
    has_plot_D <- TRUE

    plot_D <- ggplot2::ggplot(blups_assigned, ggplot2::aes(x = .data$Max_Rel, y = .data$Reliability)) +
      ggplot2::geom_point(ggplot2::aes(color = .data$status), alpha = 0.7, size = 2, stroke = 0) +
      # Fit separate regressions for Training and Target to show the distinct information flows
      ggplot2::geom_smooth(ggplot2::aes(group = .data$status, color = .data$status), method = "lm", se = FALSE, linewidth = 1) +
      ggplot2::scale_color_manual(values = c(Training = "#56B4E9", Target = "#E69F00")) +
      ggplot2::labs(
        title = "D. Relatedness vs Prediction Reliability",
        x = expression("Max Genomic Relationship to Training (" * g[max] * ")"),
        y = "Theoretical Reliability (r²)"
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(legend.position = "none")
  }

  if (has_plot_D) {
    combined_figure <- (plot_A | plot_D) / (plot_B | plot_C) +
      patchwork::plot_annotation(
        title = "Genomic Prediction Validation and Selection Differentials",
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      )
  } else {
    combined_figure <- plot_A / (plot_B | plot_C) +
      patchwork::plot_annotation(
        title = "Genomic Prediction Validation and Selection Differentials",
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      )
  }

  if (save_output) {
    pdf_file <- paste0(output_prefix, "_Figure.pdf")
    ggsave(filename = pdf_file, plot = combined_figure, device = "pdf", width = 11, height = 9, dpi = 600)
    cli::cli_alert_success("Saved enhanced manuscript figure to {.file {pdf_file}}")
  }

  res <- list(blups = blups_assigned, missing_targets = missing_targets, plots = combined_figure)
  class(res) <- c("gs_predictions", "list")
  return(res)
}

#' Print Method for Genomic Selection Predictions
#'
#' @param x An object of class `gs_predictions`.
#' @param ... Additional arguments (ignored).
#' @return Prints a summary and returns the object invisibly.
#' @export
print.gs_predictions <- function(x, ...) {
  if (!inherits(x, "gs_predictions")) stop("Object must be of class 'gs_predictions'")

  n_train <- sum(x$blups$status == "Training", na.rm = TRUE)
  n_target <- sum(x$blups$status == "Target", na.rm = TRUE)

  cli::cli_h1("Genomic Selection Predictions")
  cli::cli_text("Cohorts: {.val {n_train}} Training | {.val {n_target}} Target")

  if (length(x$missing_targets) > 0) {
    cli::cli_alert_warning("{.val {length(x$missing_targets)}} target lines were missing from the model.")
  }

  cli::cli_h2("Top 5 Candidates")
  blups_target <- x$blups[x$blups$status == "Target", ]
  blups_target <- blups_target[order(blups_target$Selection_Index, decreasing = TRUE), ]

  cols_to_print <- c("id", "GEBV_OP")
  if ("Stability_RMSD" %in% names(blups_target) && any(blups_target$Stability_RMSD > 0)) {
    cols_to_print <- c(cols_to_print, "Stability_RMSD", "Selection_Index")
  } else {
    cols_to_print <- c(cols_to_print, "PEV_OP")
  }
  cols_to_print <- c(cols_to_print, "Reliability")

  print(head(blups_target[, cols_to_print], 5))

  invisible(x)
}

#' Summary Method for Genomic Selection Predictions
#'
#' @param object An object of class `gs_predictions`.
#' @param ... Additional arguments (ignored).
#' @return Prints a detailed summary and returns the object invisibly.
#' @export
summary.gs_predictions <- function(object, ...) {
  if (!inherits(object, "gs_predictions")) stop("Object must be of class 'gs_predictions'")

  cli::cli_h1("Genomic Selection Validation Summary")

  blups <- object$blups
  mean_rel_tr <- mean(blups$Reliability[blups$status == "Training"], na.rm = TRUE)
  mean_rel_ta <- mean(blups$Reliability[blups$status == "Target"], na.rm = TRUE)

  cli::cli_text("Mean Prediction Reliability:")
  cli::cli_text("  Training: {.val {round(mean_rel_tr, 3)}}")
  cli::cli_text("  Target  : {.val {round(mean_rel_ta, 3)}}")

  if ("Stability_RMSD" %in% names(blups) && any(blups$Stability_RMSD > 0)) {
    mean_rmsd_tr <- mean(blups$Stability_RMSD[blups$status == "Training"], na.rm = TRUE)
    mean_rmsd_ta <- mean(blups$Stability_RMSD[blups$status == "Target"], na.rm = TRUE)

    cli::cli_text("")
    cli::cli_text("Mean Stability (RMSD - Lower is more stable):")
    cli::cli_text("  Training: {.val {round(mean_rmsd_tr, 3)}}")
    cli::cli_text("  Target  : {.val {round(mean_rmsd_ta, 3)}}")
  }

  cli::cli_text("")
  cli::cli_text("Use {.code plot(object)} to visualize the Breeder's Selection Space.")

  invisible(object)
}

#' Plot Method for Genomic Selection Predictions
#'
#' @param x An object of class `gs_predictions`.
#' @param ... Additional arguments (ignored).
#' @return Returns the `patchwork` combined figure.
#' @export
plot.gs_predictions <- function(x, ...) {
  if (!inherits(x, "gs_predictions")) stop("Object must be of class 'gs_predictions'")
  return(x$plots)
}
