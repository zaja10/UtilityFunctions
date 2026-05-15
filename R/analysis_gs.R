#' Evaluate Genomic Selection Predictions from ASReml FA Models
#'
#' Extracts latent scores and prediction error variances (PEVs) from a fitted
#' ASReml factor analytic (FA) model. Following the FAST framework (Smith & Cullis, 2018), 
#' this function applies a Principal Component (SVD) rotation to the raw FA parameters.
#' Overall Performance (OP) is defined by the rotated 1st Factor, while Stability (RMSD) 
#' is calculated from the higher-order rotated factors ($k \ge 2$) representing crossover GxE.
#'
#' @param model An \code{asreml} object.
#' @param target_ids Character vector. IDs of the target (unphenotyped) lines.
#' @param training_ids Character vector. IDs of the training (phenotyped) lines.
#' @param target_term Character. The base random term to extract from the ASReml model
#'   (e.g., \code{"fa(studyName, 2):vm(gkeep, G.inv)"}). 
#' @param g_inv Optional. The inverse Genomic Relationship Matrix (GRM) used in the model, 
#'   used for diagnosing missing targets. Expected to have dimnames corresponding to IDs.
#' @param output_prefix Character. Prefix for the output CSV and PDF files. Default \code{"GS_Predictions"}.
#' @param top_pct Numeric. The quantile threshold for highlighting top candidates (e.g., 0.95 for top 5%).
#' @param save_output Logical. Whether to save the CSV and PDF outputs. Default \code{TRUE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{blups}: A dataframe of scaled GEBVs, PEVs, Reliability, and Stability for all lines.
#'   \item \code{missing_targets}: A character vector of target IDs that were not found in the model.
#'   \item \code{plots}: A combined \code{patchwork} plot object.
#' }
#' 
#' @import ggplot2
#' @import patchwork
#' @import ggrepel
#' @import dplyr
#' @import cli
#' @importFrom stats quantile
#' @importFrom utils head write.csv
#' @export
evaluate_gs_predictions <- function(model, target_ids, training_ids, target_term, 
                                    g_inv = NULL, output_prefix = "GS_Predictions", 
                                    top_pct = 0.95, save_output = TRUE) {
  
  cli::cli_h1("Evaluating Genomic Selection Predictions (FAST Framework)")
  
  if (!inherits(model, "asreml")) {
    cli::cli_abort("Model must be an object of class {.cls asreml}")
  }
  
  # Ensure target_term is just the base term without a specific _Comp
  base_term <- sub("_Comp[0-9]+", "", target_term)
  cli::cli_alert_info("Analyzing FA base term: {.val {base_term}}")
  
  # -----------------------------------------------------------------------
  # 1. Identify Factors and Raw Loadings
  # -----------------------------------------------------------------------
  site_var <- sub(".*fa\\(([^,]+),.*", "\\1", base_term)
  site_var <- trimws(site_var)
  
  vparams <- model$vparameters
  fa_keys <- grep(paste0("^", site_var, ".*!fa[0-9]+$"), names(vparams), value = TRUE)
  
  if (length(fa_keys) == 0) {
    fa_keys <- grep("!fa[0-9]+$", names(vparams), value = TRUE)
  }
  if (length(fa_keys) == 0) {
    cli::cli_abort("Could not find any factor loadings (!fa) in model$vparameters.")
  }
  
  factors_detected <- unique(sub(".*!(fa[0-9]+)$", "\\1", fa_keys))
  k <- length(factors_detected)
  cli::cli_alert_info("Detected FA model with k = {k} factor(s). Applying SVD PC Rotation.")
  
  # Extract raw loadings matrix (Environments x Factors)
  env_names <- unique(sub(".*:([^!]+)!fa[0-9]+$", "\\1", fa_keys))
  if (length(env_names) == 0) env_names <- paste0("Env", seq_along(grep("!fa1$", fa_keys)))
  
  raw_loadings_list <- list()
  for (i in seq_len(k)) {
    fac_label <- paste0("fa", i)
    raw_loadings_list[[i]] <- vparams[grep(paste0("!", fac_label, "$"), names(vparams))]
  }
  Lambda_raw <- do.call(cbind, raw_loadings_list)
  rownames(Lambda_raw) <- env_names
  
  # -----------------------------------------------------------------------
  # 2. Extract Latent Scores and PEVs Safely
  # -----------------------------------------------------------------------
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
    
    if (length(idx_k) == 0) next
    
    sub_n <- names_sub[idx_k]
    sub_c <- coef_sub[idx_k, , drop = FALSE]
    
    clean_ids <- sub(paste0(".*", comp_pat, "(_|:)?"), "", sub_n)
    clean_ids <- sub("^:", "", clean_ids)
    clean_ids <- sub("vm\\([^)]+\\)_?", "", clean_ids) 
    
    if (is.null(clean_ids_master)) clean_ids_master <- clean_ids
    
    raw_scores_list[[i]] <- as.numeric(sub_c[, "solution"])
    raw_pev_list[[i]]    <- as.numeric(sub_c[, "std error"])^2 
  }
  
  if (length(raw_scores_list) == 0) {
    cli::cli_abort("Failed to process any _Comp elements. Please check the target_term format.")
  }
  
  F_raw <- do.call(cbind, raw_scores_list)
  PEV_raw <- do.call(cbind, raw_pev_list)
  rownames(F_raw) <- clean_ids_master
  
  # -----------------------------------------------------------------------
  # 3. Apply PC (SVD) Rotation for FAST Indices
  # -----------------------------------------------------------------------
  Lambda_rot <- Lambda_raw
  F_rot <- F_raw
  V <- diag(k) # Default to identity matrix for k=1
  
  if (k > 1) {
    svd_res <- svd(Lambda_raw)
    V <- svd_res$v
    Lambda_rot <- Lambda_raw %*% V
    F_rot <- F_raw %*% V
  }
  
  # Apply "mean" rotation: Ensure rotated Factor 1 represents positive Overall Performance
  if (mean(Lambda_rot[, 1]) < 0) {
    Lambda_rot[, 1] <- -Lambda_rot[, 1]
    F_rot[, 1] <- -F_rot[, 1]
    V[, 1] <- -V[, 1]
  }
  
  mean_load_rot <- colMeans(Lambda_rot)
  
  # -----------------------------------------------------------------------
  # 4. Calculate OP (Factor 1) and Stability (Factors 2:k)
  # -----------------------------------------------------------------------
  # OP is the Rotated Factor 1 scaled by the mean rotated loading
  OP <- F_rot[, 1] * mean_load_rot[1]
  Var_A_OP <- mean_load_rot[1]^2
  
  # Calculate PEV of the rotated Factor 1 OP using the SVD rotation weights
  # Var(V[1,1]*f_1 + V[2,1]*f_2 ...) = sum( V[j,1]^2 * PEV(f_j) )
  PEV_F1_rot <- rowSums(PEV_raw * (matrix(V[, 1]^2, nrow = nrow(PEV_raw), ncol = k, byrow = TRUE)))
  PEV_OP <- PEV_F1_rot * Var_A_OP
  
  blups_gs <- data.frame(
    id = clean_ids_master,
    GEBV_OP = OP,
    PEV_OP = PEV_OP,
    Reliability = 1 - (PEV_OP / Var_A_OP),
    Stability_RMSD = 0,
    stringsAsFactors = FALSE
  )
  blups_gs$Reliability <- ifelse(blups_gs$Reliability < 0, 0, blups_gs$Reliability)
  
  # Calculate RMSD for Stability from Rotated Factors 2:k
  if (k > 1) {
    cli::cli_alert_info("Calculating crossover Stability (RMSD) from rotated higher-order factors...")
    
    L_int <- Lambda_rot[, 2:k, drop = FALSE]
    F_int <- F_rot[, 2:k, drop = FALSE]
    
    # Interaction matrix for all genotypes across environments
    I_mat <- F_int %*% t(L_int)
    
    # RMSD is the Root Mean Square Deviation of these specific interactions
    blups_gs$Stability_RMSD <- sqrt(rowMeans(I_mat^2))
  }
  
  # -----------------------------------------------------------------------
  # 5. Assign Genomic Selection Status
  # -----------------------------------------------------------------------
  target_ids <- as.character(target_ids)
  training_ids <- as.character(training_ids)
  
  blups_gs$status <- NA_character_
  blups_gs$status[blups_gs$id %in% training_ids] <- "Training"
  blups_gs$status[blups_gs$id %in% target_ids]   <- "Target"
  
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
  
  # Rank targets by Smith-Cullis Index (OP - RMSD) if k > 1, else just OP
  if (k > 1) {
    blups_assigned$Selection_Index <- blups_assigned$GEBV_OP - blups_assigned$Stability_RMSD
  } else {
    blups_assigned$Selection_Index <- blups_assigned$GEBV_OP
  }
  
  blups_target <- blups_assigned[blups_assigned$status == "Target", ]
  blups_target <- blups_target[order(blups_target$Selection_Index, decreasing = TRUE), ]
  
  # Diagnostics
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
  
  # -----------------------------------------------------------------------
  # 6. Visualization Suite
  # -----------------------------------------------------------------------
  cli::cli_alert_info("Generating enhanced diagnostic plots...")
  blups_assigned$status <- factor(blups_assigned$status, levels = c("Training", "Target"))
  
  mean_train <- mean(blups_assigned$GEBV_OP[blups_assigned$status == "Training"], na.rm = TRUE)
  mean_target <- mean(blups_assigned$GEBV_OP[blups_assigned$status == "Target"], na.rm = TRUE)
  
  # Plot A: GEBV Distribution
  plot_A <- ggplot(blups_assigned, aes(x = .data$GEBV_OP, fill = .data$status)) +
    geom_density(alpha = 0.6, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Training" = "#56B4E9", "Target" = "#E69F00")) +
    geom_vline(xintercept = mean_train, linetype = "dashed", color = "#56B4E9", linewidth = 1) +
    geom_vline(xintercept = mean_target, linetype = "dashed", color = "#E69F00", linewidth = 1) +
    labs(title = "A. Distribution of Overall Performance (OP)",
         x = "Genomic Estimated Breeding Value (Rotated Factor 1)", y = "Density", fill = "Cohort") +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")
  
  # Plot B & C depend on k
  top_cutoff <- quantile(blups_assigned$GEBV_OP[blups_assigned$status == "Target"], top_pct, na.rm = TRUE)
  top_candidates <- head(blups_target, 5) # For labeling
  
  if (k > 1) {
    plot_B <- ggplot(blups_assigned, aes(x = .data$status, y = .data$Stability_RMSD, fill = .data$status)) +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.8) +
      scale_fill_manual(values = c("Training" = "#56B4E9", "Target" = "#E69F00")) +
      labs(title = paste("B. Stability (RMSD of Factors 2 to", k, ")"),
           x = "Cohort", y = "RMSD (Lower = More Stable)") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")
    
    plot_C <- ggplot(blups_assigned, aes(x = .data$Stability_RMSD, y = .data$GEBV_OP)) +
      geom_point(aes(color = .data$status), alpha = 0.7, size = 2.5, stroke = 0) +
      scale_color_manual(values = c("Training" = "#56B4E9", "Target" = "#E69F00")) +
      geom_hline(yintercept = top_cutoff, linetype = "dotted", color = "black", linewidth = 1) +
      ggrepel::geom_label_repel(data = top_candidates, aes(label = .data$id), 
                                size = 3, fontface = "bold", box.padding = 0.5, 
                                point.padding = 0.3, segment.color = "grey50", max.overlaps = Inf) +
      labs(title = "C. Breeder's Selection Space",
           x = "Crossover GxE (Stability RMSD)", y = "Genomic Merit (GEBV_OP)", color = "Cohort") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")
      
  } else {
    plot_B <- ggplot(blups_assigned, aes(x = .data$status, y = .data$Reliability, fill = .data$status)) +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.8) +
      scale_fill_manual(values = c("Training" = "#56B4E9", "Target" = "#E69F00")) +
      labs(title = "B. Prediction Reliability",
           x = "Cohort", y = "Reliability (Accuracy²)") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")
    
    plot_C <- ggplot(blups_assigned, aes(x = .data$Reliability, y = .data$GEBV_OP)) +
      geom_point(aes(color = .data$status), alpha = 0.7, size = 2.5, stroke = 0) +
      scale_color_manual(values = c("Training" = "#56B4E9", "Target" = "#E69F00")) +
      geom_hline(yintercept = top_cutoff, linetype = "dotted", color = "black", linewidth = 1) +
      ggrepel::geom_label_repel(data = top_candidates, aes(label = .data$id), 
                                size = 3, fontface = "bold", box.padding = 0.5, 
                                point.padding = 0.3, segment.color = "grey50", max.overlaps = Inf) +
      labs(title = "C. Breeder's Selection Space",
           x = "Prediction Reliability", y = "Genomic Merit (GEBV_OP)", color = "Cohort") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")
  }
  
  combined_figure <- plot_A / (plot_B | plot_C) +
    patchwork::plot_annotation(
      title = 'Genomic Prediction Validation and Selection Differentials',
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  if (save_output) {
    pdf_file <- paste0(output_prefix, "_Figure.pdf")
    ggsave(filename = pdf_file, plot = combined_figure, device = "pdf", 
           width = 11, height = 9, dpi = 600)
    cli::cli_alert_success("Saved enhanced manuscript figure to {.file {pdf_file}}")
  }
  
  return(list(
    blups = blups_assigned,
    missing_targets = missing_targets,
    plots = combined_figure
  ))
}
