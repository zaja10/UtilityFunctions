#' Master Visualization Suite for Factor Analytic Models
#'
#' The central plotting interface for the toolkit.
#' Uses **ggplot2** for rendering but Base R for data preparation to minimize dependencies.
#'
#' @param x An object of class \code{fa_model}.
#' @param type Character. Visualization type: "fast", "heatmap", "latent_reg", "biplot", "vaf", "d_opt", "diff".
#' @param factor Integer/Vector. Factors to visualize.
#' @param n_label Integer. Number of top genotypes to label.
#' @param highlight Character vector. Genotypes to highlight.
#' @param ... Additional arguments.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom stats reshape setNames aggregate
#' @export
plot.fa_model <- function(x, type = "fast", factor = NULL, n_label = 5, highlight = NULL, ...) {
  if (type == "fast") {
    .plot_fast(x, n_label, highlight)
  } else if (type == "heatmap") {
    .plot_heat(x)
  } else if (type == "latent_reg") {
    .plot_reg(x, factor, highlight, n_label)
  } else if (type == "biplot") {
    .plot_biplot(x, if (is.null(factor)) c(1, 2) else factor, highlight)
  } else if (type == "vaf") {
    .plot_vaf(x)
  } else if (type == "d_opt") {
    .plot_dopt(calculate_d_optimality(x))
  } else if (type == "diff") {
    .plot_diff(calculate_i_classes(x, if (is.null(factor)) 2 else factor), n_label, highlight)
  } else if (type == "h2") {
    stop("For reliability/heritability plots, please run 'compare_h2()' first, then plot the result.")
  } else {
    stop("Unknown type.")
  }
}

# --- Internals (Refactored to Base R Data Prep) ---

.plot_fast <- function(x, n, h) {
  if (is.null(x$fast)) stop("No FAST data available.")
  df <- x$fast

  # Label Logic (Base R)
  # Sort by OP to find Top N
  df <- df[order(df$OP, decreasing = TRUE), ]
  top_gens <- head(df$Genotype, n)

  df$Type <- "Other"
  df$Type[df$Genotype %in% top_gens] <- "Top"
  if (!is.null(h)) df$Type[df$Genotype %in% h] <- "Highlight"

  # Set Factor levels for legend order
  df$Type <- factor(df$Type, levels = c("Highlight", "Top", "Other"))

  mean_rmsd <- mean(df$RMSD, na.rm = TRUE)

  p <- ggplot(df, aes(x = RMSD, y = OP)) +
    geom_vline(xintercept = mean_rmsd, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(fill = Type, size = Type, shape = Type), color = "black") +
    scale_fill_manual(values = c("Highlight" = "#e74c3c", "Top" = "#3498db", "Other" = "grey95")) +
    scale_shape_manual(values = c("Highlight" = 23, "Top" = 21, "Other" = 21)) +
    scale_size_manual(values = c("Highlight" = 3, "Top" = 3, "Other" = 2)) +
    labs(x = "Stability (RMSD)", y = "Performance (OP)", title = "FAST Selection") +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Labels
  label_df <- df[df$Type != "Other", ]
  if (nrow(label_df) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = label_df, aes(label = Genotype), size = 3, max.overlaps = 20)
    } else {
      p <- p + geom_text(data = label_df, aes(label = Genotype), vjust = -1, size = 3)
    }
  }
  return(p)
}

.plot_heat <- function(x) {
  cor_mat <- x$matrices$Cor

  # Melt Matrix (Base R)
  # Expand grid of indices
  r_names <- rownames(cor_mat)
  c_names <- colnames(cor_mat)

  df_melt <- expand.grid(Env1 = r_names, Env2 = c_names, stringsAsFactors = TRUE)
  # Extract values using matrix indexing
  df_melt$Correlation <- as.vector(cor_mat)

  # Ensure Factor Order matches matrix order
  df_melt$Env1 <- factor(df_melt$Env1, levels = r_names)
  df_melt$Env2 <- factor(df_melt$Env2, levels = c_names)

  ggplot(df_melt, aes(x = Env2, y = Env1, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
    scale_fill_distiller(palette = "RdYlGn", limit = c(-1, 1), direction = 1) +
    theme_minimal() +
    labs(x = NULL, y = NULL, title = "Genetic Correlation Matrix") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

.plot_reg <- function(x, fac, h, n) {
  if (is.null(x$scores)) stop("No Genotype scores available.")

  k <- x$meta$k
  if (is.null(fac)) fac <- 1:k

  L <- x$loadings$rotated
  S <- x$scores$rotated
  Uc <- S %*% t(L)

  all_slopes <- data.frame()
  all_points <- data.frame()

  # Determine target genotypes
  all_gens <- rownames(S)
  top_n <- head(all_gens, n) # Assumes S is sorted/consistent with FAST
  tgt <- unique(c(top_n, h))
  gens_saf <- intersect(tgt, rownames(S))

  for (f in fac) {
    if (f > k) next

    # Calculate Y (deviations)
    Y <- if (f > 1) Uc - (S[, 1:(f - 1), drop = FALSE] %*% t(L[, 1:(f - 1), drop = FALSE])) else Uc
    xv <- L[, f]

    # Slopes Data
    slopes_sub <- data.frame(
      Genotype = gens_saf,
      Slope = S[gens_saf, f],
      Factor = paste("Factor", f),
      stringsAsFactors = FALSE
    )
    slopes_sub$ColorGroup <- "Top"
    if (!is.null(h)) slopes_sub$ColorGroup[slopes_sub$Genotype %in% h] <- "Highlight"

    # Points Data (Long Format Manual)
    y_sub <- Y[gens_saf, , drop = FALSE]

    # Vectorize creation of long format
    n_g <- length(gens_saf)
    n_e <- length(xv)

    pts_sub <- data.frame(
      Genotype = rep(gens_saf, times = n_e),
      Env = rep(rownames(L), each = n_g),
      Effect = as.vector(y_sub), # Verify order: Y is Gen x Env. as.vector fills col by col? No, R is col-major.
      # Y is (Gen x Env). as.vector(Y) flattens Gen1_E1, Gen2_E1...
      # So we need Env to vary slowly, Gen to vary fast?
      # Actually `as.vector` on matrix goes down columns (Env1, Env2...).
      # So `each=n_g` for Env is correct.
      stringsAsFactors = FALSE
    )

    # Add Loading (X-axis)
    # L is (Env x k). xv is (Env x 1).
    # We need to map Env -> Loading.
    load_map <- setNames(xv, rownames(L))
    pts_sub$Loading <- load_map[pts_sub$Env]

    pts_sub$Factor <- paste("Factor", f)
    pts_sub$ColorGroup <- "Top"
    if (!is.null(h)) pts_sub$ColorGroup[pts_sub$Genotype %in% h] <- "Highlight"

    # Add Ranges to slopes for text placement
    min_l <- min(xv, na.rm = TRUE)
    max_l <- max(xv, na.rm = TRUE)
    slopes_sub$Min <- min_l
    slopes_sub$Max <- max_l

    all_slopes <- rbind(all_slopes, slopes_sub)
    all_points <- rbind(all_points, pts_sub)
  }

  ggplot() +
    geom_point(data = all_points, aes(x = Loading, y = Effect, color = ColorGroup), alpha = 0.5) +
    geom_abline(data = all_slopes, aes(intercept = 0, slope = Slope, color = ColorGroup), alpha = 0.8) +
    facet_wrap(~Factor, scales = "free") +
    scale_color_manual(values = c("Highlight" = "red", "Top" = "navy")) +
    theme_minimal() +
    labs(title = "Latent Regression", x = "Environmental Loading", y = "GxE Deviation") +
    geom_text(
      data = all_slopes, aes(x = Max, y = Slope * Max, label = Genotype, color = ColorGroup),
      hjust = 0, vjust = 0.5, size = 3, check_overlap = TRUE
    )
}

.plot_biplot <- function(x, fac, h) {
  if (is.null(x$scores)) stop("No Genotype scores available.")

  L <- x$loadings$rotated[, fac]
  S <- x$scores$rotated[, fac]

  # Scores DF
  df_scores <- as.data.frame(S)
  colnames(df_scores) <- c("X", "Y")
  df_scores$Genotype <- rownames(S)
  df_scores$Type <- "Normal"
  if (!is.null(h)) df_scores$Type[df_scores$Genotype %in% h] <- "Highlight"

  # Loadings DF
  sf <- max(abs(S)) / max(abs(L)) * 0.8
  df_load <- as.data.frame(L * sf)
  colnames(df_load) <- c("X", "Y")
  df_load$Env <- rownames(L)

  # Hull
  hull_idx <- chull(df_scores$X, df_scores$Y)
  hull_df <- df_scores[hull_idx, ]

  p <- ggplot(df_scores, aes(x = X, y = Y)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_polygon(data = hull_df, fill = NA, color = "navy", linetype = "dashed") +
    geom_point(aes(fill = Type, size = Type, shape = Type), color = "black") +
    scale_fill_manual(values = c("Highlight" = "red", "Normal" = "grey90")) +
    scale_shape_manual(values = c("Highlight" = 23, "Normal" = 21)) +
    scale_size_manual(values = c("Highlight" = 3, "Normal" = 2)) +
    geom_segment(
      data = df_load, aes(x = 0, y = 0, xend = X, yend = Y),
      arrow = arrow(length = unit(0.2, "cm")), color = "darkgreen"
    ) +
    geom_text(data = df_load, aes(label = Env), color = "darkgreen", size = 3, hjust = -0.2) +
    labs(x = paste("Factor", fac[1]), y = paste("Factor", fac[2]), title = "GxE Biplot") +
    theme_minimal() +
    coord_fixed()

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    # Subset in base R
    high_scores <- df_scores[df_scores$Type == "Highlight", ]
    if (nrow(high_scores) > 0) {
      p <- p + ggrepel::geom_text_repel(data = high_scores, aes(label = Genotype), box.padding = 0.5)
    }
  }
  return(p)
}

.plot_vaf <- function(x) {
  df <- x$var_comp$vaf
  # Factor Reorder
  df <- df[order(df$VAF), ]
  df$Site <- factor(df$Site, levels = df$Site)

  df$IsLow <- df$VAF < 50

  ggplot(df, aes(x = VAF, y = Site)) +
    geom_col(aes(fill = IsLow)) +
    geom_vline(xintercept = 80, linetype = "dashed") +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"), guide = "none") +
    labs(title = "Site Quality (Variance Accounted For)", x = "VAF %") +
    theme_minimal()
}

.plot_dopt <- function(d, threshold = 1.0) {
  df <- d$site_impact
  df <- df[order(df$Impact_Pct), ]
  df$Site <- factor(df$Site, levels = df$Site)

  df$IsKey <- df$Impact_Pct >= threshold

  ggplot(df, aes(x = Impact_Pct, y = Site)) +
    geom_col(aes(fill = IsKey)) +
    geom_vline(xintercept = threshold, linetype = "dashed") +
    geom_text(aes(label = sprintf("%.1f%%", Impact_Pct)), hjust = -0.2, size = 3) +
    scale_fill_manual(
      values = c("TRUE" = "#27ae60", "FALSE" = "#e74c3c"),
      name = "Status", labels = c("Redundant", "Key Site")
    ) +
    labs(title = "Network Efficiency (D-Opt Contribution)", x = "% Information Loss if Removed") +
    theme_minimal()
}

.plot_diff <- function(res, n, h) {
  df <- res$gen_effects
  df <- na.omit(df)

  # Select Targets
  top_n <- head(df$Genotype, n)
  bot_n <- tail(df$Genotype, n)
  tgt <- unique(c(top_n, bot_n, h))

  sub_df <- df[df$Genotype %in% tgt, ]

  # Pivot Long manually
  # Columns: Pred_Neg, Pred_Pos
  n_rows <- nrow(sub_df)
  plot_df <- data.frame(
    Genotype = rep(sub_df$Genotype, 2),
    Class = rep(c("Pred_Neg", "Pred_Pos"), each = n_rows),
    Value = c(sub_df$Pred_Neg, sub_df$Pred_Pos),
    stringsAsFactors = FALSE
  )

  plot_df$Class <- factor(plot_df$Class, levels = c("Pred_Neg", "Pred_Pos"), labels = c("Class (-)", "Class (+)"))

  plot_df$Type <- "Normal"
  if (!is.null(h)) plot_df$Type[plot_df$Genotype %in% h] <- "Highlight"

  ggplot(plot_df, aes(x = Class, y = Value, group = Genotype, color = Type)) +
    geom_line(aes(linewidth = Type)) +
    geom_point() +
    geom_text(data = plot_df[plot_df$Class == "Class (-)", ], aes(label = Genotype), hjust = 1.1, size = 3) +
    geom_text(data = plot_df[plot_df$Class == "Class (+)", ], aes(label = Genotype), hjust = -0.1, size = 3) +
    scale_color_manual(values = c("Highlight" = "red", "Normal" = "grey50")) +
    scale_linewidth_manual(values = c("Highlight" = 1, "Normal" = 0.5)) +
    theme_minimal() +
    labs(title = "Crossover Interaction", x = NULL, y = "Genetic Value") +
    theme(legend.position = "none")
}

#' Plot Spatial Field Map
#' @export
plot_spatial <- function(input, row = "Row", col = "Column", attribute = "Yield") {
  if (inherits(input, "asreml")) {
    if (attribute == "Yield") attribute <- "residuals"
    data_name <- as.character(input$call$data)
    if (!exists(data_name)) stop("Original dataframe not found.")
    df <- get(data_name)

    if (attribute == "fitted") {
      df$Value_To_Plot <- fitted(input)
    } else {
      df$Value_To_Plot <- resid(input)
    }
    main_title <- paste("Spatial Map:", attribute)
  } else if (inherits(input, "data.frame")) {
    df <- input
    if (is.null(df[[attribute]])) stop(paste("Column", attribute, "not found."))
    df$Value_To_Plot <- df[[attribute]]
    main_title <- paste("Spatial Map:", attribute)
  } else {
    stop("Input must be an asreml model or a dataframe.")
  }

  df[[row]] <- as.numeric(as.character(df[[row]]))
  df[[col]] <- as.numeric(as.character(df[[col]]))

  ggplot(df, aes(x = .data[[col]], y = .data[[row]], fill = Value_To_Plot)) +
    geom_tile() +
    scale_fill_distiller(palette = "Spectral", na.value = "grey90", direction = -1) +
    theme_minimal() +
    labs(title = main_title, x = col, y = row) +
    coord_fixed()
}

#' Plot Heritability Comparison
#' @export
plot.h2_comparison <- function(x, ...) {
  if (!all(c("Site", "Method", "Value") %in% names(x))) stop("Invalid h2_comparison object.")

  # Base R Aggregation
  agg <- aggregate(Value ~ Site, data = x, FUN = mean)
  agg <- agg[order(agg$Value), ]
  site_order <- agg$Site

  x$Site <- factor(x$Site, levels = site_order)

  ggplot(x, aes(x = Site, y = Value, fill = Method)) +
    geom_col(position = "dodge", alpha = 0.8) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "Site Quality Comparison", y = "Reliability / Heritability")
}

#' Flexible Trial Map
#' @export
plot_trial_map <- function(data, trial_val = NULL, trial_col = "Trial",
                           row = "Row", col = "Column", val = "Yield") {
  plot_data <- if (!is.null(trial_val)) data[data[[trial_col]] == trial_val, ] else data
  title_sub <- if (!is.null(trial_val)) paste(trial_val) else "Single Site"

  if (all(c(row, col) %in% names(plot_data))) {
    ggplot(plot_data, aes(x = .data[[col]], y = .data[[row]], fill = .data[[val]])) +
      geom_tile() +
      scale_fill_distiller(palette = "Spectral", direction = -1) +
      theme_minimal() +
      labs(title = paste("Trial Map:", title_sub), x = col, y = row) +
      coord_fixed()
  } else {
    plot_data$Index <- 1:nrow(plot_data)
    ggplot(plot_data, aes(x = Index, y = .data[[val]])) +
      geom_point(aes(color = .data[[val]])) +
      geom_line(alpha = 0.3) +
      scale_color_distiller(palette = "Spectral") +
      theme_minimal() +
      labs(title = paste("Linear Plot:", title_sub))
  }
}

#' Plot MET Trends
#' @export
plot_met_trend <- function(data, x = "Year", y = "Yield", main = "Yield Trend", ...) {
  # Base R Ensure Numeric
  data$X_Num <- if (is.numeric(data[[x]])) data[[x]] else as.numeric(as.character(data[[x]]))

  ggplot(data, aes(x = factor(.data[[x]]), y = .data[[y]])) +
    geom_boxplot(fill = "lightgreen", alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.1) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "blue", linewidth = 1.2) +
    stat_summary(fun = mean, geom = "point", color = "blue", size = 3) +
    geom_smooth(aes(x = as.numeric(factor(.data[[x]])), y = .data[[y]]), method = "lm", color = "red", linetype = "dashed", se = FALSE, linewidth = 1.5) +
    theme_minimal() +
    labs(title = main, x = x, y = y)
}
