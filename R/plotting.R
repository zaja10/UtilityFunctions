#' Master Visualization Suite for Factor Analytic Models
#'
#' @description
#' The central plotting interface for the toolkit. This S3 method dispatches to
#' specific visualization engines based on the \code{type} argument, covering the
#' full analysis lifecycle: from Site Quality Control to Genotype Selection and
#' Mechanistic GxE analysis.
#'
#' It uses **ggplot2** to produce publication-quality figures.
#'
#' @param x An object of class \code{fa_model} produced by \code{fit_fa_model()}.
#' @param type A character string specifying the visualization module. Options:
#' \itemize{
#'   \item \code{"fast"} (Default): **Factor Analytic Selection Tools**. Plots Overall Performance (OP) vs Stability (RMSD).
#'   \item \code{"heatmap"}: **Genetic Correlation Matrix**. Visualizes the relationship structure between environments.
#'   \item \code{"latent_reg"}: **Latent Regression**. Visualizes the drivers of GxE by plotting genotype slopes against environmental loadings.
#'   \item \code{"biplot"}: **GxE Biplot**. A 2D projection of Genotypes (Scores) and Environments (Vectors).
#'   \item \code{"vaf"}: **Site Quality**. A bar chart of Variance Accounted For (\%) per site.
#'   \item \code{"d_opt"}: **Network Efficiency**. A bar chart of D-Optimality "Information Loss" (requires \code{get_d_optimality}).
#'   \item \code{"diff"}: **Crossover Interaction**. A slope graph showing rank changes between Interaction Classes (requires \code{get_i_classes}).
#'   \item \code{"h2"}: **Reliability/Repeatability**. A bar chart comparing Cullis/Standard metrics per site. Requires output from \code{compare_h2()}.
#' }
#' @param factor Integer or Vector. Controls which factors are visualized.
#' \itemize{
#'   \item For \code{"latent_reg"}: Can be a single integer to see specific stability drivers, or \code{NULL} to plot all factors.
#'   \item For \code{"biplot"}: A vector of length 2 (e.g., \code{c(1, 2)}) specifying the X and Y axes.
#' }
#' @param n_label Integer. The number of top-performing genotypes (sorted by OP) to automatically label. Default is 5.
#' @param highlight Character vector. Names of specific genotypes to highlight (e.g., Check varieties).
#' @param ... Additional arguments passed to ggplot2 themes or internally.
#'
#' @return A ggplot object.
#'
#' @seealso \code{\link{fit_fa_model}}, \code{\link{get_d_optimality}}, \code{\link{get_i_classes}}, \code{\link{compare_h2}}
#'
#' @import ggplot2
#' @importFrom dplyr filter mutate arrange slice_head inner_join left_join select case_when
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_vline geom_tile geom_text labs theme_minimal theme facet_wrap scale_fill_gradient2 coord_fixed scale_fill_manual scale_color_manual element_text
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

# --- Plot Internals ---

.plot_fast <- function(x, n, h) {
  if (is.null(x$fast)) stop("No FAST data available (Genotype scores missing).")

  df <- x$fast

  # Improve Labeling Logic
  df <- df %>%
    mutate(
      Type = case_when(
        Genotype %in% h ~ "Highlight",
        Genotype %in% head(df$Genotype, n) ~ "Top",
        TRUE ~ "Other"
      )
    )

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

  # Add Labels
  label_df <- df %>% filter(Type != "Other")
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

  # Melt for ggplot
  df_melt <- as.data.frame(cor_mat) %>%
    rownames_to_column("Env1") %>%
    pivot_longer(cols = -Env1, names_to = "Env2", values_to = "Correlation")

  # Order factors to match matrix order
  df_melt$Env1 <- factor(df_melt$Env1, levels = rownames(cor_mat))
  df_melt$Env2 <- factor(df_melt$Env2, levels = colnames(cor_mat))

  ggplot(df_melt, aes(x = Env2, y = Env1, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
    scale_fill_distiller(palette = "RdYlGn", limit = c(-1, 1), direction = 1) +
    theme_minimal() +
    labs(x = NULL, y = NULL, title = "Genetic Correlation Matrix") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

.plot_reg <- function(x, fac, h, n) {
  if (is.null(x$scores)) stop("No Genotype scores available for Latent Regression.")

  k <- x$meta$k
  if (is.null(fac)) fac <- 1:k

  L <- x$loadings$rotated
  S <- x$scores$rotated
  Uc <- S %*% t(L)

  # Prepare data for plotting
  plot_data <- list()

  # Determine target genotypes
  all_gens <- rownames(S)
  top_n <- head(all_gens, n)
  tgt <- unique(c(top_n, h))

  for (f in fac) {
    if (f > k) next

    # Calculate Y (deviations for f > 1)
    Y <- if (f > 1) Uc - (S[, 1:(f - 1)] %*% t(L[, 1:(f - 1)])) else Uc
    xv <- L[, f]

    # We need to construct a long dataframe: Genotype, Env_Load, Value, Factor
    # But strictly, we only plot lines for target genotypes.

    # Dataframe for segments/lines
    # For each genotype g, intercept = 0, slope = score[g, f]
    # x range is range(xv)

    gens_saf <- intersect(tgt, rownames(S))

    slopes <- data.frame(
      Genotype = gens_saf,
      Slope = S[gens_saf, f],
      Factor = paste("Factor", f)
    ) %>%
      mutate(ColorGroup = ifelse(Genotype %in% h, "Highlight", "Top"))

    # Create points for specific environments?
    # The original plot shows the actual "Y" values as points at the xv locations
    # So we need Y[g, ] vs xv

    # Melt Y[tgt, ]
    y_sub <- Y[gens_saf, , drop = FALSE]

    pts_df <- as.data.frame(y_sub) %>%
      rownames_to_column("Genotype") %>%
      pivot_longer(-Genotype, names_to = "Env", values_to = "Effect")

    # Join with loadings (xv)
    # Loadings is a matrix, convert to df
    load_df <- data.frame(Env = rownames(L), Loading = xv)

    pts_df <- pts_df %>%
      left_join(load_df, by = "Env") %>%
      mutate(Factor = paste("Factor", f)) %>%
      mutate(ColorGroup = ifelse(Genotype %in% h, "Highlight", "Top"))

    plot_data[[as.character(f)]] <- list(slopes = slopes, points = pts_df)
  }

  # Combine
  all_slopes <- do.call(rbind, lapply(plot_data, `[[`, "slopes"))
  all_points <- do.call(rbind, lapply(plot_data, `[[`, "points"))

  # Calculate line endpoints for geom_segment to mimic abline
  # x range per factor
  ranges <- all_points %>%
    group_by(Factor) %>%
    summarize(Min = min(Loading), Max = max(Loading))

  all_slopes <- all_slopes %>%
    left_join(ranges, by = "Factor")

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
  if (is.null(x$scores)) stop("No Genotype scores available for Biplot.")

  L <- x$loadings$rotated[, fac]
  S <- x$scores$rotated[, fac]

  # Convert to dataframes
  df_scores <- as.data.frame(S)
  colnames(df_scores) <- c("X", "Y")
  df_scores$Genotype <- rownames(S)
  df_scores$Type <- ifelse(df_scores$Genotype %in% h, "Highlight", "Normal")

  # Scaling
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
    # Hull
    geom_polygon(data = hull_df, fill = NA, color = "navy", linetype = "dashed") +
    # Genotypes
    geom_point(aes(fill = Type, size = Type, shape = Type), color = "black") +
    scale_fill_manual(values = c("Highlight" = "red", "Normal" = "grey90")) +
    scale_shape_manual(values = c("Highlight" = 23, "Normal" = 21)) +
    scale_size_manual(values = c("Highlight" = 3, "Normal" = 2)) +
    # Vectors
    geom_segment(
      data = df_load, aes(x = 0, y = 0, xend = X, yend = Y),
      arrow = arrow(length = unit(0.2, "cm")), color = "darkgreen"
    ) +
    geom_text(data = df_load, aes(label = Env), color = "darkgreen", size = 3, hjust = -0.2) +
    labs(x = paste("Factor", fac[1]), y = paste("Factor", fac[2]), title = "GxE Biplot") +
    theme_minimal() +
    coord_fixed()

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = filter(df_scores, Type == "Highlight"), aes(label = Genotype), box.padding = 0.5)
  }

  return(p)
}

.plot_vaf <- function(x) {
  df <- x$var_comp$vaf
  df$Site <- factor(df$Site, levels = df$Site[order(df$VAF)])

  ggplot(df, aes(x = VAF, y = Site)) +
    geom_col(aes(fill = VAF < 50)) +
    geom_vline(xintercept = 80, linetype = "dashed") +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"), guide = "none") +
    labs(title = "Site Quality (Variance Accounted For)", x = "VAF %") +
    theme_minimal()
}

.plot_dopt <- function(d, threshold = 1.0) {
  df <- d$site_impact
  df$Site <- factor(df$Site, levels = df$Site[order(df$Impact_Pct)])

  ggplot(df, aes(x = Impact_Pct, y = Site)) +
    geom_col(aes(fill = Impact_Pct >= threshold)) +
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
  tgt <- unique(c(head(df$Genotype, n), tail(df$Genotype, n), h))
  plot_df <- df %>%
    filter(Genotype %in% tgt) %>%
    pivot_longer(cols = c(Pred_Neg, Pred_Pos), names_to = "Class", values_to = "Value")

  plot_df$Class <- factor(plot_df$Class, levels = c("Pred_Neg", "Pred_Pos"), labels = c("Class (-)", "Class (+)"))
  plot_df$Type <- ifelse(plot_df$Genotype %in% h, "Highlight", "Normal")

  ggplot(plot_df, aes(x = Class, y = Value, group = Genotype, color = Type)) +
    geom_line(aes(linewidth = Type)) +
    geom_point() +
    geom_text(data = plot_df %>% filter(Class == "Class (-)"), aes(label = Genotype), hjust = 1.1, size = 3) +
    geom_text(data = plot_df %>% filter(Class == "Class (+)"), aes(label = Genotype), hjust = -0.1, size = 3) +
    scale_color_manual(values = c("Highlight" = "red", "Normal" = "grey50")) +
    scale_linewidth_manual(values = c("Highlight" = 1, "Normal" = 0.5)) +
    theme_minimal() +
    labs(title = "Crossover Interaction", x = NULL, y = "Genetic Value") +
    theme(legend.position = "none")
}


#' Plot Spatial Field Map (Robust)
#'
#' Maps raw data (from dataframe) or residuals (from model) to the physical field layout.
#'
#' @param input An `asreml` model object OR a `data.frame`.
#' @param row Column name for Field Row.
#' @param col Column name for Field Column.
#' @param attribute String. What to plot? "Yield", "residuals", "fitted".
#' @export
plot_spatial <- function(input, row = "Row", col = "Column", attribute = "Yield") {
  # Data Extraction Logic
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

  # Ensure Numeric
  df[[row]] <- as.numeric(as.character(df[[row]]))
  df[[col]] <- as.numeric(as.character(df[[col]]))

  # Determine bounds
  # We might want to fill missing spots explicitly, but geom_tile handles NA by creating holes,
  # or we can use complete().
  # Better to use complete() to show the full grid structure clearly.

  # Only complete if we have strict limits?
  # Let's trust ggplot to plot what's there, but maybe implicit missingness is hard to see.
  # Let's treat implicit missing as white/grey.

  ggplot(df, aes(x = .data[[col]], y = .data[[row]], fill = Value_To_Plot)) +
    geom_tile() +
    scale_fill_distiller(palette = "Spectral", na.value = "grey90", direction = -1) +
    theme_minimal() +
    labs(title = main_title, x = col, y = row) +
    coord_fixed()
}

#' Plot Heritability/Reliability Comparison
#'
#' @param x An object of class \code{h2_comparison}.
#' @param ... Additional arguments.
#' @export
plot.h2_comparison <- function(x, ...) {
  if (!all(c("Site", "Method", "Value") %in% names(x))) stop("Invalid h2_comparison object.")

  # Order sites by the first method?
  # Let's find mean value per site to order
  site_order <- x %>%
    group_by(Site) %>%
    summarize(MeanVal = mean(Value, na.rm = TRUE)) %>%
    arrange(MeanVal) %>%
    pull(Site)

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
    # Spatial
    ggplot(plot_data, aes(x = .data[[col]], y = .data[[row]], fill = .data[[val]])) +
      geom_tile() +
      scale_fill_distiller(palette = "Spectral", direction = -1) +
      theme_minimal() +
      labs(title = paste("Trial Map:", title_sub), x = col, y = row) +
      coord_fixed()
  } else {
    # 1D Fallback
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
  if (!all(c(x, y) %in% names(data))) stop("Variables not found.")

  # Ensure numeric X for trend line
  data$X_Num <- if (is.numeric(data[[x]])) data[[x]] else as.numeric(as.character(data[[x]]))

  # Combined Plot: Boxplot + Trend Line
  ggplot(data, aes(x = factor(.data[[x]]), y = .data[[y]])) +
    geom_boxplot(fill = "lightgreen", alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.1) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "blue", linewidth = 1.2) +
    stat_summary(fun = mean, geom = "point", color = "blue", size = 3) +
    geom_smooth(aes(x = as.numeric(factor(.data[[x]])), y = .data[[y]]), method = "lm", color = "red", linetype = "dashed", se = FALSE, linewidth = 1.5) +
    theme_minimal() +
    labs(title = main, x = x, y = y)
}
