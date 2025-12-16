#' Master Visualization Suite for Factor Analytic Models
#'
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
  } else {
    stop("Unknown type.")
  }
}

# --- Internals ---

.plot_fast <- function(x, n, h) {
  if (is.null(x$fast)) stop("No FAST data available.")
  df <- x$fast

  # Sorting and labeling logic
  df <- df[order(df$OP, decreasing = TRUE), ]
  top_gens <- head(df$Genotype, n)

  df$Type <- "Other"
  df$Type[df$Genotype %in% top_gens] <- "Top"
  if (!is.null(h)) df$Type[df$Genotype %in% h] <- "Highlight"

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
    theme_genetics() +
    theme(legend.position = "bottom")

  # Use ggrepel if available
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
  # Agnostic Labels
  grp_name <- x$meta$group

  r_names <- rownames(cor_mat)
  c_names <- colnames(cor_mat)

  df_melt <- expand.grid(Var1 = r_names, Var2 = c_names, stringsAsFactors = TRUE)
  df_melt$Correlation <- as.vector(cor_mat)

  # Enforce order
  df_melt$Var1 <- factor(df_melt$Var1, levels = r_names)
  df_melt$Var2 <- factor(df_melt$Var2, levels = c_names)

  ggplot(df_melt, aes(x = Var2, y = Var1, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
    scale_fill_distiller(palette = "RdYlGn", limit = c(-1, 1), direction = 1) +
    theme_genetics() +
    labs(x = NULL, y = NULL, title = paste("Genetic Correlation:", grp_name)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

.plot_reg <- function(x, fac, h, n) {
  if (is.null(x$scores)) stop("No Genotype scores available.")

  k <- x$meta$k
  grp_name <- x$meta$group
  if (is.null(fac)) fac <- 1:k

  L <- x$loadings$rotated
  S <- x$scores$rotated
  # Reconstruct Interaction (Genotype x Group)
  # Actually, latent regression plot is usually:
  # X-axis: Environmental Loading for Factor f
  # Y-axis: Genotype Score for Factor f * Loading ?? No.
  # Standard GGE/AMMI latent regression:
  # Y_ij = mu + ... + u_if * v_jf
  # We plot regressions of (Y_ij) against (v_jf). Slope is u_if.
  # Here Y is the interaction deviation captured by this factor.
  # Interaction_f = u_if * v_jf.
  # If we plot Interaction_f vs v_jf, the slope is u_if.

  # Interaction matrix for ALL factors
  # GxE = F L'
  # But we want specifically factor f component?
  # Component_f = F[,f] %*% t(L[,f])

  # Wait, the code provided in the prompt calculates:
  # Uc = S %*% t(L) (Total interaction + main effect if included?)
  # Then subtracts previous factors.

  Uc <- S %*% t(L)

  all_slopes <- data.frame()
  all_points <- data.frame()

  tgt <- unique(c(head(rownames(S), n), h))
  gens_saf <- intersect(tgt, rownames(S))

  for (f in fac) {
    if (f > k) next
    # Residual after removing 1..f-1
    Y <- if (f > 1) Uc - (S[, 1:(f - 1), drop = FALSE] %*% t(L[, 1:(f - 1), drop = FALSE])) else Uc
    # Actually, we just want the f-th component to show the regression?
    # Or do we want the "Total" fitted value against the loading?
    # Usually "Latent Regression" plots the interaction effect of Factor f against Loading f.
    # Effect = score * loading.
    # Points line up perfectly on lines passing through origin with slope = score.
    # This is tautological if we just plot (s*l) vs l.
    # But usually we plot the "fitted GxE" vs Loading.

    # The prompt implementation:
    # Y is the residual of Uc after removing 1..f-1.
    # For f=1, Y = Uc.
    # For f=2, Y = Uc - Comp1.

    xv <- L[, f]

    # Slopes
    slopes_sub <- data.frame(Genotype = gens_saf, Slope = S[gens_saf, f], Factor = paste("Factor", f), stringsAsFactors = FALSE)
    slopes_sub$ColorGroup <- "Top"
    if (!is.null(h)) slopes_sub$ColorGroup[slopes_sub$Genotype %in% h] <- "Highlight"

    # Points
    y_sub <- Y[gens_saf, , drop = FALSE]
    n_g <- length(gens_saf)
    n_e <- length(xv)

    pts_sub <- data.frame(
      Genotype = rep(gens_saf, times = n_e),
      Group = rep(rownames(L), each = n_g), # "Group" agnostic name
      Effect = as.vector(y_sub),
      stringsAsFactors = FALSE
    )

    load_map <- setNames(xv, rownames(L))
    pts_sub$Loading <- load_map[pts_sub$Group]
    pts_sub$Factor <- paste("Factor", f)
    pts_sub$ColorGroup <- "Top"
    if (!is.null(h)) pts_sub$ColorGroup[pts_sub$Genotype %in% h] <- "Highlight"

    slopes_sub$Min <- min(xv, na.rm = TRUE)
    slopes_sub$Max <- max(xv, na.rm = TRUE)

    all_slopes <- rbind(all_slopes, slopes_sub)
    all_points <- rbind(all_points, pts_sub)
  }

  ggplot() +
    geom_point(data = all_points, aes(x = Loading, y = Effect, color = ColorGroup), alpha = 0.5) +
    geom_abline(data = all_slopes, aes(intercept = 0, slope = Slope, color = ColorGroup), alpha = 0.8) +
    facet_wrap(~Factor, scales = "free") +
    scale_color_manual(values = c("Highlight" = "red", "Top" = "navy")) +
    theme_genetics() +
    labs(title = "Latent Regression", x = paste(grp_name, "Loading"), y = "Interaction Deviation")
}

.plot_biplot <- function(x, fac, h) {
  if (is.null(x$scores)) stop("No Genotype scores available.")

  L <- x$loadings$rotated[, fac]
  S <- x$scores$rotated[, fac]

  df_scores <- as.data.frame(S)
  colnames(df_scores) <- c("X", "Y")
  df_scores$Genotype <- rownames(S)
  df_scores$Type <- "Normal"
  if (!is.null(h)) df_scores$Type[df_scores$Genotype %in% h] <- "Highlight"

  sf <- max(abs(S)) / max(abs(L)) * 0.8
  df_load <- as.data.frame(L * sf)
  colnames(df_load) <- c("X", "Y")
  df_load$Label <- rownames(L) # Agnostic

  hull_idx <- chull(df_scores$X, df_scores$Y)

  p <- ggplot(df_scores, aes(x = X, y = Y)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    # Convex hull for context
    # geom_polygon(data = df_scores[hull_idx, ], fill = NA, color = "navy", linetype = "dashed") +
    # Removed polygon to cleaner look or use it if requested

    geom_point(aes(fill = Type, size = Type, shape = Type), color = "black") +
    scale_fill_manual(values = c("Highlight" = "red", "Normal" = "grey90")) +
    scale_shape_manual(values = c("Highlight" = 23, "Normal" = 21)) +
    scale_size_manual(values = c("Highlight" = 3, "Normal" = 2)) +
    geom_segment(
      data = df_load, aes(x = 0, y = 0, xend = X, yend = Y),
      arrow = arrow(length = unit(0.2, "cm")), color = "darkgreen"
    ) +
    geom_text(data = df_load, aes(label = Label), color = "darkgreen", size = 3, hjust = -0.2) +
    labs(x = paste("Factor", fac[1]), y = paste("Factor", fac[2]), title = "Biplot") +
    theme_genetics() +
    coord_fixed()

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    high <- df_scores[df_scores$Type == "Highlight", ]
    if (nrow(high) > 0) p <- p + ggrepel::geom_text_repel(data = high, aes(label = Genotype))
  }
  return(p)
}

.plot_vaf <- function(x) {
  df <- x$var_comp$vaf
  # Agnostic: The first column is always the grouping variable in fit_fa_model output
  grp_col <- names(df)[1]

  # Sort
  # Sort by Total VAF if available
  if ("Total_VAF" %in% names(df)) {
    df <- df[order(df$Total_VAF), ]
  } else {
    df <- df[order(df[[grp_col]]), ]
  }

  # Ensure Factor
  df$Group <- factor(df[[grp_col]], levels = df[[grp_col]])

  # Melt for stacked bar
  # We only want VAF_Fac columns
  fac_cols <- grep("VAF_Fac", names(df), value = TRUE)

  # Base R melt
  df_long <- reshape(
    df,
    direction = "long",
    varying = fac_cols,
    v.names = "VAF",
    timevar = "Factor",
    times = fac_cols,
    idvar = "Group",
    new.row.names = 1:10000 # Safety
  )

  # Clean Factor names
  df_long$Factor <- gsub("VAF_", "", df_long$Factor)

  ggplot(df_long, aes(x = VAF, y = Group, fill = Factor)) +
    geom_col() +
    geom_vline(xintercept = 80, linetype = "dashed", alpha = 0.5) +
    scale_fill_brewer(palette = "Blues") +
    labs(title = paste("Variance Accounted For by", x$meta$group), x = "% Variance", y = NULL) +
    theme_minimal()
}

.plot_dopt <- function(d) {
  df <- d$site_impact
  df <- df[order(df$Impact_Pct), ]
  df$Site <- factor(df$Site, levels = df$Site)

  ggplot(df, aes(x = Impact_Pct, y = Site)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = sprintf("%.1f%%", Impact_Pct)), hjust = -0.2, size = 3) +
    labs(title = "D-Optimality Contribution", x = "% Information Loss") +
    theme_minimal()
}

.plot_diff <- function(res, n, h) {
  df <- res$gen_effects
  df <- na.omit(df)

  top_n <- head(df$Genotype, n)
  bot_n <- tail(df$Genotype, n)
  tgt <- unique(c(top_n, bot_n, h))

  sub_df <- df[df$Genotype %in% tgt, ]

  plot_df <- data.frame(
    Genotype = rep(sub_df$Genotype, 2),
    Class = rep(c("Pred_Neg", "Pred_Pos"), each = nrow(sub_df)),
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
    scale_color_manual(values = c("Highlight" = "red", "Normal" = "grey50")) +
    scale_linewidth_manual(values = c("Highlight" = 1, "Normal" = 0.5)) +
    theme_genetics() +
    labs(title = "Crossover Interaction", x = NULL, y = "Genetic Value") +
    theme(legend.position = "none")
}

# Export these legacy/helper plotters
#' @export
plot_spatial <- function(input, row = "Row", col = "Column", attribute = "Yield") {
  if (inherits(input, "asreml")) {
    data_name <- as.character(input$call$data)
    if (!exists(data_name)) stop("Original dataframe not found.")
    df <- get(data_name)
    df$Value_To_Plot <- resid(input)
    main_title <- "Residual Map"
  } else {
    df <- input
    df$Value_To_Plot <- df[[attribute]]
    main_title <- paste("Spatial Map:", attribute)
  }

  ggplot(df, aes(x = .data[[col]], y = .data[[row]], fill = Value_To_Plot)) +
    geom_tile() +
    scale_fill_distiller(palette = "Spectral") +
    theme_genetics() +
    coord_fixed() +
    labs(title = main_title)
}

#' @export
plot_trial_map <- function(data, ...) plot_spatial(data, ...)

#' @export
plot_met_trend <- function(data, x = "Year", y = "Yield", main = "Yield Trend", ...) {
  ggplot(data, aes(x = factor(.data[[x]]), y = .data[[y]])) +
    geom_boxplot(fill = "lightgreen", alpha = 0.5) +
    geom_smooth(aes(x = as.numeric(factor(.data[[x]])), y = .data[[y]]), method = "lm", se = FALSE, color = "red") +
    theme_genetics() +
    labs(title = main, x = x, y = y)
}
