#' Default Checks
#'
#' Standard check genotypes for plots.
#' @export
DEFAULT_CHECKS <- c("17-8930", "US16-IL-063-063", "25R76", "AgriMAXX 490", "07-19334", "AgriMAXX503")

#' Master Visualization Suite for Factor Analytic Models
#'
#' @param x Object of class `fa_model`.
#' @param type String type of plot: "fast", "heatmap", "vaf", "biplot", "latent_reg", "cluster".
#' @param n_label Integer. Number of top genotypes to label/highlight (default 5).
#' @param highlight Character vector. Specific genotypes to highlight.
#' @param interactive Logical. If TRUE, wraps the ggplot in `plotly::ggplotly()`.
#' @param ... Additional arguments passed to specific plot methods.
#' @export
plot.fa_model <- function(x,
                          type = c("fast", "heatmap", "vaf", "biplot", "latent_reg", "cluster"),
                          n_label = 5,
                          highlight = NULL,
                          interactive = FALSE,
                          ...) {
  # match.arg ensures the user provides a valid option, throwing a clean error otherwise
  type <- match.arg(type)

  p <- switch(type,
    "fast"       = .plot_fast(x, n_label, highlight),
    "heatmap"    = .plot_heat(x),
    "vaf"        = .plot_vaf(x),
    "biplot"     = .plot_biplot(x, n_label, highlight),
    "latent_reg" = .plot_latent_reg(x, n_label, highlight),
    "cluster"    = .plot_cluster(x)
  )

  if (interactive && type != "cluster") {
    if (requireNamespace("plotly", quietly = TRUE)) {
      p <- plotly::ggplotly(p)
    } else {
      warning("The 'plotly' package is required for interactive plots. Returning static plot.")
    }
  }

  return(p)
}

#' Plot Trends (Phenotypic or Genetic)
#'
#' Visualizes data over time (e.g., Year).
#' \itemize{
#'   \item "phenotypic": Boxplots + Regression (Rate of Change).
#'   \item "genetic": Scatterplot of BLUPs + Check highlighting.
#' }
#'
#' @param data Dataframe.
#' @param mode String. "phenotypic" or "genetic".
#' @param x String. X-axis column (Year).
#' @param y String. Y-axis column (Yield or BLUP).
#' @param checks Character vector. Names of checks to highlight (only for genetic mode).
#' @export
plot_trend <- function(data, mode = "phenotypic", x = "Year", y = "Yield", checks = DEFAULT_CHECKS) {
  if (mode == "phenotypic") {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = factor(.data[[x]]), y = .data[[y]])) +
      ggplot2::geom_boxplot(fill = "lightgreen", alpha = 0.6, outlier.alpha = 0.5) +
      ggplot2::geom_smooth(method = "lm", ggplot2::aes(group = 1), color = "red", se = FALSE) +
      ggplot2::labs(title = paste("Phenotypic Trend:", y, "over", x), x = x, y = y) +
      theme_genetics()
    return(p)
  } else if (mode == "genetic") {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
      ggplot2::geom_point(color = "grey", alpha = 0.5) +
      ggplot2::geom_smooth(method = "lm", color = "blue", se = TRUE, alpha = 0.2) +
      ggplot2::labs(title = paste("Genetic Trend:", y, "over", x), x = x, y = "Genetic Value") +
      theme_genetics()

    if (!is.null(checks)) {
      check_data <- data[data$Genotype %in% checks, ]
      if (nrow(check_data) > 0) {
        p <- p +
          ggplot2::geom_point(data = check_data, color = "red", shape = 17, size = 3)

        # Use ggrepel to prevent labels from overlapping
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(
            data = check_data, ggplot2::aes(label = Genotype),
            color = "darkred", box.padding = 0.5
          )
        } else {
          p <- p + ggplot2::geom_text(
            data = check_data, ggplot2::aes(label = Genotype),
            vjust = -1, color = "darkred"
          )
        }
      }
    }
    return(p)
  } else {
    stop("Unknown mode. Use 'phenotypic' or 'genetic'.")
  }
}

#' Plot Spatial Map
#'
#' Visualizes a field trial as a heatmap.
#' Wraps ASReml residuals or raw data.
#'
#' @param input Either a dataframe or `asreml` object.
#' @param row Colname for Row.
#' @param col Colname for Column.
#' @param fill Colname for value to plot. If input is model, "Residuals" is used.
#' @export
plot_spatial <- function(input, row = "Row", col = "Column", fill = "Yield") {
  if (inherits(input, "asreml")) {
    if (is.null(input$data)) {
      stop("ASReml object does not contain data slot. Pass dataframe directly.")
    }
    df <- input$data
    resids <- resid(input)

    # Ensure dimensions match (sometimes missing values are dropped from residuals)
    if (length(resids) == nrow(df)) {
      df$PlotValue <- resids
    } else {
      stop("Residuals length does not match data length. Cannot map spatially.")
    }
    title <- "Spatial Residuals"
  } else {
    df <- input
    df$PlotValue <- df[[fill]]
    title <- paste("Spatial Map:", fill)
  }

  if (!all(c(row, col) %in% names(df))) {
    stop(paste("Spatial columns not found. Expected:", row, col, "Found:", paste(names(df), collapse = ", ")))
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[col]], y = .data[[row]], fill = PlotValue)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.1) +
    ggplot2::scale_fill_distiller(palette = "Spectral", na.value = "grey90") +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = col, y = row, fill = "Value") +
    theme_genetics()
}

# -------------------------------------------------------------------------
# Internal Implementations for plot.fa_model
# -------------------------------------------------------------------------

.plot_fast <- function(x, n, h) {
  if (is.null(x$fast)) stop("No FAST data available.")
  df <- x$fast
  df <- df[order(df$OP, decreasing = TRUE), ]

  top_gens <- head(df$Genotype, n)
  df$Type <- "Other"
  df$Type[df$Genotype %in% top_gens] <- "Top Candidate"
  if (!is.null(h)) df$Type[df$Genotype %in% h] <- "Highlight"

  df$Type <- factor(df$Type, levels = c("Highlight", "Top Candidate", "Other"))
  mean_rmsd <- mean(df$RMSD, na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = RMSD, y = OP)) +
    ggplot2::geom_vline(xintercept = mean_rmsd, linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(ggplot2::aes(fill = Type, size = Type, shape = Type), color = "black", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c("Highlight" = "#e74c3c", "Top Candidate" = "#3498db", "Other" = "grey90")) +
    ggplot2::scale_shape_manual(values = c("Highlight" = 23, "Top Candidate" = 21, "Other" = 21)) +
    ggplot2::scale_size_manual(values = c("Highlight" = 3.5, "Top Candidate" = 3, "Other" = 2)) +
    ggplot2::labs(x = "Crossover GxE (Stability RMSD)", y = "Overall Performance (OP)", title = "FAST Selection Space") +
    theme_genetics() +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    labels_df <- df[df$Type != "Other", ]
    p <- p + ggrepel::geom_label_repel(
      data = labels_df, ggplot2::aes(label = Genotype),
      size = 3, fontface = "bold", box.padding = 0.5
    )
  }

  return(p)
}

.plot_heat <- function(x) {
  cor_mat <- x$matrices$Cor
  if (is.null(cor_mat)) stop("No Correlation matrix available in fa_model object.")

  grp_name <- x$meta$group
  r_names <- rownames(cor_mat)
  c_names <- colnames(cor_mat)

  df_melt <- expand.grid(Var1 = r_names, Var2 = c_names, stringsAsFactors = TRUE)
  df_melt$Correlation <- as.vector(cor_mat)
  df_melt$Var1 <- factor(df_melt$Var1, levels = rev(r_names)) # Reverse to read top-to-bottom
  df_melt$Var2 <- factor(df_melt$Var2, levels = c_names)

  ggplot2::ggplot(df_melt, ggplot2::aes(x = Var2, y = Var1, fill = Correlation)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Correlation)), size = 2.5, color = "black") +
    ggplot2::scale_fill_distiller(palette = "RdYlGn", limit = c(-1, 1), direction = 1) +
    theme_genetics() +
    ggplot2::labs(x = NULL, y = NULL, title = paste("Genetic Correlation:", grp_name)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

.plot_vaf <- function(x) {
  if (is.null(x$var_comp$vaf)) stop("VAF data not available.")
  vaf_df <- x$var_comp$vaf

  vaf_df$Specific <- pmax(0, 100 - vaf_df$Total_VAF)
  vaf_df$Total_VAF <- NULL

  long_vaf <- reshape(vaf_df,
    varying = setdiff(names(vaf_df), "Group"),
    v.names = "Percentage", timevar = "Component",
    times = setdiff(names(vaf_df), "Group"), direction = "long"
  )

  long_vaf$Component <- gsub("VAF_", "", long_vaf$Component)
  group_order <- vaf_df$Group[order(vaf_df$VAF_Fac1, decreasing = TRUE)]
  long_vaf$Group <- factor(long_vaf$Group, levels = group_order)

  ggplot2::ggplot(long_vaf, ggplot2::aes(x = Group, y = Percentage, fill = Component)) +
    ggplot2::geom_col(color = "black", linewidth = 0.2) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = "Variance Accounted For (VAF) by Environment", x = "Environment", y = "Genetic Variance (%)") +
    theme_genetics() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

.plot_biplot <- function(x, n, h) {
  if (is.null(x$scores$rotated) || x$meta$k < 2) {
    stop("Biplot requires at least 2 factors and extracted scores.")
  }

  loadings <- as.data.frame(x$loadings$rotated[, 1:2])
  colnames(loadings) <- c("PC1", "PC2")
  loadings$Environment <- rownames(loadings)

  scores <- as.data.frame(x$scores$rotated[, 1:2])
  colnames(scores) <- c("PC1", "PC2")
  scores$Genotype <- rownames(scores)

  scores$Type <- "Other"
  if (!is.null(x$fast)) {
    top_gens <- head(x$fast$Genotype[order(x$fast$OP, decreasing = TRUE)], n)
    scores$Type[scores$Genotype %in% top_gens] <- "Top Candidate"
  }
  if (!is.null(h)) scores$Type[scores$Genotype %in% h] <- "Highlight"
  scores$Type <- factor(scores$Type, levels = c("Highlight", "Top Candidate", "Other"))

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
    ggplot2::geom_point(data = scores, ggplot2::aes(x = PC1, y = PC2, color = Type, size = Type, alpha = Type)) +
    ggplot2::scale_color_manual(values = c("Highlight" = "#e74c3c", "Top Candidate" = "#3498db", "Other" = "grey70")) +
    ggplot2::scale_size_manual(values = c("Highlight" = 3, "Top Candidate" = 2.5, "Other" = 1.5)) +
    ggplot2::scale_alpha_manual(values = c("Highlight" = 1, "Top Candidate" = 1, "Other" = 0.5)) +
    ggplot2::geom_segment(
      data = loadings, ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")), color = "darkblue", linewidth = 0.8
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = loadings, ggplot2::aes(x = PC1, y = PC2, label = Environment),
      color = "darkblue", size = 3, fontface = "bold"
    )
  } else {
    p <- p + ggplot2::geom_text(
      data = loadings, ggplot2::aes(x = PC1 * 1.1, y = PC2 * 1.1, label = Environment),
      color = "darkblue", size = 3
    )
  }

  p + ggplot2::labs(
    title = "Factor Analytic Biplot",
    x = paste0("Rotated Factor 1 (", round(x$meta$var_explained[1], 1), "%)"),
    y = paste0("Rotated Factor 2 (", round(x$meta$var_explained[2], 1), "%)")
  ) +
    theme_genetics() + ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
}

.plot_latent_reg <- function(x, n, h) {
  if (is.null(x$scores$rotated)) stop("Scores missing.")

  pred_mat <- x$scores$rotated %*% t(x$loadings$rotated)
  long_pred <- as.data.frame(as.table(pred_mat))
  colnames(long_pred) <- c("Genotype", "Environment", "Predicted_GEBV")

  env_loadings <- data.frame(Environment = rownames(x$loadings$rotated), Fac1 = x$loadings$rotated[, 1])
  long_pred <- merge(long_pred, env_loadings, by = "Environment")

  long_pred$Type <- "Other"
  if (!is.null(x$fast)) {
    top_gens <- head(x$fast$Genotype[order(x$fast$OP, decreasing = TRUE)], n)
    long_pred$Type[long_pred$Genotype %in% top_gens] <- "Top Candidate"
  }
  if (!is.null(h)) long_pred$Type[long_pred$Genotype %in% h] <- "Highlight"
  long_pred$Type <- factor(long_pred$Type, levels = c("Highlight", "Top Candidate", "Other"))

  ggplot2::ggplot(long_pred, ggplot2::aes(x = Fac1, y = Predicted_GEBV, group = Genotype)) +
    ggplot2::geom_line(ggplot2::aes(color = Type, alpha = Type, linewidth = Type)) +
    ggplot2::scale_color_manual(values = c("Highlight" = "#e74c3c", "Top Candidate" = "#3498db", "Other" = "grey85")) +
    ggplot2::scale_alpha_manual(values = c("Highlight" = 1, "Top Candidate" = 1, "Other" = 0.4)) +
    ggplot2::scale_linewidth_manual(values = c("Highlight" = 1.2, "Top Candidate" = 1, "Other" = 0.5)) +
    ggplot2::labs(
      title = "Latent Regression (Reaction Norms)",
      x = "Environmental Factor 1 Loading", y = "Predicted Interaction GEBV"
    ) +
    theme_genetics() +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
}

.plot_cluster <- function(x) {
  if (is.null(x$matrices$Cor)) stop("Correlation matrix missing.")
  dist_mat <- as.dist(1 - x$matrices$Cor)
  hc <- hclust(dist_mat, method = "ward.D2")

  # Base R hclust plot is very clean for this
  plot(hc, main = "Mega-Environment Clustering", xlab = "Environments", ylab = "Distance (1 - r)", sub = "", col = "darkblue")
}
