#' Default Checks
#'
#' Check varieties from the author's breeding program. Kept for convenience
#' (e.g. \code{plot_trend(..., checks = DEFAULT_CHECKS)}); no function uses
#' them by default.
#' @export
DEFAULT_CHECKS <- c("17-8930", "US16-IL-063-063", "25R76", "AgriMAXX 490", "07-19334", "AgriMAXX503")

#' Master Visualization Suite for Factor Analytic Models
#'
#' @param x Object of class \code{fa_model}.
#' @param type Plot type: \code{"fast"}, \code{"heatmap"}, \code{"vaf"},
#'   \code{"biplot"}, \code{"latent_reg"} or \code{"cluster"}.
#' @param factor Integer. Factor for the latent regression (default 1).
#' @param n_label Integer. Number of top genotypes to highlight (default 5).
#' @param highlight Character vector. Specific genotypes to highlight.
#' @param show_labels Logical. Label highlighted genotypes.
#' @param interactive Logical. Wrap the ggplot in \code{plotly::ggplotly()}.
#' @param ... Unused.
#' @return A ggplot (or plotly) object; for \code{"cluster"}, the
#'   \code{hclust} object, invisibly, after drawing the dendrogram.
#' @export
plot.fa_model <- function(x,
                          type = c("fast", "heatmap", "vaf", "biplot", "latent_reg", "cluster"),
                          factor = 1,
                          n_label = 5,
                          highlight = NULL,
                          show_labels = TRUE,
                          interactive = FALSE,
                          ...) {
  type <- match.arg(type)
  if (type == "cluster") {
    return(plot_cluster(x))
  }

  p <- switch(type,
    "fast"       = plot_fast(x, n_label, highlight, show_labels),
    "heatmap"    = plot_heat(x),
    "vaf"        = plot_vaf(x),
    "biplot"     = plot_biplot(x, n_label, highlight, show_labels),
    "latent_reg" = plot_latent_reg(x, factor, n_label, highlight)
  )

  if (interactive) {
    if (requireNamespace("plotly", quietly = TRUE)) {
      p <- plotly::ggplotly(p)
    } else {
      cli::cli_warn("Package {.pkg plotly} is required for interactive plots; returning a static plot.")
    }
  }
  p
}

#' Plot Trends (Phenotypic or Genetic)
#'
#' \itemize{
#'   \item \code{"phenotypic"}: boxplots per year with the regression line of
#'     the trait on (numeric) year, so unequally spaced years are handled.
#'   \item \code{"genetic"}: scatterplot of genetic values with regression and
#'     optional check highlighting.
#' }
#'
#' @param data Data frame.
#' @param mode \code{"phenotypic"} or \code{"genetic"}.
#' @param x X-axis column (numeric year).
#' @param y Y-axis column.
#' @param checks Character vector of checks to highlight (genetic mode).
#' @param genotype_col Genotype column used for check highlighting.
#' @return A ggplot object.
#' @export
plot_trend <- function(data, mode = c("phenotypic", "genetic"), x = "Year", y = "Yield",
                       checks = NULL, genotype_col = "Genotype") {
  mode <- match.arg(mode)
  data[[x]] <- suppressWarnings(as.numeric(as.character(data[[x]])))
  if (all(is.na(data[[x]]))) cli::cli_abort("Column {.val {x}} must be numeric (e.g. year).")

  if (mode == "phenotypic") {
    return(
      ggplot2::ggplot(data, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
        ggplot2::geom_boxplot(ggplot2::aes(group = .data[[x]]), fill = "lightgreen", alpha = 0.6, outlier.alpha = 0.5) +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE) +
        ggplot2::labs(title = paste("Phenotypic Trend:", y, "over", x), x = x, y = y) +
        theme_genetics()
    )
  }

  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
    ggplot2::geom_point(color = "grey", alpha = 0.5) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = TRUE, alpha = 0.2) +
    ggplot2::labs(title = paste("Genetic Trend:", y, "over", x), x = x, y = "Genetic Value") +
    theme_genetics()

  if (!is.null(checks)) {
    if (!genotype_col %in% names(data)) cli::cli_abort("Column {.val {genotype_col}} not found.")
    check_data <- data[data[[genotype_col]] %in% checks, , drop = FALSE]
    if (nrow(check_data) > 0) {
      p <- p +
        ggplot2::geom_point(data = check_data, color = "red", shape = 17, size = 3) +
        ggrepel::geom_text_repel(
          data = check_data, ggplot2::aes(label = .data[[genotype_col]]),
          color = "darkred", box.padding = 0.5
        )
    }
  }
  p
}

#' Plot Spatial Map
#'
#' Draws a field trial as a heatmap of raw values or of model residuals.
#'
#' @param input A data frame, or an \code{asreml} model.
#' @param row Row column name.
#' @param col Column column name.
#' @param fill Column to plot (data frame input).
#' @param data For \code{asreml} input: the data frame the model was fitted
#'   to, in the same row order. If \code{NULL}, the data named in the model
#'   call is looked up.
#' @return A ggplot object.
#' @export
plot_spatial <- function(input, row = "Row", col = "Column", fill = "Yield", data = NULL) {
  if (inherits(input, "asreml")) {
    df <- data %||% tryCatch(eval(input$call$data, parent.frame()), error = function(e) NULL)
    if (is.null(df)) cli::cli_abort("Supply {.arg data}: the data frame the model was fitted to.")
    resids <- stats::resid(input)
    if (length(resids) != nrow(df)) {
      cli::cli_abort("{length(resids)} residuals but {nrow(df)} data rows; fit with {.code na.action = na.method(y = 'include')}.")
    }
    df$PlotValue <- as.numeric(resids)
    title <- "Spatial Residuals"
  } else {
    df <- input
    if (!fill %in% names(df)) cli::cli_abort("Column {.val {fill}} not found.")
    df$PlotValue <- df[[fill]]
    title <- paste("Spatial Map:", fill)
  }

  if (!all(c(row, col) %in% names(df))) {
    cli::cli_abort("Spatial columns {.val {c(row, col)}} not found.")
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[col]], y = .data[[row]], fill = .data$PlotValue)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.1) +
    ggplot2::scale_fill_distiller(palette = "Spectral", na.value = "grey90") +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = col, y = row, fill = "Value") +
    theme_genetics()
}

# -----------------------------------------------------------------------------
# Internal implementations for plot.fa_model
# -----------------------------------------------------------------------------

type_levels <- c("Highlight", "Top Candidate", "Other")
type_colours <- c("Highlight" = "#e74c3c", "Top Candidate" = "#3498db", "Other" = "grey80")

label_types <- function(ids, fast, n, highlight) {
  type <- rep("Other", length(ids))
  if (!is.null(fast)) type[ids %in% utils::head(fast$Genotype[order(fast$OP, decreasing = TRUE)], n)] <- "Top Candidate"
  if (!is.null(highlight)) type[ids %in% highlight] <- "Highlight"
  factor(type, levels = type_levels)
}

plot_fast <- function(x, n, h, show_labels) {
  if (is.null(x$fast)) cli::cli_abort("No FAST indices available (scores missing or rotation not SVD).")
  df <- x$fast
  df$Type <- label_types(df$Genotype, df, n, h)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$RMSD, y = .data$OP)) +
    ggplot2::geom_vline(xintercept = mean(df$RMSD, na.rm = TRUE), linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(ggplot2::aes(fill = .data$Type, size = .data$Type, shape = .data$Type), color = "black", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = type_colours, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c("Highlight" = 23, "Top Candidate" = 21, "Other" = 21), drop = FALSE) +
    ggplot2::scale_size_manual(values = c("Highlight" = 3.5, "Top Candidate" = 3, "Other" = 2), drop = FALSE) +
    ggplot2::labs(x = "Crossover GxE (Stability RMSD)", y = "Overall Performance (OP)", title = "FAST Selection Space") +
    theme_genetics() +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  lab <- df[df$Type != "Other", , drop = FALSE]
  if (show_labels && nrow(lab) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = lab, ggplot2::aes(label = .data$Genotype),
      size = 3, fontface = "bold", box.padding = 0.5, show.legend = FALSE
    )
  }
  p
}

plot_heat <- function(x) {
  cor_mat <- x$matrices$Cor
  if (is.null(cor_mat)) cli::cli_abort("No correlation matrix available in the {.cls fa_model}.")
  r_names <- rownames(cor_mat)
  df <- expand.grid(Var1 = r_names, Var2 = colnames(cor_mat), stringsAsFactors = FALSE)
  df$Correlation <- as.vector(cor_mat)
  df$Var1 <- factor(df$Var1, levels = rev(r_names))
  df$Var2 <- factor(df$Var2, levels = colnames(cor_mat))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$Var2, y = .data$Var1, fill = .data$Correlation)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$Correlation)), size = 2.5) +
    ggplot2::scale_fill_distiller(palette = "RdYlGn", limits = c(-1, 1), direction = 1) +
    ggplot2::labs(x = NULL, y = NULL, title = paste("Genetic Correlation:", x$meta$group)) +
    theme_genetics() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

plot_vaf <- function(x) {
  vaf <- x$var_comp$vaf
  if (is.null(vaf) || !all(c("Group", "Total_VAF") %in% names(vaf))) {
    cli::cli_abort("VAF table with {.field Group} and {.field Total_VAF} columns is required.")
  }
  fac_cols <- grep("^VAF_", names(vaf), value = TRUE)
  long <- do.call(rbind, c(
    lapply(fac_cols, function(cn) data.frame(Group = vaf$Group, Component = sub("^VAF_", "", cn), Percentage = vaf[[cn]])),
    list(data.frame(Group = vaf$Group, Component = "Specific", Percentage = pmax(0, 100 - vaf$Total_VAF)))
  ))
  long$Group <- factor(long$Group, levels = vaf$Group[order(vaf[[fac_cols[1]]], decreasing = TRUE)])
  long$Component <- factor(long$Component, levels = c(sub("^VAF_", "", fac_cols), "Specific"))

  ggplot2::ggplot(long, ggplot2::aes(x = .data$Group, y = .data$Percentage, fill = .data$Component)) +
    ggplot2::geom_col(color = "black", linewidth = 0.2) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = "Variance Accounted For (VAF) by Environment", x = "Environment", y = "Genetic Variance (%)") +
    theme_genetics() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

plot_biplot <- function(x, n, h, show_labels) {
  if (is.null(x$scores$rotated) || x$meta$k < 2) cli::cli_abort("A biplot requires at least 2 factors and extracted scores.")

  scores <- data.frame(PC1 = x$scores$rotated[, 1], PC2 = x$scores$rotated[, 2], Genotype = rownames(x$scores$rotated))
  scores$Type <- label_types(scores$Genotype, x$fast, n, h)

  L <- x$loadings$rotated[, 1:2, drop = FALSE]
  # Loadings and scores live on different scales: stretch loadings to the
  # score cloud and say so on the plot.
  s <- 0.8 * max(abs(scores[, 1:2])) / max(abs(L))
  loadings <- data.frame(PC1 = L[, 1] * s, PC2 = L[, 2] * s, Environment = rownames(L))

  ve <- x$meta$var_explained
  ax <- function(i) if (length(ve) >= i && is.finite(ve[i])) sprintf("Rotated factor %d (%.1f%% of genetic variance)", i, ve[i]) else sprintf("Rotated factor %d", i)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
    ggplot2::geom_point(data = scores, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data$Type, size = .data$Type, alpha = .data$Type)) +
    ggplot2::scale_color_manual(values = type_colours, drop = FALSE) +
    ggplot2::scale_size_manual(values = c("Highlight" = 3, "Top Candidate" = 2.5, "Other" = 1.5), drop = FALSE) +
    ggplot2::scale_alpha_manual(values = c("Highlight" = 1, "Top Candidate" = 1, "Other" = 0.5), drop = FALSE) +
    ggplot2::geom_segment(
      data = loadings, ggplot2::aes(x = 0, y = 0, xend = .data$PC1, yend = .data$PC2),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")), color = "darkblue", linewidth = 0.8
    ) +
    ggrepel::geom_text_repel(
      data = loadings, ggplot2::aes(x = .data$PC1, y = .data$PC2, label = .data$Environment),
      color = "darkblue", size = 3, fontface = "bold"
    )

  lab <- scores[scores$Type != "Other", , drop = FALSE]
  if (show_labels && nrow(lab) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = lab, ggplot2::aes(x = .data$PC1, y = .data$PC2, label = .data$Genotype),
      size = 3, fontface = "bold", box.padding = 0.5, show.legend = FALSE
    )
  }

  p + ggplot2::labs(
    title = "Factor Analytic Biplot", x = ax(1), y = ax(2),
    caption = sprintf("Environment loadings scaled by %.3g to match genotype scores", s)
  ) +
    theme_genetics() + ggplot2::theme(legend.title = ggplot2::element_blank())
}

plot_latent_reg <- function(x, fac_num, n, h) {
  if (is.null(x$scores$rotated)) cli::cli_abort("Scores are required.")
  if (fac_num < 1 || fac_num > x$meta$k) {
    cli::cli_abort("Requested factor {fac_num} exceeds available factors (k = {x$meta$k}).")
  }

  # Common (factor-driven) genetic effect of each genotype in each environment.
  pred <- x$scores$rotated %*% t(x$loadings$rotated)
  long <- data.frame(
    Genotype = rep(rownames(pred), times = ncol(pred)),
    Environment = rep(colnames(pred), each = nrow(pred)),
    Predicted = as.vector(pred),
    stringsAsFactors = FALSE
  )
  long$Loading <- x$loadings$rotated[long$Environment, fac_num]
  long$Type <- label_types(long$Genotype, x$fast, n, h)
  long <- long[order(long$Type, decreasing = TRUE), , drop = FALSE]

  ve <- x$meta$var_explained[fac_num]
  ggplot2::ggplot(long, ggplot2::aes(x = .data$Loading, y = .data$Predicted, group = .data$Genotype)) +
    ggplot2::geom_line(ggplot2::aes(color = .data$Type, alpha = .data$Type, linewidth = .data$Type)) +
    ggplot2::scale_color_manual(values = type_colours, drop = FALSE) +
    ggplot2::scale_alpha_manual(values = c("Highlight" = 1, "Top Candidate" = 1, "Other" = 0.4), drop = FALSE) +
    ggplot2::scale_linewidth_manual(values = c("Highlight" = 1.2, "Top Candidate" = 1, "Other" = 0.5), drop = FALSE) +
    ggplot2::labs(
      title = paste("Latent Regression (Reaction Norms) - Factor", fac_num),
      x = sprintf("Factor %d loading (%.1f%% of genetic variance)", fac_num, ve),
      y = "Common genetic effect (FA component)"
    ) +
    theme_genetics() +
    ggplot2::theme(legend.title = ggplot2::element_blank())
}

plot_cluster <- function(x) {
  C <- x$matrices$Cor
  if (is.null(C)) cli::cli_abort("Correlation matrix missing.")
  # sqrt(2(1 - r)) is a Euclidean distance (1 - r is not), as Ward's method requires.
  D <- 2 * (1 - C)
  D[D < 0] <- 0
  hc <- stats::hclust(stats::as.dist(sqrt(D)), method = "ward.D2")
  graphics::plot(hc, main = "Mega-Environment Clustering", xlab = "Environments", ylab = "Distance sqrt(2(1 - r))", sub = "", col = "darkblue")
  invisible(hc)
}
