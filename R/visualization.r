#' Master Visualization Suite for Factor Analytic Models
#'
#' @description
#' The central plotting interface for the toolkit. This S3 method dispatches to
#' specific visualization engines based on the \code{type} argument, covering the
#' full analysis lifecycle: from Site Quality Control to Genotype Selection and
#' Mechanistic GxE analysis.
#'
#' It uses **Base R graphics** to ensure zero external dependencies and maximum
#' stability/speed, while producing publication-quality figures by default.
#'
#' @param x An object of class \code{fa_asreml} produced by \code{fa.asreml()}.
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
#'   \item For \code{"latent_reg"}: Can be a single integer (e.g., \code{2}) to see specific stability drivers, or \code{NULL} to plot all factors in a grid.
#'   \item For \code{"biplot"}: A vector of length 2 (e.g., \code{c(1, 2)}) specifying the X and Y axes.
#' }
#' @param n_label Integer. The number of top-performing genotypes (sorted by OP) to automatically label. Default is 5.
#' @param highlight Character vector. Names of specific genotypes to highlight (e.g., Check varieties). These will be rendered in **Red** with distinct symbols.
#' @param ... Additional graphical parameters passed to the internal base R plot functions (e.g., \code{main}, \code{cex}).
#'
#' @details
#' \strong{1. FAST Plot (\code{type = "fast"}):}
#' The primary tool for breeder decision-making.
#' \itemize{
#'   \item \strong{Y-axis (OP):} Genetic potential across the trial network. Higher is better.
#'   \item \strong{X-axis (RMSD):} Stability. Points near 0 are stable; points far right are highly reactive to GxE.
#'   \item \emph{Selection Strategy:} Select genotypes in the Top-Left quadrant (High OP, Low RMSD).
#' }
#'
#' \strong{2. Latent Regression (\code{type = "latent_reg"}):}
#' Visualizes *why* a variety is unstable.
#' \itemize{
#'   \item \strong{Factor 1:} Represents general potential. Slopes indicate responsiveness to high-yielding environments.
#'   \item \strong{Factor k > 1:} Represents specific GxE drivers (e.g., drought, heat).
#'   \item \strong{The "Peeling" Logic:} For k > 1, the function mathematically subtracts the effects of previous factors. The plot shows the \emph{deviation} from expected performance. A flat line means the variety is stable regarding that specific stressor.
#' }
#'
#' \strong{3. Biplot (\code{type = "biplot"}):}
#' A GGE-style overview.
#' \itemize{
#'   \item \strong{Vectors (Green):} Environments. The angle between vectors approximates their genetic correlation.
#'   \item \strong{Points:} Genotypes.
#'   \item \strong{Hull (Blue):} The convex hull connects the outermost genotypes. Genotypes on the hull are the "winners" in the sectors formed by perpendicular lines (not drawn to reduce clutter).
#' }
#'
#' \strong{4. Heatmap (\code{type = "heatmap"}):}
#' Displays the Genetic Correlation matrix ($C = D^{-1/2} G D^{-1/2}$).
#' \itemize{
#'   \item \strong{Red:} Negative correlation (Crossover GxE).
#'   \item \strong{Blue:} Positive correlation (Consistent ranking).
#'   \item \strong{White:} Zero correlation (Independent environments).
#'   \item Use this to identify "Mega-Environments" (blocks of blue).
#' }
#'
#' @return No return value, called for side effects (plotting).
#'
#' @seealso \code{\link{fa.asreml}}, \code{\link{get_d_optimality}}, \code{\link{get_i_classes}}, \code{\link{compare_h2}}
#'
#' @examples
#' \dontrun{
#' # Standard Selection View
#' plot(results, type = "fast", n_label = 5, highlight = c("CheckA", "CheckB"))
#'
#' # Investigate GxE Drivers (Factor 2)
#' plot(results, type = "latent_reg", factor = 2)
#'
#' # View Network Structure
#' plot(results, type = "biplot", factor = c(1, 2))
#' plot(results, type = "heatmap")
#' }
#' @export
plot.fa_asreml <- function(x, type = "fast", factor = NULL, n_label = 5, highlight = NULL, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

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
    .plot_dopt(get_d_optimality(x))
  } else if (type == "diff") {
    .plot_diff(get_i_classes(x, if (is.null(factor)) 2 else factor), n_label, highlight)
  } else if (type == "h2") {
    .plot_h2(x)
  } # Dispatch to h2 plotter
  else {
    stop("Unknown type.")
  }
}

# --- Plot Internals ---

.plot_fast <- function(x, n, h) {
  if (is.null(x$fast)) stop("No FAST data available (Genotype scores missing).")

  df <- x$fast
  top <- head(df$Genotype, n)
  col <- rep("grey60", nrow(df))
  bg <- rep("grey95", nrow(df))
  pch <- rep(21, nrow(df))
  is_h <- df$Genotype %in% h
  is_t <- df$Genotype %in% top
  bg[is_t] <- "#3498db"
  col[is_t] <- "#2980b9"
  bg[is_h] <- "#e74c3c"
  col[is_h] <- "#c0392b"
  pch[is_h] <- 23

  par(mar = c(5, 5, 4, 2))
  plot(df$RMSD, df$OP, pch = pch, bg = bg, col = col, xlab = "Stability (RMSD)", ylab = "Performance (OP)", main = "FAST Selection")
  grid()
  abline(h = 0, v = mean(df$RMSD, na.rm = T), lty = 2, col = "grey")
  lbl <- which(is_h | is_t)
  if (length(lbl) > 0) text(df$RMSD[lbl], df$OP[lbl], labels = df$Genotype[lbl], pos = 3, cex = 0.7)
}

.plot_heat <- function(x) {
  cor <- x$matrices$Cor
  p <- ncol(cor)
  rev <- cor[, p:1]
  par(mar = c(6, 6, 4, 2))
  image(1:p, 1:p, rev, axes = F, col = hcl.colors(20, "RdBu", rev = T), breaks = seq(-1, 1, l = 21), main = "Correlation")
  axis(1, at = 1:p, labels = colnames(cor), las = 2, cex.axis = 0.7)
  axis(2, at = 1:p, labels = rev(rownames(cor)), las = 1, cex.axis = 0.7)
  for (i in 1:p) for (j in 1:p) if (rev[i, j] < 0.99) text(i, j, sprintf("%.2f", rev[i, j]), cex = 0.6)
}

.plot_reg <- function(x, fac, h, n) {
  # SAFETY GATE
  if (is.null(x$scores)) stop("No Genotype scores available for Latent Regression.")

  k <- x$meta$k
  if (is.null(fac)) fac <- 1:k
  if (length(fac) > 1) par(mfrow = c(ceiling(length(fac) / 2), min(length(fac), 2)))
  L <- x$loadings$rotated
  S <- x$scores$rotated
  Uc <- S %*% t(L)
  tgt <- unique(c(head(rownames(S), n), h))

  for (f in fac) {
    if (f > k) next
    Y <- if (f > 1) Uc - (S[, 1:(f - 1)] %*% t(L[, 1:(f - 1)])) else Uc
    xv <- L[, f]
    ysub <- Y[tgt, , drop = F]

    par(mar = c(4, 4, 3, 1))
    plot(1, type = "n", xlim = range(xv), ylim = range(ysub), xlab = paste("Load", f), ylab = "Effect", main = paste("Factor", f))
    grid()
    abline(h = 0, lty = 2)
    axis(1, at = xv, labels = F, tck = -0.02, col = "green")

    for (g in tgt) {
      sl <- S[g, f]
      cc <- if (g %in% h) "red" else "navy"
      points(xv, Y[g, ], pch = 19, col = adjustcolor(cc, 0.3), cex = 0.7)
      abline(0, sl, col = cc)
      text(max(xv), sl * max(xv), labels = g, pos = 4, cex = 0.7, col = cc)
    }
  }
}

.plot_biplot <- function(x, fac, h) {
  # SAFETY GATE
  if (is.null(x$scores)) stop("No Genotype scores available for Biplot.")

  L <- x$loadings$rotated[, fac]
  S <- x$scores$rotated[, fac]
  sf <- max(abs(S)) / max(abs(L)) * 0.8
  Ls <- L * sf
  lim <- range(rbind(Ls, S))

  par(mar = c(5, 5, 4, 2))
  plot(1, type = "n", xlim = lim, ylim = lim, asp = 1, xlab = paste("F", fac[1]), ylab = paste("F", fac[2]), main = "Biplot")
  grid()
  abline(h = 0, v = 0, lty = 2, col = "grey")
  arrows(0, 0, Ls[, 1], Ls[, 2], col = "darkgreen", length = 0.1)
  text(Ls, labels = rownames(L), col = "darkgreen", pos = 4, cex = 0.7)

  bg <- rep("grey90", nrow(S))
  pch <- rep(21, nrow(S))
  is_h <- rownames(S) %in% h
  bg[is_h] <- "red"
  pch[is_h] <- 23
  points(S, pch = pch, bg = bg, cex = 0.8)
  if (any(is_h)) text(S[is_h, ], labels = rownames(S)[is_h], pos = 3, cex = 0.7)

  hpts <- chull(S)
  lines(S[c(hpts, hpts[1]), ], col = "navy", lty = 2)
}
.plot_vaf <- function(x) {
  df <- x$var_comp$vaf[order(x$var_comp$vaf$VAF), ]
  col <- ifelse(df$VAF < 50, "red", "blue")
  par(mar = c(5, 7, 4, 2))
  barplot(df$VAF, names.arg = df$Site, horiz = T, las = 1, col = col, xlab = "VAF %", main = "Site Quality")
  abline(v = 80, lty = 2)
}

.plot_dopt <- function(d, threshold = 1.0) {
  df <- d$site_impact
  # Sort ascending (so highest impact sites appear at the top of the chart)
  df <- df[order(df$Impact_Pct), ]

  # Color Logic:
  # Red (#e74c3c) for Redundant (below threshold)
  # Green (#27ae60) for High Value (above threshold)
  cols <- ifelse(df$Impact_Pct < threshold, "#e74c3c", "#27ae60")

  # Canvas Setup
  # Increase left margin (second value) to fit long site names
  par(mar = c(5, 8, 4, 2))

  # Draw Barplot
  bp <- barplot(df$Impact_Pct,
    names.arg = df$Site, horiz = TRUE,
    las = 1, col = cols, border = NA,
    xlab = "% Information Loss if Removed",
    main = "Network Efficiency (D-Opt Contribution)"
  )

  # Add Threshold Line
  abline(v = threshold, lty = 2, col = "grey40")

  # Add Value Labels at end of bars
  # xpd = TRUE ensures text isn't clipped if it goes outside the plot region
  text(
    x = df$Impact_Pct, y = bp, label = sprintf("%.1f%%", df$Impact_Pct),
    pos = 4, cex = 0.7, xpd = TRUE
  )

  # Add Legend
  legend("bottomright",
    legend = c("Key Sites (Keep)", "Redundant (Review)"),
    fill = c("#27ae60", "#e74c3c"),
    border = NA, bty = "n", cex = 0.8
  )
}

.plot_diff <- function(res, n, h) {
  df <- res$gen_effects
  df <- na.omit(df)
  tgt <- unique(c(head(df$Genotype, n), tail(df$Genotype, n), h))
  plot_df <- df[df$Genotype %in% tgt, ]

  yl <- range(c(plot_df$Pred_Neg, plot_df$Pred_Pos))
  par(mar = c(4, 4, 3, 6), xpd = F)
  plot(1, type = "n", xlim = c(0.8, 2.2), ylim = yl, xaxt = "n", ylab = "Genetic Value", xlab = "", main = "Interaction Classes")
  axis(1, at = 1:2, labels = c("Class (-)", "Class (+)"))

  for (i in 1:nrow(plot_df)) {
    g <- plot_df$Genotype[i]
    cc <- if (g %in% h) "red" else "grey50"
    segments(1, plot_df$Pred_Neg[i], 2, plot_df$Pred_Pos[i], col = cc, lwd = if (g %in% h) 2 else 1)
    text(1, plot_df$Pred_Neg[i], g, pos = 2, cex = 0.7, col = cc, xpd = T)
    text(2, plot_df$Pred_Pos[i], g, pos = 4, cex = 0.7, col = cc, xpd = T)
  }
}


#' Plot Spatial Field Map (Robust)
#'
#' Maps raw data (from dataframe) or residuals (from model) to the physical field layout.
#' Visualizes missing data with a distinct background color.
#'
#' @param input An `asreml` model object OR a `data.frame`.
#' @param row Column name for Field Row (default "Row").
#' @param col Column name for Field Column (default "Column").
#' @param attribute String. What to plot?
#'        - If input is data.frame: The column name (e.g. "Yield").
#'        - If input is asreml: "residuals" (default) or "fitted".
#' @export
plot_spatial <- function(input, row = "Row", col = "Column", attribute = "Yield") {
  if (inherits(input, "asreml")) {
    if (attribute == "Yield") attribute <- "residuals"
    data_name <- as.character(input$call$data)
    if (!exists(data_name)) stop("Original dataframe not found.")
    df <- get(data_name)

    if (attribute == "fitted") {
      vals <- fitted(input)
      main_title <- "Spatial Map: Fitted Values"
    } else {
      vals <- resid(input)
      main_title <- "Spatial Map: Residuals"
    }

    if (length(vals) == nrow(df)) {
      df$Value_To_Plot <- vals
    } else {
      warning("Residual mismatch. Padding with NA.")
      df$Value_To_Plot <- rep(NA, nrow(df))
      df$Value_To_Plot[1:length(vals)] <- vals
    }
  } else if (inherits(input, "data.frame")) {
    df <- input
    if (is.null(df[[attribute]])) stop(paste("Column", attribute, "not found."))
    vals <- df[[attribute]]
    df$Value_To_Plot <- vals
    main_title <- paste("Spatial Map:", attribute)
  } else {
    stop("Input must be an asreml model or a dataframe.")
  }

  r_vals <- as.numeric(as.character(df[[row]]))
  c_vals <- as.numeric(as.character(df[[col]]))
  n_r <- max(r_vals, na.rm = TRUE)
  n_c <- max(c_vals, na.rm = TRUE)

  field_mat <- matrix(NA, nrow = n_r, ncol = n_c)
  for (i in 1:nrow(df)) {
    if (!is.na(r_vals[i]) & !is.na(c_vals[i])) field_mat[r_vals[i], c_vals[i]] <- df$Value_To_Plot[i]
  }

  layout(matrix(1:2, ncol = 2), widths = c(4, 1))
  par(mar = c(3, 3, 3, 1))

  is_diverging <- (min(vals, na.rm = TRUE) < 0 && max(vals, na.rm = TRUE) > 0)
  cols <- if (is_diverging) hcl.colors(20, "Blue-Red 3") else hcl.colors(20, "Spectral", rev = TRUE)

  field_rev <- field_mat[n_r:1, ]
  image(1:n_c, 1:n_r, t(field_rev), axes = FALSE, col = NA, main = main_title, xlab = "", ylab = "")
  rect(0.5, 0.5, n_c + 0.5, n_r + 0.5, col = "grey90", border = NA)
  image(1:n_c, 1:n_r, t(field_rev), add = TRUE, col = cols)
  box()

  # Plot Scale Bar (FIXED LOGIC)
  par(mar = c(3, 0, 3, 3))

  # Dummy strip 1-20
  image(1, 1:20, t(as.matrix(1:20)), axes = FALSE, xlab = "", ylab = "", col = cols)

  # === FIX IS HERE ===
  data_range <- range(vals, na.rm = TRUE)
  min_v <- data_range[1]
  max_v <- data_range[2]

  # 1. Calculate the 'pretty' label values based on the data range
  pretty_vals <- pretty(data_range, n = 5)

  # Filter for values within the actual range (helpful for edge cases)
  pretty_vals <- pretty_vals[pretty_vals >= min_v & pretty_vals <= max_v]

  # 2. Map these values to the 1-20 coordinate system of the legend image
  # Formula: Position = (Value - Min) / Range * (LegendMax - 1) + 1
  if (max_v > min_v) {
    at_locs <- (pretty_vals - min_v) / (max_v - min_v) * 19 + 1
  } else {
    # Handle case where all values are the same
    at_locs <- 10 # Center
    pretty_vals <- unique(c(min_v, max_v))[1] # Use the single unique value
  }

  # 3. Draw axis using the matching locations and labels
  axis(4, at = at_locs, labels = round(pretty_vals, 2), las = 1)
  # ===================

  mtext("Grey = Missing", side = 1, line = 1, cex = 0.6)
  layout(1)
}
