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
    .plot_biplot_static(x, highlight, if (is.null(factor)) c(1, 2) else factor)
  } else if (type == "vaf") {
    .plot_vaf(x)
  } else if (type == "d_opt") {
    .plot_dopt(calculate_d_optimality(x))
  } else if (type == "diff") {
    .plot_diff_generalized(x, ...)
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

.plot_biplot_static <- function(x, highlight = NULL, fac = c(1, 2)) {
  if (is.null(x$scores$rotated)) stop("No scores.")

  L <- x$loadings$rotated[, fac]
  S <- x$scores$rotated[, fac]
  var_pct <- x$meta$var_explained[fac]

  # Smart Scaling
  max_s <- max(abs(S))
  max_l <- max(abs(L))
  scale_f <- (max_s / max_l) * 0.85

  df_scr <- data.frame(X = S[, 1], Y = S[, 2], ID = rownames(S), Type = "Genotype")
  df_lod <- data.frame(X = L[, 1] * scale_f, Y = L[, 2] * scale_f, ID = rownames(L), Type = "Site")

  # --- ANNOTATION LAYER START ---
  # Check if annotations exist in metadata
  has_anno <- !is.null(x$meta$annotation)

  if (has_anno) {
    # Merge tags: ID in df_lod matches Group in annotation
    df_lod <- merge(df_lod, x$meta$annotation, by.x = "ID", by.y = "Group", all.x = TRUE)
    # Ensure standard fallback if merge missed anything
    df_lod$Tag[is.na(df_lod$Tag)] <- "Standard"

    # Define aesthetic mapping for Sites
    site_mapping_arrow <- aes(x = 0, y = 0, xend = X, yend = Y, color = Tag, linetype = Tag)
    site_mapping_text <- aes(X, Y, label = ID, color = Tag)

    # Custom palette for sites (Distinct from Genotypes)
    # We use a viridis or brewer palette for site tags to ensure professionalism
    # But we need to handle "Standard" nicely (usually dark blue/black)
  } else {
    df_lod$Tag <- "Standard"
    site_mapping_arrow <- aes(x = 0, y = 0, xend = X, yend = Y) # Inherits fixed color
    site_mapping_text <- aes(X, Y, label = ID)
  }
  # --- ANNOTATION LAYER END ---

  # Highlight Logic for Genotypes
  df_scr$Alpha <- 0.4
  df_scr$Color <- "grey60" # Default Genotype Color

  if (!is.null(highlight)) {
    hits <- df_scr$ID %in% highlight
    df_scr$Alpha[hits] <- 1
    df_scr$Color[hits] <- "red" # Genotype Highlight Color
    df_scr <- rbind(df_scr[!hits, ], df_scr[hits, ])
  }

  p <- ggplot() +
    # 1. Genotypes (Scores) - Layered first (background)
    geom_point(data = df_scr, aes(X, Y), color = df_scr$Color, alpha = df_scr$Alpha) +

    # 2. Environments (Vectors)
    # If annotated, use dynamic color; if not, use static dark blue
    {
      if (has_anno) {
        geom_segment(
          data = df_lod, site_mapping_arrow,
          arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1
        )
      } else {
        geom_segment(
          data = df_lod, site_mapping_arrow,
          arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
          color = "#2c3e50", linewidth = 1
        )
      }
    } +

    # 3. Labels (Sites)
    {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        if (has_anno) {
          ggrepel::geom_text_repel(
            data = df_lod, site_mapping_text,
            fontface = "bold", size = 3.5, bg.color = "white", bg.r = 0.15,
            max.overlaps = 15, seed = 42
          )
        } else {
          ggrepel::geom_text_repel(
            data = df_lod, site_mapping_text,
            color = "#2c3e50", fontface = "bold", size = 3.5,
            bg.color = "white", bg.r = 0.15,
            max.overlaps = 15, seed = 42
          )
        }
      } else {
        geom_text(data = df_lod, aes(X, Y, label = ID), color = "navy")
      }
    } +

    # 4. Labels (Highlighted Genotypes)
    {
      if (!is.null(highlight)) {
        df_lbl_gen <- df_scr[df_scr$ID %in% highlight, ]
        if (nrow(df_lbl_gen) > 0) {
          if (requireNamespace("ggrepel", quietly = TRUE)) {
            ggrepel::geom_text_repel(
              data = df_lbl_gen, aes(X, Y, label = ID),
              color = "red", fontface = "bold.italic", size = 3,
              bg.color = "white", bg.r = 0.1,
              max.overlaps = 20, seed = 42
            )
          } else {
            geom_text(data = df_lbl_gen, aes(X, Y, label = ID), color = "red", vjust = -1)
          }
        }
      }
    } +

    # 4. Aesthetics
    labs(
      x = paste0("Factor ", fac[1], " (", round(var_pct[1], 1), "%)"),
      y = paste0("Factor ", fac[2], " (", round(var_pct[2], 1), "%)"),
      title = "FA Biplot",
      subtitle = if (has_anno) "Colored by Environment Condition" else NULL,
      color = "Condition", linetype = "Condition"
    ) + # Legend Titles

    coord_fixed() +
    theme_genetics() + # Use our new theme

    # 5. Scale Adjustments if Annotated
    {
      if (has_anno) scale_color_brewer(palette = "Dark2")
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

.plot_diff_generalized <- function(x, ...) {
  # Calculate classes
  classes <- calculate_i_classes(x)
  df <- classes$geno_classes

  # Plotting logic for k classes
  # Simply bar plot of counts or similar for now
  ggplot(df, aes(x = Class)) +
    geom_bar(fill = "steelblue") +
    labs(title = "Genotype Interaction Classes", x = "Sign Pattern (Factors 1..k)") +
    theme_minimal()
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

#' Plot Selection Profile (Parallel Coordinates) with Check Comparison
#'
#' Visualizes the standardized performance of top genotypes and a specific Check
#' variety across all traits relative to the population.
#'
#' @param index_df Dataframe returned by \code{calculate_mt_index}.
#' @param top_n Number of top genotypes to highlight.
#' @param check_geno Character string. Name of the Check Genotype to highlight (optional).
#' @export
plot_index_profile <- function(index_df, top_n = 5, check_geno = NULL) {
  # 1. Prepare Data
  z_cols <- grep("_OP_Z$", names(index_df), value = TRUE)
  if (length(z_cols) == 0) stop("No '_OP_Z' columns found. Did you run calculate_mt_index?")

  # Filter non-culled for background context
  if ("Status" %in% names(index_df)) {
    valid_df <- index_df[index_df$Status != "Culled", ]
  } else {
    valid_df <- index_df
  }

  # Identify Groups
  top_genos <- head(valid_df$Genotype, top_n)

  # Manual pivot to long format (base R robustness)
  df_long_list <- lapply(z_cols, function(col) {
    data.frame(
      Genotype = valid_df$Genotype,
      Trait = sub("_OP_Z", "", col),
      Z_Score = valid_df[[col]],
      stringsAsFactors = FALSE
    )
  })
  df_long <- do.call(rbind, df_long_list)

  # Assign Groups
  df_long$Group <- "Pop"
  df_long$Group[df_long$Genotype %in% top_genos] <- "Top"

  # Handle Check Genotype
  if (!is.null(check_geno)) {
    if (!check_geno %in% valid_df$Genotype) {
      warning(paste("Check genotype", check_geno, "not found in data."))
    } else {
      df_long$Group[df_long$Genotype == check_geno] <- "Check"
    }
  }

  # Order factors: Pop on bottom, Check middle, Top on top
  df_long$Group <- factor(df_long$Group, levels = c("Pop", "Check", "Top"))

  # 2. Plotting
  p <- ggplot(df_long, aes(x = Trait, y = Z_Score, group = Genotype)) +
    # A. Background Population (Faded Grey)
    geom_line(
      data = df_long[df_long$Group == "Pop", ],
      color = "grey90", linesize = 0.5, alpha = 0.6
    ) +

    # B. Check Variety (Distinct Black Dashed)
    geom_line(
      data = df_long[df_long$Group == "Check", ],
      aes(linetype = "Check"), color = "black", linesize = 1, alpha = 0.8
    ) +

    # C. Top Selections (Colored)
    geom_line(
      data = df_long[df_long$Group == "Top", ],
      aes(color = Genotype), linesize = 1.2
    ) +

    # Reference Mean Line
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +

    # Aesthetics
    scale_linetype_manual(name = "Benchmark", values = c("Check" = "longdash")) +
    theme_genetics() +
    labs(
      title = "Selection Profile: Top Candidates vs Check",
      subtitle = "Standardized Performance (0 = Population Mean)",
      y = "Standard Deviations (SD)",
      x = NULL
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_line(color = "grey90")
    )

  return(p)
}

#' Plot Predicted Response to Selection
#'
#' Visualizes the percentage change in traits if the top selection is advanced.
#'
#' @param resp_df Dataframe returned by \code{predict_response}.
#' @export
plot_response <- function(resp_df) {
  # Classify gain as Positive/Negative for color
  resp_df$Direction <- ifelse(resp_df$Pct_Change >= 0, "Positive", "Negative")

  ggplot(resp_df, aes(x = reorder(Trait, Pct_Change), y = Pct_Change, fill = Direction)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c("Positive" = "#2ecc71", "Negative" = "#e74c3c")) +
    geom_hline(yintercept = 0, color = "black") +
    theme_genetics() +
    labs(
      title = "Predicted Genetic Gain",
      subtitle = "Expected shift in population mean (Top Selection)",
      y = "% Change relative to Population Mean",
      x = NULL
    ) +
    theme(legend.position = "none")
}

#' Plot Selection Pipeline Results
#'
#' Generates visualizations for the Multi-Trait Selection Pipeline.
#'
#' @param x An object of class \code{mt_selection_results}.
#' @param type Character. "profile" for the Selection Profile (Parallel Coordinates),
#'             "gain" for the Response to Selection bar chart. Default is "profile".
#' @param ... Additional arguments.
#' @export
plot.mt_selection_results <- function(x, type = "profile", ...) {
  if (type == "profile") {
    plot_index_profile(x$index, top_n = x$params$top_n, check_geno = x$params$check_geno)
  } else if (type == "gain") {
    if (is.null(x$prediction)) stop("No prediction data available.")
    plot_response(x$prediction)
  } else {
    stop("Unknown plot type. Use 'profile' or 'gain'.")
  }
}
