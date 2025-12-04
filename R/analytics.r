# ==============================================================================
# ANALYTICS & DIAGNOSTICS ENGINE
# This file contains tools for Pre-Analysis Checks (Connectivity),
# Post-Analysis Optimization (D-Optimality, Interaction Classes),
# and Selection Decision Support (Indices).
# ==============================================================================

#' Check Trial Network Connectivity
#'
#' @description
#' A critical pre-analysis diagnostic tool. It calculates the number of common
#' genotypes between every pair of environments in a dataset.
#'
#' Factor Analytic models require genetic links (common lines) between trials to
#' estimate genetic correlations. If two sites share zero (or very few) varieties,
#' their correlation is undefined or unstable. This function flags such risks
#' **before** you attempt a computationally expensive ASReml run.
#'
#' @param data A dataframe containing the raw trial data.
#' @param genotype A character string. The column name for Genotype (default "Genotype").
#' @param trial A character string. The column name for Site/Environment (default "Site").
#' @param threshold Integer. The minimum number of shared lines required for a pair
#'        to be considered "connected" (default 10). Pairs below this are flagged.
#'
#' @return A list containing:
#' \describe{
#'   \item{matrix_count}{The Site x Site connectivity matrix (counts of shared lines).}
#'   \item{matrix_pct}{The connectivity matrix as percentages (Jaccard Index).}
#'   \item{disconnects}{A dataframe listing specific site pairs that fell below the threshold.}
#' }
#'
#' @details
#' The function also generates a heatmap plot of the connectivity matrix immediately
#' upon execution to provide visual verification of the trial network structure.
#'
#' @importFrom graphics image axis text par layout mtext
#' @importFrom grDevices hcl.colors
#' @export
check_connectivity <- function(data, genotype = "Genotype", trial = "Site", threshold = 10) {

  # 1. Incidence and Raw Counts (No changes)
  inc_table <- table(data[[genotype]], data[[trial]])
  incidence <- as.matrix(inc_table); incidence[incidence > 0] <- 1
  connect_mat <- t(incidence) %*% incidence # Intersection (A & B)

  # 2. Calculate Percentage (Jaccard Index: Intersection / Union) (No changes)
  site_totals <- diag(connect_mat)
  pct_mat <- connect_mat # Init

  p <- ncol(connect_mat)
  for(i in 1:p) {
    for(j in 1:p) {
      union_count <- site_totals[i] + site_totals[j] - connect_mat[i,j]
      pct_mat[i,j] <- round((connect_mat[i,j] / union_count) * 100, 1)
    }
  }

  # 3. Identify Issues (No changes)
  mat_check <- connect_mat; diag(mat_check) <- NA
  issues_idx <- which(mat_check < threshold, arr.ind = TRUE)

  disconnects <- NULL
  if (nrow(issues_idx) > 0) {
    disconnects <- data.frame(
      Site_A = rownames(connect_mat)[issues_idx[,1]],
      Site_B = colnames(connect_mat)[issues_idx[,2]],
      Shared_Count = connect_mat[issues_idx],
      Shared_Pct = pct_mat[issues_idx]
    )
    disconnects <- disconnects[as.character(disconnects$Site_A) < as.character(disconnects$Site_B), ]
    warning(sprintf("Found %d pairs with < %d shared lines.", nrow(disconnects), threshold))
  }

  # 4. Visualization (Heatmap with Scale Bar)

  # Use layout for legend
  layout(matrix(1:2, ncol=2), widths = c(4, 1))

  # Plot Matrix (No changes)
  par(mar = c(6, 6, 4, 1))
  mat_rev <- connect_mat[, p:1]
  cols <- hcl.colors(20, "YlGnBu", rev = TRUE) # Yellow (Low) to Blue (High)

  image(1:p, 1:p, mat_rev, axes = FALSE, col = cols,
        main = "Connectivity (Shared Count)")

  axis(1, at = 1:p, labels = colnames(connect_mat), las = 2, cex.axis = 0.7)
  axis(2, at = 1:p, labels = rev(rownames(connect_mat)), las = 1, cex.axis = 0.7)

  # Overlay Numbers (No changes)
  if(p < 20) {
    for(i in 1:p) {
      for(j in 1:p) {
        val <- mat_rev[i, j]
        # Text color logic
        txt_col <- ifelse(val < threshold, "red", "black")
        font_wt <- ifelse(val < threshold, 2, 1)
        text(i, j, labels = val, cex = 0.7, col = txt_col, font = font_wt)
      }
    }
  }

  # Plot Scale Bar (FIXED LOGIC)
  par(mar = c(6, 0, 4, 3))

  # Dummy strip 1-20
  legend_strip <- t(as.matrix(1:20))
  image(1, 1:20, legend_strip, axes = FALSE, xlab = "", ylab = "", col = cols)

  # Calculate Ticks based on Data Range
  min_v <- min(connect_mat)
  max_v <- max(connect_mat)

  # Generate "pretty" label values based on actual data
  pretty_vals <- pretty(c(min_v, max_v), n = 5)

  # === FIX IS HERE ===
  if (max_v > min_v) {
    # Map only the values that are within the min/max range of the data
    pretty_vals <- pretty_vals[pretty_vals >= min_v & pretty_vals <= max_v]

    # Calculate positions based on the actual range
    at_locs <- (pretty_vals - min_v) / (max_v - min_v) * 19 + 1
  } else {
    # If all values are the same, use only the one unique value, centered on the bar.
    # This ensures at_locs and pretty_vals have the same length (1).
    pretty_vals <- unique(c(min_v, max_v))[1]
    at_locs <- 10 # Center
  }
  # ===================

  # Draw axis
  axis(4, at = at_locs, labels = pretty_vals, las = 1, cex.axis = 0.8)
  mtext("Count", side=4, line=2, cex=0.7)

  # Reset layout
  layout(1)

  return(list(matrix_count = connect_mat, matrix_pct = pct_mat, disconnects = disconnects))
}


#' Calculate D-Optimality (Network Efficiency)
#'
#' @description
#' Quantifies the information content of the trial network based on the Rotated
#' Factor Loadings (\eqn{\Lambda}). High D-Optimality implies the sites effectively span
#' the factor space (good coverage of GxE drivers).
#'
#' It performs a "Leave-One-Out" analysis to calculate the "Impact" of each site.
#' Sites with low impact are redundant and candidates for removal to save costs.
#'
#' @param object An object of class \code{fa_asreml}.
#'
#' @return A list containing:
#' \describe{
#'   \item{total_d}{The determinant of the total information matrix (\eqn{D}).}
#'   \item{site_impact}{A dataframe ranking sites by their \% contribution to network information.}
#' }
#'
#' @details
#' The information matrix is defined as \eqn{M = \Lambda' \Lambda}. The D-criterion is \eqn{\det(M)}.
#' The impact of site \eqn{i} is calculated as:
#' \deqn{Impact_i = \frac{D_{total} - D_{-i}}{D_{total}} \times 100}
#'
#' @export
get_d_optimality <- function(object) {
  lam <- object$loadings$rotated
  M_total <- t(lam) %*% lam
  total_det <- det(M_total)

  impact_scores <- numeric(nrow(lam))
  for (i in 1:nrow(lam)) {
    lam_min <- lam[-i, , drop = FALSE]
    impact_scores[i] <- (total_det - det(t(lam_min) %*% lam_min)) / total_det * 100
  }
  site_impact <- data.frame(Site = rownames(lam), Impact_Pct = round(impact_scores, 2))
  site_impact <- site_impact[order(site_impact$Impact_Pct, decreasing = TRUE), ]
  return(list(total_d = total_det, site_impact = site_impact))
}

#' Calculate Interaction Classes (iClasses)
#'
#' @description
#' Implements the methodology of Smith et al. (2021) to detect Specific Adaptation.
#' It partitions environments into "Positive" and "Negative" classes based on the
#' sign and magnitude of their loadings for a specific Factor (usually Factor 2).
#'
#' @param object An object of class \code{fa_asreml}.
#' @param factor Integer. Which factor defines the interaction? (Default 2).
#' @param threshold Numeric. The loading magnitude required to assign a site to a class
#'        (e.g., 0.1). Sites between -0.1 and 0.1 are considered "Neutral".
#'
#' @return A list containing:
#' \describe{
#'   \item{site_classes}{Dataframe assigning each site to Class_Pos, Class_Neg, or Neutral.}
#'   \item{gen_effects}{Dataframe with predicted genetic values for each genotype in
#'   the Positive vs Negative environments, and the Differential (Crossover).}
#' }
#'
#' @details
#' Predictions are calculated as:
#' \deqn{Pred_{Class} = (Score_{F1} \times \bar{\lambda}_{F1}) + (Score_{Fk} \times \bar{\lambda}_{Class})}
#' This combines the variety's general performance (F1) with its specific response to the
#' interaction driver (Fk).
#'
#' @export
get_i_classes <- function(object, factor = 2, threshold = 0.1) {

  # SAFETY GATE
  if(is.null(object$scores)) stop("Genotype scores (BLUPs) not found. Cannot calculate iClasses.")
  if(factor > object$meta$k) stop("Factor not found.")

  lam <- object$loadings$rotated[, factor]
  sco <- object$scores$rotated[, factor]
  lam1 <- object$loadings$rotated[, 1]
  sco1 <- object$scores$rotated[, 1]

  site_class <- rep("Neutral", length(lam)); names(site_class) <- names(lam)
  site_class[lam > threshold] <- "Class_Pos"; site_class[lam < -threshold] <- "Class_Neg"

  ml_pos <- mean(lam[site_class == "Class_Pos"], na.rm = TRUE)
  ml_neg <- mean(lam[site_class == "Class_Neg"], na.rm = TRUE)
  m_perf <- mean(lam1, na.rm = TRUE)

  pred_pos <- if(is.nan(ml_pos)) rep(NA, length(sco)) else (sco1 * m_perf) + (sco * ml_pos)
  pred_neg <- if(is.nan(ml_neg)) rep(NA, length(sco)) else (sco1 * m_perf) + (sco * ml_neg)

  gen_eff <- data.frame(Genotype = names(sco), Pred_Pos = pred_pos, Pred_Neg = pred_neg,
                        Differential = pred_pos - pred_neg)
  gen_eff <- gen_eff[order(gen_eff$Differential, decreasing = TRUE, na.last = TRUE), ]

  return(list(site_classes = data.frame(Site = names(lam), Class = site_class),
              gen_effects = gen_eff, meta = list(factor = factor)))
}

#' Calculate Selection Index
#'
#' @description
#' Creates a single ranking metric by combining Overall Performance (OP) and Stability (RMSD).
#' Allows breeders to penalize instability based on their risk appetite.
#'
#' @param object An object of class \code{fa_asreml}.
#' @param weight Numeric. The economic weight applied to stability (default 1.0).
#'        \itemize{
#'          \item \code{weight = 0}: Select on Yield only (ignore stability).
#'          \item \code{weight > 1}: Heavy penalty for instability (Risk Averse).
#'        }
#'
#' @return A dataframe sorted by the calculated Index.
#' @details
#' \deqn{Index = OP - (weight \times RMSD)}
#' @export
calculate_index <- function(object, weight = 1.0) {
  if(is.null(object$fast)) stop("No FAST indices found in object.")

  df <- object$fast

  # Calculate Index
  # We subtract RMSD because high RMSD = Unstable (Bad)
  df$Index <- df$OP - (weight * df$RMSD)

  # Rank
  df$Rank <- rank(-df$Index)

  # Reorder
  df <- df[order(df$Index, decreasing = TRUE), ]

  # Return just the useful columns
  return(df[, c("Rank", "Genotype", "OP", "RMSD", "Index")])
}
