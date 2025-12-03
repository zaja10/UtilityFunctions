
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
#'   \item{matrix}{The Site x Site connectivity matrix (counts of shared lines).}
#'   \item{disconnects}{A dataframe listing specific site pairs that fell below the threshold.}
#' }
#'
#' @details
#' The function also generates a heatmap plot of the connectivity matrix immediately
#' upon execution to provide visual verification of the trial network structure.
#'
#' @importFrom graphics image axis text par
#' @importFrom grDevices hcl.colors
#' @export
check_connectivity <- function(data, genotype = "Genotype", trial = "Site", threshold = 10) {
  cat("--- Checking Network Connectivity ---\n")
  if(!all(c(genotype, trial) %in% names(data))) stop("Columns not found.")

  inc_table <- table(data[[genotype]], data[[trial]])
  incidence <- as.matrix(inc_table); incidence[incidence > 0] <- 1
  connect_mat <- t(incidence) %*% incidence

  mat_check <- connect_mat; diag(mat_check) <- NA
  issues_idx <- which(mat_check < threshold, arr.ind = TRUE)

  disconnects <- NULL
  if (nrow(issues_idx) > 0) {
    disconnects <- data.frame(
      Site_A = rownames(connect_mat)[issues_idx[,1]],
      Site_B = colnames(connect_mat)[issues_idx[,2]],
      Shared = connect_mat[issues_idx]
    )
    disconnects <- disconnects[as.character(disconnects$Site_A) < as.character(disconnects$Site_B), ]
    cat(sprintf("WARNING: Found %d pairs with low connectivity (<%d lines).\n", nrow(disconnects), threshold))
  } else {
    cat(sprintf("PASSED: All site pairs share at least %d genotypes.\n", threshold))
  }

  p <- ncol(connect_mat)
  mat_rev <- connect_mat[, p:1]
  cols <- hcl.colors(20, "RdYlBu")
  par(mar = c(5, 6, 4, 2))
  image(1:p, 1:p, mat_rev, axes = FALSE, col = cols, main = "Connectivity Matrix")
  axis(1, at = 1:p, labels = colnames(connect_mat), las = 2, cex.axis = 0.7)
  axis(2, at = 1:p, labels = rev(rownames(connect_mat)), las = 1, cex.axis = 0.7)

  return(list(matrix = connect_mat, disconnects = disconnects))
}



#' Calculate D-Optimality (Network Efficiency)
#'
#' @description
#' Quantifies the information content of the trial network based on the Rotated
#' Factor Loadings ($\Lambda$). High D-Optimality implies the sites effectively span
#' the factor space (good coverage of GxE drivers).
#'
#' It performs a "Leave-One-Out" analysis to calculate the "Impact" of each site.
#' Sites with low impact are redundant and candidates for removal to save costs.
#'
#' @param object An object of class \code{fa_asreml}.
#'
#' @return A list containing:
#' \describe{
#'   \item{total_d}{The determinant of the total information matrix ($D$).}
#'   \item{site_impact}{A dataframe ranking sites by their \% contribution to network information.}
#' }
#'
#' @details
#' The information matrix is defined as $M = \Lambda' \Lambda$. The D-criterion is $\det(M)$.
#' The impact of site $i$ is calculated as:
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
