#' Generate Synthetic Genomic Relationship Matrix (GRM)
#'
#' @description
#' Creates a Genomic Relationship Matrix (GRM) using the VanRaden (2008) Method 1.
#' \deqn{G = \frac{ZZ'}{2 \sum p_i (1-p_i)}}
#' where \eqn{Z} is the centered genotype matrix ($M - 2P$).
#'
#' @param n_gen Integer. Number of genotypes to simulate.
#' @param n_markers Integer. Number of SNP markers.
#' @param seed Integer. Seed for reproducibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{grm}{The calculated GRM (n_gen x n_gen).}
#'   \item{markers}{The raw allele dosage matrix (0, 1, 2).}
#' }
#'
#' @references
#' VanRaden, P. M. (2008). Efficient methods to compute genomic predictions.
#' \emph{Journal of dairy science}, 91(11), 4414-4423.
#'
#' @export
generate_grm <- function(n_gen = 200, n_markers = 1000, seed = 123) {
  set.seed(seed)
  maf <- runif(n_markers, 0.1, 0.5)
  M <- matrix(0, nrow = n_gen, ncol = n_markers)
  for (j in 1:n_markers) M[, j] <- rbinom(n_gen, 2, maf[j])
  gen_names <- paste0("G", sprintf("%03d", 1:n_gen))
  marker_names <- paste0("M", sprintf("%03d", 1:n_markers))
  rownames(M) <- gen_names
  colnames(M) <- marker_names
  P_freq <- colMeans(M) / 2
  Two_P <- matrix(rep(2 * P_freq, n_gen), nrow = n_gen, byrow = TRUE)
  Z <- M - Two_P
  denom <- 2 * sum(P_freq * (1 - P_freq))
  GRM <- (Z %*% t(Z)) / denom
  rownames(GRM) <- colnames(GRM) <- gen_names
  diag(GRM) <- diag(GRM) + 1e-6
  return(list(grm = GRM, markers = M))
}

#' Generate Synthetic MET Data
#'
#' Simulates a Multi-Environment Trial (MET) dataset with realistic biological
#' and spatial properties. Useful for testing analysis pipelines and training.
#'
#' @param n_sites Integer. Number of environments.
#' @param n_gen Integer. Number of genotypes.
#' @param n_reps Integer. Number of replicates per site.
#' @param sparsity Numeric (0-1). Proportion of missing plots to simulate (unbalanced data).
#' @param seed Integer. Random seed for reproducibility.
#'
#' @details
#' The simulation includes:
#' \itemize{
#'   \item \strong{Genetic Main Effects:} Base yield potential (Factor 1).
#'   \item \strong{GxE Interaction:} Differential response to stress (Factor 2).
#'   \item \strong{Spatial Trends:} A fertility gradient (sine wave + linear) across Field Rows/Cols.
#'   \item \strong{Noise:} Random residual error.
#' }
#'
#' @return A dataframe containing \code{Site}, \code{Genotype}, \code{Rep},
#'         \code{Row}, \code{Column}, and \code{Yield}.
#' @export
generate_met_data <- function(n_sites = 10, n_gen = 100, n_reps = 2, sparsity = 0.2, seed = 123) {
  set.seed(seed)
  sites <- paste0("Site_", sprintf("%02d", 1:n_sites))
  gens <- paste0("G", sprintf("%03d", 1:n_gen))

  L1 <- runif(n_sites, 5, 15)
  S1 <- rnorm(n_gen, 0, 1)
  L2 <- runif(n_sites, -8, 8)
  S2 <- rnorm(n_gen, 0, 1)
  G_mat <- outer(S1, L1) + outer(S2, L2)
  rownames(G_mat) <- gens
  colnames(G_mat) <- sites

  n_plots <- n_gen * n_reps
  n_cols_field <- 20
  n_rows_field <- ceiling(n_plots / n_cols_field)

  data_list <- list()
  for (site in sites) {
    df_site <- expand.grid(Genotype = gens, Rep = 1:n_reps)
    df_site$Site <- site
    df_site$Row <- rep(1:n_rows_field, each = n_cols_field)[1:nrow(df_site)]
    df_site$Column <- rep(1:n_cols_field, times = n_rows_field)[1:nrow(df_site)]
    g_eff <- G_mat[df_site$Genotype, site]
    spatial_trend <- (sin(df_site$Row / 5) * 2) + (df_site$Column * 0.1)
    noise <- rnorm(nrow(df_site), 0, 2.5)
    df_site$Yield <- 50 + g_eff + spatial_trend + noise
    data_list[[site]] <- df_site
  }
  data <- do.call(rbind, data_list)
  if (sparsity > 0) data <- data[-sample(1:nrow(data), floor(nrow(data) * sparsity)), ]

  data$Site <- as.factor(data$Site)
  data$Genotype <- as.factor(data$Genotype)
  data$Row <- as.factor(data$Row)
  data$Column <- as.factor(data$Column)
  return(data[, c("Site", "Genotype", "Rep", "Row", "Column", "Yield")])
}

#' Generate Genomic MET Data
#'
#' @description
#' Simulates a Multi-Environment Trial (MET) where the genetic signal is strictly
#' controlled by a user-provided Genomic Relationship Matrix (GRM). This is ideal
#' for validating `vm()` (variance model) terms in ASReml-R.
#'
#' @param n_sites Integer. Number of environments to simulate.
#' @param grm Matrix. The relationship matrix (e.g. from \code{generate_grm()}).
#'        Genotypes in the GRM define the population.
#' @param n_reps Integer. Number of replicates per site.
#' @param h2 Numeric (0 to 1). The target Heritability (Broad Sense) for the trial.
#'        Lower h2 implies higher residual noise.
#'
#' @return A dataframe with columns:
#' \code{Site, Genotype, Rep, Row, Column, Yield}.
#'
#' @export
generate_genomic_met <- function(n_sites = 10, grm, n_reps = 2, h2 = 0.5) {
  n_gen <- nrow(grm)
  gens <- rownames(grm)
  sites <- paste0("Site_", sprintf("%02d", 1:n_sites))

  L_grm <- t(chol(grm))
  f1_base <- rnorm(n_gen)
  S1_g <- as.vector(L_grm %*% f1_base)
  f2_base <- rnorm(n_gen)
  S2_g <- as.vector(L_grm %*% f2_base)
  L1 <- runif(n_sites, 5, 15)
  L2 <- runif(n_sites, -8, 8)
  G_mat <- outer(S1_g, L1) + outer(S2_g, L2)
  rownames(G_mat) <- gens
  colnames(G_mat) <- sites

  n_plots <- n_gen * n_reps
  n_cols_field <- 20
  n_rows_field <- ceiling(n_plots / n_cols_field)

  data_list <- list()
  for (site in sites) {
    df_site <- expand.grid(Genotype = gens, Rep = 1:n_reps)
    df_site$Site <- site
    df_site$Row <- rep(1:n_rows_field, each = n_cols_field)[1:nrow(df_site)]
    df_site$Column <- rep(1:n_cols_field, times = n_rows_field)[1:nrow(df_site)]
    g_eff <- G_mat[df_site$Genotype, site]
    spatial_trend <- (sin(df_site$Row / 5) * 2) + (df_site$Column * 0.1)
    var_g <- var(g_eff)
    var_e <- var_g * (1 / h2 - 1)
    noise <- rnorm(nrow(df_site), 0, sqrt(var_e))
    df_site$Yield <- 50 + g_eff + spatial_trend + noise
    data_list[[site]] <- df_site
  }
  data <- do.call(rbind, data_list)
  data$Site <- as.factor(data$Site)
  data$Genotype <- as.factor(data$Genotype)
  data$Row <- as.factor(data$Row)
  data$Column <- as.factor(data$Column)
  return(data[, c("Site", "Genotype", "Rep", "Row", "Column", "Yield")])
}

#' Predict Genetic Gain and ROI
#'
#' Simulates the expected genetic gain from selection based on heritability and selection intensity,
#' considering financial and temporal costs.
#'
#' @param h2 Numeric. Narrow-sense heritability.
#' @param sigma_p Numeric. Phenotypic standard deviation.
#' @param selection_intensity Numeric. Selection intensity (i), e.g. 1.755 for 10\%.
#' @param cycle_time Numeric. Time in years to complete one breeding cycle. Default 1.
#' @param investment cost Numeric. Cost per cycle or unit investment. Default 100000.
#'
#' @return A list containing:
#' \describe{
#'   \item{gain}{Genetic Gain (R).}
#'   \item{gain_per_year}{Genetic Gain per Year.}
#'   \item{gain_per_dollar}{Genetic Gain per Dollar investe (x 1000).}
#' }
#'
#' @export
predict_gain <- function(h2, sigma_p, selection_intensity = 1.755, cycle_time = 1, investment_cost = 100000) {
  # Breeder's Equation: R = h^2 * S, where S = i * sigma_p
  # Alternatively: R = i * h * sigma_g  where h = sqrt(h2)
  # Standard: R = i * h2 * sigma_p

  R <- selection_intensity * h2 * sigma_p

  gain_per_year <- R / cycle_time
  gain_per_dollar <- (R / investment_cost) * 1000 # Gain per $1k

  return(list(
    gain = R,
    gain_per_year = gain_per_year,
    gain_per_dollar = gain_per_dollar,
    params = list(h2 = h2, sigma_p = sigma_p, i = selection_intensity, t = cycle_time, cost = investment_cost)
  ))
}
