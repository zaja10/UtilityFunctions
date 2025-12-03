# ==============================================================================
# TITLE:        ASReml-R Factor Analytic Toolkit (fa_asreml)
# VERSION:      1.1.0
# AUTHOR:       Zachary Aldiss
# DESCRIPTION:  A robust utility suite for processing Factor Analytic (FA) models
#               fitted in ASReml-R. Implements Smith & Cullis (2018) FAST,
#               D-Optimality, Interaction Classes, and Biplots.
# DEPENDENCIES: asreml, tidyverse
# ==============================================================================

# Check dependencies
if (!requireNamespace("tidyverse", quietly = TRUE)) stop("Package 'tidyverse' required.")

library(tidyverse)

# ==============================================================================
# SECTION 1: SIMULATION TOOLS (The Sandbox)
# ==============================================================================

#' Generate Synthetic Genomic Relationship Matrix (GRM)
#' @export
generate_grm <- function(n_gen = 200, n_markers = 1000, seed = 123) {
  set.seed(seed)
  maf <- runif(n_markers, 0.1, 0.5)
  M <- matrix(0, nrow = n_gen, ncol = n_markers)
  for (j in 1:n_markers) M[, j] <- rbinom(n_gen, 2, maf[j])

  gen_names <- paste0("G", sprintf("%03d", 1:n_gen))
  rownames(M) <- colnames(M) <- gen_names # Placeholder for GRM names
  rownames(M) <- gen_names

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
  gens  <- paste0("G", sprintf("%03d", 1:n_gen))

  L1 <- runif(n_sites, 5, 15); S1 <- rnorm(n_gen, 0, 1)
  L2 <- runif(n_sites, -8, 8); S2 <- rnorm(n_gen, 0, 1)
  G_mat <- outer(S1, L1) + outer(S2, L2)
  rownames(G_mat) <- gens; colnames(G_mat) <- sites

  n_plots <- n_gen * n_reps
  n_cols_field <- 20
  n_rows_field <- ceiling(n_plots / n_cols_field)

  data_list <- list()
  for(site in sites) {
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

  data$Site <- as.factor(data$Site); data$Genotype <- as.factor(data$Genotype)
  data$Row <- as.factor(data$Row); data$Column <- as.factor(data$Column)
  return(data[, c("Site", "Genotype", "Rep", "Row", "Column", "Yield")])
}

#' Generate Genomic MET Data
#'
#' Simulates a MET where variety performance is determined by a Genomic Relationship
#' Matrix (GRM). Designed for testing \code{vm()} terms in ASReml.
#'
#' @param n_sites Integer. Number of environments.
#' @param grm Matrix. A relationship matrix (e.g., from \code{generate_grm}).
#' @param n_reps Integer. Replicates per site.
#' @param h2 Numeric. Heritability (0-1). Controls the signal-to-noise ratio.
#'
#' @return A dataframe containing spatial coordinates and yield, suitable for
#' genomic GxE analysis.
#' @export
generate_genomic_met <- function(n_sites = 10, grm, n_reps = 2, h2 = 0.5) {
  n_gen <- nrow(grm); gens <- rownames(grm)
  sites <- paste0("Site_", sprintf("%02d", 1:n_sites))

  L_grm <- t(chol(grm))
  f1_base <- rnorm(n_gen); S1_g <- as.vector(L_grm %*% f1_base)
  f2_base <- rnorm(n_gen); S2_g <- as.vector(L_grm %*% f2_base)
  L1 <- runif(n_sites, 5, 15); L2 <- runif(n_sites, -8, 8)
  G_mat <- outer(S1_g, L1) + outer(S2_g, L2)
  rownames(G_mat) <- gens; colnames(G_mat) <- sites

  n_plots <- n_gen * n_reps
  n_cols_field <- 20
  n_rows_field <- ceiling(n_plots / n_cols_field)

  data_list <- list()
  for(site in sites) {
    df_site <- expand.grid(Genotype = gens, Rep = 1:n_reps)
    df_site$Site <- site
    df_site$Row <- rep(1:n_rows_field, each = n_cols_field)[1:nrow(df_site)]
    df_site$Column <- rep(1:n_cols_field, times = n_rows_field)[1:nrow(df_site)]
    g_eff <- G_mat[df_site$Genotype, site]
    spatial_trend <- (sin(df_site$Row / 5) * 2) + (df_site$Column * 0.1)
    var_g <- var(g_eff); var_e <- var_g * (1/h2 - 1)
    noise <- rnorm(nrow(df_site), 0, sqrt(var_e))
    df_site$Yield <- 50 + g_eff + spatial_trend + noise
    data_list[[site]] <- df_site
  }
  data <- do.call(rbind, data_list)
  data$Site <- as.factor(data$Site); data$Genotype <- as.factor(data$Genotype)
  data$Row <- as.factor(data$Row); data$Column <- as.factor(data$Column)
  return(data[, c("Site", "Genotype", "Rep", "Row", "Column", "Yield")])
}

# ==============================================================================
# SECTION 2: PRE-ANALYSIS DIAGNOSTICS
# ==============================================================================

#' Check Trial Network Connectivity
#'
#' A pre-analysis diagnostic tool. Calculates the number of common genotypes
#' between every pair of environments.
#'
#' @param data Dataframe containing trial data.
#' @param genotype String. Column name for Genotype.
#' @param trial String. Column name for Site/Environment.
#' @param threshold Integer. Minimum number of shared lines required to be considered "connected".
#'
#' @return A list containing the connectivity matrix and a report of disconnected pairs.
#' @details
#' Factor Analytic models require genetic links between trials to estimate correlations.
#' If two sites share 0 varieties, their correlation is undefined. This function flags
#' those risks before you run ASReml.
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

#' Plot Spatial Field Map
#' @export
plot_spatial <- function(model, row_col = "Row", col_col = "Column", type = "residuals") {
  data_name <- as.character(model$call$data)
  if (!exists(data_name)) stop("Original dataframe not found.")
  df <- get(data_name)

  if (type == "residuals") { vals <- resid(model); title <- "Spatial Map: Residuals" }
  else { vals <- fitted(model); title <- "Spatial Map: Fitted Values" }

  if (length(vals) != nrow(df)) warning("Length mismatch. Using naive alignment.")
  df$Value_To_Plot <- vals

  r_vals <- as.numeric(as.character(df[[row_col]]))
  c_vals <- as.numeric(as.character(df[[col_col]]))
  n_r <- max(r_vals, na.rm = TRUE); n_c <- max(c_vals, na.rm = TRUE)

  field_mat <- matrix(NA, nrow = n_r, ncol = n_c)
  for(i in 1:nrow(df)) {
    if (!is.na(r_vals[i]) & !is.na(c_vals[i])) field_mat[r_vals[i], c_vals[i]] <- df$Value_To_Plot[i]
  }
  field_rev <- field_mat[n_r:1, ]

  par(mar = c(3, 3, 3, 5))
  image(1:n_c, 1:n_r, t(field_rev), axes = FALSE, xlab = "", ylab = "",
        col = hcl.colors(20, "Blue-Red 3"), main = title)
  axis(1, cex.axis = 0.7); axis(2, at = 1:n_r, labels = rev(1:n_r), cex.axis = 0.7); box()

  r_val <- range(vals, na.rm = TRUE)
  mtext(sprintf("Range: %.2f to %.2f", r_val[1], r_val[2]), side = 4, line = 1, cex = 0.8)
}

# ==============================================================================
# SECTION 3: CORE EXTRACTOR (fa.asreml)
# ==============================================================================

#' Extract and Rotate Factor Analytic Model Parameters
#'
#' The primary engine of the toolkit. Parses an ASReml-R model object to extract
#' loadings, scores, and specific variances from a Factor Analytic (FA) term.
#' It performs Singular Value Decomposition (SVD) rotation to align the solution
#' with the Principal Component axis, enabling the calculation of Smith & Cullis (2018)
#' FAST indices.
#'
#' @param model A fitted object of class \code{asreml}.
#' @param classify A character string defining the random model term. Must match the
#'        ASReml syntax used in the model call. Supports both standard terms
#'        (e.g., \code{"fa(Site, 2):Genotype"}) and genomic terms
#'        (e.g., \code{"fa(Site, 2):vm(Genotype, grm)"}).
#' @param rotate Logical. If \code{TRUE} (default), loadings and scores are rotated
#'        to the Principal Component solution. This is required for valid interpretation
#'        of Overall Performance (OP) and Stability (RMSD).
#'
#' @return An object of class \code{fa_asreml} containing:
#' \describe{
#'   \item{loadings}{List of \code{raw} and \code{rotated} site loadings.}
#'   \item{scores}{List of \code{raw} and \code{rotated} genotype scores (BLUPs).}
#'   \item{var_comp}{Specific variances (\code{psi}) and Variance Accounted For (\code{vaf}).}
#'   \item{matrices}{Estimated Genetic Covariance (\code{G}) and Correlation (\code{Cor}) matrices.}
#'   \item{fast}{Dataframe of FAST indices: Overall Performance (OP) and Stability (RMSD).}
#' }
#'
#' @references
#' Smith, A. B., & Cullis, B. R. (2018). Plant breeding selection tools built on
#' factor analytic mixed models for multi-environment trial data. \emph{Euphytica}, 214(8).
#'
#' @importFrom dplyr %>% mutate select arrange desc
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_extract
#' @importFrom tibble column_to_rownames
#' @export
fa.asreml <- function(model, classify, rotate = TRUE) {
  cat("--- Starting FA Extraction ---\n")
  vc <- summary(model)$varcomp; vc_names <- rownames(vc)
  coefs <- coef(model)$random; coef_names <- rownames(coefs)

  clean_str <- gsub("\\s+", "", classify)
  parts <- strsplit(clean_str, ":")[[1]]
  fa_part <- parts[grep("fa\\(", parts)]
  gen_part_raw <- parts[grep("fa\\(", parts, invert = TRUE)]

  site_col <- sub("fa\\(([^,]+),.*", "\\1", fa_part)
  k <- as.numeric(sub(".*,([0-9]+)\\).*", "\\1", fa_part))

  # Handle vm() terms
  if (grepl("\\(", gen_part_raw)) {
    gen_col_name <- sub("^[a-z]+\\(([^,)]+).*", "\\1", gen_part_raw)
  } else {
    gen_col_name <- gen_part_raw
  }

  cat(sprintf("-> Input: Site='%s', Term='%s', Factors=%d\n", site_col, gen_part_raw, k))

  # Dynamic Term Detection
  psi_candidates <- grep(paste0(site_col, ".*!var$"), vc_names, value = TRUE)
  if(length(psi_candidates) == 0) stop("CRITICAL: No matching specific variances found in varcomp.")
  actual_term_prefix <- strsplit(psi_candidates[1], "!")[[1]][1]
  term_regex <- gsub("\\(", "\\\\(", actual_term_prefix) %>% gsub("\\)", "\\\\)", .)

  # 1. Psi
  psi_rows <- grep(paste0("^", term_regex, "!.*!var$"), vc_names)
  psi_df <- data.frame(FullTerm = vc_names[psi_rows], Psi = vc[psi_rows, "component"]) %>%
    mutate(Site = str_extract(FullTerm, paste0("(?<=", term_regex, "!)(.*?)(?=!var)")))

  # 2. Loadings
  lambda_df <- data.frame()
  for(i in 1:k) {
    pat <- paste0("^", term_regex, "!.*!fa[_]?", i, "$")
    rows <- grep(pat, vc_names)
    if(length(rows) == 0) next
    tmp <- data.frame(FullTerm = vc_names[rows], Value = vc[rows, "component"], Factor = i) %>%
      mutate(Site = str_extract(FullTerm, paste0("(?<=", term_regex, "!)(.*?)(?=!fa[_]?", i, ")")))
    lambda_df <- bind_rows(lambda_df, tmp)
  }

  lambda_mat <- lambda_df %>% select(Site, Factor, Value) %>%
    pivot_wider(names_from = Factor, values_from = Value) %>%
    column_to_rownames("Site") %>% as.matrix()

  # 3. Scores
  scores_df <- data.frame()
  for(i in 1:k) {
    p1 <- paste0("Comp[_]?", i, ".*", gen_col_name)
    p2 <- paste0("Fac[_]?", i, ".*", gen_col_name)
    rows <- unique(c(grep(p1, coef_names), grep(p2, coef_names)))
    if(length(rows) > 0) {
      tmp <- data.frame(FullTerm = coef_names[rows], Value = coefs[rows, 1], Factor = i) %>%
        mutate(Genotype = sub(paste0(".*", gen_col_name, "[:_]"), "", FullTerm))
      scores_df <- bind_rows(scores_df, tmp)
    }
  }

  has_scores <- FALSE
  f_mat <- matrix(0, 1, k)
  if(nrow(scores_df) > 0) {
    f_mat <- scores_df %>% distinct(Genotype, Factor, .keep_all = TRUE) %>%
      select(Genotype, Factor, Value) %>% pivot_wider(names_from = Factor, values_from = Value, values_fill=0) %>%
      column_to_rownames("Genotype") %>% as.matrix()
    has_scores <- TRUE
  }

  # 4. Rotation & Calculation
  if(rotate && nrow(lambda_mat) > 1) {
    svd_res <- svd(lambda_mat); V <- svd_res$v
    lambda_rot <- lambda_mat %*% V
    if(sum(lambda_rot[,1]) < 0) { lambda_rot[,1] <- -lambda_rot[,1]; V[,1] <- -V[,1] }
    f_rot <- if(has_scores) f_mat %*% V else f_mat
    rot_mat <- V
  } else {
    lambda_rot <- lambda_mat; f_rot <- f_mat; rot_mat <- diag(k)
  }
  colnames(lambda_rot) <- colnames(f_rot) <- paste0("Fac", 1:k)

  # Genetic Matrices
  common <- intersect(rownames(lambda_rot), psi_df$Site)
  lam_ord <- lambda_rot[common, , drop=FALSE]
  psi_ord <- diag(psi_df$Psi[match(common, psi_df$Site)])
  G_est <- (lam_ord %*% t(lam_ord)) + psi_ord
  C_est <- cov2cor(G_est)

  # VAF
  vaf_pct <- (diag(lam_ord %*% t(lam_ord)) / diag(G_est)) * 100
  vaf_df <- data.frame(Site = common, VAF = vaf_pct)

  # FAST
  fast_df <- NULL
  if(has_scores) {
    OP <- f_rot[,1] * mean(lambda_rot[,1], na.rm=TRUE)
    RMSD <- rep(0, nrow(f_rot))
    if(k > 1) {
      h_pred <- f_rot[, 2:k, drop=FALSE] %*% t(lambda_rot[, 2:k, drop=FALSE])
      RMSD <- apply(h_pred, 1, function(x) sqrt(mean(x^2)))
    }
    fast_df <- data.frame(Genotype = rownames(f_rot), OP = OP, RMSD = RMSD) %>% arrange(desc(OP))
  }

  out <- list(loadings = list(raw = lambda_mat, rotated = lambda_rot),
              scores = if(has_scores) list(raw = f_mat, rotated = f_rot) else NULL,
              var_comp = list(psi = psi_df, vaf = vaf_df),
              matrices = list(G = G_est, Cor = C_est),
              fast = fast_df, rotation_matrix = rot_mat, meta = list(k = k))
  class(out) <- "fa_asreml"
  return(out)
}

# ==============================================================================
# SECTION 4: ADVANCED ANALYTICS (D-Opt, iClasses)
# ==============================================================================

#' Calculate D-Optimality (Network Efficiency)
#'
#' Quantifies the information content of the trial network based on the Factor Analytic
#' loadings ($\Lambda$). Implements the D-Optimality criterion to identify redundant sites.
#'
#' @param object An object of class \code{fa_asreml}.
#'
#' @return A list containing:
#' \describe{
#'   \item{total_d}{The total determinant of the information matrix.}
#'   \item{site_impact}{A dataframe ranking sites by their \% contribution to network information.}
#' }
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
#' Implements the methodology of Smith et al. (2021) to partition environments into
#' biological classes based on their response to specific GxE drivers (Factors).
#'
#' @param object An object of class \code{fa_asreml}.
#' @param factor Integer. Which factor defines the interaction? (Default 2).
#' @param threshold Numeric. The loading threshold to assign a site to a class
#'        (e.g., loadings > 0.1 are Positive Class).
#'
#' @return A list containing site classifications and predicted genetic values for
#'         each genotype in the Positive vs Negative classes.
#' @export
get_i_classes <- function(object, factor = 2, threshold = 0.1) {
  if(factor > object$meta$k) stop("Factor not found.")
  lam <- object$loadings$rotated[, factor]; sco <- object$scores$rotated[, factor]
  lam1 <- object$loadings$rotated[, 1]; sco1 <- object$scores$rotated[, 1]

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

# ==============================================================================
# SECTION 5: S3 METHODS & VISUALIZATION
# ==============================================================================

#' Print Method
#' @export
print.fa_asreml <- function(x, ...) {
  cat(sprintf("\n=== FACTOR ANALYTIC REPORT ===\nSites: %d, Factors: %d\n",
              nrow(x$loadings$rotated), x$meta$k))
  diag(x$matrices$Cor) <- NA
  cat(sprintf("Mean VAF: %.1f%%, Mean Cor: %.2f\n",
              mean(x$var_comp$vaf$VAF, na.rm=T), mean(x$matrices$Cor, na.rm=T)))
  if(!is.null(x$fast)) print(head(x$fast, 5))
}

#' Summary Method
#' @export
summary.fa_asreml <- function(object, ...) {
  ld <- object$loadings$rotated; G_d <- diag(object$matrices$G)
  s_stats <- data.frame(Site = rownames(ld))
  for(i in 1:object$meta$k) s_stats[[paste0("VAF_F", i)]] <- round((ld[,i]^2/G_d)*100, 1)
  s_stats$Total_VAF <- rowSums(s_stats[,-1]); s_stats$Fac1 <- round(ld[,1],3)

  g_stats <- if(!is.null(object$fast)) object$fast else NULL
  list(site_stats = s_stats[order(s_stats$Total_VAF, decreasing=T),], genotypes = g_stats)
}

#' Master Plotting Method for Factor Analytic Models
#'
#' A unified visualization interface for GxE analysis. Dispatches to specific
#' plotting routines based on the \code{type} argument.
#'
#' @param x An object of class \code{fa_asreml}.
#' @param type Character string specifying the plot type:
#' \describe{
#'   \item{\code{"fast"}}{Plots Overall Performance (OP) vs Stability (RMSD). The "Money Plot" for selection.}
#'   \item{\code{"heatmap"}}{Visualizes the Genetic Correlation Matrix between environments.}
#'   \item{\code{"latent_reg"}}{Plots Latent Regression lines (GxE Drivers). "Peels off" Factor 1 to show interaction effects.}
#'   \item{\code{"biplot"}}{Generates a GGE-style Biplot (Vectors + Points).}
#'   \item{\code{"vaf"}}{Bar chart of Variance Accounted For (\%) per site (Quality Control).}
#'   \item{\code{"d_opt"}}{Bar chart of Information Loss (D-Optimality) to identify redundant sites.}
#'   \item{\code{"diff"}}{Slope graph showing rank changes between Interaction Classes (Crossover Interaction).}
#' }
#' @param factor Integer or Vector. Which factors to plot? Used for \code{latent_reg}, \code{biplot}, and \code{diff}.
#' @param n_label Integer. Number of top-performing genotypes to automatically label.
#' @param highlight Character vector. Names of specific genotypes to highlight (e.g., check varieties).
#' @param ... Additional graphical parameters passed to base R plotting functions.
#'
#' @examples
#' \dontrun{
#'   # Selection Plot highlighting Checks
#'   plot(results, type = "fast", highlight = c("CheckA", "CheckB"))
#'
#'   # View Genetic Correlations
#'   plot(results, type = "heatmap")
#'
#'   # Biplot of Factor 1 vs Factor 2
#'   plot(results, type = "biplot")
#' }
#' @export
plot.fa_asreml <- function(x, type = "fast", factor = NULL, n_label = 5, highlight = NULL, ...) {
  old_par <- par(no.readonly = TRUE); on.exit(par(old_par))

  if (type == "fast") { .plot_fast(x, n_label, highlight) }
  else if (type == "heatmap") { .plot_heat(x) }
  else if (type == "latent_reg") { .plot_reg(x, factor, highlight, n_label) }
  else if (type == "biplot") { .plot_biplot(x, if(is.null(factor)) c(1,2) else factor, highlight) }
  else if (type == "vaf") { .plot_vaf(x) }
  else if (type == "d_opt") { .plot_dopt(get_d_optimality(x)) }
  else if (type == "diff") { .plot_diff(get_i_classes(x, if(is.null(factor)) 2 else factor), n_label, highlight) }
  else { stop("Unknown type.") }
}

# --- Plot Internals ---

.plot_fast <- function(x, n, h) {
  df <- x$fast; top <- head(df$Genotype, n)
  col <- rep("grey60", nrow(df)); bg <- rep("grey95", nrow(df)); pch <- rep(21, nrow(df))
  is_h <- df$Genotype %in% h; is_t <- df$Genotype %in% top
  bg[is_t] <- "#3498db"; col[is_t] <- "#2980b9"
  bg[is_h] <- "#e74c3c"; col[is_h] <- "#c0392b"; pch[is_h] <- 23

  par(mar = c(5, 5, 4, 2))
  plot(df$RMSD, df$OP, pch=pch, bg=bg, col=col, xlab="Stability (RMSD)", ylab="Performance (OP)", main="FAST Selection")
  grid(); abline(h=0, v=mean(df$RMSD), lty=2, col="grey")
  lbl <- which(is_h | is_t)
  if(length(lbl)>0) text(df$RMSD[lbl], df$OP[lbl], labels=df$Genotype[lbl], pos=3, cex=0.7)
}

.plot_heat <- function(x) {
  cor <- x$matrices$Cor; p <- ncol(cor); rev <- cor[, p:1]
  par(mar = c(6, 6, 4, 2))
  image(1:p, 1:p, rev, axes=F, col=hcl.colors(20, "RdBu", rev=T), breaks=seq(-1,1,l=21), main="Correlation")
  axis(1, at=1:p, labels=colnames(cor), las=2, cex.axis=0.7)
  axis(2, at=1:p, labels=rev(rownames(cor)), las=1, cex.axis=0.7)
  for(i in 1:p) for(j in 1:p) if(rev[i,j]<0.99) text(i, j, sprintf("%.2f", rev[i,j]), cex=0.6)
}

.plot_reg <- function(x, fac, h, n) {
  k <- x$meta$k; if(is.null(fac)) fac <- 1:k
  if(length(fac)>1) par(mfrow=c(ceiling(length(fac)/2), min(length(fac),2)))

  Lam <- x$loadings$rotated; Sco <- x$scores$rotated; U_com <- Sco %*% t(Lam)
  tgt <- unique(c(head(rownames(Sco), n), h))

  for(f in fac) {
    if(f>k) next
    Y <- if(f>1) U_com - (Sco[,1:(f-1)] %*% t(Lam[,1:(f-1)])) else U_com
    xv <- Lam[,f]; ysub <- Y[tgt, , drop=F]

    par(mar=c(4,4,3,1))
    plot(1, type="n", xlim=range(xv), ylim=range(ysub), xlab=paste("Load Fac", f), ylab="Effect", main=paste("Factor", f))
    grid(); abline(h=0, lty=2)
    axis(1, at=xv, labels=F, tck=-0.02, col="green")

    for(g in tgt) {
      sl <- Sco[g,f]; cc <- if(g %in% h) "red" else "navy"
      points(xv, Y[g,], pch=19, col=adjustcolor(cc,0.3), cex=0.7)
      abline(0, sl, col=cc); text(max(xv), sl*max(xv), labels=g, pos=4, cex=0.7, col=cc)
    }
  }
}

.plot_biplot <- function(x, fac, h) {
  Lam <- x$loadings$rotated[,fac]; Sco <- x$scores$rotated[,fac]
  sf <- max(abs(Sco)) / max(abs(Lam)) * 0.8; Lam_s <- Lam * sf
  lim <- range(rbind(Lam_s, Sco))

  par(mar=c(5,5,4,2))
  plot(1, type="n", xlim=lim, ylim=lim, asp=1, xlab=paste("Fac", fac[1]), ylab=paste("Fac", fac[2]), main="Biplot")
  grid(); abline(h=0, v=0, lty=2, col="grey")
  arrows(0,0, Lam_s[,1], Lam_s[,2], col="darkgreen", length=0.1)
  text(Lam_s, labels=rownames(Lam), col="darkgreen", pos=4, cex=0.7)

  bg <- rep("grey90", nrow(Sco)); col <- rep("grey60", nrow(Sco)); pch <- rep(21, nrow(Sco))
  is_h <- rownames(Sco) %in% h
  bg[is_h] <- "red"; col[is_h] <- "darkred"; pch[is_h] <- 23
  points(Sco, pch=pch, bg=bg, col=col, cex=0.8)
  if(any(is_h)) text(Sco[is_h,], labels=rownames(Sco)[is_h], pos=3, cex=0.7)

  hpts <- chull(Sco); lines(Sco[c(hpts, hpts[1]),], col="navy", lty=2)
}

.plot_vaf <- function(x) {
  df <- x$var_comp$vaf[order(x$var_comp$vaf$VAF), ]
  col <- ifelse(df$VAF < 50, "red", "blue")
  par(mar=c(5,7,4,2))
  barplot(df$VAF, names.arg=df$Site, horiz=T, las=1, col=col, xlab="VAF %", main="Site Quality")
  abline(v=80, lty=2)
}

.plot_dopt <- function(d) {
  df <- d$site_impact[order(d$site_impact$Impact_Pct), ]
  col <- ifelse(df$Impact_Pct < 1, "red", "green")
  par(mar=c(5,8,4,2))
  barplot(df$Impact_Pct, names.arg=df$Site, horiz=T, las=1, col=col, xlab="% Info Loss", main="Network Efficiency")
}

.plot_diff <- function(res, n, h) {
  df <- res$gen_effects; df <- na.omit(df)
  tgt <- unique(c(head(df$Genotype, n), tail(df$Genotype, n), h))
  plot_df <- df[df$Genotype %in% tgt, ]

  yl <- range(c(plot_df$Pred_Neg, plot_df$Pred_Pos))
  par(mar=c(4,4,3,6), xpd=F)
  plot(1, type="n", xlim=c(0.8, 2.2), ylim=yl, xaxt="n", ylab="Genetic Value", xlab="", main="Interaction Classes")
  axis(1, at=1:2, labels=c("Class (-)", "Class (+)"))

  for(i in 1:nrow(plot_df)) {
    g <- plot_df$Genotype[i]
    cc <- if(g %in% h) "red" else "grey50"
    segments(1, plot_df$Pred_Neg[i], 2, plot_df$Pred_Pos[i], col=cc, lwd=if(g %in% h) 2 else 1)
    text(1, plot_df$Pred_Neg[i], g, pos=2, cex=0.7, col=cc, xpd=T)
    text(2, plot_df$Pred_Pos[i], g, pos=4, cex=0.7, col=cc, xpd=T)
  }
}
