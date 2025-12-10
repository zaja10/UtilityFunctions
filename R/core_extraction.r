# ==============================================================================
# CORE EXTRACTION ENGINE
# ==============================================================================

#' Extract and Rotate Factor Analytic Model Parameters
#'
#' @description
#' The primary engine of the toolkit. This function parses a fitted ASReml-R model object to
#' extract loadings, scores (BLUPs), and specific variances from a Factor Analytic (FA)
#' or Reduced Rank (RR) model term.
#'
#' It performs a Singular Value Decomposition (SVD) rotation to align the solution
#' with the Principal Component axis. This rotation is required to validly calculate
#' the **Factor Analytic Selection Tools (FAST)** indices (Overall Performance and Stability)
#' as proposed by Smith & Cullis (2018).
#'
#' @param model A fitted object of class \code{asreml}.
#' @param classify A character string defining the main latent model term.
#'        Must match the ASReml syntax used in the model call.
#'        \itemize{
#'          \item Standard FA: \code{"fa(Site, 2):Genotype"}
#'          \item Genomic FA: \code{"fa(Site, 2):vm(Genotype, grm)"}
#'          \item Reduced Rank: \code{"rr(Site, 2):Genotype"}
#'        }
#' @param psi_term A character string defining the specific variance term (Optional).
#'        Required only for **Reduced Rank** models where the diagonal variance is fitted
#'        separately (e.g., \code{"diag(Site):Genotype"}). If NULL for an RR model,
#'        specific variances are assumed to be zero.
#' @param rotate Logical. If \code{TRUE} (default), loadings and scores are rotated
#'        to the Principal Component solution via SVD. If \code{FALSE}, raw parameters
#'        are returned (warning: FAST indices may be invalid without rotation).
#'
#' @return An object of class \code{fa_asreml} containing:
#' \describe{
#'   \item{loadings}{A list containing \code{raw} and \code{rotated} site loading matrices.}
#'   \item{scores}{A list containing \code{raw} and \code{rotated} genotype score matrices (BLUPs).}
#'   \item{var_comp}{A list containing the extracted specific variances (\code{psi}) and a table of Variance Accounted For (\code{vaf}).}
#'   \item{matrices}{The estimated Genetic Covariance (\code{G}) and Genetic Correlation (\code{Cor}) matrices.}
#'   \item{fast}{A dataframe containing the FAST indices: Overall Performance (\code{OP}) and Stability (\code{RMSD}).}
#'   \item{meta}{Metadata regarding the model structure (k factors, model type).}
#' }
#'
#' @details
#' \strong{Rotation Strategy:}
#' Raw FA solutions from ASReml are not unique. This function rotates the loadings (\eqn{\Lambda})
#' and scores (\eqn{f}) such that:
#' \deqn{\Lambda_{rot} = \Lambda V}
#' \deqn{f_{rot} = f V}
#' where \eqn{V} is the matrix of eigenvectors from the SVD of \eqn{\Lambda}. This ensures that
#' Factor 1 captures the maximum amount of genetic variance, allowing it to be interpreted
#' as "Overall Performance."
#'
#' \strong{FAST Indices:}
#' \itemize{
#'   \item \strong{OP (Overall Performance):} The predicted genetic value of a variety
#'   at the mean environmental loading of Factor 1.
#'   \item \strong{RMSD (Stability):} The Root Mean Square Deviation of a variety's
#'   predicted effect from the Factor 1 regression line. It quantifies the amount of
#'   GxE interaction (instability) driven by Factors 2..k.
#' }
#'
#' @references
#' Smith, A. B., & Cullis, B. R. (2018). Plant breeding selection tools built on
#' factor analytic mixed models for multi-environment trial data. \emph{Euphytica}, 214(8).
#'
#' @importFrom dplyr %>% mutate select arrange desc bind_rows distinct
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_extract
#' @importFrom tibble column_to_rownames
#' @importFrom stats coef cov2cor sd
#' @export
fa.asreml <- function(model, classify, psi_term = NULL, rotate = TRUE) {
  cat("--- Starting FA/RR Extraction ---\n")

  # 1. SETUP ------------------------------------------------------------------
  vc <- summary(model)$varcomp
  vc_names <- rownames(vc)
  coefs <- coef(model)$random
  coef_names <- rownames(coefs)

  clean_str <- gsub("\\s+", "", classify)
  parts <- strsplit(clean_str, ":")[[1]]
  latent_part <- parts[grep("(fa|rr)\\(", parts)]
  gen_part_raw <- parts[grep("(fa|rr)\\(", parts, invert = TRUE)]

  site_col <- sub("(fa|rr)\\(([^,]+),.*", "\\1", latent_part)
  k <- as.numeric(sub(".*,([0-9]+)\\).*", "\\1", latent_part))

  if (grepl("\\(", gen_part_raw)) {
    gen_col_name <- sub("^[a-z]+\\(([^,)]+).*", "\\1", gen_part_raw)
  } else {
    gen_col_name <- gen_part_raw
  }

  is_rr <- grepl("rr\\(", clean_str)
  model_type <- if (is_rr) "Reduced Rank (RR)" else "Factor Analytic (FA)"

  cat(sprintf("-> Type: %s | Site: '%s' | Factors: %d\n", model_type, site_col, k))

  # 2. IDENTIFY PREFIX --------------------------------------------------------
  actual_term_prefix <- NULL
  if (!is_rr) {
    # Standard FA: Look for !var
    psi_candidates <- grep(paste0(site_col, ".*!var$"), vc_names, value = TRUE)
    if (length(psi_candidates) == 0) stop(paste("Could not find !var for site:", site_col))
    actual_term_prefix <- strsplit(psi_candidates[1], "!")[[1]][1]
  } else {
    # Reduced Rank: Look for loadings (!rr_1)
    load_candidates <- grep(paste0(site_col, ".*!(rr|fa)[_]?1"), vc_names, value = TRUE)
    if (length(load_candidates) == 0) stop(paste("Could not find RR loadings for site:", site_col))
    actual_term_prefix <- strsplit(load_candidates[1], "!")[[1]][1]
  }
  term_regex <- gsub("\\(", "\\\\(", actual_term_prefix) %>% gsub("\\)", "\\\\)", .)

  # 3. LOADINGS ---------------------------------------------------------------
  lambda_df <- data.frame()
  for (i in 1:k) {
    pat <- paste0("^", term_regex, "!.*!(fa|rr)[_]?", i, "$")
    rows <- grep(pat, vc_names)
    if (length(rows) == 0) next
    tmp <- data.frame(FullTerm = vc_names[rows], Value = vc[rows, "component"], Factor = i) %>%
      mutate(Site = str_extract(FullTerm, paste0("(?<=", term_regex, "!)(.*?)(?=!(fa|rr)[_]?", i, ")")))
    lambda_df <- bind_rows(lambda_df, tmp)
  }
  if (nrow(lambda_df) == 0) stop("No loadings found.")
  lambda_mat <- lambda_df %>%
    select(Site, Factor, Value) %>%
    pivot_wider(names_from = Factor, values_from = Value) %>%
    column_to_rownames("Site") %>%
    as.matrix()

  # 4. PSI (SPECIFIC VARIANCES) -----------------------------------------------
  psi_df <- NULL
  if (!is_rr) {
    psi_rows <- grep(paste0("^", term_regex, "!.*!var$"), vc_names)
    psi_df <- data.frame(FullTerm = vc_names[psi_rows], Psi = vc[psi_rows, "component"]) %>%
      mutate(Site = str_extract(FullTerm, paste0("(?<=", term_regex, "!)(.*?)(?=!var)")))
  } else {
    if (!is.null(psi_term)) {
      cat("-> Hunting for specific variances...\n")
      # Identify candidate rows (exclude RR term)
      candidate_rows <- vc_names[!grepl(paste0("^", term_regex, "!"), vc_names)]
      sites_to_find <- rownames(lambda_mat)
      found_psi <- list()

      for (site in sites_to_find) {
        # Regex: Ends with (underscore or !) + SiteName + (optional !var) + End of String
        # This catches "diag(Site):Gen!Site_01" AND "diag(Site):Gen!Site_01!var"
        site_pat <- paste0("(_|!)", site, "(!var)?$")
        match <- grep(site_pat, candidate_rows, value = TRUE)

        if (length(match) > 0) {
          idx <- which(vc_names == match[1]) # Take first match
          found_psi[[site]] <- data.frame(Site = site, Psi = vc[idx, "component"])
        }
      }
      if (length(found_psi) > 0) {
        psi_df <- do.call(rbind, found_psi)
      } else {
        stop("Could not find specific variances for your sites.")
      }
    } else {
      psi_df <- data.frame(Site = rownames(lambda_mat), Psi = 0)
    }
  }

  # 5. SCORES -----------------------------------------------------------------
  scores_df <- data.frame()
  for (i in 1:k) {
    p1 <- paste0("Comp[_]?", i, ".*", gen_col_name)
    p2 <- paste0("Fac[_]?", i, ".*", gen_col_name)
    rows <- unique(c(grep(p1, coef_names), grep(p2, coef_names)))
    if (length(rows) > 0) {
      tmp <- data.frame(FullTerm = coef_names[rows], Value = coefs[rows, 1], Factor = i) %>%
        mutate(Genotype = sub(paste0(".*", gen_col_name, "[:_]"), "", FullTerm))
      scores_df <- bind_rows(scores_df, tmp)
    }
  }
  has_scores <- FALSE
  f_mat <- matrix(0, 1, k)
  if (nrow(scores_df) > 0) {
    f_mat <- scores_df %>%
      distinct(Genotype, Factor, .keep_all = T) %>%
      select(Genotype, Factor, Value) %>%
      pivot_wider(names_from = Factor, values_from = Value, values_fill = 0) %>%
      column_to_rownames("Genotype") %>%
      as.matrix()
    has_scores <- TRUE
  }

  # 6. ROTATION ---------------------------------------------------------------
  if (rotate && nrow(lambda_mat) > 1) {
    svd_res <- svd(lambda_mat)
    V <- svd_res$v
    lambda_rot <- lambda_mat %*% V
    if (sum(lambda_rot[, 1]) < 0) {
      lambda_rot[, 1] <- -lambda_rot[, 1]
      V[, 1] <- -V[, 1]
    }
    f_rot <- if (has_scores) f_mat %*% V else f_mat
    rot_mat <- V
  } else {
    lambda_rot <- lambda_mat
    f_rot <- f_mat
    rot_mat <- diag(k)
  }
  colnames(lambda_rot) <- colnames(f_rot) <- paste0("Fac", 1:k)

  # 7. RECONSTRUCTION ---------------------------------------------------------
  common_sites <- intersect(rownames(lambda_rot), psi_df$Site)
  if (length(common_sites) == 0) stop("Site mismatch between Lambda and Psi.")

  lam_ord <- lambda_rot[common_sites, , drop = FALSE]
  psi_ord <- diag(psi_df$Psi[match(common_sites, psi_df$Site)])
  G_est <- (lam_ord %*% t(lam_ord)) + psi_ord
  C_est <- tryCatch(cov2cor(G_est), error = function(e) G_est)

  vaf_pct <- ifelse(diag(G_est) > 1e-8, (diag(lam_ord %*% t(lam_ord)) / diag(G_est)) * 100, 0)
  vaf_df <- data.frame(Site = common_sites, VAF = vaf_pct)

  # 8. FAST INDICES -----------------------------------------------------------
  fast_df <- NULL
  if (has_scores) {
    OP <- f_rot[, 1] * mean(lambda_rot[, 1], na.rm = TRUE)
    RMSD <- if (k > 1) apply(f_rot[, 2:k, drop = F] %*% t(lambda_rot[, 2:k, drop = F]), 1, function(x) sqrt(mean(x^2))) else rep(0, nrow(f_rot))
    fast_df <- data.frame(Genotype = rownames(f_rot), OP = OP, RMSD = RMSD) %>% arrange(desc(OP))
  }

  out <- list(
    loadings = list(raw = lambda_mat, rotated = lambda_rot),
    scores = if (has_scores) list(raw = f_mat, rotated = f_rot) else NULL,
    var_comp = list(psi = psi_df, vaf = vaf_df),
    matrices = list(G = G_est, Cor = C_est),
    fast = fast_df, rotation_matrix = rot_mat, meta = list(k = k, type = model_type, classify = classify)
  )
  class(out) <- "fa_asreml"
  return(out)
}


# ==============================================================================
# S3 METHODS (Print & Summary)
# ==============================================================================

#' Print Method for FA Object
#' @export
print.fa_asreml <- function(x, ...) {
  k <- x$meta$k
  n_sites <- nrow(x$loadings$rotated)

  cor_mat <- x$matrices$Cor
  diag(cor_mat) <- NA
  mean_r <- mean(cor_mat, na.rm = TRUE)
  mean_vaf <- mean(x$var_comp$vaf$VAF, na.rm = TRUE)

  cat("\n=== FACTOR ANALYTIC REPORT (ASRemlFAST) ===\n")
  cat(sprintf(" Model Type: %s\n", x$meta$type))
  cat(sprintf(" Dimensions: %d Sites, %d Factors\n", n_sites, k))
  cat(sprintf(" Rotation  : SVD (PC Solution)\n"))
  cat("\n--- Diagnostics ---\n")
  cat(sprintf(" Mean Genetic VAF: %.1f%%\n", mean_vaf))
  cat(sprintf(" Mean Genetic Cor: %.2f\n", mean_r))

  if (!is.null(x$fast)) {
    cat("\n--- Top 5 Genotypes (by OP) ---\n")
    print(format(head(x$fast[, c("Genotype", "OP", "RMSD")], 5), digits = 3), row.names = FALSE)
  } else {
    cat("\n(Genotype scores not found - FAST indices unavailable)\n")
  }
  cat("===========================================\n")
}

#' Summary Method for FA Object
#' @export
summary.fa_asreml <- function(object, ...) {
  loadings <- object$loadings$rotated
  G_diag <- diag(object$matrices$G)
  k <- object$meta$k

  # Build Site Stats with breakdown per factor
  site_stats <- data.frame(Site = rownames(loadings))
  for (i in 1:k) {
    vaf_k <- (loadings[, i]^2 / G_diag) * 100
    site_stats[[paste0("VAF_Fac", i)]] <- round(vaf_k, 1)
  }
  site_stats$Total_VAF <- rowSums(site_stats[, 2:(k + 1)])
  site_stats$Fac1_Load <- round(loadings[, 1], 3)

  gen_stats <- NULL
  if (!is.null(object$fast)) {
    gen_stats <- object$fast
    gen_stats$Rank_OP <- rank(-gen_stats$OP)
    gen_stats <- gen_stats[, c("Rank_OP", "Genotype", "OP", "RMSD")]
  }

  res <- list(
    site_stats = site_stats[order(site_stats$Total_VAF, decreasing = TRUE), ],
    genotypes = gen_stats,
    correlations = round(object$matrices$Cor, 3),
    meta = object$meta
  )

  class(res) <- "summary.fa_asreml"
  return(res)
}

#' Print Method for Summary
#' @export
print.summary.fa_asreml <- function(x, ...) {
  cat("--- Site Statistics (Top 5 by Total VAF) ---\n")
  print(head(x$site_stats, 5), row.names = FALSE)

  if (!is.null(x$genotypes)) {
    cat("\n--- Genotype Selection (Top 5 by OP) ---\n")
    print(head(x$genotypes, 5), row.names = FALSE)
  }
}
