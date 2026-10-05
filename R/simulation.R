# Evaluate `code` with a fixed seed without disturbing the caller's RNG stream.
with_local_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(code)
  }
  genv <- globalenv()
  old <- if (exists(".Random.seed", envir = genv, inherits = FALSE)) get(".Random.seed", envir = genv) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = genv) else assign(".Random.seed", old, envir = genv))
  set.seed(seed)
  code
}

# Randomised complete block layout: each replicate occupies its own band of
# rows and genotypes are randomised within it.
simulate_layout <- function(gens, n_reps, n_cols) {
  n_gen <- length(gens)
  rows_per_rep <- ceiling(n_gen / n_cols)
  pos <- seq_len(n_gen) - 1
  do.call(rbind, lapply(seq_len(n_reps), function(r) {
    data.frame(
      Genotype = sample(gens), Rep = r,
      Row = (r - 1) * rows_per_rep + pos %/% n_cols + 1,
      Column = pos %% n_cols + 1,
      stringsAsFactors = FALSE
    )
  }))
}

spatial_trend <- function(row, col) sin(row / 5) * 2 + col * 0.1

#' Generate Synthetic Genomic Relationship Matrix (GRM)
#'
#' @description
#' Simulates biallelic SNP dosages and builds the VanRaden (2008) method 1 GRM
#' \deqn{G = \frac{ZZ'}{2 \sum p_i (1-p_i)}}
#' where \eqn{Z = M - 2P} is the centred dosage matrix. Centring makes
#' \eqn{G} singular, so \code{ridge} is added to the diagonal.
#'
#' @param n_gen Integer. Number of genotypes.
#' @param n_markers Integer. Number of SNP markers.
#' @param seed Integer seed (\code{NULL} to use the current RNG state). The
#'   caller's RNG state is restored afterwards.
#' @param ridge Numeric. Value added to the diagonal (default 0.01).
#' @return A list with \code{grm} (n_gen x n_gen) and \code{markers} (0/1/2 dosages).
#' @references VanRaden, P. M. (2008). Efficient methods to compute genomic
#'   predictions. \emph{Journal of Dairy Science} 91, 4414-4423.
#' @export
generate_grm <- function(n_gen = 200, n_markers = 1000, seed = 123, ridge = 0.01) {
  with_local_seed(seed, {
    maf <- stats::runif(n_markers, 0.1, 0.5)
    M <- matrix(stats::rbinom(n_gen * n_markers, 2, rep(maf, each = n_gen)), n_gen, n_markers)
  })
  gen_names <- sprintf("G%03d", seq_len(n_gen))
  dimnames(M) <- list(gen_names, sprintf("M%03d", seq_len(n_markers)))
  p <- colMeans(M) / 2
  Z <- sweep(M, 2, 2 * p)
  GRM <- tcrossprod(Z) / (2 * sum(p * (1 - p)))
  diag(GRM) <- diag(GRM) + ridge
  dimnames(GRM) <- list(gen_names, gen_names)
  list(grm = GRM, markers = M)
}

#' Generate Synthetic MET Data
#'
#' Simulates a multi-environment trial with a two-factor GxE structure, a
#' spatial trend and residual noise. Each site is laid out as a randomised
#' complete block design (one band of rows per replicate, genotypes randomised
#' within replicate), so genetic effects are not confounded with field position.
#'
#' @param n_sites Integer. Number of environments.
#' @param n_gen Integer. Number of genotypes.
#' @param n_reps Integer. Replicates per site.
#' @param sparsity Numeric in [0, 1). Proportion of plots removed at random.
#' @param seed Integer seed (\code{NULL} to use the current RNG state).
#' @param n_cols Integer. Field columns per site (default 20).
#'
#' @details Genetic values are \eqn{s_{1i} \lambda_{1j} + s_{2i} \lambda_{2j}}
#'   with \eqn{\lambda_1 \sim U(5, 15)} (main effect) and
#'   \eqn{\lambda_2 \sim U(-8, 8)} (crossover GxE); the trend is
#'   \eqn{2 \sin(row / 5) + 0.1\, col}; residual SD is 2.5.
#'
#' @return A data frame with \code{Site}, \code{Genotype}, \code{Rep},
#'   \code{Row}, \code{Column} and \code{Yield}.
#' @export
generate_met_data <- function(n_sites = 10, n_gen = 100, n_reps = 2, sparsity = 0.2, seed = 123, n_cols = 20) {
  if (sparsity < 0 || sparsity >= 1) cli::cli_abort("{.arg sparsity} must be in [0, 1).")
  sites <- sprintf("Site_%02d", seq_len(n_sites))
  gens <- sprintf("G%03d", seq_len(n_gen))

  with_local_seed(seed, {
    L1 <- stats::runif(n_sites, 5, 15)
    S1 <- stats::rnorm(n_gen)
    L2 <- stats::runif(n_sites, -8, 8)
    S2 <- stats::rnorm(n_gen)
    G_mat <- outer(S1, L1) + outer(S2, L2)
    dimnames(G_mat) <- list(gens, sites)

    data <- do.call(rbind, lapply(sites, function(site) {
      d <- simulate_layout(gens, n_reps, n_cols)
      d$Site <- site
      d$Yield <- 50 + G_mat[cbind(d$Genotype, site)] + spatial_trend(d$Row, d$Column) + stats::rnorm(nrow(d), 0, 2.5)
      d
    }))
    if (sparsity > 0) data <- data[-sample(nrow(data), floor(nrow(data) * sparsity)), , drop = FALSE]
  })
  finalise_met(data)
}

#' Generate Genomic MET Data
#'
#' @description
#' Simulates a MET whose genetic values follow a user-supplied GRM: two
#' latent factors \eqn{f_r \sim N(0, G)} combined with site loadings as in
#' \code{\link{generate_met_data}}. Layout is a randomised complete block
#' design per site.
#'
#' @param n_sites Integer. Number of environments.
#' @param grm Matrix. Relationship matrix (e.g. from \code{\link{generate_grm}}).
#' @param n_reps Integer. Replicates per site.
#' @param h2 Numeric in (0, 1). Plot-level heritability within each site,
#'   \eqn{\sigma^2_g / (\sigma^2_g + \sigma^2_e)}. The spatial trend is added on
#'   top and excluded from this ratio, because a spatial model is expected to
#'   remove it.
#' @param seed Integer seed (\code{NULL} to use the current RNG state).
#' @param n_cols Integer. Field columns per site (default 20).
#'
#' @return A data frame with \code{Site}, \code{Genotype}, \code{Rep},
#'   \code{Row}, \code{Column} and \code{Yield}.
#' @export
generate_genomic_met <- function(n_sites = 10, grm, n_reps = 2, h2 = 0.5, seed = 123, n_cols = 20) {
  if (!is.numeric(h2) || h2 <= 0 || h2 >= 1) cli::cli_abort("{.arg h2} must be in (0, 1).")
  gens <- rownames(grm)
  if (is.null(gens)) cli::cli_abort("{.arg grm} must have row names.")
  sites <- sprintf("Site_%02d", seq_len(n_sites))

  eg <- eigen(grm, symmetric = TRUE)
  root <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0)), nrow = length(eg$values))

  with_local_seed(seed, {
    S1 <- as.vector(root %*% stats::rnorm(length(gens)))
    S2 <- as.vector(root %*% stats::rnorm(length(gens)))
    G_mat <- outer(S1, stats::runif(n_sites, 5, 15)) + outer(S2, stats::runif(n_sites, -8, 8))
    dimnames(G_mat) <- list(gens, sites)

    data <- do.call(rbind, lapply(sites, function(site) {
      d <- simulate_layout(gens, n_reps, n_cols)
      d$Site <- site
      g_eff <- G_mat[cbind(d$Genotype, site)]
      var_e <- stats::var(G_mat[, site]) * (1 / h2 - 1)
      d$Yield <- 50 + g_eff + spatial_trend(d$Row, d$Column) + stats::rnorm(nrow(d), 0, sqrt(var_e))
      d
    }))
  })
  finalise_met(data)
}

finalise_met <- function(data) {
  data$Site <- factor(data$Site)
  data$Genotype <- factor(data$Genotype)
  data$Row <- factor(data$Row, levels = sort(unique(data$Row)))
  data$Column <- factor(data$Column, levels = sort(unique(data$Column)))
  data <- data[order(data$Site, data$Column, data$Row), c("Site", "Genotype", "Rep", "Row", "Column", "Yield")]
  rownames(data) <- NULL
  data
}
