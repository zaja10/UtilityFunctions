test_that("FA Extraction math is valid on synthetic data", {

  # Skip if ASReml isn't installed (e.g. on a CI server)
  skip_if_not_installed("asreml")
  library(asreml)
  library(tidyverse)

  # 1. SETUP
  set.seed(101)
  # Generate small dataset for speed
  df <- generate_met_data(n_sites = 6, n_gen = 50, n_reps = 2)

  # 2. FIT MODEL (Suppress output)
  # We use a standard FA2 model
  capture.output({
    mod <- asreml(fixed = Yield ~ Site,
                  random = ~fa(Site, 2):Genotype + Rep,
                  residual = ~idv(units),
                  data = df, trace = FALSE)
  })

  # 3. RUN EXTRACTOR
  res <- fa.asreml(mod, classify = "fa(Site, 2):Genotype")

  # --- TEST 1: Structure ---
  expect_s3_class(res, "fa_asreml")
  expect_equal(res$meta$k, 2)
  expect_equal(nrow(res$loadings$rotated), 6) # 6 sites

  # --- TEST 2: Mathematical Reconstruction (G = Lam Lam' + Psi) ---
  # Extract raw components
  lam <- res$loadings$raw
  psi <- diag(res$var_comp$psi$Psi)

  # Reconstruct G manually
  G_manual <- (lam %*% t(lam)) + psi

  # Compare to the function's calculated G
  # Tolerance allows for tiny floating point differences
  expect_equal(res$matrices$G, G_manual, tolerance = 1e-6)

  # --- TEST 3: Rotation Invariance ---
  # The G matrix should be IDENTICAL whether using Raw or Rotated loadings
  lam_rot <- res$loadings$rotated
  G_rot <- (lam_rot %*% t(lam_rot)) + psi

  expect_equal(G_manual, G_rot, tolerance = 1e-6)

  # --- TEST 4: FAST Index Logic ---
  # OP should equal Factor 1 Score * Mean Factor 1 Loading
  gen_1 <- res$fast$Genotype[1]

  score_1 <- res$scores$rotated[gen_1, 1]
  mean_load <- mean(res$loadings$rotated[, 1])
  calc_op <- score_1 * mean_load

  expect_equal(res$fast$OP[1], calc_op, tolerance = 1e-5)
})

test_that("Robust Regex handles Reduced Rank models", {
  skip_if_not_installed("asreml")

  # Simulate data
  df <- generate_met_data(n_sites = 4, n_gen = 50)

  # Fit RR model (Split terms)
  capture.output({
    mod_rr <- asreml(fixed = Yield ~ Site,
                     random = ~rr(Site, 1):Genotype + diag(Site):Genotype,
                     residual = ~idv(units),
                     data = df, trace = FALSE)
  })

  # Try extraction (Should not crash)
  expect_error(
    res <- fa.asreml(mod_rr, classify = "rr(Site, 1):Genotype", psi_term = "diag(Site):Genotype"),
    NA # NA means "Expect NO error"
  )

  # Check if Psi was actually found (not all zeros)
  expect_true(sum(res$var_comp$psi$Psi) > 0)
})
