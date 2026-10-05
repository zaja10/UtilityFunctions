test_that("model ranking puts converged fits first, then AIC", {
  res <- data.frame(Model = c("a", "b", "c"), AIC = c(10, 5, 7), Converged = c(TRUE, FALSE, TRUE))
  expect_equal(select_best_model(res), "c")
  expect_equal(UtilityFunctions:::rank_model_table(res)$Model, c("c", "a", "b"))
  none <- data.frame(Model = c("a", "b"), AIC = c(3, 2), Converged = c(FALSE, FALSE))
  expect_warning(best <- select_best_model(none), "No model converged")
  expect_equal(best, "b")
})

test_that("mean_pairwise_vd equals the brute-force mean squared SED", {
  set.seed(1)
  A <- matrix(rnorm(25), 5)
  C <- crossprod(A)
  pairs <- combn(5, 2)
  brute <- mean(apply(pairs, 2, function(p) C[p[1], p[1]] + C[p[2], p[2]] - 2 * C[p[1], p[2]]))
  expect_equal(UtilityFunctions:::mean_pairwise_vd(C), brute)
  # Squaring the mean SED (old approach) is biased low (Jensen).
  sed <- sqrt(apply(pairs, 2, function(p) C[p[1], p[1]] + C[p[2], p[2]] - 2 * C[p[1], p[2]]))
  expect_lt(mean(sed)^2, brute)
})

test_that("harmonic_reps uses observed replication", {
  expect_equal(UtilityFunctions:::harmonic_reps(c("a", "a", "b", "b", "c")), 3 / (1 / 2 + 1 / 2 + 1))
})

test_that("compare_spatial_models rejects different fixed effects", {
  models <- list(m1 = list(fixed = y ~ 1, random = ~g), m2 = list(fixed = y ~ x, random = ~g))
  expect_error(compare_spatial_models(models, data.frame()), "different fixed effects")
  expect_error(compare_spatial_models(list(list(fixed = y ~ 1)), data.frame()), "named list")
})

test_that("extract_variance_components validates input", {
  expect_error(extract_variance_components(list()), "must be an asreml object")
})

test_that("failed fits are reported, not hidden", {
  skip_if_no_asreml()
  d <- data.frame(Yield = rnorm(10), Genotype = factor(rep(c("A", "B"), 5)))
  res <- compare_spatial_models(list(m1 = list(fixed = Yield ~ 1, random = ~Genotype, residual = ~ ar1(nope):ar1(nada))), d)
  expect_equal(res$Model, "m1")
  expect_true(is.infinite(res$AIC))
  expect_false(res$Converged)
  expect_false(is.na(res$Error))
})

test_that("extract_variance_components matches the varcomp table", {
  skip_if_no_asreml()
  s <- sim_fa_data(n_gen = 40, n_sites = 1)
  d <- s$data
  m <- suppressWarnings(asreml::asreml(Yield ~ 1, random = ~Genotype, residual = ~ ar1(Column):ar1(Row), data = d, trace = FALSE))
  vc <- summary(m)$varcomp
  out <- extract_variance_components(m, "Genotype", n_reps = 2)
  expect_equal(out$Vg, vc["Genotype", "component"])
  expect_equal(out$Ve, vc["Column:Row!R", "component"])
  expect_equal(out$Vr, 0)
  expect_equal(out$h2, out$Vg / (out$Vg + out$Ve / 2))
  p <- predict(m, classify = "Genotype", only = "Genotype", vcov = TRUE)
  expect_equal(out$vBLUP, UtilityFunctions:::mean_pairwise_vd(as.matrix(p$vcov)))
  expect_true(out$Gen_Reliability > 0 && out$Gen_Reliability < 1)
  expect_equal(out$Acc, sqrt(out$Gen_Reliability))
})

test_that("user-named matrices resolve from the caller's environment", {
  skip_if_no_asreml()
  s <- sim_fa_data(n_gen = 30, n_sites = 1)
  run_in_function <- function(dat, rel) {
    my_G <- rel # only visible in this function's frame
    compare_spatial_models(
      list(id = list(fixed = Yield ~ 1, random = ~ vm(Genotype, my_G), residual = ~units)),
      dat
    )
  }
  res <- run_in_function(s$data, s$G)
  expect_true(is.na(res$Error))
  expect_true(is.finite(res$AIC))
})

test_that("fit_bivariate_models returns numeric correlations", {
  skip_if_no_asreml()
  s <- sim_fa_data(n_gen = 40, n_sites = 1)
  d1 <- s$data
  long <- rbind(
    data.frame(Trait = "Yield", Genotype = d1$Genotype, Value_Scaled = as.numeric(scale(d1$Yield))),
    data.frame(Trait = "NDVI", Genotype = d1$Genotype, Value_Scaled = as.numeric(scale(d1$Yield + rnorm(nrow(d1), 0, 3))))
  )
  out <- suppressWarnings(suppressMessages(fit_bivariate_models(
    long, "Yield", "NDVI",
    rand_f = ~ us(Trait):vm(Genotype, grm), res_f = ~ dsum(~ units | Trait),
    genotype_term = "Trait:vm(Genotype, grm)", spatial_term = NULL, grm = s$G
  )))
  expect_true(is.na(out$Error))
  expect_type(out$rG, "double")
  expect_true(out$rG > 0 && out$rG <= 1)
  expect_true(out$h2 > 0 && out$h2 < 1)
  expect_false(is.na(out$vBLUP))
})

test_that("generalised reliability does not depend on ASReml's parameterisation", {
  skip_if_no_asreml()
  s <- sim_fa_data(n_gen = 40, n_sites = 1)
  d <- s$data
  ma <- suppressWarnings(asreml::asreml(Yield ~ 1, random = ~Genotype, data = d, trace = FALSE))
  mb <- suppressWarnings(asreml::asreml(Yield ~ 1, random = ~Genotype, residual = ~ idv(units), data = d, trace = FALSE))
  va <- extract_variance_components(ma, "Genotype", n_reps = 2)
  vb <- extract_variance_components(mb, "Genotype", n_reps = 2) # units!R = 1, variance in units!units
  expect_equal(va$Ve, vb$Ve, tolerance = 1e-4)
  expect_equal(va$Vg, vb$Vg, tolerance = 1e-4)
  expect_equal(va$Gen_Reliability, vb$Gen_Reliability, tolerance = 1e-4)
  expect_equal(va$h2, vb$h2, tolerance = 1e-4)
})
