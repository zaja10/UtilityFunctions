fit_fa_model <- function() {
  s <- sim_fa_data(n_gen = 60, n_sites = 6, seed = 11)
  d <- s$data
  G <- s$G
  m <- suppressWarnings(asreml::asreml(
    Yield ~ Site,
    random = ~ fa(Site, 2):vm(Genotype, G),
    residual = ~ dsum(~ units | Site),
    data = d, trace = FALSE, maxit = 40
  ))
  list(m = m, d = d, G = G)
}

test_that("fa.asreml reconstructs G and rotates correctly", {
  skip_if_no_asreml()
  f <- fit_fa_model()
  fa <- suppressMessages(fa.asreml(f$m, grm = f$G))
  vc <- summary(f$m)$varcomp
  term <- "fa(Site, 2):vm(Genotype, G)"
  sites <- levels(f$d$Site)

  expect_equal(fa$meta$term, term)
  expect_equal(fa$meta$k, 2L)
  lam <- cbind(vc[paste0(term, "!", sites, "!fa1"), "component"], vc[paste0(term, "!", sites, "!fa2"), "component"])
  psi <- vc[paste0(term, "!", sites, "!var"), "component"]
  expect_equal(unname(fa$loadings$raw), lam)
  expect_equal(unname(fa$matrices$G), lam %*% t(lam) + diag(psi))
  # rotation is orthogonal and orients factor 1 positively
  expect_equal(fa$loadings$rotated %*% t(fa$loadings$rotated), fa$loadings$raw %*% t(fa$loadings$raw), ignore_attr = TRUE)
  expect_gt(mean(fa$loadings$rotated[, 1]), 0)
  expect_equal(sum(fa$meta$var_explained), 100 * sum(lam^2) / sum(diag(fa$matrices$G)))

  # scores: raw scores reproduce the common GxE effects after rotation
  expect_equal(fa$scores$rotated %*% t(fa$loadings$rotated), fa$scores$raw %*% t(fa$loadings$raw), ignore_attr = TRUE)
  expect_equal(nrow(fa$scores$raw), nlevels(f$d$Genotype))

  # reliability r2_g uses g_ii; no clipping should be needed for a sound fit
  ids <- fa$fast$Genotype
  manual <- 1 - fa$scores$pev_rotated[ids, 1] / diag(f$G)[ids]
  expect_true(all(manual > 0 & manual < 1))
  expect_equal(fa$fast$Reliability, unname(manual))
})

test_that("FA scores, OP and reliability do not depend on ASReml's parameterisation", {
  skip_if_no_asreml()
  s <- sim_fa_data(n_gen = 40, n_sites = 4, seed = 3)
  d <- s$data
  G <- s$G
  # Same model: implicit residual (sigma parameterisation, sigma2 != 1) vs explicit idv(units)
  ma <- suppressWarnings(asreml::asreml(Yield ~ Site, random = ~ fa(Site, 2):vm(Genotype, G), data = d, trace = FALSE, maxit = 50))
  mb <- suppressWarnings(asreml::asreml(Yield ~ Site, random = ~ fa(Site, 2):vm(Genotype, G), residual = ~ idv(units), data = d, trace = FALSE, maxit = 50))
  expect_gt(ma$sigma2, 10)
  expect_equal(mb$sigma2, 1)
  fa_a <- suppressMessages(fa.asreml(ma, grm = G))
  fa_b <- suppressMessages(fa.asreml(mb, grm = G))
  expect_equal(fa_a$matrices$G, fa_b$matrices$G, tolerance = 5e-3)
  fb <- fa_b$fast[match(fa_a$fast$Genotype, fa_b$fast$Genotype), ]
  expect_equal(fa_a$fast$OP, fb$OP, tolerance = 5e-3)
  expect_equal(fa_a$fast$RMSD, fb$RMSD, tolerance = 5e-3)
  expect_equal(fa_a$fast$Reliability, fb$Reliability, tolerance = 5e-3)
  expect_true(all(fa_a$fast$Reliability > 0 & fa_a$fast$Reliability < 1))
})

test_that("fa.asreml options: no rotation, iClasses, summaries", {
  skip_if_no_asreml()
  f <- fit_fa_model()
  none <- suppressMessages(fa.asreml(f$m, rotate = FALSE))
  expect_null(none$fast)
  expect_error(calculate_i_classes(none), "SVD")

  fa <- suppressMessages(fa.asreml(f$m))
  ic <- calculate_i_classes(fa)
  expect_setequal(unique(ic$environments$iClass), unique(ifelse(fa$loadings$rotated[, 2] >= 0, "p", "n")))
  cl <- ic$environments$iClass[1]
  in_c <- ic$environments$iClass == cl
  lbar <- colMeans(fa$loadings$rotated[in_c, 1:2, drop = FALSE])
  g1 <- ic$genotypes[ic$genotypes$iClass == cl, ]
  expect_equal(g1$iClassOP, as.vector(fa$scores$rotated[g1$Genotype, 1:2] %*% lbar))

  s <- summary(fa)
  expect_s3_class(s, "summary.fa_model")
  expect_equal(nrow(s$top_genotypes), 5)
  expect_s3_class(plot(fa, type = "biplot"), "ggplot")
})

test_that("evaluate_gs_predictions builds on fa.asreml and writes nothing by default", {
  skip_if_no_asreml()
  f <- fit_fa_model()
  withr::local_dir(withr::local_tempdir())
  ids <- levels(f$d$Genotype)
  gs <- suppressMessages(evaluate_gs_predictions(
    f$m,
    target_ids = ids[1:20], training_ids = ids[21:60],
    target_term = "fa(Site, 2):vm(Genotype, G)", grm = f$G
  ))
  expect_s3_class(gs, "gs_predictions")
  expect_length(list.files(), 0)
  expect_equal(sort(gs$blups$id), sort(ids))
  fast <- gs$fa$fast
  expect_equal(gs$blups$Reliability[match(fast$Genotype, gs$blups$id)], fast$Reliability)
  expect_equal(gs$blups$Selection_Index, gs$blups$GEBV_OP - gs$blups$Stability_RMSD)
  expect_true("Max_Rel" %in% names(gs$blups))
  expect_s3_class(gs$plots, "patchwork")
  expect_s3_class(summary(gs), "summary.gs_predictions")
  expect_warning(
    suppressMessages(evaluate_gs_predictions(f$m, ids[1:20], ids[15:60], "fa(Site, 2):vm(Genotype, G)")),
    "both training and target"
  )
})
