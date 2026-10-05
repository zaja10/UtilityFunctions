test_that("calculate_vegetation_indices validates input", {
  expect_error(calculate_vegetation_indices("not a dataframe"), "spectra_df must be a data frame")
  expect_warning(out <- calculate_vegetation_indices(data.frame(Z = 1)), "No indices")
  expect_equal(out, data.frame(Z = 1))
})

test_that("indices match the Awesome Spectral Indices formulas", {
  b <- data.frame(N = 0.8, R = 0.1, G = 0.2, B = 0.05)
  res <- suppressMessages(calculate_vegetation_indices(b))
  expect_equal(res$NDVI, (0.8 - 0.1) / (0.8 + 0.1))
  expect_equal(res$GNDVI, (0.8 - 0.2) / (0.8 + 0.2))
  expect_equal(res$EVI, 2.5 * (0.8 - 0.1) / (0.8 + 6 * 0.1 - 7.5 * 0.05 + 1))
  # GARI (Gitelson 1996), lmb = 1: (N - (G - (B - R))) / (N + (G - (B - R)))
  expect_equal(res$GARI, (0.8 - (0.2 - (0.05 - 0.1))) / (0.8 + (0.2 - (0.05 - 0.1))))
  # OCVI uses cexp = 1.16 by default
  expect_equal(res$OCVI, (0.8 / 0.2) * (0.1 / 0.2)^1.16)
})

test_that("unrelated data columns cannot override index constants", {
  b <- data.frame(N = 0.8, R = 0.1, B = 0.05)
  a <- suppressMessages(calculate_vegetation_indices(b))$EVI
  b$L <- 100
  b$g <- -1
  expect_equal(suppressMessages(calculate_vegetation_indices(b))$EVI, a)
})

test_that("process_longitudinal_indices stacks timepoints", {
  df <- data.frame(
    PlotID = c("P1", "P2"),
    Blue.1 = c(0.1, 0.11), Green.1 = c(0.2, 0.22), Red.1 = c(0.1, 0.12), RedEdge.1 = c(0.3, 0.33), NIR.1 = c(0.8, 0.82),
    Blue.2 = c(0.15, 0.16), Green.2 = c(0.25, 0.26), Red.2 = c(0.15, 0.18), RedEdge.2 = c(0.35, 0.38), NIR.2 = c(0.85, 0.88)
  )
  res <- process_longitudinal_indices(df, "PlotID", names(df)[-1])
  expect_equal(nrow(res), 4)
  expect_equal(sort(unique(res$Day)), c(1, 2))
  expect_equal(res$NDVI[res$PlotID == "P1" & res$Day == 2], 0.7)
})

test_that("process_longitudinal_indices warns about and skips incomplete timepoints", {
  df <- data.frame(id = 1, Blue.1 = 0.1, Blue.2 = 0.15)
  expect_warning(res <- process_longitudinal_indices(df, "id", c("Blue.1", "Blue.2")), "Skipped 2")
  expect_null(res)
  expect_error(process_longitudinal_indices(df, "id", "Nope.1"), "not found")
})

test_that("evaluate_vi_loco reports pooled and within-column accuracy", {
  set.seed(3)
  n_col <- 8
  d <- data.frame(colNumber = rep(seq_len(n_col), each = 25))
  d$V1 <- rnorm(nrow(d))
  d$V2 <- rnorm(nrow(d))
  col_eff <- rnorm(n_col, 0, 3)[d$colNumber]
  d$Yield <- 2 * d$V1 - d$V2 + rnorm(nrow(d), 0, 0.5) + col_eff
  r <- suppressMessages(evaluate_vi_loco(d, "Yield", c("V1", "V2")))
  expect_named(r, c("Acc", "Acc_within", "Importance", "Col_Acc", "Predictions"))
  expect_gt(r$Acc_within, r$Acc) # column effects dilute pooled accuracy here
  expect_equal(r$Importance$Predictor[1], "V1")
  expect_true(all(r$Importance$VIF < 2))
  expect_equal(nrow(r$Col_Acc), n_col)
})

test_that("evaluate_vi_loco warns on collinear predictors and skips rank-deficient folds", {
  set.seed(4)
  d <- data.frame(colNumber = rep(1:4, each = 10), V1 = rnorm(40))
  d$V2 <- d$V1 + rnorm(40, 0, 0.01)
  d$Yield <- d$V1 + rnorm(40)
  expect_warning(suppressMessages(evaluate_vi_loco(d, "Yield", c("V1", "V2"))), "collinear")
})
