mock_fa <- function(vg) {
  G <- diag(vg, nrow = length(vg))
  rownames(G) <- colnames(G) <- paste0("Site", seq_along(vg))
  structure(list(matrices = list(G = G)), class = "fa_model")
}

test_that("predicted gain applies the breeder's equation with accuracy", {
  i_val <- dnorm(qnorm(0.95)) / 0.05
  g <- calculate_genetic_gain(mock_fa(c(100, 25)), selection_percent = 5, accuracy = 0.5, cycle_length = 2)
  expect_equal(g$Sigma_G, c(10, 5))
  expect_equal(g$Intensity[1], i_val)
  expect_equal(g$Predicted_Gain, i_val * 0.5 * c(10, 5) / 2)
})

test_that("predicted gain warns that r = 1 is an upper bound", {
  expect_warning(g <- calculate_genetic_gain(mock_fa(100)), "upper bound")
  expect_equal(g$Accuracy, 1)
})

test_that("predicted gain validates input", {
  expect_error(calculate_genetic_gain(mock_fa(1), selection_percent = 0), "selection_percent")
  expect_error(calculate_genetic_gain(mock_fa(1), selection_percent = 100), "selection_percent")
  expect_error(calculate_genetic_gain(mock_fa(1), accuracy = 1.2), "accuracy")
})

test_that("realized gain accepts factor years and adjusts for check trend", {
  set.seed(2)
  yrs <- c(2015, 2017, 2023)
  d <- expand.grid(Year = yrs, Genotype = c(paste0("L", 1:20), "CHK1", "CHK2"))
  env <- 0.3 * (d$Year - 2015)
  gen <- ifelse(grepl("CHK", d$Genotype), 0, 0.5 * (d$Year - 2015))
  d$Yield <- 50 + env + gen + rnorm(nrow(d), 0, 0.01)
  d$Year <- factor(d$Year)
  raw <- calculate_genetic_gain(d)
  expect_equal(raw$Slope, (0.5 * 20 / 22 + 0.3), tolerance = 0.01)
  adj <- calculate_genetic_gain(d, check_list = c("CHK1", "CHK2"))
  expect_equal(adj$Slope, 0.5, tolerance = 0.01)
  expect_equal(adj$Env_Slope, 0.3, tolerance = 0.01)
})

vp <- list(
  `studyName:gdrop` = c(`studyName:gdrop!studyName_TrialB` = 0.3, `studyName:gdrop!studyName_TrialA` = 0.2),
  `studyName:vm(gkeep, G.inv)` = c(`studyName:vm(gkeep, G.inv)!studyName_TrialA` = 0.5, `studyName:vm(gkeep, G.inv)!studyName_TrialB` = 0.4),
  `studyName_TrialA!R` = 0.1,
  `studyName_TrialB!R` = 0.15,
  `at(studyName, 'TrialA'):units` = 0.05,
  `at(studyName, 'TrialB'):units` = 0.08
)
h2_fun <- function(...) UtilityFunctions:::met_h2_from_vparameters(vp, "studyName:gdrop", "studyName:vm(gkeep, G.inv)", "studyName", ...)

test_that("MET heritability matches components by trial name, not position", {
  res <- h2_fun()
  expect_equal(res$Trial, c("TrialA", "TrialB"))
  expect_equal(res$Vg_add, c(0.5, 0.4))
  expect_equal(res$Vg_na, c(0.2, 0.3)) # gdrop listed in the opposite order
  expect_equal(res$V_spatial, c(0.1, 0.15))
  expect_equal(res$V_units, c(0.05, 0.08))
})

test_that("disjoint line sets never sum genetic variances", {
  res <- h2_fun(disjoint = TRUE)
  expect_equal(res$h2, c(0.5 / 0.65, 0.4 / 0.63))
  expect_equal(res$H2_ungenotyped, c(0.2 / 0.35, 0.3 / 0.53))
  res2 <- h2_fun(disjoint = FALSE)
  expect_equal(res2$H2, c(0.7 / 0.85, 0.7 / 0.93))
  expect_equal(res2$h2, c(0.5 / 0.85, 0.4 / 0.93))
})

test_that("calculate_met_h2 validates input", {
  expect_error(calculate_met_h2("not a model"), "Input must be a valid asreml model.")
})
