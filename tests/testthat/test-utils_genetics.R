test_that("match_germplasm matches names and synonyms", {
  valid_ids <- c("LineA", "LineB", "LineC")
  expect_equal(match_germplasm("LineA", NA, valid_ids), "LineA")
  expect_equal(match_germplasm("Unknown", "Alias1, LineB", valid_ids), "LineB")
  expect_true(is.na(match_germplasm("Unknown", "Alias1, Alias2", valid_ids)))
})

test_that("match_germplasm is vectorised", {
  out <- match_germplasm(c("LineA", "x", "y", "z"), c(NA, "LineC", "", NA), c("LineA", "LineC"))
  expect_equal(out, c("LineA", "LineC", NA, NA))
})

test_that("prepare_trial_grm pads and scales correctly", {
  master_grm <- matrix(c(2, 1, 1, 2), nrow = 2, dimnames = list(c("LineA", "LineB"), c("LineA", "LineB")))
  trial_data <- data.frame(Genotype = factor(c("LineA", "LineC")))
  trial_grm <- prepare_trial_grm(master_grm, trial_data, "Genotype")

  expect_equal(rownames(trial_grm), c("LineA", "LineC"))
  expect_equal(trial_grm["LineA", "LineA"], 1 + 1e-4) # rescaled to mean diag 1
  expect_equal(trial_grm["LineC", "LineC"], 1 + 1e-4)
  expect_equal(trial_grm["LineA", "LineC"], 0)
  expect_error(prepare_trial_grm(master_grm, data.frame(Genotype = "Q"), "Genotype"), "No trial lines")
})
