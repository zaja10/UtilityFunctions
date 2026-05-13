test_that("match_germplasm matches correctly", {
  valid_ids <- c("LineA", "LineB", "LineC")
  
  # Primary match
  expect_equal(match_germplasm("LineA", NA, valid_ids), "LineA")
  
  # Synonym match
  expect_equal(match_germplasm("Unknown", "Alias1, LineB", valid_ids), "LineB")
  
  # No match
  expect_true(is.na(match_germplasm("Unknown", "Alias1, Alias2", valid_ids)))
})

test_that("prepare_trial_grm pads and scales correctly", {
  master_grm <- matrix(c(1, 0.5, 0.5, 1), nrow=2, dimnames=list(c("LineA", "LineB"), c("LineA", "LineB")))
  trial_data <- data.frame(Genotype = factor(c("LineA", "LineC"))) # LineC is ungenotyped
  
  trial_grm <- prepare_trial_grm(master_grm, trial_data, "Genotype")
  
  # Should be 2x2: LineA and LineC
  expect_equal(nrow(trial_grm), 2)
  expect_equal(rownames(trial_grm), c("LineA", "LineC"))
  
  # LineC should be padded with 1 + ridge on diagonal
  expect_true(trial_grm["LineC", "LineC"] > 1)
  expect_equal(trial_grm["LineA", "LineC"], 0)
})
