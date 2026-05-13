test_that("compare_spatial_models handles basic errors", {
  # Without asreml installed, or without real data, it should return a data.frame with AIC=Inf
  # We'll just test the error handling pathway.
  mock_data <- data.frame(Yield = rnorm(10), Genotype = factor(rep(c("A","B"), 5)))
  
  mock_models <- list(
    m1 = list(fixed = Yield ~ 1, random = ~ Genotype, residual = ~ 1)
  )
  
  # This should run through the tryCatch error block because we don't have a real asreml license/setup here or asreml will fail to converge
  res <- compare_spatial_models(mock_models, mock_data)
  expect_true(is.data.frame(res))
  expect_equal(res$Model, "m1")
  expect_true(is.infinite(res$AIC))
  expect_false(res$Converged)
})

test_that("extract_variance_components handles invalid input", {
  expect_error(extract_variance_components(list()), "must be an asreml object")
})
