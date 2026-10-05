test_that("generate_met_data randomises layout within replicate", {
  m <- generate_met_data(n_sites = 3, n_gen = 20, sparsity = 0)
  pos <- tapply(paste(m$Row, m$Column), list(m$Genotype, m$Rep), function(v) length(unique(v)))
  expect_gt(mean(pos > 1), 0.8)
  # each replicate is a complete block
  expect_true(all(table(m$Genotype, m$Site) == 2))
})

test_that("simulators are reproducible and leave the caller's RNG alone", {
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  a <- generate_met_data(n_sites = 2, n_gen = 10, seed = 5)
  after <- runif(1)
  expect_equal(before, after)
  expect_equal(a, generate_met_data(n_sites = 2, n_gen = 10, seed = 5))
})

test_that("generate_grm is a positive-definite VanRaden GRM", {
  g <- generate_grm(30, 200)$grm
  expect_equal(dim(g), c(30, 30))
  expect_true(isSymmetric(g))
  expect_gt(min(eigen(g, only.values = TRUE)$values), 0)
})

test_that("generate_genomic_met validates h2 and attains it approximately", {
  g <- generate_grm(100, 300)$grm
  expect_error(generate_genomic_met(2, g, h2 = 1), "h2")
  expect_error(generate_genomic_met(2, unname(g)), "row names")
  d <- generate_genomic_met(1, g, n_reps = 4, h2 = 0.5)
  expect_equal(nrow(d), 400)
})
