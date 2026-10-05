mock_fa_model <- function() {
  set.seed(1)
  scores <- matrix(rnorm(20), nrow = 10, dimnames = list(paste0("G", 1:10), c("Fac1", "Fac2")))
  loadings <- matrix(rnorm(10), nrow = 5, dimnames = list(paste0("E", 1:5), c("Fac1", "Fac2")))
  G <- loadings %*% t(loadings) + diag(0.5, 5)
  structure(list(
    scores = list(rotated = scores),
    loadings = list(rotated = loadings),
    matrices = list(G = G, Cor = cov2cor(G)),
    meta = list(k = 2, rotation = "svd", var_explained = c(60, 20), group = "Site"),
    var_comp = list(vaf = data.frame(Group = paste0("E", 1:5), VAF_Fac1 = runif(5, 20, 50), VAF_Fac2 = runif(5, 10, 30), Total_VAF = 70)),
    fast = data.frame(Genotype = paste0("G", 1:10), OP = rnorm(10), RMSD = runif(10))
  ), class = "fa_model")
}

test_that("plot.fa_model returns ggplots for every ggplot type", {
  x <- mock_fa_model()
  for (type in c("fast", "heatmap", "vaf", "biplot", "latent_reg")) {
    expect_s3_class(plot(x, type = type, highlight = "G3"), "ggplot")
  }
  expect_error(plot(x, type = "latent_reg", factor = 3), "exceeds")
})

test_that("cluster plot uses a Euclidean distance and returns the tree", {
  pdf(NULL)
  on.exit(dev.off())
  hc <- plot(mock_fa_model(), type = "cluster")
  expect_s3_class(hc, "hclust")
  expect_true(all(hc$height >= 0))
})

test_that("plot_trend treats years as numeric", {
  d <- data.frame(Year = rep(c(2015, 2017, 2023), each = 5), Yield = rnorm(15), Genotype = "a")
  p <- plot_trend(d)
  expect_s3_class(p, "ggplot")
  expect_true(is.numeric(p$data$Year))
  expect_s3_class(plot_trend(d, mode = "genetic", checks = "a"), "ggplot")
})

test_that("plot_spatial validates columns", {
  d <- data.frame(Row = 1:2, Column = 1:2, Yield = 1:2)
  expect_s3_class(plot_spatial(d), "ggplot")
  expect_error(plot_spatial(d, fill = "nope"), "not found")
})
