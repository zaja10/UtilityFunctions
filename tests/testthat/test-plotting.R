test_that("plot.fa_model returns a ggplot object", {
    # Mock an fa_model object structure
    # We don't need a real fit, just the structure required by the plot function logic

    # Mock scores and loadings for plotting
    scores <- matrix(rnorm(20), nrow = 10, dimnames = list(paste0("G", 1:10), c("Factor1", "Factor2")))
    loadings <- matrix(rnorm(10), nrow = 5, dimnames = list(paste0("E", 1:5), c("Factor1", "Factor2")))

    mock_model <- list(
        scores = list(rotated = scores),
        loadings = list(rotated = loadings),
        matrices = list(Cor = cor(t(loadings))), # Mock correlation matrix
        meta = list(k = 2),
        var_comp = list(vaf = data.frame(Site = paste0("E", 1:5), VAF_Fac1 = runif(5, 20, 50), VAF_Fac2 = runif(5, 20, 40))),
        fast = data.frame(Genotype = paste0("G", 1:10), OP = rnorm(10), RMSD = runif(10))
    )
    class(mock_model) <- "fa_model"

    # Test FAST plot
    p1 <- plot(mock_model, type = "fast")
    expect_s3_class(p1, "ggplot")

    # Test Heatmap
    p2 <- plot(mock_model, type = "heatmap")
    expect_s3_class(p2, "ggplot")

    # Test Biplot
    expect_error(plot(mock_model, type = "biplot"))

    # Test VAF
    expect_error(plot(mock_model, type = "vaf"))
})


