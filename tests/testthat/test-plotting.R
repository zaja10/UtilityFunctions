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
    p3 <- plot(mock_model, type = "biplot")
    expect_s3_class(p3, "ggplot")

    # Test VAF
    p4 <- plot(mock_model, type = "vaf")
    expect_s3_class(p4, "ggplot")
})

test_that("plot.padded_trial returns a ggplot object", {
    df <- expand.grid(Row = 1:5, Col = 1:5)
    df$Yield <- rnorm(25)
    attr(df, "coords") <- c(row = "Row", col = "Col")
    attr(df, "is_met") <- FALSE
    class(df) <- c("padded_trial", "data.frame")

    p <- plot(df)
    expect_s3_class(p, "ggplot")
})
