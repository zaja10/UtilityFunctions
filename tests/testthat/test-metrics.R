test_that("predict_gain calculates correctly with mock fa_model", {
    # Mock fa_model object
    # We need object$matrices$G
    # G matrix diagonals are variances

    # Example: 2 sites
    # Site 1 Vg = 100 -> sigma_g = 10
    # Site 2 Vg = 25  -> sigma_g = 5

    G <- matrix(0, nrow = 2, ncol = 2)
    diag(G) <- c(100, 25)
    rownames(G) <- colnames(G) <- c("SiteA", "SiteB")

    mock_model <- list(
        matrices = list(G = G)
    )
    class(mock_model) <- "fa_model"

    # Calculation:
    # i for 5% selection:
    p <- 0.05
    i_val <- dnorm(qnorm(1 - p)) / p
    # qnorm(0.95) ~= 1.645
    # dnorm(1.645) ~= 0.103
    # i ~= 2.063

    # Site 1: i * 10 = 20.63
    # Site 2: i * 5 = 10.31

    gain_df <- predict_gain(mock_model, selection_percent = 5)

    expect_s3_class(gain_df, "data.frame")
    expect_equal(nrow(gain_df), 2)
    expect_equal(gain_df$Sigma_G, c(10, 5))

    # Check i value consistent with R's calculation
    expect_equal(gain_df$Intensity_i[1], round(i_val, 2))

    # Check Potential Gain
    expect_equal(gain_df$Potential_Gain[1], round(i_val * 10, 3))
    expect_equal(gain_df$Potential_Gain[2], round(i_val * 5, 3))
})

test_that("predict_gain validates input", {
    mock_model <- list(matrices = list(G = matrix(1)))
    class(mock_model) <- "fa_model"

    expect_error(predict_gain(mock_model, selection_percent = 0))
    expect_error(predict_gain(mock_model, selection_percent = 100))
})
