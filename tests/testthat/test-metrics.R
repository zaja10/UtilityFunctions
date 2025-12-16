test_that("predict_gain calculates correctly", {
    # R = i * r * sigma_g
    # i = 1.755 (10%)
    # r = 0.8 (accuracy) => reliability = 0.64
    # sigma_g = 10

    # Using reliability input
    gain <- predict_gain(selection_intensity = 1.755, reliability = 0.64, genetic_variance = 100)
    # sigma_g = sqrt(100) = 10
    # acc = sqrt(0.64) = 0.8
    # R = 1.755 * 0.8 * 10 = 14.04

    expect_equal(gain, 14.04, tolerance = 1e-3)
})

test_that("predict_gain handles vector inputs", {
    gain <- predict_gain(1, c(0.25, 1.0), c(100, 100))
    # 1: 1 * 0.5 * 10 = 5
    # 2: 1 * 1.0 * 10 = 10
    expect_equal(gain, c(5, 10))
})
