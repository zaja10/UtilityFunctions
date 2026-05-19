test_that("calculate_genetic_gain calculates correctly with mock fa_model", {
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

    gain_df <- calculate_genetic_gain(mock_model, selection_percent = 5)

    expect_s3_class(gain_df, "data.frame")
    expect_equal(nrow(gain_df), 2)
    # The new function returns standard deviations in Sigma_G column
    expect_equal(gain_df$Sigma_G, c(10, 5))

    # Check i value consistent with R's calculation
    expect_equal(gain_df$Intensity[1], i_val, tolerance = 1e-4)

    # Check Potential Gain
    expect_equal(gain_df$Predicted_Gain[1], i_val * 10, tolerance = 1e-3)
    expect_equal(gain_df$Predicted_Gain[2], i_val * 5, tolerance = 1e-3)
})

test_that("predict_gain validates input", {
    mock_model <- list(matrices = list(G = matrix(1)))
    class(mock_model) <- "fa_model"

    expect_error(calculate_genetic_gain(mock_model, selection_percent = 0))
    expect_error(calculate_genetic_gain(mock_model, selection_percent = 100))
})

test_that("calculate_met_h2 validates input", {
    expect_error(calculate_met_h2("not a model"), "Input must be a valid asreml model.")
})

test_that("calculate_met_h2 works with mock asreml model", {
    skip_if_not_installed("asreml")

    mock_model <- list()
    class(mock_model) <- "asreml"

    mock_summary <- function(object, ...) {
        list(vparameters = list(
            `studyName:gdrop` = c(
                `studyName:gdrop!studyName_TrialA` = 0.2,
                `studyName:gdrop!studyName_TrialB` = 0.3
            ),
            `studyName:vm(gkeep, G.inv)` = c(
                `studyName:vm(gkeep, G.inv)!studyName_TrialA` = 0.5,
                `studyName:vm(gkeep, G.inv)!studyName_TrialB` = 0.4
            ),
            `studyName_TrialA!R` = 0.1,
            `studyName_TrialB!R` = 0.15,
            `at(studyName, 'TrialA'):units` = 0.05,
            `at(studyName, 'TrialB'):units` = 0.08
        ))
    }

    with_mocked_bindings(
        summary.asreml = mock_summary,
        .package = "asreml",
        code = {
            res <- calculate_met_h2(mock_model)
            expect_s3_class(res, "data.frame")
            expect_equal(nrow(res), 2)
            expect_equal(res$Trial, c("TrialA", "TrialB"))
            expect_equal(res$Vg_add, c(0.5, 0.4))
            expect_equal(res$Vg_na, c(0.2, 0.3))
            expect_equal(res$Vg_total, c(0.7, 0.7))
            expect_equal(res$V_spatial, c(0.1, 0.15))
            expect_equal(res$V_units, c(0.05, 0.08))
            expect_equal(res$H2[1], 0.7 / 0.85)
            expect_equal(res$H2[2], 0.7 / 0.93)
            expect_equal(res$h2[1], 0.5 / 0.65)
            expect_equal(res$h2[2], 0.4 / 0.63)
        }
    )
})

