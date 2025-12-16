library(testthat)
library(UtilityFunctions)

test_that("pad_trial_layout fills missing coordinates", {
    # Create sparse data
    df <- data.frame(
        Row = c(1, 1, 2),
        Column = c(1, 2, 2),
        Yield = c(10, 20, 30)
    )
    # Grid should be 2x2. Missing (2,1).
    padded <- pad_trial_layout(df, "Row", "Column")

    expect_equal(nrow(padded), 4)
    expect_true(is.na(padded$Yield[padded$Row == 2 & padded$Column == 1]))
    expect_equal(padded$Yield[padded$Row == 1 & padded$Column == 2], 20)
})

test_that("pad_trial_layout handles MET groups", {
    df_gap <- data.frame(
        Trial = c("T1", "T1", "T2"),
        Row = c(1, 3, 1),
        Column = c(1, 1, 1),
        Yield = c(10, 20, 30)
    )

    padded <- pad_trial_layout(df_gap, group = "Trial")
    expect_equal(nrow(padded), 4)

    # Check T1 missing row (2,1)
    t1_pad <- padded[padded$Trial == "T1" & padded$Row == 2, ]
    expect_equal(nrow(t1_pad), 1)
    expect_true(is.na(t1_pad$Yield))
})

test_that("calculate_connectivity returns correct counts", {
    df <- data.frame(
        Year = c("Y1", "Y1", "Y2", "Y2"),
        Genotype = c("G1", "G2", "G2", "G3")
    )

    mat <- calculate_connectivity(df, "Year", "Year", "Genotype")

    expect_equal(mat["Y1", "Y1"], 2)
    expect_equal(mat["Y1", "Y2"], 1)

    jacc <- calculate_connectivity(df, "Year", "Year", "Genotype", method = "jaccard")
    expect_equal(jacc["Y1", "Y2"], 1 / 3, tolerance = 1e-4)
})

test_that("convert_buac_to_tha calculates correctly", {
    y_wheat <- convert_buac_to_tha(100, crop = "wheat")
    expect_equal(y_wheat, 6.725, tolerance = 1e-2)
})

test_that("plot functions run without error", {
    set.seed(123)
    df <- data.frame(
        Year = rep(2020:2022, each = 20),
        Genotype = sample(LETTERS[1:5], 60, replace = TRUE),
        Yield = rnorm(60, 5, 2),
        Row = rep(1:5, 12),
        Column = rep(1:4, 15),
        Trial = rep(c("T1", "T2", "T3"), each = 20)
    )

    pdf(NULL)
    expect_silent(plot_connectivity(df, x = "Year", trace = "Genotype"))
    expect_silent(plot_met_trend(df, x = "Year", y = "Yield"))
    expect_silent(plot_trial_map(df, trial_val = "T1"))
    dev.off()
})
