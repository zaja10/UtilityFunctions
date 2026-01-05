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

    padded <- pad_trial_layout(df_gap, group_cols = "Trial")
    expect_equal(nrow(padded), 4)

    # Check T1 missing row (2,1)
    t1_pad <- padded[padded$Trial == "T1" & padded$Row == 2, ]
    expect_equal(nrow(t1_pad), 1)
    expect_true(is.na(t1_pad$Yield))
})

# test_that("convert_buac_to_tha calculates correctly", { ... }) DELETED

test_that("plot functions run without error", {
    set.seed(123)
    df <- data.frame(
        Year = rep(c("2020", "2021"), each = 10),
        Yield = rnorm(20),
        Trial = rep(c("T1", "T2"), each = 10),
        Row = rep(1:2, 10),
        Column = rep(1:5, 4),
        Long = rnorm(20),
        Lat = rnorm(20),
        Genotype = rep(c("A", "B"), 10)
    )


    p2 <- plot_trend(df, mode = "phenotypic", x = "Year", y = "Yield")
    expect_s3_class(p2, "ggplot")

    p3 <- plot_spatial(df[df$Trial == "T1", ])
    expect_s3_class(p3, "ggplot")
})
