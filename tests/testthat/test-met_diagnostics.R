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
    df <- data.frame(
        Trial = c("T1", "T1", "T2"),
        Row = c(1, 2, 1),
        Column = c(1, 1, 1), # T1 is 2x1 grid, T2 is 1x1
        Yield = c(10, 20, 30)
    )
    # Modify T1 to have missing spot: Row 1,2 exist. Suppose max is 3?
    # Let's make T1 miss (2,1).
    df <- data.frame(
        Trial = c("T1", "T1", "T2"),
        Row = c(1, 3, 1),
        Column = c(1, 1, 1),
        Yield = c(10, 20, 30)
    ) # T1 ranges 1-3. Missing 2.

    padded <- pad_trial_layout(df, group_cols = "Trial")

    # T1 should have 3 rows now (1, 2(NA), 3). T2 should have 1 row. Total 4.
    expect_equal(nrow(padded), 4)

    # Check T1 missing row
    t1_pad <- padded[padded$Trial == "T1" & padded$Row == 2, ]
    expect_equal(nrow(t1_pad), 1)
    expect_true(is.na(t1_pad$Yield))
})


test_that("check_met_connectivity detects islands", {
    # Connected set: T1-T2 share G2
    # Disconnected: T3 has G3 only
    df <- data.frame(
        Trial = c("T1", "T1", "T2", "T3"),
        Genotype = c("G1", "G2", "G2", "G3")
    )

    res <- check_met_connectivity(df, trial_col = "Trial", gen_col = "Genotype")

    expect_equal(res$n_clusters, 2)
    expect_true("T3" %in% res$islands)
    expect_false("T1" %in% res$islands)
})

test_that("plot functions run without error", {
    # Synthetic MET Data
    set.seed(123)
    df <- data.frame(
        Year = rep(2020:2022, each = 20),
        Genotype = sample(LETTERS[1:5], 60, replace = TRUE),
        Yield = rnorm(60, 5, 2),
        Row = rep(1:5, 12),
        Column = rep(1:4, 15),
        Trial = rep(c("T1", "T2", "T3"), each = 20)
    )

    # Test Plotting (Run for side effects)
    pdf(NULL) # Sink plot output
    # expect_silent(plot_connectivity(df, x = "Year", trace = "Genotype")) # Deprecated
    # expect_silent(plot_connectivity(df, x = "Year", trace = "Genotype")) # Deprecated
    expect_silent(plot_trend(df, mode = "phenotypic", x = "Year", y = "Yield"))
    expect_silent(plot_spatial(df[df$Trial == "T1", ]))
    dev.off()
})
