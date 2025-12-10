source("../../R/met_diagnostics.r")

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

test_that("get_connectivity returns correct counts", {
    # 2 Years, 3 Genotypes
    # Y1: G1, G2
    # Y2: G2, G3
    df <- data.frame(
        Year = c("Y1", "Y1", "Y2", "Y2"),
        Genotype = c("G1", "G2", "G2", "G3")
    )

    mat <- get_connectivity(df, "Year", "Year", "Genotype")

    # Diagonal: Total Genotypes per Year
    expect_equal(mat["Y1", "Y1"], 2)
    expect_equal(mat["Y2", "Y2"], 2)

    # Off-diagonal: Shared (G2) -> 1
    expect_equal(mat["Y1", "Y2"], 1)

    # Test Jaccard
    # Y1 has 2, Y2 has 2. Shared 1. Union = 3. J = 1/3 = 0.333
    jacc <- get_connectivity(df, "Year", "Year", "Genotype", method = "jaccard")
    expect_equal(jacc["Y1", "Y2"], 1 / 3, tolerance = 1e-4)
})

test_that("convert_buac_to_tha calculates correctly", {
    # Wheat: 100 bu/ac -> ~6.725 t/ha
    y_wheat <- convert_buac_to_tha(100, crop = "wheat")
    expect_equal(y_wheat, 6.725, tolerance = 1e-2)

    # Corn: 100 bu/ac -> ~6.277 t/ha
    y_corn <- convert_buac_to_tha(100, crop = "corn")
    expect_equal(y_corn, 6.277, tolerance = 1e-2)

    # Custom: 50lbs -> 50 * 0.00112... * 100 = 5.604
    y_cust <- convert_buac_to_tha(100, lbs_per_bu = 50)
    expect_equal(y_cust, 5.604, tolerance = 1e-2)
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
    expect_silent(plot_connectivity(df, x = "Year", trace = "Genotype"))
    expect_silent(plot_met_trend(df, x = "Year", y = "Yield"))
    expect_silent(plot_trial_map(df, trial_val = "T1"))
    dev.off()
})
