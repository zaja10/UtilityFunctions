test_that("calculate_vegetation_indices input validation", {
    expect_error(calculate_vegetation_indices("not a dataframe"), "spectra_df must be a data frame")
})

test_that("calculate_vegetation_indices works with valid input", {
    # Create a minimal dataset conducive to calculating NDVI (N, R) and GNDVI (N, G)
    df <- data.frame(id = 1, N = 0.8, R = 0.1, G = 0.2)

    # Suppress CLI messages during test
    suppressMessages({
        res <- calculate_vegetation_indices(df)
    })

    expect_true(is.data.frame(res))
    # NDVI = (N - R) / (N + R) = 0.7 / 0.9 = 0.7777...
    expect_equal(res$NDVI[1], (0.8 - 0.1) / (0.8 + 0.1))
    # GNDVI = (N - G) / (N + G) = 0.6 / 1.0 = 0.6
    expect_equal(res$GNDVI[1], (0.8 - 0.2) / (0.8 + 0.2))

    # Ensure original columns are kept
    expect_true(all(c("N", "R", "G") %in% names(res)))
})

test_that("process_longitudinal_indices works with mocked data", {
    # Mock data: 2 timepoints, full suite of bands
    # Traits: Blue, Green, Red, RedEdge, NIR
    # Timepoints: 1, 2

    df <- data.frame(
        PlotID = c("P1", "P2"),
        Blue.1 = c(0.1, 0.11), Green.1 = c(0.2, 0.22), Red.1 = c(0.1, 0.12), RedEdge.1 = c(0.3, 0.33), NIR.1 = c(0.8, 0.82),
        Blue.2 = c(0.15, 0.16), Green.2 = c(0.25, 0.26), Red.2 = c(0.15, 0.18), RedEdge.2 = c(0.35, 0.38), NIR.2 = c(0.85, 0.88)
    )

    meta_cols <- "PlotID"
    spectral_cols <- names(df)[-1] # All except PlotID

    suppressMessages({
        res <- process_longitudinal_indices(df, meta_cols, spectral_cols)
    })

    # Expecting: 2 plots * 2 timepoints = 4 rows
    expect_equal(nrow(res), 4)
    expect_true("Day" %in% names(res))
    expect_equal(sort(unique(res$Day)), c(1, 2))

    # Check if NDVI is calculated
    expect_true("NDVI" %in% names(res))

    # Check value for P1 at Time 2
    # P1, Time 2: N=0.85, R=0.15 -> NDVI = 0.7 / 1.0 = 0.7
    p1_t2 <- subset(res, PlotID == "P1" & Day == 2)
    expect_equal(p1_t2$NDVI, 0.7)
})

test_that("process_longitudinal_indices handles missing bands gracefully", {
    # Data only has Blue for Time 1 -> Required bands (default map) missing.

    df <- data.frame(
        id = 1,
        Blue.1 = 0.1,
        Blue.2 = 0.15
    )

    # Should silently return NULL if all days are missing bands
    suppressMessages({
        res <- process_longitudinal_indices(df, "id", c("Blue.1", "Blue.2"))
    })
    expect_null(res)
})
