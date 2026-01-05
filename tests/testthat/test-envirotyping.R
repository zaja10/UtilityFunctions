# Tests for Envirotyping Helpers

test_that("calculate_days_since_start works correctly", {
    # Planting: 2023-Oct-04 (DOY 277)
    planting_date <- "2023-October-04"

    # DOY 277 (Oct 4) -> 0 days
    # DOY 278 (Oct 5) -> 1 day
    # DOY 1 (Jan 1, 2024) -> Should clearly be after Oct 4 -> Year + 1 logic

    doys <- c(277, 278, 1) # Oct 4, Oct 5, Jan 1

    res <- calculate_days_since_start(doy_vector = doys, planting_date_string = planting_date)

    expect_equal(res[1], 0)
    expect_equal(res[2], 1)

    # Oct 4 to Dec 31 2023 is (365 - 277) = 88 days?
    # Oct 4 is day 0.
    # Jan 1 2024 is...
    # difftime logic handles leap years etc.
    # 2024 is leap year, but that doesn't affect Jan 1 date necessarily

    expected_jan1 <- as.numeric(difftime(as.Date("2024-01-01"), as.Date("2023-10-04"), units = "days"))
    expect_equal(res[3], expected_jan1)
})

test_that("calculate_days_since_start throws error on bad date format", {
    expect_error(calculate_days_since_start(100, "2023-10-04"), "Planting date string format is incorrect")
})
