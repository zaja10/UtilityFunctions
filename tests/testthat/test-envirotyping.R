test_that("calculate_days_since_start spans New Year", {
  res <- calculate_days_since_start(c(277, 278, 1), "2023-October-04")
  expect_equal(res, c(0, 1, as.numeric(as.Date("2024-01-01") - as.Date("2023-10-04"))))
})

test_that("calculate_days_since_start accepts ISO dates and Date objects", {
  expect_equal(calculate_days_since_start(278, "2023-10-04"), 1)
  expect_equal(calculate_days_since_start(278, as.Date("2023-10-04")), 1)
})

test_that("calculate_days_since_start does not depend on the locale", {
  withr::with_locale(c(LC_TIME = "fr_FR.UTF-8"), {
    expect_equal(calculate_days_since_start(278, "2023-October-04"), 1)
  })
})

test_that("calculate_days_since_start rejects bad input", {
  expect_error(calculate_days_since_start(100, "04/10/2023"), "format is incorrect")
  expect_error(calculate_days_since_start(400, "2023-10-04"), "1..366")
})

weather <- function(days, et0 = TRUE) {
  w <- data.frame(Date = days, Location = "L1", MinT = 10, MaxT = 32, Precip = 2)
  if (et0) w$ET0 <- 5
  w
}
meta <- data.frame(TrialID = "T1", Location = "L1", PlantingDate = "2024-05-01", FloweringDate = "2024-05-11", HarvestDate = "2024-05-20")

test_that("process_weather_indices computes window covariates", {
  w <- weather(seq(as.Date("2024-05-01"), as.Date("2024-05-20"), by = "day"))
  m <- process_weather_indices(w, meta, t_base = 0, heat_threshold = 30)
  expect_equal(rownames(m), "T1")
  expect_equal(m[, "Vegetative_GDD"], 10 * 21)
  expect_equal(m[, "GrainFill_GDD"], 10 * 21)
  expect_equal(m[, "Vegetative_HeatStress"], 10)
  expect_equal(m[, "GrainFill_WaterBalance"], 10 * (2 - 5)) # negative = deficit
})

test_that("process_weather_indices omits water balance without ET0 and NA's gappy windows", {
  w <- weather(seq(as.Date("2024-05-01"), as.Date("2024-05-20"), by = "day"), et0 = FALSE)
  m <- process_weather_indices(w, meta)
  expect_false(any(grepl("WaterBalance", colnames(m))))
  gappy <- w[-(12:20), ]
  expect_warning(m2 <- process_weather_indices(gappy, meta), "coverage")
  expect_true(all(is.na(m2[, grep("^GrainFill", colnames(m2))])))
  expect_false(anyNA(m2[, grep("^Vegetative", colnames(m2))]))
})
