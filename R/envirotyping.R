#' Process Weather Indices for Envirotyping
#'
#' Aggregates daily weather into two developmental windows, Vegetative
#' (planting to the day before flowering) and GrainFill (flowering to harvest),
#' producing environmental covariates (ECs) for factorial regression.
#'
#' @param weather_data Daily weather with columns \code{Date},
#'   \code{Location}, \code{MinT}, \code{MaxT}, \code{Precip}; optional
#'   \code{ET0}.
#' @param trial_metadata Trial dates with columns \code{TrialID},
#'   \code{Location}, \code{PlantingDate}, \code{FloweringDate},
#'   \code{HarvestDate}.
#' @param t_base Numeric. Base temperature for growing degree days (default 0).
#' @param heat_threshold Numeric. Daily maximum above which a day counts as
#'   heat stress (default 30).
#' @param min_coverage Numeric in (0, 1]. Minimum fraction of days in a window
#'   with weather records; windows below it return \code{NA} (default 0.9).
#'
#' @return A numeric matrix (rows = trials) with, per window, \code{_GDD}
#'   (sum of daily \eqn{\max(0, (T_{max}+T_{min})/2 - T_{base})}),
#'   \code{_HeatStress} (days above \code{heat_threshold}), \code{_Precip},
#'   and, when \code{ET0} is present, \code{_WaterBalance}
#'   (\eqn{\sum(P - ET_0)}; negative values indicate a deficit).
#' @export
process_weather_indices <- function(weather_data, trial_metadata, t_base = 0, heat_threshold = 30, min_coverage = 0.9) {
  req_weather <- c("Date", "Location", "MinT", "MaxT", "Precip")
  req_meta <- c("TrialID", "Location", "PlantingDate", "FloweringDate", "HarvestDate")
  miss_w <- setdiff(req_weather, names(weather_data))
  miss_m <- setdiff(req_meta, names(trial_metadata))
  if (length(miss_w) > 0) cli::cli_abort("Missing column{?s} in {.arg weather_data}: {.val {miss_w}}.")
  if (length(miss_m) > 0) cli::cli_abort("Missing column{?s} in {.arg trial_metadata}: {.val {miss_m}}.")

  weather_data$Date <- as.Date(weather_data$Date)
  has_et0 <- "ET0" %in% names(weather_data)
  stats_names <- c("GDD", "HeatStress", "Precip", if (has_et0) "WaterBalance")

  window_stats <- function(w, start, end) {
    n_days <- as.numeric(end - start) + 1
    dat <- w[w$Date >= start & w$Date <= end, , drop = FALSE]
    if (n_days <= 0 || nrow(dat) / n_days < min_coverage) {
      return(stats::setNames(rep(NA_real_, length(stats_names)), stats_names))
    }
    out <- c(
      GDD = sum(pmax(0, (dat$MaxT + dat$MinT) / 2 - t_base), na.rm = TRUE),
      HeatStress = sum(dat$MaxT > heat_threshold, na.rm = TRUE),
      Precip = sum(dat$Precip, na.rm = TRUE)
    )
    if (has_et0) out["WaterBalance"] <- sum(dat$Precip - dat$ET0, na.rm = TRUE)
    out
  }

  trials <- unique(as.character(trial_metadata$TrialID))
  low_cov <- character(0)
  rows <- lapply(trials, function(tid) {
    meta <- trial_metadata[trial_metadata$TrialID == tid, , drop = FALSE][1, ]
    p_date <- as.Date(meta$PlantingDate)
    f_date <- as.Date(meta$FloweringDate)
    h_date <- as.Date(meta$HarvestDate)
    if (anyNA(c(p_date, f_date, h_date)) || !(p_date < f_date && f_date <= h_date)) {
      cli::cli_abort("Trial {.val {tid}}: dates must satisfy planting < flowering <= harvest.")
    }
    w <- weather_data[weather_data$Location == meta$Location, , drop = FALSE]
    veg <- window_stats(w, p_date, f_date - 1)
    gf <- window_stats(w, f_date, h_date)
    if (anyNA(c(veg, gf))) low_cov <<- c(low_cov, tid)
    c(stats::setNames(veg, paste0("Vegetative_", names(veg))), stats::setNames(gf, paste0("GrainFill_", names(gf))))
  })

  if (length(low_cov) > 0) {
    cli::cli_warn("Weather coverage below {100 * min_coverage}% for {length(low_cov)} trial{?s} (ECs set to NA): {.val {low_cov}}.")
  }
  mat <- do.call(rbind, rows)
  rownames(mat) <- trials
  mat
}

#' Calculate Days Since Start
#'
#' Converts day-of-year (DOY) values to days elapsed since planting. A DOY
#' smaller than the planting DOY is taken to fall in the following year, so
#' seasons may span one New Year.
#'
#' @param doy_vector Numeric vector of DOY values (1-366).
#' @param planting_date_string Planting date: a \code{Date}, an ISO string
#'   (\code{"2023-10-04"}) or \code{"YYYY-Month-DD"} with an English month
#'   name (\code{"2023-October-04"}). Parsing does not depend on the locale.
#' @return Numeric vector of days elapsed since planting.
#' @export
calculate_days_since_start <- function(doy_vector, planting_date_string) {
  start_date <- parse_planting_date(planting_date_string)
  if (any(doy_vector < 1 | doy_vector > 366, na.rm = TRUE)) cli::cli_abort("DOY values must be in 1..366.")

  start_year <- as.integer(format(start_date, "%Y"))
  start_doy <- as.integer(format(start_date, "%j"))
  year <- ifelse(doy_vector < start_doy, start_year + 1L, start_year)
  dates <- as.Date(paste0(year, "-01-01")) + (doy_vector - 1)
  as.numeric(dates - start_date)
}

parse_planting_date <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }
  x <- as.character(x)
  iso <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  if (!is.na(iso) && grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", x)) {
    return(iso)
  }
  m <- regmatches(x, regexec("^(\\d{4})-([A-Za-z]+)-(\\d{1,2})$", x))[[1]]
  if (length(m) == 4) {
    mon <- match(tolower(m[3]), tolower(month.name))
    if (is.na(mon)) mon <- match(tolower(m[3]), tolower(month.abb))
    if (!is.na(mon)) {
      d <- as.Date(sprintf("%s-%02d-%02d", m[2], mon, as.integer(m[4])), format = "%Y-%m-%d")
      if (!is.na(d)) {
        return(d)
      }
    }
  }
  cli::cli_abort("Planting date string format is incorrect. Use 'YYYY-MM-DD' or 'YYYY-Month-DD', e.g. '2023-October-04'.")
}
