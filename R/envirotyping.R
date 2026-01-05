#' Process Weather Indices for Envirotyping
#'
#' Aggregates daily weather data into biological windows (Vegetative, Flowering, Grain Fill)
#' to create environmental covariates (ECs) suitable for Factorial Regression models.
#'
#' @param weather_data A dataframe containing daily weather. Columns must include:
#'        "Date" (Date object), "Location", "MinT", "MaxT", "Precip".
#'        Optional: "Radiation", "ET0".
#' @param trial_metadata A dataframe containing trial dates. Columns must include:
#'        "TrialID", "Location", "PlantingDate", "FloweringDate", "HarvestDate".
#' @param t_base Numeric. Base temperature for GDD calculation (default 0).
#' @param heat_threshold Numeric. Threshold for heat stress calculation (default 30).
#'
#' @return A matrix of environmental covariates (Rows = Trials, Cols = Indices).
#' @importFrom dplyr %>% filter summarize mutate select
#' @export
process_weather_indices <- function(weather_data, trial_metadata, t_base = 0, heat_threshold = 30) {
    # Validation
    req_weather <- c("Date", "Location", "MinT", "MaxT", "Precip")
    if (!all(req_weather %in% names(weather_data))) stop("Missing columns in weather_data.")

    req_meta <- c("TrialID", "Location", "PlantingDate", "FloweringDate", "HarvestDate")
    if (!all(req_meta %in% names(trial_metadata))) stop("Missing columns in trial_metadata.")

    trials <- unique(trial_metadata$TrialID)
    ec_list <- list()

    for (tid in trials) {
        meta <- trial_metadata[trial_metadata$TrialID == tid, ]
        if (nrow(meta) > 1) meta <- meta[1, ]

        loc <- meta$Location
        p_date <- as.Date(meta$PlantingDate)
        f_date <- as.Date(meta$FloweringDate)
        h_date <- as.Date(meta$HarvestDate)

        # Slicing Weather Data
        w_sub <- weather_data[weather_data$Location == loc, ]
        if (nrow(w_sub) == 0) {
            warning(paste("No weather data for location:", loc))
            next
        }

        # Define Windows
        windows <- list(
            Vegetative = w_sub %>% filter(Date >= p_date & Date < f_date),
            GrainFill = w_sub %>% filter(Date >= f_date & Date <= h_date)
        )

        # Calculate Indices per Window
        trial_ecs <- list(TrialID = tid)

        for (win_name in names(windows)) {
            dat <- windows[[win_name]]

            if (nrow(dat) > 0) {
                # 1. Thermal Time (GDD)
                # Standard: Clamp daily GDD to 0 so cold days don't subtract from total
                daily_gdd <- ((dat$MaxT + dat$MinT) / 2) - t_base
                gdd <- sum(pmax(0, daily_gdd), na.rm = TRUE)
                trial_ecs[[paste0(win_name, "_GDD")]] <- max(0, gdd)

                # 2. Heat Stress
                hs <- sum(dat$MaxT > heat_threshold, na.rm = TRUE)
                trial_ecs[[paste0(win_name, "_HeatStress")]] <- hs

                # 3. Water Deficit (Precip - ET0)
                # If ET0 missing, use simplified proxy or just Precip
                et0 <- if ("ET0" %in% names(dat)) dat$ET0 else 0
                wd <- sum(dat$Precip - et0, na.rm = TRUE)
                trial_ecs[[paste0(win_name, "_WaterDeficit")]] <- wd

                # 4. Total Precip
                trial_ecs[[paste0(win_name, "_Precip")]] <- sum(dat$Precip, na.rm = TRUE)
            } else {
                trial_ecs[[paste0(win_name, "_GDD")]] <- NA
                trial_ecs[[paste0(win_name, "_HeatStress")]] <- NA
                trial_ecs[[paste0(win_name, "_WaterDeficit")]] <- NA
                trial_ecs[[paste0(win_name, "_Precip")]] <- NA
            }
        }
        ec_list[[tid]] <- as.data.frame(trial_ecs)
    }

    ec_df <- do.call(rbind, ec_list)
    rownames(ec_df) <- ec_df$TrialID

    # Return as numeric matrix (excluding ID) for modeling
    mat <- as.matrix(ec_df[, -1])
    return(mat)
}

#' Calculate Days Since Start
#'
#' Converts a vector of Day of Year (DOY) values into "Days Since Start" based on a planting date string.
#' Handles year transitions (e.g., planting in Oct, sampling in Jan).
#'
#' @param doy_vector Numeric vector of DOY values.
#' @param planting_date_string Character string format "YYYY-Month-DD" (e.g., "2023-October-04").
#' @return Numeric vector of days elapsed since planting.
#' @export
calculate_days_since_start <- function(doy_vector, planting_date_string) {
    # --- Helper function to convert DOY to Date ---
    doy_to_date <- function(doy, year) {
        date_time_obj <- strptime(paste(year, doy), format = "%Y %j")
        return(as.Date(date_time_obj))
    }

    # 1. Parse the planting date string
    #    Note: "%Y-%B-%d" matches "2023-October-04"
    start_date <- as.Date(planting_date_string, format = "%Y-%B-%d")
    if (is.na(start_date)) {
        stop("Planting date string format is incorrect. Please use 'YYYY-Month-DD', e.g., '2023-October-04'")
    }

    # 2. Get the planting year and planting day-of-year
    start_year <- as.numeric(format(start_date, "%Y"))
    start_doy <- as.numeric(format(start_date, "%j"))

    # 3. Determine the correct year for each doy in the vector
    #    If a doy is LESS THAN the start doy, it's from the next year.
    year_vector <- ifelse(doy_vector < start_doy,
        start_year + 1,
        start_year
    )

    # 4. Convert all (doy, year) pairs into real dates
    date_vector <- doy_to_date(doy = doy_vector, year = year_vector)

    # 5. Calculate the difference in days from the start date
    days_diff <- difftime(date_vector, start_date, units = "days")

    # Return as a simple number
    return(as.numeric(days_diff))
}
