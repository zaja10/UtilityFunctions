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
                gdd <- sum(((dat$MaxT + dat$MinT) / 2) - t_base, na.rm = TRUE)
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

