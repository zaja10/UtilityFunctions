#' Calculate Realized Genetic Gain
#'
#' Estimates the rate of genetic gain by regressing breeding values (OP) against year of release.
#' Excludes long-term checks from the regression to avoid bias.
#'
#' @param data A dataframe containing genotype performance.
#' @param year_col String. Column for Year of Release (numeric).
#' @param value_col String. Column for Breeding Value / OP.
#' @param check_list Character vector. Names of checks to exclude from slope calculation.
#'
#' @return A list containing the regression model, gain percentage, and summary stats.
#' @export
calculate_realized_gain <- function(data, year_col = "Year", value_col = "OP", check_list = NULL) {
    df <- data[!is.na(data[[year_col]]) & !is.na(data[[value_col]]), ]

    # Filter breeding population (exclude checks)
    breed_pop <- df
    if (!is.null(check_list)) {
        breed_pop <- df[!df$Genotype %in% check_list, ]
    }

    if (nrow(breed_pop) < 5) stop("Not enough data points for regression.")

    form <- as.formula(paste(value_col, "~", year_col))
    mod <- lm(form, data = breed_pop)

    slope <- coef(mod)[2]
    intercept <- coef(mod)[1]

    # Calculate gain relative to 1990 baseline (or min year)
    base_year <- min(breed_pop[[year_col]])
    base_val <- predict(mod, newdata = data.frame(setNames(list(base_year), year_col)))

    pct_gain <- (slope / base_val) * 100

    return(list(
        model = mod,
        slope = slope,
        r_squared = summary(mod)$r.squared,
        pct_gain_per_year = pct_gain,
        base_year = base_year
    ))
}

#' Plot Genetic Trend (Era Plot)
#'
#' Visualizes the genetic trend over time, highlighting checks and the regression line.
#'
#' @param data A dataframe containing genotype performance.
#' @param year_col String. Column for Year of Release.
#' @param value_col String. Column for Breeding Value / OP.
#' @param check_list Character vector. Names of checks to highlight.
#' @param ... Additional arguments to plot.
#' @export
plot_genetic_trend <- function(data, year_col = "Year", value_col = "OP", check_list = NULL, ...) {
    df <- data[!is.na(data[[year_col]]) & !is.na(data[[value_col]]), ]

    # Plot Base
    plot(df[[year_col]], df[[value_col]],
        type = "n",
        xlab = "Year of Release", ylab = "Genetic Value (BLUP)",
        main = "Genetic Trend (Era Plot)", ...
    )
    grid()

    # Add Regression (Breeding Pop only)
    breed_pop <- df
    if (!is.null(check_list)) {
        breed_pop <- df[!df$Genotype %in% check_list, ]
    }

    if (nrow(breed_pop) > 5) {
        mod <- lm(as.formula(paste(value_col, "~", year_col)), data = breed_pop)
        abline(mod, col = "blue", lwd = 2)

        # Add Gain Text
        slope <- coef(mod)[2]
        base_val <- predict(mod, newdata = data.frame(setNames(list(min(breed_pop[[year_col]])), year_col)))
        pct <- round((slope / base_val) * 100, 2)
        legend("topleft", legend = paste("Gain:", pct, "% / yr"), bty = "n", text.col = "blue")
    }

    # Points
    # Breeding Pop
    points(breed_pop[[year_col]], breed_pop[[value_col]], pch = 19, col = alpha("grey", 0.6))

    # Checks
    if (!is.null(check_list)) {
        checks <- df[df$Genotype %in% check_list, ]
        if (nrow(checks) > 0) {
            points(checks[[year_col]], checks[[value_col]], pch = 17, col = "red", cex = 1.2)
            text(checks[[year_col]], checks[[value_col]], labels = checks$Genotype, pos = 3, cex = 0.7, col = "red")
        }
    }
}

# Helper for manual alpha if scales not loaded
alpha <- function(col, alpha_val) {
    rgb_vals <- col2rgb(col)
    rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], max = 255, alpha = alpha_val * 255)
}
