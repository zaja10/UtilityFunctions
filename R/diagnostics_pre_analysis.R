#' Pad Trial Layout (Robust Integer Merge)
#'
#' Ensures a trial dataset has a complete grid of Row x Column positions by padding
#' missing coordinates with NA values. Uses strict integer matching to prevent duplication.
#'
#' @param data A data.frame containing the trial data.
#' @param row String. Column name for Row positions (Default: "Row").
#' @param col String. Column name for Column positions (Default: "Column").
#' @param group String. Column name for Trial/Environment to pad independently.
#'
#' @return A data.frame with the original data plus rows for missing physical positions.
#' @export
pad_trial_layout <- function(data, row = "Row", col = "Column", group = NULL) {
    # 1. Recursive handling for groups
    if (!is.null(group)) {
        if (!group %in% names(data)) stop(paste("Group column", group, "not found."))

        # Filter NA groups
        valid_data <- data[!is.na(data[[group]]), ]

        # Split and apply
        out_list <- split(valid_data, valid_data[[group]])
        padded_list <- lapply(names(out_list), function(grp_name) {
            sub_df <- out_list[[grp_name]]
            # Recurse for this specific trial
            padded <- pad_trial_layout(sub_df, row = row, col = col, group = NULL)

            # Restore group label
            val <- sub_df[[group]][1]
            padded[[group]] <- val
            return(padded)
        })

        out_df <- do.call(rbind, padded_list)
        rownames(out_df) <- NULL
        return(out_df)
    }

    # 2. Base Case: Pad Single Experiment
    # Create SAFE integer coordinates for merging
    data$..Row_Int.. <- tryCatch(as.integer(as.character(data[[row]])), warning = function(w) NA)
    data$..Col_Int.. <- tryCatch(as.integer(as.character(data[[col]])), warning = function(w) NA)

    # If conversion fails (e.g. Row is "Plot1"), abort padding to be safe
    if (all(is.na(data$..Row_Int..)) || all(is.na(data$..Col_Int..))) {
        warning(paste("Could not convert", row, "or", col, "to integers. Skipping padding."))
        data$..Row_Int.. <- NULL
        data$..Col_Int.. <- NULL
        return(data)
    }

    # 3. Create Grid based on Integers
    min_r <- min(data$..Row_Int.., na.rm = TRUE)
    max_r <- max(data$..Row_Int.., na.rm = TRUE)
    min_c <- min(data$..Col_Int.., na.rm = TRUE)
    max_c <- max(data$..Col_Int.., na.rm = TRUE)

    grid <- expand.grid(
        ..Row_Int.. = seq(min_r, max_r),
        ..Col_Int.. = seq(min_c, max_c)
    )

    # 4. Strict Merge on Integers
    data$..present.. <- TRUE
    out <- merge(grid, data, by = c("..Row_Int..", "..Col_Int.."), all = TRUE)

    # 5. Post-Processing
    padded_rows <- is.na(out$..present..)

    if (any(padded_rows)) {
        out[padded_rows, row] <- as.character(out[padded_rows, "..Row_Int.."])
        out[padded_rows, col] <- as.character(out[padded_rows, "..Col_Int.."])

        # Convert back to factor if original was factor
        if (is.factor(data[[row]])) out[[row]] <- as.factor(out[[row]])
        if (is.factor(data[[col]])) out[[col]] <- as.factor(out[[col]])

        # Fill other constant columns (like Location, Year)
        for (nm in names(out)) {
            if (nm %in% c(row, col, "..Row_Int..", "..Col_Int..", "..present..")) next

            vals <- out[!padded_rows, nm]
            u_vals <- unique(vals[!is.na(vals)])

            if (length(u_vals) == 1) {
                out[padded_rows, nm] <- u_vals[1]
            }
        }
    }

    # Remove temporary columns
    out$..Row_Int.. <- NULL
    out$..Col_Int.. <- NULL
    out$..present.. <- NULL

    # Sort
    out <- out[order(out[[col]], out[[row]]), ]

    return(out)
}

#' Find Missing Row/Column Combinations in MET Data
#'
#' Identifies missing plots required to complete a rectangular grid for spatial analysis.
#' Useful for diagnosing ASREml "residual model implies X" errors.
#'
#' @param data Dataframe containing the trial data.
#' @param experiment String. Column name for the Experiment/Trial ID.
#' @param row String. Column name for Row.
#' @param col String. Column name for Column.
#'
#' @return A dataframe listing the Experiment, Row, and Column of missing plots.
#' @export
find_missing_plots <- function(data, experiment = "Experiment", row = "Row", col = "Column") {
    if (!all(c(experiment, row, col) %in% names(data))) {
        stop("Specified columns not found in dataframe.")
    }

    # Ensure strict types for merging
    data[[experiment]] <- as.factor(data[[experiment]])

    # Iterate through each experiment
    missing_list <- lapply(levels(data[[experiment]]), function(exp_id) {
        # Subset data for this experiment
        sub_df <- data[data[[experiment]] == exp_id, ]

        if (nrow(sub_df) == 0) {
            return(NULL)
        }

        r_vals <- as.numeric(as.character(sub_df[[row]]))
        c_vals <- as.numeric(as.character(sub_df[[col]]))

        min_r <- min(r_vals, na.rm = TRUE)
        max_r <- max(r_vals, na.rm = TRUE)
        min_c <- min(c_vals, na.rm = TRUE)
        max_c <- max(c_vals, na.rm = TRUE)

        expected_grid <- expand.grid(
            Row_Num = seq(min_r, max_r),
            Col_Num = seq(min_c, max_c)
        )

        sub_df$key <- paste(r_vals, c_vals, sep = "_")
        expected_grid$key <- paste(expected_grid$Row_Num, expected_grid$Col_Num, sep = "_")

        missing_keys <- setdiff(expected_grid$key, sub_df$key)

        if (length(missing_keys) > 0) {
            coords <- do.call(rbind, strsplit(missing_keys, "_"))
            return(data.frame(
                Experiment = exp_id,
                Row = as.numeric(coords[, 1]),
                Column = as.numeric(coords[, 2]),
                Status = "MISSING",
                stringsAsFactors = FALSE
            ))
        } else {
            return(NULL)
        }
    })

    out_df <- do.call(rbind, missing_list)

    if (is.null(out_df)) {
        message("Success: All experiments have complete rectangular grids.")
        return(invisible(NULL))
    } else {
        message(paste("Found", nrow(out_df), "missing plots across", length(unique(out_df$Experiment)), "experiments."))
        return(out_df)
    }
}

#' Convert Yield Units
#'
#' Converts yield from Bushels per Acre (bu/ac) to Tonnes per Hectare (t/ha).
#'
#' @param yield Numeric vector of yield values.
#' @param crop String. Crop type to determine bushel weight. Options: "wheat", "soybeans", "corn", "barley", "oats".
#' @param lbs_per_bu Numeric. Custom weight per bushel (overrides crop).
#'
#' @return Numeric vector of yield in t/ha.
#' @export
convert_buac_to_tha <- function(yield, crop = "wheat", lbs_per_bu = NULL) {
    weights <- c(
        "wheat" = 60, "soybeans" = 60, "soy" = 60,
        "corn" = 56, "maize" = 56, "barley" = 48, "oats" = 32
    )

    if (!is.null(lbs_per_bu)) {
        weight <- lbs_per_bu
    } else {
        crop <- tolower(crop)
        if (crop %in% names(weights)) {
            weight <- weights[[crop]]
        } else {
            warning("Crop not recognized. Defaulting to 60 lbs/bu (Wheat). Use 'lbs_per_bu' to specify.")
            weight <- 60
        }
    }

    # Factor = (weight_lbs * 0.453592 kg/lb) / (0.404686 ha/ac * 1000 kg/t)
    factor <- weight * 0.00112085116
    return(yield * factor)
}

#' Flexible Trial Map
#'
#' Plots the spatial layout of a specific trial. Adapts to missing spatial coordinates
#' by falling back to a 1D visualization if necessary.
#'
#' @param data A data.frame.
#' @param trial_val Value. The specific trial ID to filter.
#' @param trial_col String. Column name for Trial ID.
#' @param row String. Column name for Row.
#' @param col String. Column name for Column.
#' @param val String. Variable to fill heatmap (e.g. "Yield").
#'
#' @importFrom grDevices hcl.colors
#' @importFrom graphics image layout par box axis points grid plot
#' @export
plot_trial_map <- function(data, trial_val = NULL, trial_col = "Trial",
                           row = "Row", col = "Column", val = "Yield") {
    # Filter Data
    if (!is.null(trial_val) && trial_col %in% names(data)) {
        plot_data <- data[data[[trial_col]] == trial_val, ]
        title_sub <- paste(trial_val)
    } else {
        plot_data <- data
        title_sub <- "Single Site"
    }

    # Check if Row/Col exist
    has_spatial <- all(c(row, col) %in% names(plot_data))

    # Setup Palette
    vals <- plot_data[[val]]
    if (is.null(vals)) stop(paste("Variable", val, "not found."))

    cols <- hcl.colors(20, "Spectral", rev = TRUE)

    if (has_spatial) {
        # Padded Map
        padded <- pad_trial_layout(plot_data, row, col)

        r_vals <- as.numeric(padded[[row]])
        c_vals <- as.numeric(padded[[col]])

        mat <- matrix(NA, nrow = max(r_vals, na.rm = T), ncol = max(c_vals, na.rm = T))
        for (i in 1:nrow(padded)) {
            if (!is.na(r_vals[i]) & !is.na(c_vals[i])) {
                mat[r_vals[i], c_vals[i]] <- padded[[val]][i]
            }
        }

        # Image Plot
        layout(matrix(1:2, ncol = 2), widths = c(4, 1))
        par(mar = c(4, 4, 3, 1))
        image(1:ncol(mat), 1:nrow(mat), t(mat),
            col = cols,
            main = paste("Trial Map:", title_sub),
            xlab = col, ylab = row
        )
        box()

        # Legend
        par(mar = c(4, 1, 3, 2))
        image(1, seq(min(vals, na.rm = T), max(vals, na.rm = T), length.out = 20),
            t(as.matrix(1:20)),
            axes = F, xlab = "", ylab = "", col = cols
        )
        axis(4, las = 1)
        layout(1)
    } else {
        # 1D Fallback Plan
        warning("Spatial columns missing. Plotting linear sequence.")
        plot(vals, type = "n", main = paste("Linear Plot:", title_sub), ylab = val, xlab = "Index")
        points(vals, pch = 21, bg = cols[cut(vals, 20)], cex = 1.5)
        grid()
    }
}

#' Plot MET Trends
#'
#' Visualizes performance trends over time using a combined Boxplot (distribution)
#' and Linear Regression (rate of change) approach.
#'
#' @param data A data.frame.
#' @param x String. Categorical variable for X-axis (e.g. "Year"). Must be convertible to numeric for trend analysis.
#' @param y String. Response variable (e.g. "Yield").
#' @param main String. Plot title prefix.
#' @param ... Additional arguments passed to `plot` or `boxplot`.
#'
#' @importFrom stats lm coef
#' @importFrom graphics boxplot abline legend
#' @export
plot_met_trend <- function(data, x = "Year", y = "Yield", main = "Yield Trend", ...) {
    if (!all(c(x, y) %in% names(data))) stop("Variables not found in data.")

    raw_x <- data[[x]]
    if (is.numeric(raw_x)) {
        x_num <- raw_x
    } else {
        x_num <- suppressWarnings(as.numeric(as.character(raw_x)))
        if (all(is.na(x_num))) stop(paste("Column", x, "must be numeric (e.g. Year) for trend analysis."))
    }
    x_fac <- as.factor(raw_x)

    # Fits
    tmp_df <- data.frame(Y = data[[y]], X = x_num)
    model <- lm(Y ~ X, data = tmp_df)
    slope <- coef(model)["X"]
    r_sq <- summary(model)$r.squared

    # Visuals
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    layout(matrix(1:2, nrow = 2))

    # Boxplot
    par(mar = c(3, 5, 3, 2))
    boxplot(data[[y]] ~ x_fac,
        col = "lightgreen", border = "darkgreen",
        main = paste(main, "(Distribution)"),
        ylab = y, las = 1, ...
    )

    # Line Plot
    par(mar = c(5, 5, 2, 2))
    means <- tapply(tmp_df$Y, tmp_df$X, mean, na.rm = TRUE)
    years_plot <- as.numeric(names(means))

    plot(years_plot, means,
        type = "b", pch = 19, lwd = 2, col = "blue",
        main = paste(main, "(Rate of Change)"),
        xlab = x, ylab = paste("Mean", y),
        xaxt = "n", las = 1
    )
    axis(1, at = years_plot, labels = years_plot)
    grid()
    abline(model, col = "red", lwd = 2)

    legend("topleft",
        legend = c(paste("Rate:", round(slope, 4), "unit/yr"), paste("R2:", round(r_sq, 3))),
        bty = "n", cex = 0.9, text.col = "black"
    )

    layout(1)
}
