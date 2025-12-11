#' MET Dataset Diagnostics & Visualization
#'
#' @description
#' A suite of tools for inspecting large Multi-Environment Trial (MET) datasets
#' prior to complex analysis. Includes checking connectivity, temporal trends,
#' and spatial integrity.
#'
#' @name met_diagnostics
NULL

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
    # This strips away factor levels, whitespace, etc.
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
    # all=TRUE keeps grid rows (padding) AND data rows
    data$..present.. <- TRUE
    out <- merge(grid, data, by = c("..Row_Int..", "..Col_Int.."), all = TRUE)

    # 5. Post-Processing
    padded_rows <- is.na(out$..present..)

    # If the original columns (Row/Column) are NA in the new rows, fill them
    # This puts the factor labels back (e.g., Integer 1 becomes Factor "1")
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

#' Calculate Connectivity Matrix
#'
#' Computes the connectivity (shared levels of a trace factor) between levels of
#' grouping factors (e.g. Years or Trials). Can return raw counts or normalized
#' similarity indices pairwise.
#'
#' @param data A data.frame.
#' @param x_fac String. First grouping factor (e.g., "Year").
#' @param y_fac String. Second grouping factor (default same as x_fac for symmetric matrix).
#' @param trace String. The linking factor (e.g., "Genotype").
#' @param method String. Metric to calculate:
#'   \itemize{
#'     \item \code{"count"} (Default): Raw number of shared genotypes.
#'     \item \code{"jaccard"}: Intersection / Union. Good for varying sizes.
#'     \item \code{"prop_min"}: Intersection / size of smaller group.
#'   }
#'
#' @return A matrix of connectivity scores.
#' @export
get_connectivity <- function(data, x_fac = "Year", y_fac = NULL, trace = "Genotype", method = "count") {
    if (is.null(y_fac)) y_fac <- x_fac

    # Filter complete cases
    df <- data[!is.na(data[[x_fac]]) & !is.na(data[[trace]]), ]
    if (y_fac != x_fac) df <- df[!is.na(df[[y_fac]]), ]

    f1 <- as.factor(df[[x_fac]])
    f2 <- as.factor(df[[y_fac]])
    tr <- as.factor(df[[trace]])

    # Create incidence matrices
    # M1: Rows = Trace, Cols = X_fac
    M1 <- table(tr, f1)
    M1[M1 > 0] <- 1

    if (x_fac == y_fac) {
        # Symmetric
        conn <- crossprod(M1) # Intersection

        if (method != "count") {
            # Sizes of each group (Diagonal of intersection matrix)
            sizes <- diag(conn)

            if (method == "jaccard") {
                # Union = SizeA + SizeB - Intersection
                denom <- outer(sizes, sizes, "+") - conn
                conn <- conn / denom
            } else if (method == "prop_min") {
                # Denom = min(SizeA, SizeB)
                # pmin is not outer-aware, use simple logic
                s_mat_row <- matrix(sizes, nrow = length(sizes), ncol = length(sizes))
                s_mat_col <- t(s_mat_row)
                denom <- pmin(s_mat_row, s_mat_col)
                conn <- conn / denom
            }
        }
    } else {
        # Asymmetric
        M2 <- table(tr, f2)
        M2[M2 > 0] <- 1

        # Align rows
        common <- intersect(rownames(M1), rownames(M2))
        subM1 <- M1[common, , drop = FALSE]
        subM2 <- M2[common, , drop = FALSE]
        conn <- crossprod(subM1, subM2)

        if (method != "count") {
            s1 <- colSums(subM1)
            s2 <- colSums(subM2)

            if (method == "jaccard") {
                denom <- outer(s1, s2, "+") - conn
                conn <- conn / denom
            } else if (method == "prop_min") {
                s1_mat <- matrix(s1, nrow = length(s1), ncol = length(s2))
                s2_mat <- matrix(s2, nrow = length(s1), ncol = length(s2), byrow = TRUE)
                denom <- pmin(s1_mat, s2_mat)
                conn <- conn / denom
            }
        }
    }

    # Clean up potential NaNs (e.g. division by zero if empty groups)
    conn[is.na(conn)] <- 0

    return(as.matrix(conn))
}

#' Visualise MET Connectivity
#'
#' Visualizes the connectivity structure of a MET dataset using a heatmap.
#' Identify if trials or years are connected by common germplasm.
#'
#' @param data A data.frame.
#' @param x String. Factor for X-axis (e.g. "Year" or "Trial").
#' @param y String. Factor for Y-axis (Default `NULL` uses x for symmetric).
#' @param trace String. Variable establishing connection (e.g. "Genotype").
#' @param method String. "count", "jaccard", or "prop_min".
#' @param order_by String. Column name to use for sorting the axes (e.g. "Year").
#'        If NULL, defaults to "Year" if present. Set to FALSE to disable sorting.
#' @param ... Additional arguments passed to `image()`.
#'
#' @export
plot_connectivity <- function(data, x = "Year", y = NULL, trace = "Genotype", method = "count", order_by = NULL, ...) {
    mat <- get_connectivity(data, x, y, trace, method)

    # Determine sorting column
    sort_col <- NULL
    if (is.null(order_by)) {
        if ("Year" %in% names(data)) sort_col <- "Year"
    } else if (!isFALSE(order_by)) {
        if (order_by %in% names(data)) {
            sort_col <- order_by
        } else {
            warning(paste("Ordering column", order_by, "not found. Skipping sort."))
        }
    }

    if (!is.null(sort_col)) {
        # Helper to get order
        get_order <- function(fac_name, current_levels) {
            if (fac_name == sort_col) {
                return(current_levels)
            } # Already sorted usually, or numeric

            # Find mapping of Factor -> SortCol
            meta <- unique(data[, c(fac_name, sort_col)])
            meta <- meta[meta[[fac_name]] %in% current_levels, ]

            # Aggregate if multiple sort-values per factor level (use min)
            ord_df <- aggregate(as.formula(paste(sort_col, "~", fac_name)), data = meta, min)

            # Return levels sorted by the sort_col
            ordering <- ord_df[order(ord_df[[sort_col]], ord_df[[fac_name]]), fac_name]
            return(as.character(ordering))
        }

        # Sort Rows (X axis in get_connectivity)
        row_ord <- get_order(x, rownames(mat))

        # Sort Cols (Y axis in get_connectivity)
        y_name <- if (is.null(y)) x else y
        col_ord <- get_order(y_name, colnames(mat))

        # Intersect to be safe (in case get_order missed something)
        row_ord <- intersect(row_ord, rownames(mat))
        col_ord <- intersect(col_ord, colnames(mat))

        mat <- mat[row_ord, col_ord, drop = FALSE]
    }

    # Plotting
    # Image puts (0,0) at bottom left.
    # We want a standard matrix view (Row 1 at top). Reorder rows.
    mat_rev <- mat[nrow(mat):1, , drop = FALSE]

    # Color scale: White (0) -> Blue (Max)
    cols <- colorRampPalette(c("white", "aliceblue", "#3498db", "#2c3e50"))(20)

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mar = c(5, 5, 4, 3))
    image(1:ncol(mat_rev), 1:nrow(mat_rev), t(mat_rev),
        axes = FALSE, col = cols, xlab = ifelse(is.null(y) || x == y, x, y), ylab = x,
        main = paste("Connectivity by", trace)
    )

    # Axis labels
    axis(1, at = 1:ncol(mat_rev), labels = colnames(mat_rev), las = 2, cex.axis = 0.7)
    axis(2, at = 1:nrow(mat_rev), labels = rownames(mat_rev), las = 2, cex.axis = 0.7)
    box()

    # Add counts for small matrices
    if (nrow(mat) < 20 && ncol(mat) < 20) {
        grid_x <- slice.index(mat_rev, 2)
        grid_y <- slice.index(mat_rev, 1)
        text(grid_x, grid_y, labels = mat_rev, cex = 0.7, col = ifelse(mat_rev > max(mat) / 2, "white", "black"))
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
#' @export
plot_met_trend <- function(data, x = "Year", y = "Yield", main = "Yield Trend", ...) {
    if (!all(c(x, y) %in% names(data))) stop("Variables not found in data.")

    # Prepare Data Types
    raw_x <- data[[x]]

    # create numeric version for regression
    if (is.numeric(raw_x)) {
        x_num <- raw_x
    } else {
        # Try converting character/factor to numeric (assuming "2020", "2021" etc.)
        x_num <- suppressWarnings(as.numeric(as.character(raw_x)))
        if (all(is.na(x_num))) stop(paste("Column", x, "must be numeric (e.g. Year) for trend analysis."))
    }

    # create factor version for boxplot grouping
    x_fac <- as.factor(raw_x)

    # --- Regression Analysis ---
    # Fit model: y ~ x_numeric
    # Use a temp dataframe to avoid environment issues with lm formulas
    tmp_df <- data.frame(Y = data[[y]], X = x_num)
    model <- lm(Y ~ X, data = tmp_df)

    slope <- coef(model)["X"]
    r_sq <- summary(model)$r.squared

    # --- Visualisation ---
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    # 2 Rows (1 Column)
    layout(matrix(1:2, nrow = 2))

    # Plot 1: Boxplot (Distribution)
    par(mar = c(3, 5, 3, 2))
    boxplot(data[[y]] ~ x_fac,
        col = "lightgreen",
        border = "darkgreen",
        main = paste(main, "(Distribution)"),
        ylab = y,
        las = 1,
        ...
    )

    # Plot 2: Line Plot (Means + Trend)
    par(mar = c(5, 5, 2, 2))

    # Calculate Means using the numeric x values to correctly space them on x-axis
    # Note: tapply sorts by the factor levels. We need to match means to specific numeric years.
    means <- tapply(tmp_df$Y, tmp_df$X, mean, na.rm = TRUE)
    years_plot <- as.numeric(names(means))

    plot(years_plot, means,
        type = "b",
        pch = 19, lwd = 2, col = "blue",
        main = paste(main, "(Rate of Change)"),
        xlab = x, ylab = paste("Mean", y),
        xaxt = "n", # Custom axis
        las = 1
    )

    # Smart Axis
    axis(1, at = years_plot, labels = years_plot)
    grid()

    # Add Regression Line
    abline(model, col = "red", lwd = 2)

    # Stats Legend
    legend_text <- c(
        paste("Rate:", round(slope, 4), "unit/yr"),
        paste("R2:", round(r_sq, 3))
    )

    legend("topleft",
        legend = legend_text, bty = "n", cex = 0.9,
        text.col = "black"
    )

    layout(1)
}

#' Flexible Trial Map
#'
#' Plots the spatial layout of a specific trial. Adapts to missing spatial coordinates
#' by falling back to a 1D visualization if necessary.
#'
#' @param data A data.frame.
#' @param trial_val Value. The specific trial ID to filter (if data contains multiple).
#' @param trial_col String. Column name for Trial ID.
#' @param row String. Column name for Row.
#' @param col String. Column name for Column.
#' @param val String. Variable to fill heatmap (e.g. "Yield").
#'
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

        # Cast to Matrix
        # We need to ensure we fill the matrix correctly.
        # Rows are Y axis, Cols are X axis.
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

#' Convert Yield Units
#'
#' Converts yield from Bushels per Acre (bu/ac) to Tonnes per Hectare (t/ha).
#'
#' @param yield Numeric vector of yield values.
#' @param crop String. Crop type to determine bushel weight. Options:
#'   "wheat" (60 lbs), "soybeans" (60 lbs), "corn" (56 lbs),
#'   "barley" (48 lbs), "oats" (32 lbs).
#' @param lbs_per_bu Numeric. Custom weight per bushel (overrides crop).
#'
#' @return Numeric vector of yield in t/ha.
#' @export
convert_buac_to_tha <- function(yield, crop = "wheat", lbs_per_bu = NULL) {
    # Standard weights (lbs)
    weights <- c(
        "wheat" = 60,
        "soybeans" = 60,
        "soy" = 60,
        "corn" = 56,
        "maize" = 56,
        "barley" = 48,
        "oats" = 32
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

    # Conversion logic
    # 1 bu/ac = [weight] lbs/ac
    # 1 lb = 0.45359237 kg
    # 1 ac = 0.404685642 ha
    # 1 tonne = 1000 kg
    # factor = (weight * 0.45359237) / (0.404685642 * 1000)
    #        = weight * 0.00112085116

    factor <- weight * 0.00112085116
    return(yield * factor)
}

#' Check Experimental Design Factors by Study
#'
#' Scans a MET dataset to determine which design terms (Row, Column, Block, Replicate)
#' are valid for inclusion in a model for each specific study.
#'
#' @param data The dataframe containing the MET data.
#' @param study_col String. Name of the column defining the Study/Experiment (default "studyName").
#' @param row_col String. Name of the Row column (default "rowNumber").
#' @param col_col String. Name of the Column column (default "colNumber").
#' @param block_col String. Name of the Block column (default "blockNumber").
#' @param rep_col String. Name of the Replicate column (default "replicate").
#'
#' @return A dataframe summarizing which terms are available for each study.
#' @export
check_design_terms <- function(data,
                               study_col = "studyName",
                               row_col = "rowNumber",
                               col_col = "colNumber",
                               block_col = "blockNumber",
                               rep_col = "replicate") {
    # Ensure the study column exists
    if (!study_col %in% names(data)) {
        stop(paste("Study column", study_col, "not found in dataframe."))
    }

    # Split data by Study
    # usage of factor ensures we don't drop empty levels if they exist,
    # but here we usually only care about present data.
    study_list <- split(data, data[[study_col]])

    # Helper function to check if a column has > 1 unique value (ignoring NAs)
    has_var <- function(df, col_name) {
        if (!col_name %in% names(df)) {
            return(FALSE)
        }

        vals <- df[[col_name]]
        vals <- vals[!is.na(vals)]

        # Needs at least 2 levels to be a valid factor in a model
        # (e.g., if Block is always "1", you can't fit it as a random effect)
        return(length(unique(vals)) > 1)
    }

    # Iterate through studies
    results <- lapply(names(study_list), function(nm) {
        sub_df <- study_list[[nm]]

        data.frame(
            Study = nm,
            Has_Row = has_var(sub_df, row_col),
            Has_Col = has_var(sub_df, col_col),
            Has_Block = has_var(sub_df, block_col),
            Has_Rep = has_var(sub_df, rep_col),
            stringsAsFactors = FALSE
        )
    })

    # Combine results
    out_df <- do.call(rbind, results)
    rownames(out_df) <- NULL

    return(out_df)
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

        # Skip empty levels
        if (nrow(sub_df) == 0) {
            return(NULL)
        }

        # Get actual numeric coordinates (handle factors safely)
        r_vals <- as.numeric(as.character(sub_df[[row]]))
        c_vals <- as.numeric(as.character(sub_df[[col]]))

        # 1. Define the Expected Grid (Min to Max)
        # ASREml assumes the grid is defined by the range of the factors
        min_r <- min(r_vals, na.rm = TRUE)
        max_r <- max(r_vals, na.rm = TRUE)
        min_c <- min(c_vals, na.rm = TRUE)
        max_c <- max(c_vals, na.rm = TRUE)

        expected_grid <- expand.grid(
            Row_Num = seq(min_r, max_r),
            Col_Num = seq(min_c, max_c)
        )

        # Create a composite key "Row_Col" for easy comparison
        # We use the numeric values to match the grid
        sub_df$key <- paste(r_vals, c_vals, sep = "_")
        expected_grid$key <- paste(expected_grid$Row_Num, expected_grid$Col_Num, sep = "_")

        # 2. Find Missing Keys
        missing_keys <- setdiff(expected_grid$key, sub_df$key)

        if (length(missing_keys) > 0) {
            # Split key back into Row/Col
            coords <- do.call(rbind, strsplit(missing_keys, "_"))

            return(data.frame(
                Experiment = exp_id,
                Row = as.numeric(coords[, 1]),
                Column = as.numeric(coords[, 2]),
                Status = "MISSING",
                stringsAsFactors = FALSE
            ))
        } else {
            return(NULL) # No missing plots
        }
    })

    # Combine results
    out_df <- do.call(rbind, missing_list)

    if (is.null(out_df)) {
        message("Success: All experiments have complete rectangular grids.")
        return(invisible(NULL))
    } else {
        message(paste("Found", nrow(out_df), "missing plots across", length(unique(out_df$Experiment)), "experiments."))
        return(out_df)
    }
}
