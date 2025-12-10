#' MET Dataset Diagnostics & Visualization
#'
#' @description
#' A suite of tools for inspecting large Multi-Environment Trial (MET) datasets
#' prior to complex analysis. Includes checking connectivity, temporal trends,
#' and spatial integrity.
#'
#' @name met_diagnostics
NULL

#' Pad Trial Layout for Spatial Analysis
#'
#' Ensures a trial dataset has a complete grid of Row x Column positions by padding
#' missing coordinates with NA values. Useful for checking experimental design integrity.
#'
#' @param data A data.frame containing the trial data.
#' @param row String. Column name for Row positions (Default: "Row").
#' @param col String. Column name for Column positions (Default: "Column").
#'
#' @return A data.frame with the original data plus rows for missing physical positions.
#' @export
pad_trial_layout <- function(data, row = "Row", col = "Column") {
    if (!all(c(row, col) %in% names(data))) {
        warning("Row or Column columns not found. Returning original data.")
        return(data)
    }

    # Ensure coordinates are numeric
    r_vals <- tryCatch(as.numeric(as.character(data[[row]])), warning = function(w) NULL)
    c_vals <- tryCatch(as.numeric(as.character(data[[col]])), warning = function(w) NULL)

    if (is.null(r_vals) || is.null(c_vals)) {
        warning("Row/Column coordinates must be coercible to numeric. Returning original.")
        return(data)
    }

    min_r <- min(r_vals, na.rm = TRUE)
    max_r <- max(r_vals, na.rm = TRUE)
    min_c <- min(c_vals, na.rm = TRUE)
    max_c <- max(c_vals, na.rm = TRUE)

    # Create complete grid
    grid <- expand.grid(
        Row_Temp = seq(min_r, max_r),
        Col_Temp = seq(min_c, max_c)
    )
    names(grid) <- c(row, col)

    # Ensure types match for join
    grid[[row]] <- as(grid[[row]], class(data[[row]]))
    grid[[col]] <- as(grid[[col]], class(data[[col]]))

    # Merge
    # Use base R merge (equivalent to left_join on the grid)
    out <- merge(grid, data, by = c(row, col), all.x = TRUE)
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
#' @param ... Additional arguments passed to `image()`.
#'
#' @export
plot_connectivity <- function(data, x = "Year", y = NULL, trace = "Genotype", method = "count", ...) {
    mat <- get_connectivity(data, x, y, trace, method)

    # Sorting Logic: Try to sort axes by Year if they aren't "Year" itself
    if ("Year" %in% names(data)) {
        # Helper to get order
        get_order <- function(fac_name, current_levels) {
            if (fac_name == "Year") {
                return(current_levels)
            } # Already sorted usually, or numeric

            # Find mapping of Factor -> Year (use median or min year)
            # We need to make sure we only use the data that went into the matrix?
            # Actually, just using global data is fine for finding the year.
            meta <- unique(data[, c(fac_name, "Year")])
            meta <- meta[meta[[fac_name]] %in% current_levels, ]

            # Aggregate if multiple years per factor level (unlikely for Trial, possible for others)
            # Use min year to sort
            ord_df <- aggregate(as.formula(paste("Year ~", fac_name)), data = meta, min)
            ordering <- ord_df[order(ord_df$Year, ord_df[[fac_name]]), fac_name]
            return(as.character(ordering))
        }

        # Sort Rows (X axis in get_connectivity, usually Y axis in image)
        # Note: get_connectivity returns Rows=X_fac, Cols=Y_fac if symmetric?
        # Wait, check get_connectivity:
        # M1 Rows=Trace, Cols=X. M2 Rows=Trace, Cols=Y.
        # crossprod(M1, M2) -> Rows=X, Cols=Y.

        # So Rows of Mat correspond to 'x' arg.
        row_ord <- get_order(x, rownames(mat))

        # Sort Cols (Y axis in get_connectivity)
        y_name <- if (is.null(y)) x else y
        col_ord <- get_order(y_name, colnames(mat))

        # Intersect to be safe (in case get_order missed something, rare)
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
#' Visualizes performance trends over time or across environments using boxplots.
#'
#' @param data A data.frame.
#' @param x String. Categorical variable for X-axis (e.g. "Year").
#' @param y String. Response variable (e.g. "Yield").
#' @param main String. Plot title.
#' @param ... Additional arguments passed to `boxplot`.
#'
#' @export
plot_met_trend <- function(data, x = "Year", y = "Yield", main = "Yield Trend", ...) {
    if (!all(c(x, y) %in% names(data))) stop("Variables not found in data.")

    form <- as.formula(paste(y, "~", x))

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mar = c(5, 5, 4, 2))
    boxplot(form,
        data = data,
        col = "lightblue",
        border = "navy",
        main = main,
        ylab = y,
        xlab = x,
        las = 2, # Vertical labels if many years
        ...
    )

    # Add trend line of means
    means <- tapply(data[[y]], data[[x]], mean, na.rm = TRUE)
    lines(seq_along(means), means, col = "red", lwd = 2, type = "o", pch = 19)
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
