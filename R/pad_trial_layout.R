#' Pad Trial Layout for Spatial Analysis (MET-Ready)
#'
#' @description Prepares a dataset for spatial analysis by ensuring a complete
#' Row x Column grid. Missing plots are filled with NAs.
#'
#' For Multi-Environment Trials (METs), specify `group_cols` to pad each
#' environment independently.
#'
#' @param data A data.frame containing the trial data.
#' @param row_col Character string. Name of the Row column.
#' @param col_col Character string. Name of the Column column.
#' @param group_cols Character vector. Columns identifying unique trials (e.g. c("Site", "Year")).
#' The function will split data by these columns, pad each subset, and combine them.
#'
#' @return An object of class `padded_trial` (inherits from data.frame).
#' @export
pad_trial_layout <- function(data, row_col = "Row", col_col = "Column", group_cols = NULL) {
    # Silence R CMD check notes for .data in tidyverse
    utils::globalVariables(c(".data"))

    # --- 1. Input Validation ---
    req_cols <- c(row_col, col_col, group_cols)
    if (!all(req_cols %in% names(data))) {
        missing <- req_cols[!req_cols %in% names(data)]
        stop(paste("Columns not found in data:", paste(missing, collapse = ", ")))
    }

    # Ensure Row/Col are numeric
    if (is.factor(data[[row_col]])) data[[row_col]] <- as.numeric(as.character(data[[row_col]]))
    if (is.factor(data[[col_col]])) data[[col_col]] <- as.numeric(as.character(data[[col_col]]))

    # --- 2. Helper Function: Pad Single Trial ---
    pad_single <- function(sub_data) {
        # Determine grid extents for THIS specific trial
        r_range <- range(sub_data[[row_col]], na.rm = TRUE)
        c_range <- range(sub_data[[col_col]], na.rm = TRUE)

        # Create perfect grid
        full_grid <- expand.grid(
            Row = seq(r_range[1], r_range[2]),
            Col = seq(c_range[1], c_range[2])
        )
        names(full_grid) <- c(row_col, col_col)

        # Preserve grouping info (e.g. if this is Site A, put "Site A" in the new rows)
        if (!is.null(group_cols)) {
            for (grp in group_cols) {
                # We take the first value because they are all identical in this subset
                full_grid[[grp]] <- sub_data[[grp]][1]
            }
        }

        # Merge to keep data
        merged <- merge(full_grid, sub_data, by = c(row_col, col_col, group_cols), all.x = TRUE)
        return(merged)
    }

    # --- 3. Execution Logic (Single vs MET) ---

    if (is.null(group_cols)) {
        # Scenario A: Single Trial (No grouping)
        padded_data <- pad_single(data)
        n_added <- nrow(padded_data) - nrow(data)
    } else {
        # Scenario B: MET (Split - Apply - Combine)

        # Create a splitting factor based on all group columns
        split_fac <- interaction(data[group_cols], drop = TRUE)
        data_split <- split(data, split_fac)

        # Apply padding to each list element
        padded_list <- lapply(data_split, pad_single)

        # Combine back into one dataframe
        padded_data <- do.call(rbind, padded_list)
        rownames(padded_data) <- NULL # Clean row names

        n_added <- nrow(padded_data) - nrow(data)
    }

    # --- 4. Diagnostics & Attributes ---
    sparsity <- n_added / nrow(padded_data)

    attr(padded_data, "sparsity") <- sparsity
    attr(padded_data, "n_added") <- n_added
    attr(padded_data, "is_met") <- !is.null(group_cols)
    attr(padded_data, "coords") <- c(row = row_col, col = col_col)

    class(padded_data) <- c("padded_trial", "data.frame")

    return(padded_data)
}

#' Plot Trial Layout
#'
#' Usage: plot(padded_object, fill_col = "Yield")
#'
#' @param x Object of class `padded_trial`.
#' @param fill_col Character. Column to map to fill color.
#' @param ... Unused.
#' @export
plot.padded_trial <- function(x, fill_col = NULL, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plotting.")

    row_name <- attr(x, "coords")["row"]
    col_name <- attr(x, "coords")["col"]

    # Helper to find first non-structural column to detect missingness
    candidate_cols <- names(x)[!names(x) %in% c(row_name, col_name)]
    check_col <- candidate_cols[1]

    x$Status <- ifelse(is.na(x[[check_col]]),
        "Missing (Padded)", "Data Present"
    )

    p <- ggplot2::ggplot(x, ggplot2::aes(x = .data[[col_name]], y = .data[[row_name]])) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Trial Layout Integrity") +
        ggplot2::coord_fixed()

    # If it's a MET, we MUST facet by the grouping columns
    if (attr(x, "is_met")) {
        p <- p + ggplot2::facet_wrap(~., scales = "free")
    }

    if (!is.null(fill_col)) {
        p <- p + ggplot2::geom_tile(ggplot2::aes(fill = .data[[fill_col]]), color = "grey90")
    } else {
        p <- p + ggplot2::geom_tile(ggplot2::aes(fill = Status), color = "white") +
            ggplot2::scale_fill_manual(values = c(
                "Data Present" = "steelblue",
                "Missing (Padded)" = "firebrick"
            ))
    }

    return(p)
}
