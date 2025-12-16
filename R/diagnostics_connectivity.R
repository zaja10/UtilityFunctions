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
calculate_connectivity <- function(data, x_fac = "Year", y_fac = NULL, trace = "Genotype", method = "count") {
    if (is.null(y_fac)) y_fac <- x_fac

    # Filter complete cases
    df <- data[!is.na(data[[x_fac]]) & !is.na(data[[trace]]), ]
    if (y_fac != x_fac) df <- df[!is.na(df[[y_fac]]), ]

    f1 <- as.factor(df[[x_fac]])
    f2 <- as.factor(df[[y_fac]])
    tr <- as.factor(df[[trace]])

    # Create incidence matrices
    M1 <- table(tr, f1)
    M1[M1 > 0] <- 1

    if (x_fac == y_fac) {
        # Symmetric
        conn <- crossprod(M1)

        if (method != "count") {
            sizes <- diag(conn)
            if (method == "jaccard") {
                denom <- outer(sizes, sizes, "+") - conn
                conn <- conn / denom
            } else if (method == "prop_min") {
                s_mat <- matrix(sizes, nrow = length(sizes), ncol = length(sizes))
                denom <- pmin(s_mat, t(s_mat))
                conn <- conn / denom
            }
        }
    } else {
        # Asymmetric
        M2 <- table(tr, f2)
        M2[M2 > 0] <- 1

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
#' @param ... Additional arguments passed to `image()`.
#'
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
#' @param ... Additional arguments passed to ggplot.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 labs theme_minimal theme element_text geom_text
#' @export
plot_connectivity <- function(data, x = "Year", y = NULL, trace = "Genotype", method = "count", order_by = NULL, ...) {
    mat <- calculate_connectivity(data, x, y, trace, method)

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
        get_order <- function(fac_name, current_levels) {
            if (fac_name == sort_col) {
                return(current_levels)
            }
            # Only aggregate if column exists
            if (!sort_col %in% names(data) || !fac_name %in% names(data)) {
                return(current_levels)
            }

            meta <- unique(data[, c(fac_name, sort_col)])
            meta <- meta[meta[[fac_name]] %in% current_levels, ]
            # Robust aggregation
            ord_df <- aggregate(as.formula(paste(sort_col, "~", fac_name)), data = meta, min)
            ordering <- ord_df[order(ord_df[[sort_col]], ord_df[[fac_name]]), fac_name]
            return(as.character(ordering))
        }

        row_ord <- intersect(get_order(x, rownames(mat)), rownames(mat))
        y_name <- if (is.null(y)) x else y
        # Check if y column exists in data before trying to order by it (might be same as x)
        col_ord <- intersect(get_order(y_name, colnames(mat)), colnames(mat))

        mat <- mat[row_ord, col_ord, drop = FALSE]
    }

    # Reshape for ggplot (Base R melt)
    df_long <- as.data.frame(as.table(mat))
    names(df_long) <- c("X", "Y", "Value")

    # Enforce Factor ordering from matrix
    df_long$X <- factor(df_long$X, levels = rownames(mat))
    df_long$Y <- factor(df_long$Y, levels = colnames(mat))

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data$Y, y = .data$X, fill = .data$Value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_fill_gradient2(low = "white", high = "steelblue", mid = "aliceblue", midpoint = max(mat) / 10) +
        ggplot2::labs(
            title = paste("Connectivity by", trace),
            subtitle = paste("Method:", method),
            x = ifelse(is.null(y) || x == y, x, y),
            y = x,
            fill = "Score"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = ggplot2::element_text(size = 8)
        )

    # Add text labels if small
    if (nrow(mat) < 20 && ncol(mat) < 20) {
        # Contrast color logic
        p <- p + ggplot2::geom_text(ggplot2::aes(label = round(Value, 1)),
            color = ifelse(df_long$Value > max(df_long$Value) / 2, "white", "black"),
            size = 3
        )
    }

    return(p)
}

#' Check Trial Network Connectivity (Wrapper)
#'
#' Diagnostic wrapper that calculates connectivity and produces a visual plot.
#'
#' @param data A dataframe containing the raw trial data.
#' @param genotype A character string. The column name for Genotype (default "Genotype").
#' @param trial A character string. The column name for Site/Environment (default "Site").
#' @param threshold Integer. The minimum number of shared lines required.
#'
#' @return A list containing the count matrix, pct matrix, and disconnects table.
#' @export
check_connectivity <- function(data, genotype = "Genotype", trial = "Site", threshold = 10) {
    # Calculate counts
    connect_mat <- calculate_connectivity(data, x_fac = trial, trace = genotype, method = "count")

    # Calculate pct
    pct_mat <- calculate_connectivity(data, x_fac = trial, trace = genotype, method = "jaccard") * 100

    # Identify disconnects
    mat_check <- connect_mat
    diag(mat_check) <- NA
    issues_idx <- which(mat_check < threshold, arr.ind = TRUE)

    disconnects <- NULL
    if (nrow(issues_idx) > 0) {
        disconnects <- data.frame(
            Site_A = rownames(connect_mat)[issues_idx[, 1]],
            Site_B = colnames(connect_mat)[issues_idx[, 2]],
            Shared_Count = connect_mat[issues_idx],
            Shared_Pct = round(pct_mat[issues_idx], 1)
        )
        disconnects <- disconnects[as.character(disconnects$Site_A) < as.character(disconnects$Site_B), ]
        warning(sprintf("Found %d pairs with < %d shared lines.", nrow(disconnects), threshold))
    }

    # Plot
    plot_connectivity(data, x = trial, trace = genotype, method = "count")

    return(list(matrix_count = connect_mat, matrix_pct = pct_mat, disconnects = disconnects))
}
