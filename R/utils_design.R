#' Diagnose Trial Design Structure
#'
#' @description
#' Analyzes the experimental design of a dataset to check for connectivity, balance,
#' and aliasing between structural factors (e.g., Reps nested within Sites).
#'
#' @param data A dataframe containing the trial data.
#' @param genotype Character. Column name for Genotype.
#' @param trial Character. Column name for Trial/Environment.
#' @param rep Character. Column name for Replicate/Block (Optional).
#'
#' @return A `design_diagnosis` object (S3 class) containing:
#' \describe{
#'   \item{structure}{Counts of levels for key factors.}
#'   \item{connectivity}{Connectivity matrix (shared genotypes) between sites.}
#'   \item{aliasing}{Check if Reps are nested within Sites (Recycled vs Unique coding).}
#'   \item{balance}{Summary of missing plots.}
#' }
#'
#' @export
diagnose_design <- function(data, genotype, trial, rep = NULL) {
    if (!all(c(genotype, trial) %in% names(data))) {
        stop("Specified columns not found in data.")
    }

    # 1. Structural Counts
    n_gen <- length(unique(data[[genotype]]))
    n_site <- length(unique(data[[trial]]))
    counts <- list(Genotypes = n_gen, Sites = n_site)

    # 2. Connectivity (Mini-Implementation)
    # Check disjoint sets
    inc_table <- table(data[[genotype]], data[[trial]])
    inc_mat <- as.matrix(inc_table)
    inc_mat[inc_mat > 0] <- 1
    con_mat <- t(inc_mat) %*% inc_mat

    # 3. Aliasing / Nesting Check
    nesting_status <- "Unknown"
    if (!is.null(rep) && rep %in% names(data)) {
        # Check if Rep levels are recycled across sites
        tbl <- table(data[[trial]], data[[rep]])
        # If Rep 1 exists in Site A and Site B, it is likely "nested" physically but coded same
        is_recycled <- all(colSums(tbl > 0) > 1)
        nesting_status <- if (is_recycled) "Recycled (Nested)" else "Unique (Crossed/Implicitly Nested)"
    }

    # 4. Construct Object
    obj <- list(
        data_summary = counts,
        connectivity = con_mat,
        nesting = nesting_status,
        call = list(genotype = genotype, trial = trial, rep = rep)
    )

    class(obj) <- "design_diagnosis"
    return(obj)
}

#' Print Design Diagnosis
#'
#' @param x A design_diagnosis object.
#' @param ... Additional arguments.
#' @export
print.design_diagnosis <- function(x, ...) {
    cat("=== Design Diagnosis ===\n")
    cat(sprintf("Genotypes: %d | Sites: %d\n", x$data_summary$Genotypes, x$data_summary$Sites))
    cat(sprintf("Rep Coding: %s\n", x$nesting))

    # Connectivity Warning
    disconnected <- which(x$connectivity == 0, arr.ind = TRUE)
    # Remove self-loops and duplicates
    disconnected <- disconnected[disconnected[, 1] < disconnected[, 2], , drop = FALSE]

    if (nrow(disconnected) > 0) {
        cat("\n[WARNING]: Use check_connectivity() to inspect disjoint sites!\n")
    } else {
        cat("Connectivity: Full network graph connected.\n")
    }
}

#' Plot Design Tableau
#'
#' @description
#' Visualizes the randomization design using the first two variables in the plot formula
#' (e.g., Site and Rep) to show coverage/balance.
#'
#' @param x A design_tableau object.
#' @param ... Additional arguments.
#' @export
plot.design_tableau <- function(x, ...) {
    # Implementation placeholder - will use graphics::image or similar
    # to visualize the contingency table of the design factors.

    pf_vars <- all.vars(x$formulas$random)
    if (length(pf_vars) < 2) {
        message("Need at least 2 random effect factors to plot design.")
        return(invisible(NULL))
    }

    v1 <- pf_vars[1]
    v2 <- pf_vars[2]

    tab <- table(x$data[[v1]], x$data[[v2]])

    image(tab,
        main = paste("Design Structure:", v1, "vs", v2),
        axes = FALSE, col = c("white", "steelblue")
    )
    axis(1, at = seq(0, 1, length.out = nrow(tab)), labels = rownames(tab), las = 2)
    axis(2, at = seq(0, 1, length.out = ncol(tab)), labels = colnames(tab), las = 1)
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
            meta <- unique(data[, c(fac_name, sort_col)])
            meta <- meta[meta[[fac_name]] %in% current_levels, ]
            ord_df <- aggregate(as.formula(paste(sort_col, "~", fac_name)), data = meta, min)
            ordering <- ord_df[order(ord_df[[sort_col]], ord_df[[fac_name]]), fac_name]
            return(as.character(ordering))
        }

        row_ord <- intersect(get_order(x, rownames(mat)), rownames(mat))
        y_name <- if (is.null(y)) x else y
        col_ord <- intersect(get_order(y_name, colnames(mat)), colnames(mat))
        mat <- mat[row_ord, col_ord, drop = FALSE]
    }

    # Plotting
    mat_rev <- mat[nrow(mat):1, , drop = FALSE]
    cols <- colorRampPalette(c("white", "aliceblue", "#3498db", "#2c3e50"))(20)

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mar = c(5, 5, 4, 3))
    image(1:ncol(mat_rev), 1:nrow(mat_rev), t(mat_rev),
        axes = FALSE, col = cols, xlab = ifelse(is.null(y) || x == y, x, y), ylab = x,
        main = paste("Connectivity by", trace)
    )

    axis(1, at = 1:ncol(mat_rev), labels = colnames(mat_rev), las = 2, cex.axis = 0.7)
    axis(2, at = 1:nrow(mat_rev), labels = rownames(mat_rev), las = 2, cex.axis = 0.7)
    box()

    if (nrow(mat) < 20 && ncol(mat) < 20) {
        grid_x <- slice.index(mat_rev, 2)
        grid_y <- slice.index(mat_rev, 1)
        text(grid_x, grid_y, labels = mat_rev, cex = 0.7, col = ifelse(mat_rev > max(mat) / 2, "white", "black"))
    }
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

