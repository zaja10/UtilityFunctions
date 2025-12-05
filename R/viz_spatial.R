#' Plot Spatial Field Map (Robust)
#'
#' Maps raw data (from dataframe) or residuals (from model) to the physical field layout.
#' Visualizes missing data with a distinct background color.
#'
#' @param input An `asreml` model object OR a `data.frame`.
#' @param row Column name for Field Row (default "Row").
#' @param col Column name for Field Column (default "Column").
#' @param attribute String. What to plot?
#'        - If input is data.frame: The column name (e.g. "Yield").
#'        - If input is asreml: "residuals" (default) or "fitted".
#' @export
plot_spatial <- function(input, row = "Row", col = "Column", attribute = "Yield") {
    if (inherits(input, "asreml")) {
        if (attribute == "Yield") attribute <- "residuals"
        data_name <- as.character(input$call$data)
        if (!exists(data_name)) stop("Original dataframe not found.")
        df <- get(data_name)

        if (attribute == "fitted") {
            vals <- fitted(input)
            main_title <- "Spatial Map: Fitted Values"
        } else {
            vals <- resid(input)
            main_title <- "Spatial Map: Residuals"
        }

        if (length(vals) == nrow(df)) {
            df$Value_To_Plot <- vals
        } else {
            warning("Residual mismatch. Padding with NA.")
            df$Value_To_Plot <- rep(NA, nrow(df))
            df$Value_To_Plot[1:length(vals)] <- vals
        }
    } else if (inherits(input, "data.frame")) {
        df <- input
        if (is.null(df[[attribute]])) stop(paste("Column", attribute, "not found."))
        vals <- df[[attribute]]
        df$Value_To_Plot <- vals
        main_title <- paste("Spatial Map:", attribute)
    } else {
        stop("Input must be an asreml model or a dataframe.")
    }

    r_vals <- as.numeric(as.character(df[[row]]))
    c_vals <- as.numeric(as.character(df[[col]]))
    n_r <- max(r_vals, na.rm = TRUE)
    n_c <- max(c_vals, na.rm = TRUE)

    field_mat <- matrix(NA, nrow = n_r, ncol = n_c)
    for (i in 1:nrow(df)) {
        if (!is.na(r_vals[i]) & !is.na(c_vals[i])) field_mat[r_vals[i], c_vals[i]] <- df$Value_To_Plot[i]
    }

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    layout(matrix(1:2, ncol = 2), widths = c(4, 1))
    par(mar = c(3, 3, 3, 1))

    is_diverging <- (min(vals, na.rm = TRUE) < 0 && max(vals, na.rm = TRUE) > 0)
    cols <- if (is_diverging) hcl.colors(20, "Blue-Red 3") else hcl.colors(20, "Spectral", rev = TRUE)

    field_rev <- field_mat[n_r:1, ]
    image(1:n_c, 1:n_r, t(field_rev), axes = FALSE, col = NA, main = main_title, xlab = "", ylab = "")
    rect(0.5, 0.5, n_c + 0.5, n_r + 0.5, col = "grey90", border = NA)
    image(1:n_c, 1:n_r, t(field_rev), add = TRUE, col = cols)
    box()

    # Plot Scale Bar
    par(mar = c(3, 0, 3, 3))
    image(1, 1:20, t(as.matrix(1:20)), axes = FALSE, xlab = "", ylab = "", col = cols)

    data_range <- range(vals, na.rm = TRUE)
    min_v <- data_range[1]
    max_v <- data_range[2]

    pretty_vals <- pretty(data_range, n = 5)
    pretty_vals <- pretty_vals[pretty_vals >= min_v & pretty_vals <= max_v]

    if (max_v > min_v) {
        at_locs <- (pretty_vals - min_v) / (max_v - min_v) * 19 + 1
    } else {
        at_locs <- 10
        pretty_vals <- unique(c(min_v, max_v))[1]
    }

    axis(4, at = at_locs, labels = round(pretty_vals, 2), las = 1)
    mtext("Grey = Missing", side = 1, line = 1, cex = 0.6)
    layout(1)
}
