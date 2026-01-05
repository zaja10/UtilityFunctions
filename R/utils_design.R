#' Diagnose Trial Design Structure
#'
#' Check using `diagnose_design(data, "Geno", "env")`.
#' @export
diagnose_design <- function(data, genotype, trial, rep = NULL) {
    if (!all(c(genotype, trial) %in% names(data))) stop("Columns missing")

    n_gen <- length(unique(data[[genotype]]))
    n_site <- length(unique(data[[trial]]))

    # Connectivity
    inc <- table(data[[genotype]], data[[trial]])
    inc[inc > 0] <- 1
    # Check disjoint

    obj <- list(
        counts = list(G = n_gen, S = n_site),
        connectivity = t(inc) %*% inc
    )
    class(obj) <- "design_diagnosis"
    obj
}

#' Convert Factor to Numeric safely
#'
#' @param df Vector or column to convert.
#' @return Numeric vector.
#' @export
makeNm <- function(df) {
    as.numeric(as.character(df))
}

#' Find Missing Plots
#'
#' Identifies gaps in rectangular grids (Row x Col).
#' @export
find_missing_plots <- function(data, experiment = "Experiment", row = "Row", col = "Column") {
    if (!all(c(experiment, row, col) %in% names(data))) stop("Cols missing")

    out <- do.call(rbind, lapply(unique(data[[experiment]]), function(id) {
        sub <- data[data[[experiment]] == id, ]
        if (nrow(sub) == 0) {
            return(NULL)
        }

        r <- as.numeric(sub[[row]])
        c <- as.numeric(sub[[col]])

        grid <- expand.grid(R = seq(min(r), max(r)), C = seq(min(c), max(c)))
        grid$Key <- paste(grid$R, grid$C)
        sub$Key <- paste(r, c)

        miss <- setdiff(grid$Key, sub$Key)
        if (length(miss) > 0) {
            coords <- do.call(rbind, strsplit(miss, " "))
            data.frame(Experiment = id, Row = as.numeric(coords[, 1]), Column = as.numeric(coords[, 2]))
        } else {
            NULL
        }
    }))

    if (is.null(out)) message("Complete grid.")
    out
}
