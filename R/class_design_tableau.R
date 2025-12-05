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
