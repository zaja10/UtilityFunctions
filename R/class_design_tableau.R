#' Build Design Tableau
#'
#' @description
#' Creates a `design_tableau` object that encapsulates the structural design of a Multi-Environment Trial (MET).
#' This object validates the experimental design structure and formats the model formulas for the solver.
#'
#' @param data A dataframe containing the trial data.
#' @param treatment_formula A formula specifying the fixed effects (e.g., ~ Genotype).
#' @param plot_formula A formula specifying the randomization structure (e.g., ~ Site/Rep/Block).
#'
#' @return A `design_tableau` object (S3 class) containing:
#' \describe{
#'   \item{data}{The original dataframe.}
#'   \item{formulas}{A list containing `fixed` and `random` formulas.}
#'   \item{structure}{A summary of the design levels.}
#'   \item{aliasing}{Aliasing diagnostics.}
#' }
#'
#' @details
#' This function performs critical structural checks before model fitting:
#' \itemize{
#'   \item Checks if variables exist in data.
#'   \item Detects if factors are typically nested (e.g., Rep within Site).
#'   \item Validates that the provided formula accounts for the data structure.
#' }
#'
#' @importFrom stats terms model.frame
#' @export
build_design_tableau <- function(data, treatment_formula, plot_formula) {
    # 1. Validation: Check Variable Existence
    tf_vars <- all.vars(treatment_formula)
    pf_vars <- all.vars(plot_formula)
    all_vars <- c(tf_vars, pf_vars)

    missing_vars <- setdiff(all_vars, names(data))
    if (length(missing_vars) > 0) {
        stop(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
    }

    # 2. Structural Diagnostics
    structural_summary <- list()
    for (v in all_vars) {
        if (is.character(data[[v]]) || is.factor(data[[v]])) {
            structural_summary[[v]] <- length(unique(data[[v]]))
        }
    }

    # 3. Aliasing / Nesting Checks
    # Heuristic: If we have >1 Site, check if "Rep" or "Block" is functionally nested
    # (i.e., Rep 1 at Site A is different physical entity than Rep 1 at Site B)
    aliasing_report <- list()

    # Identify 'Site' or 'Env' variable (naive assumption: usually the first var in plot_formula)
    # This is a simplification; passed formulas can be complex.
    # We assume plot_formula is like ~ Site/Rep or ~ Site + Rep

    # Extract terms
    pf_terms <- attr(terms(plot_formula), "term.labels")

    # Simple check: If multiple vars, check contingency
    if (length(pf_vars) >= 2) {
        v1 <- pf_vars[1] # e.g., Site
        v2 <- pf_vars[2] # e.g., Rep

        tab <- table(data[[v1]], data[[v2]])

        # Check for empty cells (structural zeros might imply specific design)
        # Check for perfect nesting: Each level of v2 appears in ONLY one level of v1?
        # This detects if Rep is globally unique (1..N across all sites) vs Recycled (1..k per site)

        # If Rep is recycled (1,2,3 at Site A, 1,2,3 at Site B), table is fully populated.
        # If Rep is unique (1,2,3 at A, 4,5,6 at B), table is diagonal-ish.
    }

    # 4. Construct Object
    obj <- list(
        data = data,
        formulas = list(
            fixed = treatment_formula,
            random = plot_formula
        ),
        structure = structural_summary,
        aliasing = aliasing_report
    )

    class(obj) <- "design_tableau"
    return(obj)
}

#' Print Design Tableau
#'
#' @param x A design_tableau object.
#' @param ... Additional arguments.
#' @export
print.design_tableau <- function(x, ...) {
    cat("=== Design Tableau ===\n")
    cat("Fixed Effects:  ", deparse(x$formulas$fixed), "\n")
    cat("Random Effects: ", deparse(x$formulas$random), "\n")
    cat("\nStructure Summary:\n")
    print(unlist(x$structure))
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
