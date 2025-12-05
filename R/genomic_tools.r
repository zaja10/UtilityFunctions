#' Validate Genomic Relationship Matrix (GRM) Scaling
#'
#' @description
#' Checks if a GRM is properly scaled according to the VanRaden method.
#' For valid heritability and reliability (Oakey) estimates, the diagonal of the GRM
#' should average approximately 1 (or 1 + f, where f is inbreeding).
#' Unscaled GRMs (e.g. centering only) can result in inflated or deflated variance estimates.
#'
#' @param grm A numeric matrix representing the GRM.
#' @param tolerance Numeric. Allowed deviation from 1. Default 0.2.
#'
#' @return Logical TRUE if valid, FALSE if warning issued.
#' @export
validate_grm <- function(grm, tolerance = 0.2) {
    cat("--- GRM Validation ---\n")

    if (!is.matrix(grm)) {
        stop("GRM must be a matrix.")
    }

    # 1. Check Diagonal
    diag_vals <- diag(grm)
    mean_diag <- mean(diag_vals, na.rm = TRUE)

    cat(sprintf("-> Mean Diagonal: %.4f\n", mean_diag))

    if (abs(mean_diag - 1) > tolerance) {
        cat("[WARNING] GRM diagonal deviates significantly from 1.\n")
        cat("          Variance estimates (Vg) will be inversely scaled.\n")
        cat("          Heritability (Oakey) estimates may be incorrect.\n")
        cat("          Recommendation: Rescale GRM using VanRaden Method 1.\n")
        return(invisible(FALSE))
    }

    # 2. Check Symmetry
    if (!isSymmetric(grm, tol = 1e-8)) {
        warning("GRM is not symmetric.")
        return(invisible(FALSE))
    }

    cat("-> GRM appears correctly scaled and valid.\n")
    return(invisible(TRUE))
}
