#' Prepare GRM for ASReml (Inverse & Sparse)
#'
#' Checks a Genomic Relationship Matrix (GRM) for singularity, bends it if necessary,
#' inverts it, and converts it to the sparse 'giv' format required by ASReml.
#'
#' @param G A symmetric numeric matrix (GRM).
#' @param blend Numeric. Amount of identity matrix to blend (0-1) to ensure invertibility.
#'        Default 0.02 (2% bending).
#' @return A sparse dataframe (Row, Col, Value) with attribute 'INVERSE' = TRUE.
#' @importFrom MASS ginv
#' @importFrom Matrix forceSymmetric dsCMatrix
#' @importFrom cli cli_alert_warning cli_alert_info cli_alert_success
#' @export
prepare_asreml_grm <- function(G, blend = 0.02) {
    if (!isSymmetric(G)) {
        cli::cli_warn("GRM is not symmetric. Forcing symmetry.")
        G <- (G + t(G)) / 2
    }

    # Check dimensions
    n <- nrow(G)
    cli::cli_alert_info("Processing {n} x {n} GRM...")

    # Bending (Blending with I) to ensure positive definite
    # G_star = (1 - w) * G + w * I
    if (blend > 0) {
        cli::cli_alert_info("Applying {blend*100}% bending to ensure invertibility.")
        I <- diag(n)
        G <- ((1 - blend) * G) + (blend * I)
    }

    # Invert
    # Try standard solve first (faster), fall back to generalized inverse
    G_inv <- tryCatch(solve(G), error = function(e) {
        cli::cli_warn("Standard inversion failed. Using Generalized Inverse (slower).")
        MASS::ginv(G)
    })

    # Ensure Dimnames match original
    rownames(G_inv) <- rownames(G)
    colnames(G_inv) <- colnames(G)

    # Convert to Sparse Triplets (Row, Col, Value) for ASReml
    # Use upper triangle only
    G_inv[lower.tri(G_inv)] <- 0

    # Base R sparse conversion (to avoid heavy Matrix dep if desired, but Matrix is Imports)
    # Efficient extraction using Matrix package
    # Convert to TsparseMatrix (triplets) explicitly to safely access slots
    # forceSymmetric ensures it's treated as symmetric, casting to TsparseMatrix gives i, j, x
    s_obj <- as(Matrix::forceSymmetric(G_inv), "TsparseMatrix")

    # Extract triplets
    # Slots @i and @j are 0-based in Matrix packages
    out <- data.frame(
        Row = s_obj@i + 1,
        Col = s_obj@j + 1,
        Value = s_obj@x
    )

    # Return the inverse matrix object itself with 'INVERSE' attribute is most modern ASReml-R way
    # But often users need the .giv dataframe.
    # We will return the dataframe but attach the matrix as an attribute if needed?
    # Or just return the dataframe which is what ASReml needs for 'vm'.

    # Actually, usually for vm(Name, giv), 'giv' is the sparse matrix.
    # We'll return the dataframe with the INVERSE attribute as per request.

    attr(out, "INVERSE") <- TRUE

    cli::cli_alert_success("GRM inverted and prepared.")
    return(out)
}
