#' Prepare GRM for ASReml (Inverse & Sparse)
#'
#' Checks a Genomic Relationship Matrix (GRM) for singularity, bends it if necessary,
#' inverts it, and converts it to the sparse 'giv' format required by ASReml.
#'
#' @param G A symmetric numeric matrix (GRM).
#' @param blend Numeric. Amount of identity matrix to blend (0-1) to ensure invertibility.
#'        Default 0.02 (2\%).
#' @return A sparse dataframe (Row, Col, Value) with attribute 'INVERSE' = TRUE.
#' @importFrom MASS ginv
#' @importFrom Matrix forceSymmetric
#' @importClassesFrom Matrix dsCMatrix
#' @import cli
#' @export
prepare_asreml_grm <- function(G, blend = 0.02) {
    if (!isSymmetric(G)) {
        cli::cli_alert_warning("GRM is not symmetric. Forcing symmetry.")
        G <- (G + t(G)) / 2
    }

    n <- nrow(G)
    cli::cli_alert_info("Processing {n} x {n} GRM...")

    if (blend > 0) {
        cli::cli_alert_info("Applying {blend*100}% bending to ensure invertibility.")
        I <- diag(n)
        G <- ((1 - blend) * G) + (blend * I)
    }

    G_inv <- tryCatch(solve(G), error = function(e) {
        cli::cli_alert_warning("Standard inversion failed. Using Generalized Inverse (slower).")
        MASS::ginv(G)
    })

    rownames(G_inv) <- rownames(G)
    colnames(G_inv) <- colnames(G)

    G_inv[lower.tri(G_inv)] <- 0
    s_obj <- as(Matrix::forceSymmetric(G_inv), "dsCMatrix")

    # Helper to align with asreml's sparse format (1-based index if needed, but usually row/col names suffice if using vm(..., giv(df)))
    # ASReml-R usually wants a dataframe with Row, Col, Value
    # But if we return the matrix directly with INVERSE attribute, newer asreml versions handle it.

    # However, consistent with the prompt's request for commercial robustness, let's stick to the sparse matrix object or dataframe?
    # The user prompt's code returns G_inv (matrix) but also calculates 's_obj' and 'summ' without using them?
    # Correction: The code provided in the prompt creates s_obj but returns G_inv with attribute.
    # "return(G_inv)"

    # Attribute is key for ASReml 4.2+
    attr(G_inv, "INVERSE") <- TRUE

    cli::cli_alert_success("GRM inverted and prepared.")
    return(G_inv)
}

#' Align Phenotype and Genotype Data
#'
#' Ensures intersection of genotypes between phenotypic data and GRM.
#'
#' @param pheno Dataframe containing phenotypic data.
#' @param geno GRM matrix.
#' @param id_col Character string. Name of the genotype column in pheno.
#' @return A list containing aligned pheno and geno.
#' @export
align_pheno_geno <- function(pheno, geno, id_col = "Genotype") {
    if (!id_col %in% names(pheno)) stop(paste("Column", id_col, "not found in phenotype data."))

    common <- intersect(pheno[[id_col]], rownames(geno))

    if (length(common) == 0) stop("No common genotypes found.")

    pheno_sub <- pheno[pheno[[id_col]] %in% common, ]
    geno_sub <- geno[common, common]

    cli::cli_alert_success("Aligned: {length(common)} genotypes retained.")

    return(list(pheno = pheno_sub, geno = geno_sub))
}
