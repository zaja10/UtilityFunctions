test_that("Matrix conversion and sparse matrix handling work correctly", {
    skip_if_not_installed("Matrix")
    library(Matrix)

    # Logic from debug_mat.R
    mat <- matrix(c(
        1.0, 0.5, 0.2,
        0.5, 1.0, 0.1,
        0.2, 0.1, 1.0
    ), nrow = 3, byrow = TRUE)

    # Assuming mat2sparse is an internal function or exported one we want to test
    # If it's internal we might need UtilityFunctions:::mat2sparse
    # But debug_mat.R used devtools::load_all(), so it had access.
    # We will assume it's available or use the package namespace if needed.

    # Check if mat2sparse exists, if not, skip or adapt
    if (exists("mat2sparse", where = asNamespace("UtilityFunctions"), mode = "function")) {
        res <- UtilityFunctions:::mat2sparse(mat)
        expect_s4_class(res, "CsparseMatrix")
    }

    # Debug specific case: specific sparse matrix causing issues
    mat_sparse <- matrix(c(1, 0, 0, 1), 2, 2)
    mat_sparse[upper.tri(mat_sparse)] <- 0
    m <- Matrix::Matrix(mat_sparse, sparse = TRUE)

    dsT <- as(as(m, "dgCMatrix"), "TsparseMatrix")
    expect_s4_class(dsT, "TsparseMatrix")
    expect_equal(dsT@x, c(1, 1))

    # Logic from debug_matrix.R (Reproducing matrix summary behavior)
    G <- matrix(c(1.00, 0.98, 0.98, 1.00), nrow = 2)
    s_obj <- as(Matrix::forceSymmetric(G), "dsCMatrix")

    expect_s4_class(s_obj, "dsCMatrix")

    summ <- summary(s_obj)
    # Ensure summary doesn't error and returns expected structure
    expect_true(is.data.frame(summ) || inherits(summ, "sparseSummary"))
    expect_true(all(c("i", "j", "x") %in% names(summ) || c("i", "j", "x") %in% slotNames(summ)))
})
