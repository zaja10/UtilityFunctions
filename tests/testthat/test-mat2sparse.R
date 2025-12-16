test_that("mat2sparse correctly converts symmetric matrix to triplet format", {
    # Create small symmetric matrix
    mat <- matrix(c(
        1.0, 0.5, 0.2,
        0.5, 1.0, 0.1,
        0.2, 0.1, 1.0
    ), nrow = 3, byrow = TRUE)

    res <- mat2sparse(mat)

    expect_equal(nrow(res), 6) # Lower triangle including diagonal
    expect_equal(names(res), c("Row", "Col", "Value"))

    # Check content
    # Row 1: (1,1) -> 1.0
    r1 <- res[res$Row == 1 & res$Col == 1, ]
    expect_equal(r1$Value, 1.0)

    # Row 2: (2,1) -> 0.5, (2,2) -> 1.0
    r2 <- res[res$Row == 2 & res$Col == 1, ]
    expect_equal(r2$Value, 0.5)

    # Check Strict Lower Triangle condition
    expect_true(all(res$Row >= res$Col))

    # Verify memory optimization (sparse) logic
    # Set a zero and see if it's potentially handled?
    # Current logic keeps explicit zeros if they are in the lower triangle unless thresholded BEFORE calling.
    # But mat2sparse logic calls Matrix(..., sparse=TRUE) which drops zeros.

    mat_sparse <- matrix(c(1, 0, 0, 1), 2, 2)
    res_sparse <- mat2sparse(mat_sparse)

    # Should have 2 rows (1,1) and (2,2). (2,1) is zero and should be dropped by Matrix conversion
    expect_equal(nrow(res_sparse), 2)
    expect_equal(res_sparse$Value, c(1, 1))
})
