test_that("theme_genetics returns a valid theme", {
    t <- theme_genetics()
    expect_s3_class(t, "theme")
    expect_s3_class(t, "gg")
})

test_that("prepare_asreml_grm handles singular matrices", {
    # 1. Create Singular Matrix (Perfect correlation)
    G <- matrix(c(1, 1, 1, 1), nrow = 2)
    rownames(G) <- colnames(G) <- c("G1", "G2")

    # 2. Run with bending
    # 98% G + 2% I
    # I = diag(1, 1)
    # G_bend = 0.98*[1 1; 1 1] + 0.02*[1 0; 0 1]
    #        = [0.98 0.98; 0.98 0.98] + [0.02 0; 0 0.02]
    #        = [1.00 0.98; 0.98 1.00]
    # Determinant = 1 - 0.98^2 > 0. Valid.

    G_inv <- expect_warning(prepare_asreml_grm(G, blend = 0.02), regexp = NA) # No warning for singularity

    expect_true(attr(G_inv, "INVERSE"))
    # Updated: Now returns Matrix object, not dataframe
    expect_true(is.matrix(G_inv) || inherits(G_inv, "Matrix"))

    # Check if symmetric
    # For triplet, we might just check uniqueness of pairs if we enforced upper tri?
    # The function forces upper tri zeros for the sparse, but the dataframe has Row/Col.
    # Let's check dimensions suitable for ASReml (Row >= Col usually preferred or vice versa)
    # Function enforces upper triangle for matrix, then conversion.
    # summary() of symmetric matrix usually gives unique triplets?
})

test_that("prepare_asreml_grm handles non-symmetric input", {
    G <- matrix(c(1, 0.5, 0.6, 1), nrow = 2)
    rownames(G) <- colnames(G) <- c("A", "B")

    # Should warn
    expect_warning(prepare_asreml_grm(G, blend = 0), "Forcing symmetry")
})
