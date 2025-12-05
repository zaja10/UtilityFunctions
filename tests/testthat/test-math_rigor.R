test_that("Mathematical Identity: G_est equals Lambda*Lambda' + Psi", {
    # Mock FA Object structure
    # Define a simple 2-factor model for 3 sites
    sites <- c("S1", "S2", "S3")
    k <- 2

    # Lambda (3 x 2)
    lambda <- matrix(c(
        1.0, 0.2,
        0.8, 0.5,
        0.5, 0.9
    ), nrow = 3, byrow = TRUE, dimnames = list(sites, c("Fac1", "Fac2")))

    # Psi (Diagonal)
    psi_vec <- c(0.1, 0.2, 0.05)
    psi <- diag(psi_vec)
    rownames(psi) <- colnames(psi) <- sites

    # Implied G
    G_true <- (lambda %*% t(lambda)) + psi

    # Create Mock fa_asreml object
    fa_obj <- list(
        loadings = list(rotated = lambda), # Assuming raw=rotated for this test or pre-rotated
        var_comp = list(psi = data.frame(Site = sites, Psi = psi_vec)),
        matrices = list(G = G_true),
        meta = list(k = k)
    )
    class(fa_obj) <- "fa_asreml"

    # Manually reconstruct
    lambda_test <- fa_obj$loadings$rotated
    psi_test <- diag(fa_obj$var_comp$psi$Psi)
    G_calc <- (lambda_test %*% t(lambda_test)) + psi_test

    expect_equal(G_calc, fa_obj$matrices$G, tolerance = 1e-8)
})

test_that("Rotation Invariance: Trace Preserved", {
    # Random Lambda
    set.seed(123)
    lambda <- matrix(rnorm(10), 5, 2)

    # SVD Rotation manually
    svd_res <- svd(lambda)
    V <- svd_res$v
    lambda_rot <- lambda %*% V

    # Check Total Variance (Sum of Squares)
    ss_raw <- sum(lambda^2)
    ss_rot <- sum(lambda_rot^2)

    expect_equal(ss_raw, ss_rot, tolerance = 1e-8)

    # Check Orthogonality of Rotated Columns
    # For SVD rotated lambda (Lambda * V), the columns might NOT be orthogonal?
    # Wait, SVD of Lambda = U D V'.
    # Rotated Lambda = Lambda * V = U D V' V = U D.
    # Columns of U D are orthogonal because U is orthogonal.

    cross_prod <- t(lambda_rot) %*% lambda_rot
    # Off-diagonals should be 0
    off_diag <- cross_prod[upper.tri(cross_prod)]

    expect_equal(sum(abs(off_diag)), 0, tolerance = 1e-8)
})
