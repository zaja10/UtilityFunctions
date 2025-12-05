test_that("FA Object has correct structure", {
  # Mock object
  fa_obj <- list(
    loadings = list(rotated = matrix(1:4, 2)),
    matrices = list(G = matrix(1:4, 2), Cor = matrix(1:4, 2)),
    meta = list(k = 2, type = "FA")
  )
  class(fa_obj) <- "fa_asreml"

  expect_s3_class(fa_obj, "fa_asreml")
  expect_true(is.list(fa_obj$loadings))
  expect_true(is.list(fa_obj$matrices))
})

test_that("Theory: Rotated Loadings are Orthogonal", {
  # This corresponds to Smith & Cullis (2018) requirement for PC rotation

  # Mock a Lambda that needs rotation
  # L = U D V'
  # We construct a known L
  U <- matrix(c(0.6, 0.8, -0.8, 0.6), 2, 2) # Orthogonal
  D <- diag(c(10, 5))
  V <- matrix(c(1, 0, 0, 1), 2, 2)
  Lambda <- U %*% D %*% t(V) # Simple

  # Apply SVD Rotation manually as done in core_extraction.r
  svd_res <- svd(Lambda)
  V_calc <- svd_res$v
  Lambda_rot <- Lambda %*% V_calc

  # Check Orthogonality: t(L_rot) %*% L_rot should be Diagonal
  inner <- t(Lambda_rot) %*% Lambda_rot
  off_diags <- inner[upper.tri(inner)]

  expect_equal(max(abs(off_diags)), 0, tolerance = 1e-8)
})
