test_that("theme_genetics returns a valid theme", {
  t <- theme_genetics()
  expect_s3_class(t, "theme")
})

test_that("makeNm uses labels, not factor codes", {
  expect_equal(makeNm(c("1", "2.5", "3")), c(1, 2.5, 3))
  expect_equal(makeNm(factor(c("10", "20"))), c(10, 20))
  expect_warning(makeNm(c("A", "1")), "NAs introduced by coercion")
})

dense_from_sparse <- function(sp) {
  n <- length(attr(sp, "rowNames"))
  M <- matrix(0, n, n)
  M[sp[, 1:2]] <- sp[, 3]
  M[sp[, 2:1]] <- sp[, 3]
  M
}

test_that("prepare_asreml_grm returns the sparse inverse of the blended GRM", {
  G <- matrix(c(1, 1, 1, 1), nrow = 2, dimnames = list(c("G1", "G2"), c("G1", "G2")))
  expect_no_warning(G_inv <- prepare_asreml_grm(G, blend = 0.02))
  expect_true(attr(G_inv, "INVERSE"))
  expect_equal(attr(G_inv, "rowNames"), c("G1", "G2"))
  expect_equal(colnames(G_inv), c("row", "column", "value"))
  expect_true(all(G_inv[, "row"] >= G_inv[, "column"]))
  Gb <- 0.98 * G + 0.02 * diag(2)
  expect_equal(dense_from_sparse(G_inv), unname(solve(Gb)), tolerance = 1e-10)
})

test_that("prepare_asreml_grm symmetrises and rejects non-PD input", {
  G <- matrix(c(1, 0.5, 0.6, 1), nrow = 2, dimnames = list(c("A", "B"), c("A", "B")))
  expect_message(prepare_asreml_grm(G, blend = 0), "GRM is not symmetric")
  singular <- matrix(1, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(prepare_asreml_grm(singular, blend = 0), "not positive definite")
  expect_error(prepare_asreml_grm(unname(G)), "row names")
})
