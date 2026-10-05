vc_two_terms <- c(
  "fa(Site, 2):vm(Genotype, Ginv)!S1!var" = 0.1, "fa(Site, 2):vm(Genotype, Ginv)!S1!fa1" = 0.8,
  "fa(Site, 2):gdrop!S1!var" = 0.2, "fa(Site, 2):gdrop!S1!fa1" = 0.3,
  "units!R" = 1
)

test_that("an ambiguous FA prefix is refused instead of merging terms", {
  rt <- UtilityFunctions:::resolve_fa_term
  expect_error(rt(names(vc_two_terms), NULL), "Several FA/RR terms")
  expect_error(rt(names(vc_two_terms), "fa(Site, 2)"), "matches 2")
  expect_equal(rt(names(vc_two_terms), "fa(Site,2):vm(Genotype, Ginv)"), "fa(Site, 2):vm(Genotype, Ginv)")
})

test_that("parameters are taken from the exact term only", {
  tab <- UtilityFunctions:::fa_parameter_table(vc_two_terms, "fa(Site, 2):vm(Genotype, Ginv)")
  expect_equal(tab$Value[tab$Param == "fa1"], 0.8) # not 0.8 + 0.3
  expect_equal(nrow(tab), 2)
})

test_that("FA terms are parsed", {
  s <- UtilityFunctions:::parse_fa_term("fa(studyName, 3):vm(gkeep, G.inv)")
  expect_equal(
    s[c("type", "site_var", "k", "fa_part", "rest", "geno_var")],
    list(type = "fa", site_var = "studyName", k = 3L, fa_part = "fa(studyName, 3)", rest = "vm(gkeep, G.inv)", geno_var = "gkeep")
  )
})

test_that("score prefixes are stripped exactly, keeping spaces in IDs", {
  sp <- UtilityFunctions:::strip_prefix
  x <- c("fa(Site, 2)_Comp1:vm(G, A)_AgriMAXX 490", "fa(Site, 2)_Comp10:vm(G, A)_X", "fa(Site, 2)_S1:vm(G, A)_Y")
  expect_equal(sp(x, "fa(Site, 2)_Comp1:vm(G, A)_"), c("AgriMAXX 490", NA, NA))
  expect_equal(sp("fa(Site,2)_Comp1:vm(G,A)_Z", "fa(Site, 2)_Comp1:vm(G, A)_"), "Z")
})

test_that("rotation arguments are normalised", {
  nr <- UtilityFunctions:::normalise_rotation
  expect_equal(nr(TRUE), "svd")
  expect_equal(nr(FALSE), "none")
  expect_message(expect_equal(nr("varimax"), "svd"), "deprecated")
  expect_error(nr("promax"))
})

test_that("FAST indices refuse unrotated solutions and use g_ii for reliability", {
  L <- matrix(c(1, 2, 3, 0.5, -0.5, 0), 3, dimnames = list(paste0("E", 1:3), c("Fac1", "Fac2")))
  Fs <- matrix(c(1, -1, 0.5, 0.2), 2, dimnames = list(c("a", "b"), c("Fac1", "Fac2")))
  P <- matrix(c(0.5, 0.5, 0.1, 0.1), 2, dimnames = dimnames(Fs))
  obj <- list(loadings = list(rotated = L), scores = list(rotated = Fs, pev_rotated = P), meta = list(rotation = "svd"))
  f <- calculate_fast_indices(obj, g_diag = c(a = 2, b = 1))
  f <- f[order(f$Genotype), ]
  expect_equal(f$OP, c(1, -1) * 2)
  expect_equal(f$Reliability, c(1 - 0.5 / 2, 1 - 0.5 / 1))
  I <- Fs[, 2, drop = FALSE] %*% t(L[, 2, drop = FALSE])
  expect_equal(f$RMSD, sqrt(rowMeans(I^2)), ignore_attr = TRUE)
  obj$meta$rotation <- "none"
  expect_error(calculate_fast_indices(obj), "SVD-rotated")
})

test_that("rank_genotypes applies the stability weight", {
  d <- data.frame(Genotype = c("a", "b"), OP = c(10, 9), RMSD = c(3, 0))
  expect_equal(rank_genotypes(d, 0)$Genotype, c("a", "b"))
  expect_equal(rank_genotypes(d, 1)$Genotype, c("b", "a"))
  expect_equal(rank_genotypes(d)$Rank, 1:2)
})
