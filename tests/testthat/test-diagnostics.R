test_that("missing genotypes do not bridge disconnected trials", {
  d <- data.frame(Trial = c("T1", "T1", "T2", "T2"), Genotype = c("A", NA, "B", NA))
  r <- suppressMessages(check_met_connectivity(d))
  expect_equal(r$n_clusters, 2)
  expect_equal(length(r$islands), 1)
})

test_that("connectivity reports shared counts and weak links", {
  d <- data.frame(
    Trial = rep(c("T1", "T2", "T3"), each = 6),
    Genotype = c(paste0("G", 1:6), paste0("G", 1:6), paste0("G", c(6, 7:11)))
  )
  r <- suppressMessages(check_met_connectivity(d, warn_below = 2))
  expect_equal(r$n_clusters, 1)
  expect_equal(r$shared["T1", "T2"], 6)
  expect_equal(r$shared["T1", "T3"], 1)
  expect_equal(r$weak, "T3")
  r2 <- suppressMessages(check_met_connectivity(d, min_shared = 2))
  expect_equal(r2$islands, "T3")
  expect_error(suppressMessages(check_met_connectivity(d, min_shared = 2, stop_if_disconnected = TRUE)), "Disconnected")
})

test_that("diagnose_design summarises trials", {
  d <- data.frame(T = rep(c("a", "b"), each = 4), G = c("x", "x", "y", "z", "x", "y", "y", "w"), R = rep(1:2, 4))
  dd <- diagnose_design(d, "G", "T", rep = "R")
  expect_s3_class(dd, "design_diagnosis")
  expect_equal(dd$counts$G, 4)
  expect_equal(dd$per_trial$Genotypes, c(3, 3))
  expect_equal(dd$connectivity["a", "b"], 2)
  expect_output(suppressMessages(print(dd)), "Prop_Replicated")
  expect_error(diagnose_design(d, "nope", "T"), "not found")
})

test_that("find_missing_plots names output after the input columns", {
  d <- data.frame(Exp = "E1", R = c(1, 1, 2), C = c(1, 2, 1))
  out <- find_missing_plots(d, "Exp", "R", "C")
  expect_equal(names(out), c("Exp", "R", "C"))
  expect_equal(unlist(out[1, c("R", "C")]), c(R = 2, C = 2))
  expect_null(suppressMessages(find_missing_plots(d[1:2, ], "Exp", "R", "C")))
})
