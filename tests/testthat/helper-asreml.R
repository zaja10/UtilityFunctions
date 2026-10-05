skip_if_no_asreml <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("asreml")
  ok <- tryCatch(
    {
      suppressMessages(loadNamespace("asreml"))
      TRUE
    },
    error = function(e) FALSE
  )
  if (!ok) testthat::skip("asreml could not be loaded (licence?)")
}

# Small genomic two-factor MET with known structure.
sim_fa_data <- function(n_gen = 60, n_sites = 6, seed = 1) {
  g <- generate_grm(n_gen, 400, seed = seed)$grm
  d <- generate_genomic_met(n_sites, g, n_reps = 2, h2 = 0.6, seed = seed)
  d$Genotype <- factor(as.character(d$Genotype), levels = rownames(g))
  list(data = d, G = g)
}
