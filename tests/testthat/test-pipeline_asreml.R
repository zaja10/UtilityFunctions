make_pipeline_data <- function() {
  g <- generate_grm(50, 300, seed = 21)$grm
  d <- generate_genomic_met(2, g, n_reps = 2, h2 = 0.6, seed = 21)
  ungeno <- rownames(g)[46:50]
  data.frame(
    studyName = as.character(d$Site),
    colNumber = as.integer(as.character(d$Column)),
    rowNumber = as.integer(as.character(d$Row)),
    gkeep = ifelse(d$Genotype %in% ungeno, NA, as.character(d$Genotype)),
    gdrop = ifelse(d$Genotype %in% ungeno, as.character(d$Genotype), NA),
    grain_yield = d$Yield,
    stringsAsFactors = FALSE
  ) -> out
  out$grain_yield[out$studyName == "Site_01"][7] <- out$grain_yield[out$studyName == "Site_01"][7] + 60
  list(data = out, G = g, Ginv = prepare_asreml_grm(g[1:45, 1:45]))
}

test_that("the univariate -> review -> MET pipeline runs end to end", {
  skip_if_no_asreml()
  p <- make_pipeline_data()
  out_dir <- withr::local_tempdir()
  tables <- file.path(out_dir, "tables")

  st <- suppressWarnings(suppressMessages(run_univariate_pipeline(p$data, p$Ginv, out_dir)))
  expect_equal(st$Status, c("done", "done"))
  expect_true(all(file.exists(file.path(out_dir, "diagnostics", paste0(c("Site_01", "Site_02"), "_Diagnostics.pdf")))))
  sel <- read.csv(file.path(tables, "Site_01_Spatial_Model_Selection.csv"))
  expect_true(all(c("Model", "AIC", "Converged", "Comparable") %in% names(sel)))
  met <- read.csv(file.path(tables, "Site_01_Final_Metrics.csv"))
  expect_true(all(c("h2", "vBLUP", "Gen_Reliability", "n_reps") %in% names(met)))
  expect_equal(met$n_reps, 2)
  expect_true(file.exists(file.path(tables, "Site_01_Flagged_Outliers.csv")))

  # Resume: completed trials are skipped
  st2 <- suppressMessages(run_univariate_pipeline(p$data, p$Ginv, out_dir))
  expect_equal(st2$Status, c("skipped", "skipped"))

  # Review: remove outliers in Site_01 only
  tmpl_path <- file.path(out_dir, "template.xlsx")
  tmpl <- suppressMessages(generate_diagnostic_template(tables, p$data, tmpl_path))
  expect_true("residual_terms" %in% names(tmpl))
  tmpl$remove_outliers[tmpl$studyName == "Site_01"] <- TRUE
  writexl::write_xlsx(tmpl, tmpl_path)

  st3 <- suppressWarnings(suppressMessages(run_modified_trials(p$data, tmpl_path, p$Ginv, out_dir)))
  expect_equal(st3$Status[st3$Trial == "Site_01"], "done")
  expect_equal(st3$Status[st3$Trial == "Site_02"], "unchanged")
  fs <- as.data.frame(readxl::read_xlsx(file.path(out_dir, "Final_Model_Summary.xlsx")))
  s1 <- fs[fs$studyName == "Site_01", ]
  expect_equal(s1$final_model, "MOD")
  expect_match(s1$selection_basis, "Reviewer")
  expect_true(s1$outliers_removed)
  expect_equal(fs$selection_basis[fs$studyName == "Site_02"], "AIC")
  sel1 <- read.csv(file.path(tables, "Site_01_Spatial_Model_Selection.csv"))
  expect_false(sel1$Comparable[sel1$Model == "MOD"])

  met_model <- suppressWarnings(suppressMessages(run_met_pipeline(
    p$data, p$Ginv, file.path(out_dir, "Final_Model_Summary.xlsx"),
    tables_dir = tables, save_dir = out_dir, trace = FALSE
  )))
  expect_s3_class(met_model, "asreml")
  expect_true(file.exists(file.path(out_dir, "MET_Step1_Diag_Model.rds")))
  # missing responses = padded grid gaps + the masked outlier(s)
  gaps <- sum(sapply(split(p$data, p$data$studyName), function(t) {
    length(unique(t$colNumber)) * length(unique(t$rowNumber)) - nrow(t)
  }))
  n_out <- nrow(read.csv(file.path(tables, "Site_01_Flagged_Outliers.csv")))
  expect_equal(sum(is.na(met_model$mf$Yield)), gaps + n_out)

  h <- calculate_met_h2(met_model, gkeep_term = "studyName:vm(gkeep, grm)")
  expect_equal(h$Trial, c("Site_01", "Site_02"))
  expect_true(all(h$h2 > 0 & h$h2 < 1))
})

test_that("extract_asreml_outliers aligns residuals to rows", {
  skip_if_no_asreml()
  p <- make_pipeline_data()
  d <- suppressMessages(UtilityFunctions:::prepare_trial_data(p$data[p$data$studyName == "Site_01", ], "grain_yield", p$Ginv))
  models <- UtilityFunctions:::build_model_list(d)
  fit <- UtilityFunctions:::fit_asreml(models$AxA$fixed, models$AxA$random, models$AxA$residual, data = d, env = UtilityFunctions:::model_env(p$Ginv))
  grm <- p$Ginv
  out <- extract_asreml_outliers(fit, d, threshold = 4)
  expect_gte(nrow(out), 1)
  # the spiked plot (largest raw value) is among the flagged rows
  expect_true(which.max(d$grain_yield) %in% as.integer(rownames(out)))
})
