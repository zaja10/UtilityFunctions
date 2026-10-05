test_that("choose_final_model: reviewer MOD wins only when not AIC-comparable and converged", {
  base <- data.frame(Model = c("RCD", "AxA"), AIC = c(100, 90), Converged = TRUE, Comparable = TRUE)
  ch <- UtilityFunctions:::choose_final_model
  mod_nc <- rbind(base, data.frame(Model = "MOD", AIC = 500, Converged = TRUE, Comparable = FALSE))
  expect_equal(ch(mod_nc)$model, "MOD")
  mod_nc$Converged[3] <- FALSE
  expect_equal(ch(mod_nc)$model, "AxA")
  mod_c <- rbind(base, data.frame(Model = "MOD", AIC = 95, Converged = TRUE, Comparable = TRUE))
  expect_equal(ch(mod_c)$model, "AxA")
  mod_c$AIC[3] <- 80
  expect_equal(ch(mod_c)$model, "MOD")
  # Older selection files without a Comparable column
  expect_equal(ch(base[, 1:3])$model, "AxA")
})

test_that("formula_key ignores whitespace and term order", {
  fk <- UtilityFunctions:::formula_key
  expect_equal(fk("~ a + vm(g, grm)"), fk("~vm(g,grm)+a"))
  expect_false(fk("Yield ~ 1") == fk("Yield ~ lin(rowNumber)"))
})

test_that("numeric_factor orders numerically", {
  f <- UtilityFunctions:::numeric_factor(c("10", "2", "1"))
  expect_equal(levels(f), c("1", "2", "10"))
})

test_that("prepare_trial_data handles dense GRMs and aborts when nothing matches", {
  g <- generate_grm(5, 50)$grm
  d <- data.frame(
    studyName = "T1", colNumber = c(1, 1, 2), rowNumber = c(1, 2, 1),
    gkeep = c("G001", "G002", NA), gdrop = c(NA, NA, "X1"), grain_yield = c(1, 2, 3)
  )
  out <- suppressMessages(UtilityFunctions:::prepare_trial_data(d, "grain_yield", g))
  expect_equal(sum(!is.na(out$gkeep)), 2)
  d$gkeep <- c("nope", "nada", NA)
  expect_error(suppressMessages(UtilityFunctions:::prepare_trial_data(d, "grain_yield", g)), "None of the")
})

test_that("MET equations follow per-trial decisions and quote names safely", {
  fs <- data.frame(
    studyName = c("A", "B'2"),
    fixed_terms = c("Yield ~ 1", "Yield ~ lin(rowNumber)"),
    random_terms = c("~ vm(gkeep, grm) + colNumber + rowNumber + units", "~ vm(gkeep, grm) + rowNumber"),
    residual_terms = c("~ ar1(colNumber):ar1(rowNumber)", "~ colNumber:rowNumber")
  )
  eq <- UtilityFunctions:::build_met_equations(fs)
  expect_match(eq$fixed, 'at(studyName, "B\'2"):rowNum_num', fixed = TRUE)
  expect_match(eq$random, 'at(studyName, "A"):units', fixed = TRUE)
  expect_false(grepl('"B\'2"):colNumber', eq$random, fixed = TRUE))
  expect_match(eq$residual, 'ar1(colNumber):ar1(rowNumber) | studyName, levels = c("A")', fixed = TRUE)
  expect_match(eq$residual, 'id(colNumber):id(rowNumber) | studyName, levels = c("B\'2")', fixed = TRUE)
  expect_no_error(stats::as.formula(eq$fixed))
  expect_no_error(stats::as.formula(eq$residual))
})

test_that("pipeline validates required columns", {
  expect_error(UtilityFunctions:::validate_pipeline_data(data.frame(a = 1), "grain_yield"), "missing required")
})

test_that("residual_variance handles scale-only and idv-type residuals", {
  rv <- UtilityFunctions:::residual_variance
  expect_equal(rv(c(g = 1, "colNumber:rowNumber!R" = 2.5, "colNumber:rowNumber!colNumber!cor" = 0.4, "colNumber:rowNumber!rowNumber!cor" = 0.2)), 2.5)
  expect_equal(rv(c(g = 1, "units!units" = 3, "units!R" = 1)), 3)
  expect_error(rv(c("Trait_A!R" = 1, "Trait_B!R" = 2)), "residual sections")
})
