# UtilityFunctions 1.2.0

A correctness and robustness release. Several results from 1.1.0 were
statistically wrong; re-run analyses that relied on the items marked **(results change)**.

## Packaging

* The package installs and passes `R CMD check` without ASReml-R:
  `asreml` and `ASExtras4` moved to Suggests, `readxl`/`writexl` are declared,
  and unused `MASS`, `Matrix` and `stringr` dependencies were dropped.
* Whole-namespace imports and the `globalVariables()` list were removed; data
  masking now uses `.data`/`.env` explicitly.
* Every asreml fit, `update()` and `predict()` runs in a dedicated
  environment, so user-named matrices (e.g. `G.inv`) and data frames resolve
  from inside package functions.
* Tests: the suite was rewritten (unit tests plus asreml integration tests that
  skip when asreml is unavailable). CI runs on Linux and macOS without asreml.

## Statistical fixes

* `fa.asreml()` **(results change)**: parameters are matched to the exact FA
  term; two FA terms sharing a prefix (e.g. `fa(Site, 2):vm(gkeep, G)` and
  `fa(Site, 2):gdrop`) are no longer silently summed. Scores are matched by exact
  coefficient prefix (`Comp1` no longer matches `Comp10`; IDs keep their spaces).
* Reliability, \(r^2_g\), is `1 - PEV / g_ii` using the new `grm` argument,
  instead of assuming unit score variance **(results change for genomic models)**.
  Negative values are now reported, not silently clipped.
* FA scores and their PEVs are rescaled by ASReml's `sigma2`. Under the sigma
  parameterisation (no explicit `residual`), ASReml reports score coefficients
  multiplied by \(\sqrt{\sigma^2}\) while the loadings are not, so OP and RMSD
  were inflated by \(\sqrt{\sigma^2}\) and reliabilities were wrong (often
  negative) **(results change for models without an explicit residual)**.
  Rankings at a fixed stability weight were unaffected. Results are now
  identical whichever way the residual is specified.
* `extract_variance_components()` detects the residual variance of
  `idv()`-type residuals (`units!R` fixed at 1, variance in `units!units`)
  instead of returning 1.
* `var_explained` is a percentage of total genetic variance (including
  specific variances) **(results change)**.
* `rotate = "varimax"`/`"mean"` were the same SVD rotation; both are now
  deprecated aliases of `"svd"`. FAST indices are not computed for unrotated
  solutions.
* `calculate_i_classes()` implements Smith et al. (2021) interaction classes
  (environment classes from loading signs, within-class genotype performance);
  the unused `alpha` argument was removed. **(output format changes)**
* `evaluate_gs_predictions()` now uses `fa.asreml()` (one FAST
  implementation), takes `weight_stability`, no longer writes files by default
  (`save_output = FALSE`) and resolves training/target overlaps as training.
* The Cullis et al. (2006) H² is now reported as a generalised reliability,
  \(\bar r^2_g\), in column `Gen_Reliability` (formerly `GH_c`), with
  `Acc` = \(\sqrt{\bar r^2_g}\). It uses the exact mean squared SED (`vBLUP`)
  from the prediction error variance matrix instead of the squared mean SED,
  which was biased upwards **(results change)**. Entry-mean h² uses the observed harmonic-mean
  replication instead of a hard-coded 2.
* `run_final_analysis()` reads the nugget from the fitted model, not the model
  name, so reviewer ("MOD") models with `units` no longer drop it
  **(results change)**.
* Model selection ranks converged models first; `mkConv()` warns when it gives up.
* `run_modified_trials()`: MOD models that remove outliers or change fixed
  effects are flagged `Comparable = FALSE` and adopted as reviewer decisions
  rather than ranked by an invalid AIC comparison **(results change)**.
  Unchanged template rows no longer refit a duplicate model.
* `run_met_pipeline()` uses each trial's selected residual structure (AR1×AR1
  or independent), quotes trial names safely, accepts `...` for
  `workspace`, and raises errors instead of returning `NULL`.
* Pipeline: trials are padded to complete row × column grids (required by
  AR1×AR1); dense GRMs are supported (previously every genotype became `NA`).
* `calculate_genetic_gain()` takes `accuracy` and `cycle_length` (the breeder's
  equation, warning when accuracy is assumed to be 1), an explicit
  `varcomp_term`, numeric-or-factor years, and check-adjusted realized gain.
* `calculate_met_h2()` matches components by trial name and, by default
  (`disjoint = TRUE`), no longer adds the variances of the disjoint genotyped
  and ungenotyped line sets **(results change)**.
* `check_met_connectivity()` ignores missing genotypes (which previously
  connected unrelated trials) and reports shared-genotype counts and weakly
  linked trials.
* Spectral indices: GARI formula corrected; `epsilon`, `fdelta` and `c`
  defaults aligned with Awesome Spectral Indices **(results change for GARI,
  EBI, SEVI, OCVI)**; data columns can no longer override index constants.
* `evaluate_vi_loco()` adds within-column accuracy, standardised coefficients
  and VIFs, and skips rank-deficient folds instead of predicting from them.
* `generate_met_data()`/`generate_genomic_met()` randomise genotypes within
  replicate (positions were previously identical across sites) and no longer
  reset the global random seed.
* Ward clustering of environments uses the Euclidean distance √(2(1 − r)).

## Other changes

* New: `compare_spatial_models()`, `select_best_model()`,
  `extract_variance_components()`, `prepare_trial_grm()`,
  `prepare_asreml_grm()` (sparse ASReml inverse; dense inverses can crash
  ASReml-R 4.2), and a `print()` method for `diagnose_design()`.
* `summary()` methods return summary objects with `print()` methods.
* `process_weather_indices()` reports `WaterBalance` (P − ET0; negative =
  deficit) instead of the misnamed `WaterDeficit`, only when ET0 is present,
  and returns `NA` for windows with poor weather coverage.
* `calculate_days_since_start()` is locale-independent and accepts ISO dates.
* `match_germplasm()` is vectorised.
* `plot_trend()` treats years as numeric; `DEFAULT_CHECKS` is no longer a default.
* `fit_bivariate_models()` returns numeric columns.
