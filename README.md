# UtilityFunctions

**Mixed-model and factor analytic tools for plant breeding**

`UtilityFunctions` post-processes multi-environment trials (METs) fitted with
ASReml-R: spatial model selection, factor analytic (FA) selection tools
(FAST), genomic prediction diagnostics, design and connectivity checks,
vegetation indices and envirotyping.

## Installation

```r
# install.packages("remotes")
remotes::install_github("zaja10/UtilityFunctions")
```

ASReml-R (licensed) is a *suggested* dependency: the package installs and the
non-model functions work without it, and model functions abort with a clear
message when it is missing. `ASExtras4` is optional and only adds variograms to
the pipeline's diagnostic PDFs.

## Workflows

### 1. Spatial model selection

All candidate models must share the same response and fixed effects so that
their REML AICs are comparable; this is checked. Converged fits are ranked
first.

```r
library(UtilityFunctions)

models <- list(
  AxAe = list(fixed = Yield ~ 1, random = ~ Genotype + units, residual = ~ ar1(Col):ar1(Row)),
  AxA  = list(fixed = Yield ~ 1, random = ~ Genotype,         residual = ~ ar1(Col):ar1(Row)),
  RCD  = list(fixed = Yield ~ 1, random = ~ Genotype,         residual = ~ Col:Row)
)
comp <- compare_spatial_models(models, df)
best <- select_best_model(comp)

# Variance components, entry-mean h2, generalised reliability (Cullis H2) and accuracy
extract_variance_components(best_fit, genotype_term = "Genotype", n_reps = 2)
```

### 2. Univariate → review → MET pipeline

```r
Ginv <- prepare_asreml_grm(G)                       # sparse inverse for vm()
run_univariate_pipeline(pheno, Ginv, "outputs")      # per-trial spatial models, outliers, diagnostics
generate_diagnostic_template("outputs/tables", pheno, "outputs/review.xlsx")
# ... edit review.xlsx (terms, remove_outliers) ...
run_modified_trials(pheno, "outputs/review.xlsx", Ginv, "outputs")
met <- run_met_pipeline(pheno, Ginv, "outputs/Final_Model_Summary.xlsx", workspace = "8gb")
```

Each step returns a per-trial status table and warns about failures rather
than hiding them. A reviewer modification that removes outliers or changes the
fixed effects is not AIC-comparable with the original models, so it is
recorded as an explicit reviewer decision.

### 3. Factor analytic selection tools (FAST)

```r
model <- asreml::asreml(yield ~ env, random = ~ fa(env, 2):vm(gen, G), data = met)
fa <- fa.asreml(model, grm = G)        # G (not its inverse) gives correct reliabilities (r2_g)
summary(fa)
rank_genotypes(fa, weight_stability = 1)
calculate_i_classes(fa)
plot(fa, type = "fast")                # also "heatmap", "vaf", "biplot", "latent_reg", "cluster"
```

### 4. Genomic selection validation

```r
gs <- evaluate_gs_predictions(model, target_ids, training_ids,
                              target_term = "fa(env, 2):vm(gen, G.inv)", grm = G)
gs; summary(gs); plot(gs)
```

### 5. High-throughput phenotyping

```r
vis <- calculate_vegetation_indices(bands)            # Awesome Spectral Indices formulas
loco <- evaluate_vi_loco(df, "Yield", c("NDVI", "NDRE"))
loco$Acc_within                                       # within-column accuracy
```

## Plot types for `fa_model`

| Type | Shows |
|------|-------|
| `"fast"` | Overall performance (OP) vs stability (RMSD) |
| `"heatmap"` | Genetic correlations between environments |
| `"vaf"` | Variance accounted for per environment |
| `"biplot"` | Rotated scores and (scaled) loadings |
| `"latent_reg"` | Reaction norms on a factor's loadings |
| `"cluster"` | Ward clustering of environments on \(\sqrt{2(1-r)}\) |

## References

- Smith, A. B., & Cullis, B. R. (2018). Plant breeding selection tools built on factor analytic mixed models for multi-environment trial data. *Euphytica*, 214, 143.
- Smith, A., Norman, A., Kuchel, H., & Cullis, B. (2021). Plant variety selection using interaction classes derived from factor analytic linear mixed models. *Frontiers in Plant Science*, 12, 737462.
- Cullis, B. R., Smith, A. B., & Coombes, N. E. (2006). On the design of early generation variety trials with correlated data. *JABES*, 11, 381–393.
- Montero, D., et al. (2023). A standardized catalogue of spectral indices to advance the use of remote sensing in Earth system research. *Scientific Data*, 10, 197.
- Gevartosky, R., Munaro, L. B., & Rutkoski, J. E. (2026). Phenotype imputation using high-throughput phenotyping produces a new secondary trait for further selection modeling. *The Plant Phenome Journal*, 9, e70081.
