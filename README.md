# UtilityFunctions

**Comprehensive Mixed Modeling & Factor Analytic Tools for Plant Breeding**

`UtilityFunctions` provides a specialized suite of tools for robust post-processing of Multi-Environment Trials (MET) using ASReml-R, handling Factor Analytic (FA) models, spatial models, bivariate traits, and genomic/phenomic diagnostics.

## Key Features

1. **Spatial & Bivariate Modeling**: Streamlined wrappers (`compare_spatial_models`, `fit_bivariate_models`) to automate multi-trait and spatial model fitting.
2. **FA Model Diagnostics**: A universal parser (`fa.asreml`) to extract, rotate, and summarize Factor Analytic/Reduced Rank parameters using the "Smith & Cullis" FAST framework (OP & RMSD).
3. **Outlier Detection**: Automated Average Outlier Measure (AOM) calculation using `extract_asreml_outliers()`.
4. **Spectral Analysis**: Evaluate Vegetation Indices using LOCO (Leave-One-Column-Out) cross-validation via `evaluate_vi_loco()`.
5. **Genomic Utilities**: Subsetting and padding Genomic Relationship Matrices (GRM) safely for ASReml via `prepare_trial_grm()`.

## Installation

```r
# Install from GitHub
devtools::install_github("zaja10/UtilityFunctions")
```

*Requires `asreml` for model fitting workflows.*

## Basic Workflow Examples

### 1. Spatial Model Selection

Evaluate multiple spatial structures for a single trial and extract variance components.

```r
library(UtilityFunctions)

# Define models
models <- list(
  RCD = list(fixed = Yield ~ 1, random = ~ Genotype, residual = ~ ar1(Col):ar1(Row)),
  Units = list(fixed = Yield ~ 1, random = ~ Genotype, residual = ~ dsum(~ ar1(Col):ar1(Row) | units))
)

# Compare AIC/LogLik
comp_df <- compare_spatial_models(models, df)

# Extract h2, PEV, and Accuracy
var_df <- extract_variance_components(best_model, genotype_term = "vm(Genotype, G.inv)", error_term = "Col:Row!R")
```

### 2. Factor Analytic Models (FAST)

Fit an FA model and extract rotated solutions for selection.

```r
# Fit Model
model <- asreml(fixed = yield ~ env,
                random = ~ fa(env, 2):gen,
                data = dasilva.maize)

# Extract & Rotate
results <- fa.asreml(model, classify = "fa(env, 2):gen", rotate = "mean")

# Clean Summary (VAF, Loadings, OP Index)
summary(results)

# Extract matrices directly
loadings <- extract_fa_loadings(results)$Loadings
```

### 3. High-Throughput Phenotyping (HTP)

Evaluate Vegetation Indices (VIs) predicting Yield using LOCO-CV.

```r
vi_results <- evaluate_vi_loco(df, target_trait = "Yield", vi_names = c("NDVI", "NDRE"))
print(vi_results$Col_Acc) # Column-wise accuracy
```

## Visualization Options

The `plot()` function supports various types to inspect the `fa_model`:

| Type | Description |
|------|-------------|
| `"fast"` | **Selection**: Overall Performance (OP) vs Stability (RMSD) |
| `"heatmap"` | **Connectivity**: Genetic Correlation Matrix between sites |
| `"latent_reg"` | **GxE**: Latent regression showing stability drivers |
| `"vaf"` | **Quality**: Variance Accounted For (%) per site |

## References

Smith, A. B., & Cullis, B. R. (2018). Plant breeding selection tools built on factor analytic mixed models for multi-environment trial data. _Euphytica_, 214(8).
