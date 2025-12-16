# UtilityFunctions

**Factor Analytic Model Extraction and Selection Tools**

`UtilityFunctions` provides a specialized suite for post-processing Factor Analytic (FA) models in plant breeding. It decouples the analysis integration, allowing you to fit models with `ASReml-R` and use this package for:

1.  **Structure**: `diagnose_design` validates experimental design.
2.  **Extraction**: `fa.asreml` extracts and rotates parameters from fitted models.
3.  **Selection**: "FAST" indices for variety advancement (OP & RMSD).
4.  **Prediction**: Rank crosses using gametic variance ($UC = \mu + i\sigma_g$).

## Installation

```r
# Install from GitHub
devtools::install_github("zaja10/UtilityFunctions")
```

_Requires `asreml` for model fitting._

## Basic Workflow

### 1. Structural Validation

Ensure your data structure (Rep nesting, Aliasing) is correct.

```r
library(UtilityFunctions)
library(agridat)

diag_report <- diagnose_design(dasilva.maize, genotype="gen", trial="env", rep="rep")
print(diag_report)
```

### 2. Fit Model & Extract (FAST)

Fit an FA2 model using ASReml, then extract rotated solution.

```r
# 1. Fit Model (Standard ASReml)
model <- asreml(fixed = yield ~ env,
                random = ~ fa(env, 2):gen,
                data = dasilva.maize)

# 2. Extract & Rotate
# 2. Extract & Rotate
results <- fa.asreml(model, classify = "fa(env, 2):gen", rotation = "varimax")

# 3. Access BLUEs & BLUPs (New)
blues <- results$blues                   # Centered Fixed Effects
blups <- results$scores$blups_in_met     # Reconstructed Site Effects

# 3. Visualization: OP vs RMSD
plot(results, type = "fast")
```

### 3. Cross Prediction

Prioritize crosses based on the Usefulness Criterion ($UC = \mu + i\sigma_g$).

```r
preds <- predict_cross_utility(parents, markers, effects, map)
head(preds)
```

## References

Smith, A. B., & Cullis, B. R. (2018). Plant breeding selection tools built on factor analytic mixed models for multi-environment trial data. _Euphytica_, 214(8).
