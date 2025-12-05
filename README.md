# GenomicFlow

**Genomic Selection and Factor Analytic Tools for Plant Breeding**

`GenomicFlow` (formerly `UtilityFunctions`) is a comprehensive suite for analyzing Multi-Environment Trials (MET) and predicting cross utility in breeding programs.

It unifies the analysis pipeline:

1.  **Structure**: The `Design Tableau` validates your trial design before analysis.
2.  **Solver**: `fit_met_model` wraps `ASReml-R` to fit Factor Analytic (FA) models and automatically rotates parameters.
3.  **Selection**: Implements "FAST" (Factor Analytic Selection Tools) for variety advancement (OP & RMSD).
4.  **Prediction**: Uses a C++ backend to predict the utility of crosses (Mean + Variance) from genomic data.

## Installation

```r
# Install from GitHub
devtools::install_github("zaja10/GenomicFlow")
```

_Requires `asreml` for model fitting._

## Basic Workflow

### 1. Structural Validation (Design Tableau)

Ensure your data structure (Rep nesting, Aliasing) is correct.

```r
library(GenomicFlow)
library(agridat)

tableu <- build_design_tableau(dasilva.maize, ~ gen, ~ env/rep)
print(tableu)
```

### 2. Smart Selection (FAST)

Fit an FA2 model and select stable, high-performing varieties.

```r
# Fit Model (Wrapper for asreml)
model <- fit_met_model(tableu, k=2, genotype="gen", site="env", rotation="varimax")

# Visualization: OP vs RMSD
plot(model, type = "fast")
```

### 3. Cross Prediction

Prioritize crosses based on the Usefulness Criterion ($UC = \mu + i\sigma_g$).

```r
preds <- predict_cross_utility(parents, markers, effects, map)
head(preds)
```

## References

Smith, A. B., & Cullis, B. R. (2018). Plant breeding selection tools built on factor analytic mixed models for multi-environment trial data. _Euphytica_, 214(8).
