# **UtilityFunctions** 

**Factor Analytic Model Extraction and Selection Tools for Multi-Environment Trials (METs)**

The `UtilityFunctions` package provides a specialized suite of tools for processing and interpreting Factor Analytic (FA) and Reduced Rank (RR) models fitted using the `ASReml-R` software.

Its primary goal is to implement the Factor Analytic Selection Tools (FAST) methodology proposed by Smith & Cullis (2018), calculating clear, actionable metrics and visualizations for genotype selection and program optimization.

`fa.asreml()` is the core function to extract, rotate, and standardize FA/RR model components (loadings, scores, specific variances).

**FAST Indices**: Calculates Overall Performance (OP) and Stability (RMSD) metrics essential for variety selection.

**Visualization**: A single `plot()` function dispatches seven different types of publication-ready diagnostic and selection plots.

**Network Diagnostics**: Tools to check trial connectivity and assess network efficiency using D-Optimality.

**Advanced Analytics**: Calculation of Interaction Classes (iClasses) to detect specific adaptation.

**General Utilities**: Functions for simulating data, calculating broad-sense heritability, and calculating spectral vegetation indices.

This package is designed for a development environment. You can install the latest version directly from GitHub using `remotes`. You must have a working installation of `ASReml-R` and a valid license.

# Basic Example

**Step 1: Simulate a dataset for demonstration**

```{r}
df <- generate_met_data(n_sites = 8, n_gen = 120)

# Fit FA2 Model
model_fit <- asreml(fixed = Yield ~ Site,
                    random = ~fa(Site, 2):Genotype + Rep,
                    residual = ~idv(units), 
                    data = df, 
                    trace = FALSE)


```

**Step 2: Extract & Calculate FAST Indices**

Use fa.asreml() to extract the parameters, perform the necessary SVD rotation, and calculate the FAST selection indices. The 'classify' argument must exactly match the FA random term in your model.

```{r}
fa_results <- fa.asreml(model_fit, classify = "fa(Site, 2):Genotype")

# Print the top 5 genotypes based on Overall Performance (OP)
print(fa_results) 

# Access the full selection index dataframe (OP vs. RMSD)
head(fa_results$fast)

```

# Factor Analytic Selection Tools (FAST) 

The selection metrics are derived from the rotated Factor Scores and Loadings:

| Metric | Calculation | Interpretation |
|----|----|----|
| OP (Overall Performance) | Factor 1 Score X Mean Factor 1 Loading | Predicted genetic value across the network. Higher is generally better |
| RMSD (Stability) | Root Mean Square Deviation from Factor 1 Regression | Quantifies non-linear GxE interaction (instability). For product development lower is generally better. |

You can use the `calculate_index()` function to combine these metrics based on a custom weight:

$$
Index = OP - (\text{weight} \times RMSD)
$$

Where the `weight` parameter (default = 1) determines the penalty applied to instability. A higher weight favors more stable genotypes.

```{r}
selection_ranking <- calculate_index(fa_results, weight = 2)

head(selection_ranking)
```

## Visualization Suite 

The package includes `plot.fa_asreml` method with multiple `type` options to visualize different aspects of the analysis.

1.  The Selection View (type = "fast") This is the primary tool for variety selection. It plots Performance (OP) vs. Stability (RMSD).

```{r}
# Label the top 5 performers and highlight specific checks
plot(fa_results, type = "fast", 
     n_label = 5, 
     highlight = c("G001", "G002"), 
     main = "FAST Selection: Performance vs. Stability")

```

2.  Genetic Correlation Heatmap (type = "heatmap") Visualize the genetic relationship structure between environments to identify mega-environments.

```{r}
plot(fa_results, type = "heatmap")
```

3.  Latent Regression (type = "latent_reg") A mechanistic view that visualizes why a variety is unstable by showing its specific response (slope) to the Factor 2 GxE driver.

```{r}
#View specific GxE drivers (FA2)
plot(fa_results, type = "latent_reg", factor = 2)
```

## Diagnostics and Network Optimization

| Function | Purpose | Plot Option |
|:----------------|:--------------------------|:-----------------------|
| check_connectivity() | Verifies common lines between all trial pairs before fitting the model to ensure genetic correlations are estimable. | Automatic Heatmap |
| get_d_optimality() | Assesses the unique information content of each site in the network, helping to identify redundant sites. | plot(..., type = "d_opt") |
| get_i_classes() | Partitions environments into "Positive" and "Negative" classes based on a GxE driver (Factor k\>1) to detect crossover interactions. | plot(..., type = "diff") |
| plot_spatial() | Maps raw data or residuals onto the physical field layout (Row/Column coordinates) for spatial diagnostic checks. | N/A (Dedicated function) |

**Examples**

```{r}
# Check trial connectivity
connectivity_report <- check_connectivity(df, genotype = "Genotype", 
                                          trial = "Site", threshold = 10)

# View network efficiency
plot(fa_results, type = "d_opt")

# Calculate and plot Interaction Classes
iclass <- get_i_classes(fa_results, factor = 2)

plot(fa_results, type = "diff", 
     factor = 2, n_label = 5)

# Plot residuals spatially (requires Row/Column in model data)
plot_spatial(model_fit, row = "Row", col = "Column", 
             attribute = "residuals")

```

# Utility Functions 

The package also includes several utility functions for related tasks:

-   `generate_grm()` Generates a synthetic Genomic Relationship Matrix (GRM) for simulation purposes.

-   `generate_genomic_met()` Creates MET data where phenotypes are driven by a GRM.

-   `calculate_h2_from_list()` Calculates broad-sense heritability (H2) for a list of ASReml models.

-   `calculate_vegetation_indices()` Calculates over 100 spectral vegetation indices (like NDVI, EVI) from standard band names (N, R, G, B, RE1, etc.)

**Reference**

Smith, A. B., & Cullis, B. R. (2018). Plant breeding selection tools built on factor analytic mixed models for multi-environment trial data. Euphytica, 214(8).
