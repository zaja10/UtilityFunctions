# test_mating_planner.R

# We rely on devtools::load_all() running before this script to load the package environment

cat("Loaded UtilityFunctions.\n")

# Create dummy data
set.seed(123)
n_markers <- 100
n_geno <- 50
M <- matrix(sample(c(0, 2), n_geno * n_markers, replace = TRUE), n_geno, n_markers)
rownames(M) <- paste0("G", 1:n_geno)
eff <- rnorm(n_markers)
map <- seq(0, 10, length.out = n_markers)

# Run Optimizer
parents <- rownames(M)[1:10]
cat("Running optimization...\n")

tryCatch(
    {
        res <- optimize_mating_list(
            parents = parents,
            marker_effects = eff,
            map = map,
            markers = M,
            n_crosses = 5,
            iterations = 10
        )
        print(head(res))
        cat("\nSUCCESS: Optimization completed.\n")
    },
    error = function(e) {
        cat("\nFAILURE: ", e$message, "\n")
    }
)
