#' Validate MET Data Structure
#'
#' Performs robust quality checks on a Multi-Environment Trial dataset.
#'
#' @import cli
#' @export
validate_met_data <- function(data, trait, genotype, env, check_numeric = TRUE) {
    # 1. Column Existence
    req_cols <- c(trait, genotype, env)
    missing <- setdiff(req_cols, names(data))

    if (length(missing) > 0) {
        cli::cli_abort(c(
            "x" = "Required columns are missing from the dataset.",
            "i" = "Missing: {.val {missing}}"
        ))
    }

    # 2. Data Types
    if (check_numeric && !is.numeric(data[[trait]])) {
        cli::cli_abort(c(
            "x" = "Trait column {.val {trait}} must be numeric.",
            "i" = "Found class: {.cls {class(data[[trait])}}"
        ))
    }

    # 3. Missingness checks
    n_na_gen <- sum(is.na(data[[genotype]]))
    if (n_na_gen > 0) {
        cli::cli_abort(c(
            "x" = "Genotype column {.val {genotype}} contains {n_na_gen} missing values (NA).",
            "i" = "Please filter or impute missing identities before analysis."
        ))
    }

    n_na_env <- sum(is.na(data[[env]]))
    if (n_na_env > 0) {
        cli::cli_abort(c(
            "x" = "Environment column {.val {env}} contains {n_na_env} missing values (NA)."
        ))
    }

    # 4. Level Checks
    n_g <- length(unique(data[[genotype]]))
    n_e <- length(unique(data[[env]]))

    if (n_g < 2) {
        cli::cli_warn(c(
            "!" = "Data contains only 1 Genotype.",
            "i" = "Analysis usually requires genetic variation."
        ))
    }

    if (n_e < 1) {
        cli::cli_abort("Data must contain at least 1 Environment.")
    }

    # Success Message
    cli::cli_alert_success("Data Validated: {n_g} Genotypes across {n_e} Environments.")
    invisible(TRUE)
}
