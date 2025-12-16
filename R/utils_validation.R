#' Validate MET Data Structure
#'
#' Performs robust quality checks on a Multi-Environment Trial dataset before analysis.
#' Checks for missing values in key columns, data types, and structural integrity.
#'
#' @param data A dataframe.
#' @param trait Character. Name of the response variable.
#' @param genotype Character. Name of the genotype column.
#' @param env Character. Name of the environment/site column.
#' @param check_numeric Logical. Check if trait is strictly numeric? Default TRUE.
#'
#' @return Invisible TRUE if valid. Stops with informative error otherwise.
#' @export
validate_met_data <- function(data, trait, genotype, env, check_numeric = TRUE) {
    # 1. Column Existence
    missing_cols <- setdiff(c(trait, genotype, env), names(data))
    if (length(missing_cols) > 0) {
        stop(paste("Columns not found in data:", paste(missing_cols, collapse = ", ")))
    }

    # 2. Data Types
    if (check_numeric) {
        if (!is.numeric(data[[trait]])) {
            stop(paste0("Trait column '", trait, "' must be numeric. Found: ", class(data[[trait]])[1]))
        }
    }

    # 3. Missingness in Structure
    if (any(is.na(data[[genotype]]))) {
        n_na <- sum(is.na(data[[genotype]]))
        stop(paste("Genotype column contains", n_na, "missing values (NA). Please filter or impute."))
    }

    if (any(is.na(data[[env]]))) {
        n_na <- sum(is.na(data[[env]]))
        stop(paste("Site/Env column contains", n_na, "missing values (NA)."))
    }

    # 4. Level Checks
    n_gen <- length(unique(data[[genotype]]))
    n_env <- length(unique(data[[env]]))

    if (n_gen < 2) warning("Data contains only 1 Genotype. Analysis may fail.")
    if (n_env < 1) stop("Data must contain at least 1 Environment.")

    message(sprintf("✔ Data Validated: %d Genotypes across %d Environments.", n_gen, n_env))
    invisible(TRUE)
}
