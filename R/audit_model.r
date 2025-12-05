#' Audit ASReml Model for Design Validation
#'
#' @description
#' Performs a "Design Tableau" audit on a fitted ASReml model.
#' It checks if the model structure respects the physical trial design (Row/Column).
#' Failing to account for spatial variation in field trials is a common source of error.
#'
#' @param model A fitted \code{asreml} object.
#'
#' @return Prints diagnostic messages. Returns TRUE if pass, FALSE if warnings found.
#' @export
audit_model_terms <- function(model) {
    chk_pass <- TRUE
    cat("--- Model Design Audit ---\n")

    # 1. Recover Data
    data_name <- as.character(model$call$data)
    if (!exists(data_name)) {
        warning("Original dataframe not found in environment. Skipping data checks.")
        df <- NULL
    } else {
        df <- get(data_name)
    }

    # 2. Check for Spatial Coordinates
    has_row <- !is.null(df) && any(grepl("^row$", names(df), ignore.case = TRUE))
    has_col <- !is.null(df) && any(grepl("^col(umn)?$", names(df), ignore.case = TRUE))

    if (has_row && has_col) {
        cat("-> Detected Spatial Data (Row/Col columns found).\n")

        # helper: convert call to string
        res_str <- paste(deparse(model$call$residual), collapse = "")
        rand_str <- paste(deparse(model$call$random), collapse = "")

        # 3. Check for Spatial Terms in Residual or Random
        # Looking for ar1(Row), idv(Row), at(Site):ar1(Row), etc.
        # We look for "Row" or "Col" inside the model formulas

        terms_found <- c()
        if (grepl("Row", res_str, ignore.case = TRUE)) terms_found <- c(terms_found, "Row (Residual)")
        if (grepl("Col", res_str, ignore.case = TRUE)) terms_found <- c(terms_found, "Col (Residual)")
        if (grepl("Row", rand_str, ignore.case = TRUE)) terms_found <- c(terms_found, "Row (Random)")
        if (grepl("Col", rand_str, ignore.case = TRUE)) terms_found <- c(terms_found, "Col (Random)")

        if (length(terms_found) == 0) {
            cat("[WARNING] Spatial columns exist but no Row/Col terms found in model.\n")
            cat("          Consider adding ar1(Row):ar1(Col) or similar spatial corrections.\n")
            chk_pass <- FALSE
        } else {
            cat(sprintf("-> Spatial terms detected: %s\n", paste(terms_found, collapse = ", ")))
        }
    } else {
        cat("-> No obvious spatial columns (Row/Col) detected in data.\n")
    }

    # 4. Check for 'units' (Nugget)
    # Often important if ar1 is used
    res_str <- paste(deparse(model$call$residual), collapse = "")
    if (grepl("units", res_str)) {
        cat("-> 'units' term detected in residual (Nugget effect).\n")
    }

    if (chk_pass) {
        cat("-> Audit Passed: Design terms appear consistent.\n")
    } else {
        cat("-> Audit Finished with Warnings.\n")
    }

    return(invisible(chk_pass))
}
