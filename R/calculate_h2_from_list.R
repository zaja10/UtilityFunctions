#' Calculate Heritability from a List of ASReml Models
#'
#' This function iterates through a list of fitted asreml model objects and
#' calculates broad-sense heritability (H2) based on the Cullis method
#' (1 - mean(PEV) / Vg), using the variance-covariance matrix of predictions.
#' Note: This function uses `predict(..., vcov = TRUE)` and may fail with a
#' "long vector not supported" error if the genetic factor has many levels,
#' leading to a very large prediction variance-covariance matrix. A safer
#' alternative uses standard errors from predict.
#'
#' @param model_list A named list containing fitted `asreml` objects.
#' @param id A character string specifying the name of the factor representing
#'   the genetic random effect (e.g., "Genotype", "Variety"). Defaults to "Genotype".
#'
#' @return A named list or vector containing the calculated heritability (H2)
#'   for each valid and converged model in the input list. Returns NA for
#'   models that are invalid, did not converge, or encountered errors during
#'   calculation.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'fm' is a list like fm <- list(fm0=model0, fm1=model1, ...)
#' heritabilities <- calculate_h2_from_list(fm, id = "Variety")
#' print(heritabilities)
#' }
calculate_h2_from_list <- function(model_list, id = "Genotype") {

  # Check if input is a list
  if (!is.list(model_list)) {
    stop("Input 'model_list' must be a list.")
  }

  # Initialize an empty list to store the results
  h2_results <- list()

  # Loop through each model in the list 'model_list'
  for (model_name in names(model_list)) {

    print(paste("Processing model:", model_name))
    asr <- model_list[[model_name]]

        # --- Basic Checks ---
    # Check if the element is actually an asreml object and converged
    if (is.null(asr) || !inherits(asr, "asreml")) {
      print(paste("Skipping", model_name, "- Not a valid asreml object."))
      h2_results[[model_name]] <- NA # Store NA for invalid objects
      next # Go to the next iteration
    }
    if (!asr$converge) {
      print(paste("Skipping", model_name, "- Model did not converge."))
      h2_results[[model_name]] <- NA # Store NA for non-converged models
      next # Go to the next iteration
    }

    # --- Calculation within tryCatch ---
    h2_calc <- tryCatch({

      # Extract response variable name (safer way)
      response <- NA
      # Try different ways to access the response variable name based on asreml versions/call structure
      if (!is.null(asr$call$fixed) && length(as.list(asr$call$fixed)) >= 2) {
        response <- as.character(as.list(asr$call$fixed)[[2]])
      } else if (!is.null(asr$call[[2]]) && length(asr$call[[2]]) >= 2) {
        # Fallback for older style or direct formula input potentially
        response <- as.character(asr$call[[2]][[2]])
      }
      if(is.na(response)){
        warning(paste("Could not determine response variable name for", model_name))
        response <- "h2" # Default name if extraction fails
      }

      # Run predict with vcov=TRUE (this is where the "long vector" error might occur)
      # Note: 'only=id' is removed as it wasn't in the final user script and can sometimes cause issues
      mypred <- predict(asr, classify = id, maxit = 1, vcov = TRUE) #Added average, sed

      # Get the prediction vcov matrix
      if (!"vcov" %in% names(mypred)) {
        stop("vcov component missing from predict output.")
      }
      my.vcov <- mypred$vcov

      # Find the genetic variance component index
      summary_asr <- summary(asr, coef = FALSE) # Get summary efficiently
      var_comp_table <- summary_asr$varcomp

      # Improved pattern to match id!id or vm(id, ...) more reliably
      pattern <- paste0("^(", id, ")!|^(vm\\(", id, "[ ,])") # Matches id! or vm(id, or vm(id<space>
      which.vc <- grep(pattern, rownames(var_comp_table))

      # Check if exactly one component was found
      if (length(which.vc) != 1) {
        # Try a simpler grep if the first failed, maybe simpler naming used
        which.vc <- grep(paste0("^",id,"$"), rownames(var_comp_table)) # Match exact term name
        if (length(which.vc) != 1) {
          stop(paste0("Could not uniquely identify variance component for '", id,
                      "'. Found: ", length(which.vc), " matches. Check component names: ",
                      paste(rownames(var_comp_table), collapse=", ")))
        }
      }
      Vg <- var_comp_table[which.vc, "component"]

      # Check if Vg is valid
      if (is.na(Vg) || Vg <= 1e-8) {
        stop(paste("Genetic variance component for '", id, "' is zero, negative, or NA (", Vg, ")."))
      }

      # Calculate h2 using diag() - this might fail for large matrices
      h2_value <- 1 - sum(diag(my.vcov)) / (Vg * nrow(my.vcov))

      # Ensure h2 is within bounds [0, 1]
      h2_value <- max(0, min(1, h2_value))

      # Name the result
      names(h2_value) <- response

      # Return the calculated value
      h2_value

    }, error = function(e) {
      # If any error occurred in the try block
      print(paste("Error processing", model_name, ":", e$message))
      return(NA) # Return NA if calculation fails
    })

    # Store the result (either the h2 value or NA)
    h2_results[[model_name]] <- h2_calc
  }

  # --- Combine results ---
  # Unlist to return a named vector
  h2_vector <- unlist(h2_results)

  # Print summary of NAs at the end
  print("------------------------------------")
  print("Heritability calculation summary:")
  print(paste("Number of models processed:", length(model_list)))
  print(paste("Number of successful calculations:", sum(!is.na(h2_vector))))
  print(paste("Number of failed calculations (NA):", sum(is.na(h2_vector))))
  print("------------------------------------")

  return(h2_vector)
}

# --- How to use the function ---
# Assuming 'fm' is your list of models:
# h2_values <- calculate_h2_from_list(fm, id = "Genotype")
# print(h2_values)
