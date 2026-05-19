#' Calculate Genetic Gain
#'
#' Estimates genetic gain using either "Realized" (Historical Regression) or
#' "Predicted" (Breeder's Equation) methods.
#'
#' @param input Either a dataframe (for "Realized") or an `fa_model`/`asreml` object (for "Predicted").
#' @param ... arguments handling specifics.
#' \describe{
#'   \item{Realized}{Requires `year_col` (Year) and `value_col` (Yield/BV).}
#'   \item{Predicted}{Requires `selection_percent` (default 5).}
#' }
#'
#' @return
#' For Realized: A list with slope, R2, and pct_gain.
#' For Predicted: A dataframe of potential gain per trait/site.
#' @importFrom stats lm coef dnorm qnorm as.formula
#' @export
calculate_genetic_gain <- function(input, ...) {
    args <- list(...)

    # MODE 1: PREDICTED (Model Object)
    if (inherits(input, c("asreml", "fa_model"))) {
        sel_pct <- if ("selection_percent" %in% names(args)) args$selection_percent else 5
        if (sel_pct <= 0 || sel_pct >= 100) stop("Selection % must be 0-100")

        # Intensity
        p <- sel_pct / 100
        i_val <- dnorm(qnorm(1 - p)) / p

        # Sigma G
        if (inherits(input, "fa_model")) {
            Vg <- diag(input$matrices$G)
            lbl <- names(Vg)
        } else {
            # Simple scalar model assumption
            summ <- summary(input)$varcomp
            Vg <- summ[grep("Genotype", rownames(summ))[1], "component"]
            lbl <- "Global"
        }

        param_gain <- i_val * sqrt(Vg) # Assuming r=1

        return(data.frame(
            Group = lbl,
            Sigma_G = sqrt(Vg),
            Intensity = i_val,
            Predicted_Gain = param_gain
        ))
    }

    # MODE 2: REALIZED (Dataframe)
    else if (is.data.frame(input)) {
        year_col <- if ("year_col" %in% names(args)) args$year_col else "Year"
        val_col <- if ("value_col" %in% names(args)) args$value_col else "Yield"
        checks <- if ("check_list" %in% names(args)) args$check_list else NULL

        if (!all(c(year_col, val_col) %in% names(input))) stop("Columns missing.")

        df <- input
        if (!is.null(checks)) df <- df[!df$Genotype %in% checks, ]

        form <- as.formula(paste(val_col, "~", year_col))
        mod <- lm(form, data = df)

        slope <- coef(mod)[2]
        base_yr <- min(df[[year_col]], na.rm = TRUE)
        base_val <- predict(mod, newdata = data.frame(setNames(list(base_yr), year_col)))

        return(list(
            Method = "Realized (Regression)",
            Slope = slope,
            Pct_Gain_Per_Year = (slope / base_val) * 100,
            R2 = summary(mod)$r.squared
        ))
    } else {
        stop("Input must be a dataframe or asreml/fa_model object.")
    }
}

#' Calculate Multi-Environment Heritability
#'
#' Calculates plot-basis heritability (broad-sense H2 and narrow-sense h2)
#' for multi-environment trials (MET) from a fitted ASReml model with diagonal
#' or factor analytic covariance structures.
#'
#' @param model A fitted \code{asreml} object.
#' @param gdrop_term Character string. The name of the non-additive/ungenotyped genetic term in vparameters (default "studyName:gdrop").
#' @param gkeep_term Character string. The name of the additive/genotyped genetic term in vparameters (default "studyName:vm(gkeep, G.inv)").
#' @param env_var Character string. The environment factor name prefix (default "studyName").
#'
#' @return A data frame containing:
#' \describe{
#'   \item{Trial}{The name of the trial/environment.}
#'   \item{Vg_add}{Additive genetic variance component.}
#'   \item{Vg_na}{Non-additive genetic variance component.}
#'   \item{Vg_total}{Total genetic variance (Vg_add + Vg_na).}
#'   \item{V_spatial}{Spatial/residual variance component.}
#'   \item{V_units}{Units/nugget variance component.}
#'   \item{H2}{Broad-sense heritability on a plot-basis.}
#'   \item{h2}{Narrow-sense heritability on a plot-basis.}
#' }
#' @export
calculate_met_h2 <- function(model, gdrop_term = "studyName:gdrop", gkeep_term = "studyName:vm(gkeep, G.inv)", env_var = "studyName") {
    check_asreml_availability()

    if (!inherits(model, "asreml")) {
        stop("Input must be a valid asreml model.")
    }

    # 1. Extract vparameters
    summ_list <- summary(model, vparameters = TRUE)
    if (!"vparameters" %in% names(summ_list)) {
        stop("Could not extract vparameters from the model summary.")
    }
    varp <- summ_list$vparameters

    # 2. Extract Genetic Components
    if (!gdrop_term %in% names(varp)) {
        stop(paste("gdrop term", gdrop_term, "not found in vparameters. Available names:", paste(names(varp), collapse = ", ")))
    }
    if (!gkeep_term %in% names(varp)) {
        stop(paste("gkeep term", gkeep_term, "not found in vparameters. Available names:", paste(names(varp), collapse = ", ")))
    }

    vna_raw <- varp[[gdrop_term]]
    va_raw <- varp[[gkeep_term]]

    clean_trial <- function(n) {
        n <- sub(".*!", "", n)
        n <- sub(paste0("^", env_var, "_"), "", n)
        return(n)
    }

    trial_keys <- clean_trial(names(vna_raw))

    # 3. Build Results Data Frame
    results <- data.frame(
        Trial = trial_keys,
        Vg_add = as.numeric(va_raw),
        Vg_na = as.numeric(vna_raw),
        V_spatial = NA_real_,
        V_units = 0, # Default to 0 as some trials lack nugget
        stringsAsFactors = FALSE
    )
    results$Vg_total <- results$Vg_add + results$Vg_na

    # 4. Map Residuals and Units
    all_param_names <- names(varp)

    for (i in 1:nrow(results)) {
        tr <- results$Trial[i]

        # A. Extract Spatial Variance
        if (tr %in% all_param_names) {
            results$V_spatial[i] <- varp[[tr]][1]
        } else {
            spatial_pattern <- paste0(env_var, "_", tr, "!R")
            match_idx <- grep(spatial_pattern, all_param_names, fixed = TRUE)
            if (length(match_idx) > 0) {
                results$V_spatial[i] <- varp[[match_idx[1]]][1]
            }
        }

        # B. Extract Units (Nugget)
        unit_pattern <- paste0("'", tr, "'):units")
        match_unit <- grep(unit_pattern, all_param_names, fixed = TRUE)
        if (length(match_unit) > 0) {
            results$V_units[i] <- varp[[match_unit[1]]][1]
        }
    }

    # 5. Calculate Heritability (Plot-Basis)
    results$H2 <- results$Vg_total / (results$Vg_total + results$V_spatial + results$V_units)
    results$h2 <- results$Vg_add / (results$Vg_add + results$V_spatial + results$V_units)

    return(results)
}
