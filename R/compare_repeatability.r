#' Compare Site-Specific Repeatability Metrics
#'
#' @description
#' Calculates the repeatability ($R$) for each environment in the trial.
#' Repeatability represents the upper limit of heritability and is defined as:
#' \deqn{R = \frac{\sigma^2_g}{\sigma^2_g + \sigma^2_e}}
#' where $\sigma^2_g$ is the genetic variance (heterogeneous in FA models) and
#' $\sigma^2_e$ is the residual variance (extracted from the model).
#'
#' @param model A fitted \code{asreml} object.
#' @param fa_object An object of class \code{fa_asreml} produced by \code{fa.asreml()}.
#'
#' @return A data.frame summarizing Repeatability, Genetic Variance (Vg), and Residual Variance (Ve) per Site.
#'
#' @details
#' This function attempts to automatically extract site-specific residual variances.
#' It supports common ASReml residual specifications:
#' \itemize{
#'   \item Homogeneous: \code{rcov = ~ units} or \code{idv(units)}
#'   \item Heterogeneous: \code{rcov = ~ at(Site):units} or \code{diag(Site):units}
#'   \item Spatial AR1: \code{rcov = ~ ar1(Row):ar1(Col)} (Uses the nugget/variance component)
#' }
#'
#' @export
compare_repeatability <- function(model, fa_object) {
    if (!inherits(fa_object, "fa_asreml")) stop("fa_object must be of class 'fa_asreml'.")

    # 1. SETUP: Extract Site-Specific Genetic Variances
    G_mat <- fa_object$matrices$G
    site_stats <- data.frame(Site = rownames(G_mat), Vg = diag(G_mat))
    sites <- as.character(site_stats$Site)

    # 2. EXTRACT RESIDUAL VARIANCES
    vc <- summary(model)$varcomp
    vc_names <- rownames(vc)

    # Initialize Ve (Error Variance) vector
    Ve_vec <- setNames(rep(NA, length(sites)), sites)

    cat("-> Extracting residual variances...\n")

    # Strategy 1: Look for explicit 'units' or 'R' terms linked to sites
    # e.g. "Site_Loc1!units!R", "at(Site, Loc1):units!R"

    for (site in sites) {
        # Regex to find this site's error component
        # Matches:
        # 1. "at(Site, SiteName):units!R" (Heterogeneous IDV)
        # 2. "SiteName!units!R" (Old school)
        # 3. "SiteName:units!R"
        # 4. "ids(Site, SiteName)!units!R"

        # Escape site name for regex
        site_esc <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", site)

        # Heuristic Patterns
        pats <- c(
            paste0("at\\(.*", site_esc, ".*\\).*(units|R)"), # at(Site, X):units
            paste0("^", site_esc, "([:!_]).*(units|R)"), # Site:units
            paste0("units.*", site_esc), # units:at(Site, X)
            "units!R", # Homogeneous fallback (checked last)
            "R!variance" # General fallback
        )

        match_found <- FALSE

        for (pat in pats) {
            if (pat == "units!R" || pat == "R!variance") {
                # Only use global fallback if NO specific matches found for ANY site yet?
                # Or if we are sure it's homogeneous.
                # Let's check specific site patterns first.
                next
            }

            matches <- grep(pat, vc_names, value = TRUE)
            if (length(matches) > 0) {
                # Take the row with "variance" or "R" usually
                # summary(model)$varcomp usually has component in column 1
                # Use the first match
                Ve_vec[site] <- vc[matches[1], "component"]
                match_found <- TRUE
                break
            }
        }

        # Fallback to Global Homogeneous if not found
        if (!match_found) {
            global_matches <- grep("(units!R|R!variance|^units$)", vc_names, value = TRUE)
            if (length(global_matches) > 0) {
                Ve_vec[site] <- vc[global_matches[1], "component"]
                # Note: This assigns the same Ve to this site.
            }
        }
    }

    # 3. CALCULATE REPEATABILITY
    results_list <- list()

    for (site in sites) {
        Vg <- site_stats$Vg[site_stats$Site == site]
        Ve <- Ve_vec[site]

        if (is.na(Ve)) {
            warning(sprintf("Could not determine residual variance for site %s. R set to NA.", site))
            R_val <- NA
        } else {
            R_val <- Vg / (Vg + Ve)
        }

        results_list[[site]] <- data.frame(
            Site = site,
            Vg = Vg,
            Ve = Ve,
            Repeatability = R_val
        )
    }

    final_df <- do.call(rbind, results_list)
    rownames(final_df) <- NULL

    # Sort by Repeatability descending
    final_df <- final_df[order(-final_df$Repeatability), ]

    return(final_df)
}
