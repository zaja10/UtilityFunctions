#' Factor Analytic Selection Tools (FAST)
#'
#' @description
#' Functions to calculate Smith & Cullis (2018) selection indices from a Factor Analytic model.
#'
#' @param model An object of class \code{fa_asreml}.
#'
#' @return A \code{fast_selection} object containing the selection table and metadata.
#' @name fast_tools
NULL

#' Calculate Overall Performance (OP) and RMSD
#'
#' @describeIn fast_tools Calculates OP and RMSD and returns a selection object.
#' @export
calculate_op <- function(model) {
    if (!inherits(model, "fa_asreml")) stop("Model must be of class fa_asreml")

    # Extract Rotated Parameters
    # Ensure these are ROTATED. fit_met_model ensures this by default.
    lam <- model$loadings$rotated
    sco <- model$scores$rotated

    if (is.null(sco)) stop("No genotype scores available in model.")

    k <- model$meta$k

    # OP = Factor 1 Score * Mean(Factor 1 Loading)
    # This relies on Factor 1 capturing "Overall Performance" (via Rotation)
    mean_lam1 <- mean(lam[, 1], na.rm = TRUE)
    op <- sco[, 1] * mean_lam1

    # RMSD calculation
    # Deviations from Factor 1 regression
    if (k > 1) {
        # Predicted Genetic Value excluding Fac1
        # G_dev = Fac2*Lam2 + ... + Fack*Lamk
        # We want the magnitude (RMSE) of this deviation part across sites?
        # Smith & Cullis define RMSD_i = sqrt( sum( (u_ij - u_i_OP)^2 ) / p ) ?
        # Approximation: RMSD is the variability explained by Factors 2..k

        # Effectively: variance of the effects of Factors 2..k across sites
        # For genotype i:
        # effects_i = scores[i, 2:k] %*% t(loadings[, 2:k])  (Vector of length p sites)
        # RMSD_i = sqrt( mean( effects_i^2 ) )

        effects_dev <- sco[, 2:k, drop = FALSE] %*% t(lam[, 2:k, drop = FALSE])
        rmsd <- apply(effects_dev, 1, function(x) sqrt(mean(x^2)))
    } else {
        rmsd <- rep(0, nrow(sco))
    }

    # Create Table
    df <- data.frame(
        Genotype = rownames(sco),
        OP = op,
        RMSD = rmsd
    )

    # Sort by OP
    df <- df[order(df$OP, decreasing = TRUE), ]

    obj <- list(
        selection = df,
        source_model = model
    )
    class(obj) <- "fast_selection"

    return(obj)
}

#' Calculate RMSD (Alias/Specific Wrapper)
#'
#' @describeIn fast_tools Calculates RMSD specifically (currently integrated in calculate_op).
#' @export
calculate_rmsd <- function(model) {
    # Wraps calculate_op but emphasizes RMSD
    res <- calculate_op(model)
    # Could re-sort by RMSD if desired?
    # For now, just return the same object
    return(res)
}

#' Plot FAST Selection
#'
#' @param x A \code{fast_selection} object.
#' @param type Character. "op_rmsd" (default) or "latent_reg".
#' @param ... Additional arguments.
#' @export
plot.fast_selection <- function(x, type = "op_rmsd", ...) {
    df <- x$selection

    if (type == "op_rmsd") {
        plot(df$RMSD, df$OP,
            xlab = "Stability (RMSD) [Lower is Stable]",
            ylab = "Overall Performance (OP) [Higher is Better]",
            main = "FAST Selection: OP vs RMSD",
            pch = 19, col = "steelblue", ...
        )
        grid()
        abline(h = mean(df$OP), col = "red", lty = 2)
        abline(v = mean(df$RMSD), col = "red", lty = 2)
    } else if (type == "latent_reg") {
        # Requires access to source model (which we stored)
        # Plot Genotype Score Fac 1 (Slope) vs Fac 1 Loading?
        # Actually Latent Regression usually plots Site Loading 1 vs Site Loading 2 (Biplot)
        # OR Genotype Effect across sites (Reaction Norms) vs Fac 1 Loading.

        # Smith & Cullis Latent Reg: E[u_ij] = f_i1 * lam_j1.
        # Plot u_ij against lam_j1. Slope is f_i1.
        # We can visualize this for top genotypes.
        message("Latent Regression plot not yet fully implemented.")
    }
}
