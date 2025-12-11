#' Factor Analytic Selection Tools (FAST)
#'
#' Functions to calculate Smith & Cullis (2018) selection indices from a Factor Analytic model.
#'
#' @name fast_tools
#' @return A \code{fast_selection} object or dataframes of indices.
NULL

#' Calculate Overall Performance (OP) and RMSD
#'
#' Calculates OP (Factor 1 performance) and RMSD (deviation/instability) and returns a selection object.
#'
#' @param model An object of class \code{fa_asreml}.
#' @return A \code{fast_selection} object.
#' @export
calculate_op <- function(model) {
    if (!inherits(model, "fa_asreml")) stop("Model must be of class fa_asreml")

    lam <- model$loadings$rotated
    sco <- model$scores$rotated

    if (is.null(sco)) stop("No genotype scores available in model.")

    k <- model$meta$k

    # OP = Factor 1 Score * Mean(Factor 1 Loading)
    mean_lam1 <- mean(lam[, 1], na.rm = TRUE)
    op <- sco[, 1] * mean_lam1

    # RMSD calculation
    if (k > 1) {
        effects_dev <- sco[, 2:k, drop = FALSE] %*% t(lam[, 2:k, drop = FALSE])
        rmsd <- apply(effects_dev, 1, function(x) sqrt(mean(x^2)))
    } else {
        rmsd <- rep(0, nrow(sco))
    }

    df <- data.frame(
        Genotype = rownames(sco),
        OP = op,
        RMSD = rmsd
    )
    df <- df[order(df$OP, decreasing = TRUE), ]

    obj <- list(
        selection = df,
        source_model = model
    )
    class(obj) <- "fast_selection"
    return(obj)
}

#' Calculate RMSD (Wrapper)
#'
#' @param model An object of class \code{fa_asreml}.
#' @return A \code{fast_selection} object (same as calculate_op).
#' @export
calculate_rmsd <- function(model) {
    return(calculate_op(model))
}

#' Calculate D-Optimality (Network Efficiency)
#'
#' Quantifies the information content of the trial network.
#'
#' @param object An object of class \code{fa_asreml}.
#' @return A list containing total_d and site_impact.
#' @export
calculate_d_optimality <- function(object) {
    lam <- object$loadings$rotated
    M_total <- t(lam) %*% lam
    total_det <- det(M_total)

    impact_scores <- numeric(nrow(lam))
    for (i in seq_len(nrow(lam))) {
        lam_min <- lam[-i, , drop = FALSE]
        impact_scores[i] <- (total_det - det(t(lam_min) %*% lam_min)) / total_det * 100
    }
    site_impact <- data.frame(Site = rownames(lam), Impact_Pct = round(impact_scores, 2))
    site_impact <- site_impact[order(site_impact$Impact_Pct, decreasing = TRUE), ]
    return(list(total_d = total_det, site_impact = site_impact))
}

#' Calculate Interaction Classes (iClasses)
#'
#' Implements Smith et al. (2021) to detect Specific Adaptation.
#'
#' @param object An object of class \code{fa_asreml}.
#' @param factor Integer. Which factor defines the interaction? (Default 2).
#' @param threshold Numeric. Loading magnitude for class assignment.
#'
#' @return A list with site classes and genetic effects.
#' @export
calculate_i_classes <- function(object, factor = 2, threshold = 0.1) {
    if (is.null(object$scores)) stop("Genotype scores (BLUPs) not found.")
    if (factor > object$meta$k) stop("Factor not found.")

    lam <- object$loadings$rotated[, factor]
    sco <- object$scores$rotated[, factor]
    lam1 <- object$loadings$rotated[, 1]
    sco1 <- object$scores$rotated[, 1]

    site_class <- rep("Neutral", length(lam))
    names(site_class) <- names(lam)
    site_class[lam > threshold] <- "Class_Pos"
    site_class[lam < -threshold] <- "Class_Neg"

    ml_pos <- mean(lam[site_class == "Class_Pos"], na.rm = TRUE)
    ml_neg <- mean(lam[site_class == "Class_Neg"], na.rm = TRUE)
    m_perf <- mean(lam1, na.rm = TRUE)

    pred_pos <- if (is.nan(ml_pos)) rep(NA, length(sco)) else (sco1 * m_perf) + (sco * ml_pos)
    pred_neg <- if (is.nan(ml_neg)) rep(NA, length(sco)) else (sco1 * m_perf) + (sco * ml_neg)

    gen_eff <- data.frame(
        Genotype = names(sco), Pred_Pos = pred_pos, Pred_Neg = pred_neg,
        Differential = pred_pos - pred_neg
    )
    gen_eff <- gen_eff[order(gen_eff$Differential, decreasing = TRUE, na.last = TRUE), ]

    return(list(
        site_classes = data.frame(Site = names(lam), Class = site_class),
        gen_effects = gen_eff, meta = list(factor = factor)
    ))
}

#' Calculate Selection Index
#'
#' Combines Overall Performance (OP) and specific Stability penalty (RMSD).
#'
#' @param object An object of class \code{fa_asreml}.
#' @param weight Numeric. Risk aversion weight applied to RMSD (default 1.0).
#' @return A dataframe sorted by the calculated Index.
#' @export
calculate_index <- function(object, weight = 1.0) {
    if (is.null(object$fast)) stop("No FAST indices found in object.")
    df <- object$fast
    df$Index <- df$OP - (weight * df$RMSD)
    df$Rank <- rank(-df$Index)
    df <- df[order(df$Index, decreasing = TRUE), ]
    return(df[, c("Rank", "Genotype", "OP", "RMSD", "Index")])
}

#' Plot FAST Selection
#'
#' @param x A \code{fast_selection} object.
#' @param type Character. "op_rmsd" (default).
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
    }
}
