#' Master Visualization Suite for GenomicFlow
#'
#' @description
#' The central plotting interface for the toolkit. Dispatches to specific visualization engines.
#'
#' @param x An object of class \code{fa_asreml}.
#' @param type A character string specifying the visualization module.
#' @param factor Integer or Vector. Controls which factors are visualized.
#' @param n_label Integer. Number of top genotypes to label.
#' @param highlight Character vector. Names of specific genotypes to highlight.
#' @param ... Additional graphical parameters.
#'
#' @export
plot.fa_asreml <- function(x, type = "fast", factor = NULL, n_label = 5, highlight = NULL, ...) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    if (type == "fast") {
        if (!is.null(x$fast)) {
            # Delegate to the FAST object plot method
            op_rmsd_obj <- list(selection = x$fast, source_model = x)
            class(op_rmsd_obj) <- "fast_selection"
            plot(op_rmsd_obj, type = "op_rmsd", ...)
        } else {
            stop("FAST indices not available.")
        }
    } else if (type == "latent_reg") {
        # Delegate to FAST object latent_reg
        # Note: My internal FAST plot implementation for latent_reg is a "placeholder" warning
        # So maybe I should keep the robust implementation HERE for now, OR move it to class_fast_selection.R
        # Given I implemented a placeholder there, I will keep the robust one here as legacy .plot_reg
        # UNLESS I update class_fast_selection.R.
        # For now, I'll use the local internal one to ensure functionality works.
        .plot_reg(x, factor, highlight, n_label)
    } else if (type == "heatmap") {
        .plot_heat(x)
    } else if (type == "biplot") {
        .plot_biplot(x, if (is.null(factor)) c(1, 2) else factor, highlight)
    } else if (type == "vaf") {
        .plot_vaf(x)
    } else if (type == "d_opt") {
        if (exists("get_d_optimality")) {
            .plot_dopt(get_d_optimality(x))
        } else {
            warning("get_d_optimality function not found.")
        }
    } else if (type == "diff") {
        if (exists("get_i_classes")) {
            .plot_diff(get_i_classes(x, if (is.null(factor)) 2 else factor), n_label, highlight)
        } else {
            warning("get_i_classes function not found.")
        }
    } else {
        stop("Unknown type.")
    }
}

# --- Plot Internals ---

.plot_heat <- function(x) {
    cor <- x$matrices$Cor
    p <- ncol(cor)
    rev <- cor[, p:1]
    par(mar = c(6, 6, 4, 2))
    image(1:p, 1:p, rev, axes = F, col = hcl.colors(20, "RdBu", rev = T), breaks = seq(-1, 1, l = 21), main = "Correlation")
    axis(1, at = 1:p, labels = colnames(cor), las = 2, cex.axis = 0.7)
    axis(2, at = 1:p, labels = rev(rownames(cor)), las = 1, cex.axis = 0.7)
    for (i in 1:p) for (j in 1:p) if (rev[i, j] < 0.99) text(i, j, sprintf("%.2f", rev[i, j]), cex = 0.6)
}

.plot_reg <- function(x, fac, h, n) {
    if (is.null(x$scores)) stop("No Genotype scores available for Latent Regression.")

    k <- x$meta$k
    if (is.null(fac)) fac <- 1:k
    if (length(fac) > 1) par(mfrow = c(ceiling(length(fac) / 2), min(length(fac), 2)))
    L <- x$loadings$rotated
    S <- x$scores$rotated
    Uc <- S %*% t(L)
    if (!is.null(x$fast)) {
        tgt <- unique(c(head(x$fast$Genotype, n), h))
    } else {
        tgt <- head(rownames(S), n)
    }

    for (f in fac) {
        if (f > k) next
        Y <- if (f > 1) Uc - (S[, 1:(f - 1)] %*% t(L[, 1:(f - 1)])) else Uc
        xv <- L[, f]
        ysub <- Y[tgt, , drop = F]

        par(mar = c(4, 4, 3, 1))
        plot(1, type = "n", xlim = range(xv), ylim = range(ysub), xlab = paste("Load", f), ylab = "Effect", main = paste("Factor", f))
        grid()
        abline(h = 0, lty = 2)
        axis(1, at = xv, labels = F, tck = -0.02, col = "green")

        for (g in tgt) {
            sl <- S[g, f]
            cc <- if (g %in% h) "red" else "navy"
            points(xv, Y[g, ], pch = 19, col = adjustcolor(cc, 0.3), cex = 0.7)
            abline(0, sl, col = cc)
            text(max(xv), sl * max(xv), labels = g, pos = 4, cex = 0.7, col = cc)
        }
    }
}

.plot_biplot <- function(x, fac, h) {
    if (is.null(x$scores)) stop("No Genotype scores available for Biplot.")

    L <- x$loadings$rotated[, fac]
    S <- x$scores$rotated[, fac]
    sf <- max(abs(S)) / max(abs(L)) * 0.8
    Ls <- L * sf
    lim <- range(rbind(Ls, S))

    par(mar = c(5, 5, 4, 2))
    plot(1, type = "n", xlim = lim, ylim = lim, asp = 1, xlab = paste("F", fac[1]), ylab = paste("F", fac[2]), main = "Biplot")
    grid()
    abline(h = 0, v = 0, lty = 2, col = "grey")
    arrows(0, 0, Ls[, 1], Ls[, 2], col = "darkgreen", length = 0.1)
    text(Ls, labels = rownames(L), col = "darkgreen", pos = 4, cex = 0.7)

    bg <- rep("grey90", nrow(S))
    pch <- rep(21, nrow(S))
    is_h <- rownames(S) %in% h
    bg[is_h] <- "red"
    pch[is_h] <- 23
    points(S, pch = pch, bg = bg, cex = 0.8)
    if (any(is_h)) text(S[is_h, ], labels = rownames(S)[is_h], pos = 3, cex = 0.7)

    hpts <- chull(S)
    lines(S[c(hpts, hpts[1]), ], col = "navy", lty = 2)
}

.plot_vaf <- function(x) {
    df <- x$var_comp$vaf[order(x$var_comp$vaf$VAF), ]
    col <- ifelse(df$VAF < 50, "red", "blue")
    par(mar = c(5, 7, 4, 2))
    barplot(df$VAF, names.arg = df$Site, horiz = T, las = 1, col = col, xlab = "VAF %", main = "Site Quality")
    abline(v = 80, lty = 2)
}

.plot_dopt <- function(d, threshold = 1.0) {
    df <- d$site_impact
    df <- df[order(df$Impact_Pct), ]
    cols <- ifelse(df$Impact_Pct < threshold, "#e74c3c", "#27ae60")
    par(mar = c(5, 8, 4, 2))
    bp <- barplot(df$Impact_Pct,
        names.arg = df$Site, horiz = TRUE,
        las = 1, col = cols, border = NA,
        xlab = "% Information Loss if Removed",
        main = "Network Efficiency (D-Opt Contribution)"
    )
    abline(v = threshold, lty = 2, col = "grey40")
    text(
        x = df$Impact_Pct, y = bp, label = sprintf("%.1f%%", df$Impact_Pct),
        pos = 4, cex = 0.7, xpd = TRUE
    )
}

.plot_diff <- function(res, n, h) {
    df <- res$gen_effects
    df <- na.omit(df)
    tgt <- unique(c(head(df$Genotype, n), tail(df$Genotype, n), h))
    plot_df <- df[df$Genotype %in% tgt, ]
    yl <- range(c(plot_df$Pred_Neg, plot_df$Pred_Pos))
    par(mar = c(4, 4, 3, 6), xpd = F)
    plot(1, type = "n", xlim = c(0.8, 2.2), ylim = yl, xaxt = "n", ylab = "Genetic Value", xlab = "", main = "Interaction Classes")
    axis(1, at = 1:2, labels = c("Class (-)", "Class (+)"))
    for (i in 1:nrow(plot_df)) {
        g <- plot_df$Genotype[i]
        cc <- if (g %in% h) "red" else "grey50"
        segments(1, plot_df$Pred_Neg[i], 2, plot_df$Pred_Pos[i], col = cc, lwd = if (g %in% h) 2 else 1)
        text(1, plot_df$Pred_Neg[i], g, pos = 2, cex = 0.7, col = cc, xpd = T)
        text(2, plot_df$Pred_Pos[i], g, pos = 4, cex = 0.7, col = cc, xpd = T)
    }
}
