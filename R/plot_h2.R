
#' Plot Heritability/Reliability Comparison
#'
#' @description
#' Visualizes the output of \code{compare_h2} using a grouped bar chart.
#' Compares different metrics (Cullis, Oakey, Standard) across sites.
#'
#' @param x An object of class \code{h2_comparison} produced by \code{compare_h2()}.
#' @param ... Additional arguments passed to \code{barplot}.
#'
#' @importFrom graphics barplot legend par text abline box
#' @importFrom grDevices hcl.colors
#' @export
plot.h2_comparison <- function(x, ...) {
  # x has columns: Site, Method, Value, Vg
  
  # Check structure
  if(!all(c("Site", "Method", "Value") %in% names(x))) stop("Invalid h2_comparison object.")
  
  # Pivot to wide format for barplot matrix: Rows = Methods, Cols = Sites
  # Use base R reshape
  wide <- reshape(x, idvar = "Site", timevar = "Method", v.names = "Value", direction = "wide")
  
  # Clean column names (remove "Value.")
  names(wide) <- gsub("Value\\.", "", names(wide))
  
  # Matrix for barplot
  methods <- unique(x$Method)
  sites <- unique(x$Site)
  
  # Ensure matching order
  mat <- as.matrix(wide[, methods, drop=FALSE])
  rownames(mat) <- wide$Site
  
  # Sort by first method values (descending)
  ord <- order(mat[,1], decreasing = FALSE) # Ascending for horiz plot (top is highest)
  mat <- mat[ord, , drop=FALSE]
  
  # Plotting
  old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
  
  # Colors
  cols <- hcl.colors(length(methods), "Zissou 1")
  
  # Layout
  par(mar = c(5, 8, 4, 10), xpd = TRUE) # Right margin for legend
  
  bp <- barplot(t(mat), 
          beside = TRUE, 
          horiz = TRUE, 
          col = cols, 
          border = NA,
          las = 1,
          xlab = "Reliability / Heritability",
          main = "Site Quality Comparison",
          xlim = c(0, 1),
          ...)
  
  # Gridlines
  abline(v = seq(0, 1, 0.2), col = "white", lwd = 0.5)
  box()
  
  # Legend
  legend("right", 
         inset = c(-0.35, 0),
         legend = methods, 
         fill = cols, 
         border = NA, 
         bty = "n", 
         title = "Metric")
         
  # Add Threshold
  abline(v = 0.05, col="red", lty=2) # Random threshold often used? Or 0.8?
  # Let's not add arbitrary lines
}
