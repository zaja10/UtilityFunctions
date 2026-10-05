#' Standard Genetics Theme
#'
#' A ggplot2 theme optimized for genetics and breeding plots.
#' Features clean white background, grid lines for readability, and readable fonts.
#'
#' @param base_size Base font size (default 12).
#' @param base_family Base font family.
#' @return A ggplot2 theme object.
#' @export
theme_genetics <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Text
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1), color = "grey40", margin = ggplot2::margin(b = 10)),
      axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9)),
      axis.text = ggplot2::element_text(color = "grey30", size = ggplot2::rel(0.8)),

      # Background & Grid
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      panel.grid.minor = ggplot2::element_blank(),

      # Legend
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.8)),

      # Facets
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9))
    )
}
