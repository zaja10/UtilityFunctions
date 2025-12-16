#' Standard Genetics Theme
#'
#' A ggplot2 theme optimized for genetics and breeding plots.
#' Features clean white background, grid lines for readability, and readable fonts.
#'
#' @param base_size Base font size (default 12).
#' @param base_family Base font family.
#' @return A ggplot2 theme object.
#' @import ggplot2
#' @export
theme_genetics <- function(base_size = 12, base_family = "") {
    theme_minimal(base_size = base_size, base_family = base_family) +
        theme(
            # Text
            plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0),
            plot.subtitle = element_text(size = rel(1), color = "grey40", margin = margin(b = 10)),
            axis.title = element_text(face = "bold", size = rel(0.9)),
            axis.text = element_text(color = "grey30", size = rel(0.8)),

            # Background & Grid
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            panel.grid.major = element_line(color = "grey92"),
            panel.grid.minor = element_blank(),

            # Legend
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = rel(0.8)),

            # Facets
            strip.background = element_rect(fill = "grey95", color = NA),
            strip.text = element_text(face = "bold", size = rel(0.9))
        )
}
