#' Commercial Publication Theme for ggplot2
#'
#' A clean, serif-based theme suitable for scientific journals and commercial reports.
#'
#' @param base_size Numeric. Base font size. Default 12.
#' @param base_family Character. Font family. Default "serif".
#' @importFrom ggplot2 theme_classic theme element_text element_line element_rect rel margin
#' @export
theme_genetics <- function(base_size = 12, base_family = "serif") {
    ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
        ggplot2::theme(
            # Text
            plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = ggplot2::rel(1.2)),
            plot.subtitle = ggplot2::element_text(color = "grey40", hjust = 0, margin = ggplot2::margin(b = 10)),
            plot.caption = ggplot2::element_text(color = "grey60", hjust = 1),

            # Axes
            axis.line = ggplot2::element_line(color = "black", linewidth = 0.8),
            axis.ticks = ggplot2::element_line(color = "black"),
            axis.text = ggplot2::element_text(color = "black"),

            # Legend
            legend.position = "bottom",
            legend.title = ggplot2::element_text(face = "bold"),
            legend.background = ggplot2::element_rect(fill = NA, color = NA),

            # Facets
            strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
            strip.text = ggplot2::element_text(face = "bold")
        )
}
