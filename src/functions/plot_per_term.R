

#' Plot a dotplot from GSEA results with a custom representative term for group of gene sets
#'
#' This function generates a dot plot visualizing normalized enrichment scores (NES)
#' and adjusted p-values for a specific GSEA term (e.g., a biological process or pathway).
#' The plot displays each subterm (Description_str) on the y-axis, colored by the
#' direction of enrichment (positive or negative NES), and sized by the adjusted p-value.
#'
#' @param term Character string specifying the custom representative term for group of gene sets to plot.
#'   Must match a value in the `New_representative_term` column of `data`.
#' @param data A data frame containing GSEA results with the following required columns:
#'   \describe{
#'     \item{New_representative_term}{Grouping term name for sets of GSEA results.}
#'     \item{Description_str}{Text description of each GSEA pathway or subterm.}
#'     \item{NES}{Normalized enrichment score for each subterm.}
#'     \item{p.adjust}{Adjusted p-value for each subterm.}
#'   }
#' @param x_lim Numeric value defining the symmetric x-axis limit for NES values.
#'   Typically set as `max(abs(range(data$NES)))` to ensure balanced scaling.
#'
#' @details
#' The function filters the input data to a single GSEA term and plots NES
#' against the subterm (`Description_str`). Points are colored by enrichment direction
#' (`positive` or `negative`) and scaled by adjusted p-value. A vertical line at `x = 0`
#' marks the boundary between positive and negative enrichment.
#'
#' The plot aesthetics (color palette, text sizes, aspect ratio, and background)
#' are optimized for inclusion in multi-panel figures.
#'
#' @return
#' A `ggplot` object representing the GSEA dot plot for the specified term.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' gsea_data <- readr::read_csv("data/VectorDb_HKEC__vs__Naive.csv")
#' x_lim <- max(abs(range(gsea_data$NES)))
#' p <- plot_per_term("Immunity", gsea_data, x_lim)
#' print(p)
#' }
#'
#' @seealso
#' [ggplot2::ggplot()], [cowplot::align_patches()], [forcats::fct_reorder()]
#'
#' @export
#' 
plot_per_term <- function(term, data, x_lim) {
  df_sub <- data %>%
    filter(New_representative_term == term) %>%
    mutate(
      Description_str = fct_reorder(Description_str, NES),
      direction = factor(
        ifelse(NES < 0, "negative", "positive"),
        levels = c("negative", "positive")
      )
    )
  wrapped_term <- str_wrap(term, width = 26)
  
  p <- ggplot(df_sub, aes(x=NES, y=Description_str, size=p.adjust, color=direction)) +
    geom_point(show.legend = FALSE) +
    scale_color_manual(values=c("negative" = "#377eb8", "positive" = "#e41a1c"),
                       breaks = c("negative", "positive"),
                       drop = FALSE) +
    scale_size_continuous(
      range = c(7, 1),
      breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
      labels = c("≤ 0.01", "0.02", "0.03", "0.04", "0.05"),
      limits = c(0, 0.05),
      name = "Adjusted p-value"
    ) +
    geom_vline(xintercept=0, linetype='solid', color='gray') +
    labs(
      title = wrapped_term,
      x = "NES",
      y = "",
      color = "Direction"
    ) +
    theme(
      plot.title = element_text(hjust=0, face='bold', color = "black", size = 14),
      axis.text.x = element_text(angle = 0, hjust = 1, margin = margin(t = 10)),
      panel.background = element_rect(fill = "#f0f0f0", color = "grey80"),
      text = element_text(size=14, color = "black"),
      panel.grid.major.y = element_line(color = "grey30", linetype = "dotted"),
      panel.grid.major.x = element_line(color = "white", linetype = "solid", size = 0.5),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
    ) +
    coord_cartesian(clip = "off", xlim = c(-x_lim, x_lim)) +
    theme(aspect.ratio = 0.5 * nrow(df_sub) / 6)
  
  return(p)
}