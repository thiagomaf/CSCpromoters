#' plot_promoter_maps
#'
#' @inheritParams get_promoters
#' @param .promoter_sizes {TODO}
#'
#' @export
#'
plot_promoter_maps <- function(.promoter_sizes, .annotations) {
  .promoter_sizes %>%
    get_plotdata_locus(.annotations) %>%
    ggplot2::ggplot(ggplot2::aes(y = locus_tag)) +
    gggenes::theme_genes() +
    gggenes::geom_gene_arrow(ggplot2::aes(
      xmin    = begin,
      xmax    = end,
      forward = strand,
      fill    = "reference locus"
    )) +
    gggenes::geom_gene_arrow(ggplot2::aes(
      xmin    = closest_begin,
      xmax    = closest_end,
      forward = closest_strand,
      fill    = "closest locus"
    )) +
    ggplot2::geom_segment(
      mapping   = ggplot2::aes(
        x     = x,
        xend  = xend,
        y     = locus_tag,
        yend  = locus_tag,
        color = "promoter"
      ),
      linewidth = 1.5,
      data      = . %>% get_plotdata_promoter()
    ) +
    ggplot2::facet_wrap(~ locus_tag, scales = "free", ncol = 1) +
    ggplot2::scale_fill_brewer(palette = "Set3")
}