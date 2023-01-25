#' plot_promoter_maps
#'
#' @inheritParams get_promoters
#' @param .promoter_sizes {TODO}
#' @param .relative_x     Whether to transform X-axis scale do indicate distance 
#'                        to each locus start coordinate.
#'
#' @export
#'
plot_promoter_maps <- function(
    .promoter_sizes,
    .annotations,
    .relative_x = TRUE
  ) {
  .plot_data <- .promoter_sizes %>%
    get_plotdata_locus(.annotations) 
  
  if(.relative_x) {
    .plot_data <- .plot_data %>%
      dplyr::mutate(dplyr::across(
        .cols = dplyr::all_of(
          c("begin", "end", "closest_begin", "closest_end")
        ),
        .fns  = function(.col) {
          dplyr::case_when(
            strand ==  1 ~ .col - begin,
            strand == -1 ~ .col - end,
          )
        }
      ))
  }
  
  .plot_data %>%
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