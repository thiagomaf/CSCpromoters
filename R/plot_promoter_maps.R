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
    .locus_var  = "locus_tag",
    .strand_var = "strand",
    .start_var  = "begin",
    .end_var    = "end",
    .relative_x = TRUE
  ) {
  .plot_data <- .promoter_sizes %>%
    get_plotdata_locus(.annotations) 
  
  if(.relative_x) {
    .plot_data <- .plot_data %>%
      dplyr::mutate(dplyr::across(
        .cols = dplyr::all_of(
          c(
            .start_var,
            .end_var,
            paste0("closest_", .start_var),
            paste0("closest_", .end_var)
          )
        ),
        .fns  = function(.col) {
          dplyr::case_when(
            get(.strand_var) ==  1 ~ .col - get(.start_var),
            get(.strand_var) == -1 ~ .col - get(.end_var),
          )
        }
      ))
  }
  
  .x    <- "x"
  .xend <- "xend"
  
  .plot_data %>%
    ggplot2::ggplot(ggplot2::aes(y = get(.locus_var))) +
    gggenes::theme_genes() +
    gggenes::geom_gene_arrow(ggplot2::aes(
      xmin    = get(.start_var),
      xmax    = get(.end_var),
      forward = get(.strand_var),
      fill    = "reference locus"
    )) +
    gggenes::geom_gene_arrow(ggplot2::aes(
      xmin    = get(paste0("closest_", .start_var)),
      xmax    = get(paste0("closest_", .end_var)),
      forward = get(paste0("closest_", .strand_var)),
      fill    = "closest locus"
    )) +
    ggplot2::geom_segment(
      mapping   = ggplot2::aes(
        x     = get(.x),
        xend  = get(.xend),
        y     = get(.locus_var),
        yend  = get(.locus_var),
        color = "promoter"
      ),
      linewidth = 1.5,
      data      = . %>% get_plotdata_promoter()
    ) +
    ggplot2::facet_wrap(
      stats::as.formula(paste0("~ ", .locus_var)), scales = "free", ncol = 1
    ) +
    ggplot2::scale_fill_brewer(palette = "Set3")
}
