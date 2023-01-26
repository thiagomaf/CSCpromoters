#' @title get_plotdata_promoter
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' @inheritParams get_promoters
#' @param .promoter_sizes {TODO}
#'
#' @export
#'
get_plotdata_promoter <- function(
    .promoter_sizes,
    .locus_var   = "locus_tag",
    .chr_var     = "chr",
    .size_var    = "promoter_size",
    .strand_var  = "strand",
    .start_var   = "begin",
    .end_var     = "end"
) {
  .promoter_sizes %>% 
    dplyr::select(dplyr::all_of(
      c(.chr_var, .locus_var, .strand_var, .start_var, .end_var, .size_var)
    )) %>% 
    dplyr::mutate(x = dplyr::case_when(                        # 'x' on the fly?
      get(.strand_var) ==  1 ~ get(.start_var),
      get(.strand_var) == -1 ~ get(.end_var)
    )) %>%
    dplyr::mutate(xend = dplyr::case_when(                  # 'xend' on the fly?
      get(.strand_var) ==  1 ~ get(.start_var) - get(.size_var),
      get(.strand_var) == -1 ~ get(.end_var) + get(.size_var)
    ))
}