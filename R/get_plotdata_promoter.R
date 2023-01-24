#' get_plotdata_promoter
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
    .strand_var  = "strand",
    .start_var   = "begin",
    .end_var     = "end"
) {
  .promoter_sizes %>% 
    dplyr::select(dplyr::all_of(
      c(.chr_var, .locus_var, .strand_var, .start_var, .end_var, "promoter_size")
    )) %>% 
    dplyr::mutate(x = dplyr::case_when(
      strand ==  1 ~ begin, # on the fly
      strand == -1 ~ end    # on the fly
    )) %>%
    dplyr::mutate(xend = dplyr::case_when(
      strand ==  1 ~ begin - promoter_size, # on the fly
      strand == -1 ~ end + promoter_size    # on the fly
    ))
}