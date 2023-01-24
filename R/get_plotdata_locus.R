#' get_plotdata_locus
#'
#' @inheritParams get_promoters
#' @param .promoter_sizes {TODO}
#'
#' @export
#'
get_plotdata_locus <- function(
    .promoter_sizes,
    .annotations,
    .locus_var   = "locus_tag",
    .chr_var     = "chr",
    .closest_var = "closest_locus",
    .strand_var  = "strand",
    .start_var   = "begin",
    .end_var     = "end"
  ) {
  .promoter_sizes %>%
    dplyr::select(dplyr::all_of(
      c(.chr_var, .locus_var, .closest_var, "dist", "promoter_size")
    )) %>% 
    dplyr::left_join(.annotations, by = c(.chr_var, .locus_var)) %>% 
    dplyr::left_join(
      .annotations %>% 
        dplyr::rename(c(
          "closest_locus"  = all_of(.locus_var),
          "closest_strand" = all_of(.strand_var),
          "closest_begin"  = all_of(.start_var),
          "closest_end"    = all_of(.end_var)
        )),
      by = c(.chr_var, .closest_var)
    )
}