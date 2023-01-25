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
    .dist_var    = "dist",
    .size_var    = "promoter_size",
    .strand_var  = "strand",
    .start_var   = "begin",
    .end_var     = "end"
) {
  .promoter_sizes %>%
    dplyr::select(dplyr::all_of(
      c(.chr_var, .locus_var, .closest_var, .dist_var, .size_var)
    )) %>% 
    dplyr::left_join(.annotations, by = c(.chr_var, .locus_var)) %>% 
    dplyr::left_join(
      .annotations %>%
        dplyr::rename(
          !!.closest_var := dplyr::all_of(.locus_var)
        ) %>% 
        dplyr::rename(
          !!paste0("closest_", .strand_var) := dplyr::all_of(.strand_var)
        ) %>% 
        dplyr::rename(
          !!paste0("closest_", .start_var) := dplyr::all_of(.start_var)
        ) %>% 
        dplyr::rename(
          !!paste0("closest_", .end_var) := dplyr::all_of(.end_var)
        ),
      by = c(.chr_var, .closest_var)
    )
}