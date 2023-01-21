#' Get the distance from the gene-of-interest transcription start site to the 
#' closest upstream gene.
#'
#' @inheritParams get_promoters
#'
#' @return A data frame with the closest upstream gene ID and distance from the 
#'         gene-of-interest transcription start site to the closest upstream 
#'         gene.
#' @export
#' 
get_promoter_distances <- function(
    .annotations,
    .pb        = NULL,
    .pb_format = ":what - [:bar] :percent (:spin)",
    .locus_var = "locus_tag",
    .chr_var   = "chr"
) {
  # Initiate (or not) the progress bar - can also receive external object
  if(is.null(.pb)) {
    .pb <- progress::progress_bar$new(
      format = .pb_format,
      total  = .annotations %>% 
        dplyr::select(dplyr::all_of(c(.locus_var, .chr_var))) %>% 
        unique() %>% 
        nrow()
    )
  }
  
  # Start the progress bar
  .pb$tick(0)
  
  # THIS MUST BE PARALLELIZED
  .distances <- .annotations %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(.locus_var, .chr_var)))) %>%
    dplyr::group_modify(.f = function(.each_annotation, .keys) {
      # Assign the current locus and chromosome names to objects ("character")
      .curr_locus <- .keys %>% dplyr::pull(.locus_var)
      .curr_chr   <- .keys %>% dplyr::pull(.chr_var)
      
      # Subset relevant loci - a.k.a. "loci of interest"
      .loci_OI <- .annotations %>%
        tidyr::pivot_longer(cols = c("begin", "end")) %>%
        # get all loci in the current chromosome
        #subset(chr == .keys$chr) %>%
        subset(get(.chr_var) == .curr_chr) %>%
        # get all loci but the current locus
        #subset(locus_tag != .keys$locus_tag) %>%
        subset(get(.locus_var) != .curr_locus) %>%
        # get all loci upstream the current locus
        subset(
          dplyr::case_when(
            .each_annotation$strand ==  1 ~ value < .each_annotation$begin,
            .each_annotation$strand == -1 ~ value > .each_annotation$end,
            TRUE               ~ NA
          )
        )
      
      # Initiate the output object
      .each_output <- data.frame(
        closest_locus = NA_character_, 
        dist          = NA
      )
      
      # If there are "loci of interest", updates the output object with the 
      # upstream distance to the closest gene (still, in a given chromosome). 
      # The if() below handles e.g. the first locus in each chromosome which 
      # doesn't have other loci upstream.
      if(nrow(.loci_OI) != 0) {
        .each_output <- .loci_OI %>%
          # get the reference coordinate to the upstream locus closest to the 
          # current locus
          subset(dplyr::case_when(
            .each_annotation$strand ==  1 ~ value == max(value, na.rm = TRUE),
            .each_annotation$strand == -1 ~ value == min(value, na.rm = TRUE)
          )) %>%
          # Calculate the distance to the closest upstream locus we got above
          dplyr::mutate(dist = dplyr::case_when(
            .each_annotation$strand ==  1 ~ .each_annotation$begin - value,
            .each_annotation$strand == -1 ~ value - .each_annotation$end
          )) %>%
          dplyr::select(dplyr::all_of(c(.locus_var, "dist"))) %>%
          dplyr::rename(closest_locus = dplyr::all_of(.locus_var))
      }
      
      # progress bar oogie boogie
      .pb$tick(
        tokens = list(what = paste(.curr_chr, .curr_locus, sep = " - "))
      )
      
      .each_output
    })
  
  # End progress bar
  rm(.pb)
  
  .distances
}