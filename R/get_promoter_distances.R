#' @title Get the distance from the gene-of-interest transcription start site to the 
#'        closest upstream gene.
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' @inheritParams get_promoters
#' @inheritDotParams get_promoters
#'
#' @return A data frame with the closest upstream gene ID and distance from the 
#'         gene-of-interest transcription start site to the closest upstream 
#'         gene.
#' @export
#' 
get_promoter_distances <- function(
    .annotations,
    .pb          = NULL,
    .pb_format   = ":what - [:bar] :percent (:spin)",
    .locus_var   = "locus_tag",
    .chr_var     = "chr",
    .closest_var = "closest_locus",
    .dist_var    = "dist",
    .strand_var  = "strand",
    .start_var   = "begin",
    .end_var     = "end",
    ...
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
  
  # Query which loci to calculate distances from
  .to_keep <- .annotations %>% attr(which = "keep")
  if(is.null(.to_keep)) {
    .to_keep <- 1:nrow(.annotations)
  }
  
  # THIS MUST BE PARALLELIZED
  # Calculate distances
  .distances <- .annotations %>%
    dplyr::slice(.to_keep) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(.locus_var, .chr_var)))) %>%
    dplyr::group_modify(.f = function(.each_annotation, .keys) {
      # Assign the current locus and chromosome names to objects ("character")
      .curr_locus <- .keys %>% dplyr::pull(.locus_var)
      .curr_chr   <- .keys %>% dplyr::pull(.chr_var)
      
      .curr_strand <- .each_annotation %>% dplyr::pull(.strand_var)
      .curr_begin  <- .each_annotation %>% dplyr::pull(.start_var)
      .curr_end    <- .each_annotation %>% dplyr::pull(.end_var)
      
      # Subset relevant loci - a.k.a. "loci of interest"
      .loci_OI <- .annotations %>%
        tidyr::pivot_longer(cols = c(.start_var, .end_var)) %>%
        # get all loci in the current chromosome
        subset(get(.chr_var) == .curr_chr) %>%
        # get all loci but the current locus
        subset(get(.locus_var) != .curr_locus) %>%
        # get all loci upstream the current locus
        subset(
          dplyr::case_when(
            .curr_strand ==  1 ~ value < .curr_begin,
            .curr_strand == -1 ~ value > .curr_end,
            TRUE               ~ NA
          )
        )
      
      # Initiate the output object
      .each_output <- data.frame(
        closest_locus = NA_character_,
        dist          = NA
      ) %>%
        dplyr::rename(!!.closest_var := dplyr::all_of("closest_locus")) %>%
        dplyr::rename(!!.dist_var    := dplyr::all_of("dist"))

      # If there are "loci of interest", updates the output object with the 
      # upstream distance to the closest gene (still, in a given chromosome). 
      # The if() below handles e.g. the first locus in each chromosome which 
      # doesn't have other loci upstream.
      if(nrow(.loci_OI) != 0) {
        .each_output <- .loci_OI %>%
          # get the reference coordinate to the upstream locus closest to the 
          # current locus
          subset(dplyr::case_when(
            .curr_strand ==  1 ~ value == max(value, na.rm = TRUE),
            .curr_strand == -1 ~ value == min(value, na.rm = TRUE)
          )) %>%
          # Calculate the distance to the closest upstream locus we got above
          dplyr::mutate(!!.dist_var := dplyr::case_when(
            .curr_strand ==  1 ~ .curr_begin - value,
            .curr_strand == -1 ~ value - .curr_end
          )) %>%
          dplyr::mutate(!!.dist_var := as.double(get(.dist_var))) %>%
          dplyr::select(dplyr::all_of(c(.locus_var, .dist_var))) %>%
          dplyr::rename(!!.closest_var := dplyr::all_of(.locus_var))
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