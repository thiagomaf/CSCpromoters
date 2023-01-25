#' Set the size of promoters based on the distance to the closest upstream gene 
#' and the user-defined minimum and maximum length.
#'
#' @inheritParams get_promoters
#' @inheritDotParams get_promoters
#'
#' @return A data frame with filtered genes that have promoter lengths longer 
#'         than the user-defined minimum length, and a new column containing the 
#'         length of the final promoter sizes.
#' @export
#'
set_promoter_sizes <- function(
    .distances,
    .min_size = 100,
    .max_size = 2000,
    .dist_var = "dist",
    .size_var = "promoter_size",
    ...
) {
  .value_var <- "dist"
  
  .distances <- .distances  %>%
    dplyr::ungroup() %>% # breaks the magig if removed!
    subset(get(.dist_var) >= .min_size) %>%
    dplyr::mutate(
      #promoter_size = dplyr::case_when(
      !!.size_var := dplyr::case_when(
        get(.dist_var) > .max_size ~ .max_size,
        #TRUE ~ dist
        TRUE ~ get(.dist_var)
      )
    )
  
  .distances
}