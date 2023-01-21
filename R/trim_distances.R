#' Trim the length of promoters based on the distance to the closest upstream gene and the user-defined minimum and maximum length.
#'
#' @param .distances 
#' @param .min_size 
#' @param .max_size 
#'
#' @return A data frame with filtered genes that have promoter lengths longer than the user-defined minimum length, and a new column containing the length of the final promoter sizes.
#' @export
#'
#' @examples
trim_distances <- function(.distances, .min_size = 100, .max_size = 2000) {
  .distances <- .distances  %>%
    dplyr::ungroup() %>% # breaks the magig if removed!
    subset(dist >= .min_size) %>%
    dplyr::mutate(
      promoter_size = dplyr::case_when(
        dist > .max_size ~ .max_size,
        TRUE ~ dist
      )
    )
  
  .distances # absolutelly useless return but makes it clearer what is happening
}