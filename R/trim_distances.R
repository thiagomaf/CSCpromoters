#' Trim the length of promoters based on the distance to the closest upstream gene and the user-defined minimum and maximum length.
#'
#' @inheritParams get_promoters
#' @inheritDotParams get_promoters
#'
#' @return A data frame with filtered genes that have promoter lengths longer than the user-defined minimum length, and a new column containing the length of the final promoter sizes.
#' @export
#'
trim_distances <- function(.distances, .min_size = 100, .max_size = 2000, ...) {
  .value_var <- "dist"
  
  .distances <- .distances  %>%
    dplyr::ungroup() %>% # breaks the magig if removed!
    #subset(dist >= .min_size) %>%
    subset(get(.value_var) >= .min_size) %>%
    dplyr::mutate(
      promoter_size = dplyr::case_when(
        #dist > .max_size ~ .max_size,
        get(.value_var) > .max_size ~ .max_size,
        TRUE ~ dist
      )
    )
  
  .distances # absolutelly useless return but makes it clearer what is happening
}