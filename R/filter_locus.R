# This is one of the few cases in which we will have more than one function per
# file. In this case, the `filter_locus()` function uses two auxiliary functions
# `filter_locus.numeric()` and `filter_locus.character()` for our convenience. 
# Each function will have its own `roxygen2` comments.


#' @title Filter annotation table by loci names or loci indexes
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' @inheritParams get_promoters
#' @inheritDotParams get_promoters
#'
#' @return An annotations table object with the same class as `.annotations`
#' 
#' @export
filter_locus <- function(
    .annotations,
    .keep       = NULL,
    .locus_var  = "locus_tag",
    ...
) {
  # Check if .keep argument was set.
  if(is.null(.keep)) {
    stop(
      "CSCpromoters: .keep argument is NULL. Use `help(filter_locus)` for help."
    )
  }

  # Use the correct auxiliary function depending on the class of .keep argument.
  switch (
    class(.keep),
    numeric   = .annotations %>% filter_locus.numeric(.keep),
    integer   = .annotations %>% filter_locus.numeric(.keep),
    character = .annotations %>% filter_locus.character(.keep),
    # If the argument .keep is not character or numeric, complain!
    stop(
      paste(
        "CSCpromoters: .keep argument must either be a character vector of",
        "loci names or an interger vector of indexes. Use `help(filter_locus)`",
        "for help."
      )
    )
  )
}



#' Filter annotation table by loci indexes
#'
#' @inheritParams filter_locus
#' @param .keep A `numeric` vector of indexes (of the `.annotations` table).
#'
#' @export
filter_locus.numeric <- function(.annotations, .keep = NULL) {
  # .annotations %>%
  #   dplyr::slice(.keep)
  
  .annotations %>% 
    magrittr::set_attr(which = "keep", value = .keep)
}



#' Filter annotation table by loci names
#'
#' @inheritParams filter_locus
#' @param .keep A `character` vector of loci names (present in the 
#' `.annotations` table).
#'              
#' @export
filter_locus.character <- function(
    .annotations,
    .keep,
    .locus_var  = "locus_tag"
) {
  # .annotations %>%
  #   subset(get(.locus_var) %in% .keep)
  
  .keep_indexes <- .annotations %>%
    dplyr::pull(.locus_var) %in% .keep %>%
    which()
  
  .annotations %>% 
    magrittr::set_attr(
      which = "keep", 
      value = .keep_indexes
    )
}