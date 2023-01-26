#' @title Export extracted promoter sequence results in a R Data File.
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' @param .sequences   {TODO}
#' @param .output_file {TODO}
#' @param ...          Future-proof stuff.
#'
#' @return A R Data File `rda`.
#' @export
#'
write_rda <- function(.sequences, .output_file, ...) {
  save(.sequences, file = .output_file)
}