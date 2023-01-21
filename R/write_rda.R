#' Export extracted promoter sequence results in a R Data File.
#' 
#' @inheritParams get_promoters
#'
#' @return A R Data File `rda`.
#' @export
#'
write_rda <- function(.sequences, .output_file) {
  save(.sequences, file = .output_file)
}