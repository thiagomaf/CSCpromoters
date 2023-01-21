#' Export extracted promoter sequence results in a R Data File.
#'
#' @param .sequences 
#' @param .output_file 
#'
#' @return A R Data File `rda`.
#' @export
#'
#' @examples
write_rda <- function(.sequences, .output_file) {
  save(.sequences, file = .output_file)
}