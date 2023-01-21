#' Export promoter sequences in a fasta file.
#'
#' @param .sequences 
#' @param .output_file 
#'
#' @return A fasta file.
#' @export
#'
#' @examples
write_fasta <- function(.sequences, .output_file) {
  Biostrings::writeXStringSet(.sequences, filepath = .output_file)
}