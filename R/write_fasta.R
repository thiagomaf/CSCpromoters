#' Export promoter sequences in a fasta file.
#' 
#' @param .sequences   {TODO}
#' @param .output_file {TODO}
#'
#' @return A fasta file.
#' @export
#'
write_fasta <- function(.sequences, .output_file) {
  Biostrings::writeXStringSet(.sequences, filepath = .output_file)
}