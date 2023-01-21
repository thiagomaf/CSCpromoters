#' Export promoter sequences in a fasta file.
#' 
#' @inheritParams get_promoters
#'
#' @return A fasta file.
#' @export
#'
write_fasta <- function(.sequences, .output_file) {
  Biostrings::writeXStringSet(.sequences, filepath = .output_file)
}