#' @title Export promoter sequences in a fasta file.
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' @param .sequences   {TODO}
#' @param .output_file {TODO}
#' @param ...          Future-proof stuff.
#'
#' @return A fasta file.
#' @export
#'
write_fasta <- function(.sequences, .output_file, ...) {
  Biostrings::writeXStringSet(.sequences, filepath = .output_file)
}