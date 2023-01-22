#' Extract annotation table from TxDB object
#'
#' @inheritParams get_promoters

#' @export
get_txdb_annotation <- function(
    .txdb,
    .locus_var  = "locus_tag",
    .chr_var    = "chr",
    .strand_var = "strand",
    .start_var  = "begin",
    .end_var    = "end"
) {
  .txdb_dump <- GenomicFeatures::as.list(.txdb)
  
  .txdb_dump$transcripts %>%
    data.table::as.data.table() %>% 
    dplyr::mutate(tx_name = stringr::str_remove(.$tx_name, ".*\\|")) %>%
    dplyr::mutate(tx_strand = dplyr::case_when(
      tx_strand == "+" ~ 1,
      tx_strand == "-" ~ -1,
    )) %>%
    dplyr::select(
      dplyr::all_of(c("tx_name", "tx_chrom", "tx_strand", "tx_start", "tx_end"))
    ) %>%
    magrittr::set_colnames(c(
      "tx_name"   = .locus_var,
      "tx_chrom"  = .chr_var,
      "tx_strand" = .strand_var,
      "tx_start"  = .start_var,
      "tx_end"    = .end_var
    ))
}