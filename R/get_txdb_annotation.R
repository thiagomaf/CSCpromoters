#' @title Extract annotation table from TxDB object
#'
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' @inheritParams get_promoters
#' @inheritDotParams get_promoters
#' 
#' @export
#' 
get_txdb_annotation <- function(
    .txdb,
    .locus_var  = "locus_tag",
    .chr_var    = "chr",
    .strand_var = "strand",
    .start_var  = "begin",
    .end_var    = "end",
    ...
) {
  .tx_name   <- "tx_name"
  .tx_chrom  <- "tx_chrom"
  .tx_strand <- "tx_strand"
  .tx_start  <- "tx_start"
  .tx_end    <- "tx_end"
  
  .txdb_dump <- GenomicFeatures::as.list(.txdb)
  
  .txdb_dump$transcripts %>%
    data.table::as.data.table() %>% 
    dplyr::mutate(
      !!.tx_name   := stringr::str_remove(get(.tx_name), ".*\\|")
    ) %>%
    dplyr::mutate(!!.tx_strand := dplyr::case_when(
      get(.tx_strand) == "+" ~  1,
      get(.tx_strand) == "-" ~ -1,
    )) %>%
    dplyr::select(
      dplyr::all_of(c(.tx_name, .tx_chrom, .tx_strand, .tx_start, .tx_end))
    ) %>%
    dplyr::rename(!!.locus_var  := dplyr::all_of(.tx_name)) %>% 
    dplyr::rename(!!.chr_var    := dplyr::all_of(.tx_chrom)) %>% 
    dplyr::rename(!!.strand_var := dplyr::all_of(.tx_strand)) %>% 
    dplyr::rename(!!.start_var  := dplyr::all_of(.tx_start)) %>% 
    dplyr::rename(!!.end_var    := dplyr::all_of(.tx_end))
}