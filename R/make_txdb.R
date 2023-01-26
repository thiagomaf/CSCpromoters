#' @title Make a TxDb object from transcript annotations available as a GFF3 or 
#'        GTF file.
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#'
#' @param .filename    {TODO}
#' @param .data_source {TODO}
#' @param .organism    {TODO}
#' @param ...          Future-proof stuff.
#'
#' @return A TxDb object.
#' @export
#'
make_txdb <- function(
    .filename    = NULL,
    .data_source = NULL,
    .organism    = NULL,
    ...
  ) {
  if(is.null(.filename))    stop("CSCpromoters: .filename is NULL")
  if(is.null(.data_source)) stop("CSCpromoters: .data_source is NULL")
  if(is.null(.organism))    stop("CSCpromoters: .organism is NULL")
  
  GenomicFeatures::makeTxDbFromGFF(
    file       = .filename,
    dataSource = .data_source,
    organism   = .organism
  )
}