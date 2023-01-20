#' make a TxDb object from transcript annotations available as a GFF3 or GTF file.
#'
#' @param .filename 
#' @param .data_source 
#' @param .organism 
#'
#' @return A TxDb object.
#' @export
#'
#' @examples
make_txdb <- function(.filename = NULL, .data_source = NULL, .organism = NULL) {
  if(is.null(.filename))    stop("CSCpromoters: .filename is NULL")
  if(is.null(.data_source)) stop("CSCpromoters: .data_source is NULL")
  if(is.null(.organism))    stop("CSCpromoters: .organism is NULL")
  
  GenomicFeatures::makeTxDbFromGFF(
    file       = .filename,
    dataSource = .data_source,
    organism   = .organism
  )
}