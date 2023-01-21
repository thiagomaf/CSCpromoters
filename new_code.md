CSCpromoters - new code
================

# SETUP

## Load libraries

``` r
library(magrittr)
```

## Define constants

``` r
pb_format <- ":what - [:bar] :percent (:spin)"

folder1 <- paste0(
  "../../../",
  "Data/FASTA sequences/Genomes/Hordeum vulgare - Golden Promise/"
)

folder2 <- paste0(
  "../../../",
  "rcs-gedo2-team_coldstorage/LAB_Share/RNAseq/barley/genome_GP/"
)
```

# LOAD DATA

## H. vulgare cv. Golden Promise annotations

``` r
annotations <- paste0(
  folder1,
  "Annotation_Golden_Promise_v1r1_Apollo_300620_mRNA.fasta"
) %>% 
  CSCfasta::load_fasta() %>%
  CSCfasta::get_fasta_annotation()
```

## H. vulgare cv. Golden Promise TxDB from GFF file

``` r
txdb <- paste0(
  folder1,
  "Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.gff3"
) %>% 
  CSCpromoters::make_txdb(.data_source = "Hv - Golden Promise", .organism = "Hordeum vulgare")
```

    ## Import genomic features from the file as a GRanges object ... OK
    ## Prepare the 'metadata' data frame ... OK
    ## Make the TxDb object ... OK

``` r
# For some reason this cannot be properly loaded from an .RData or .rda file,
# must be run on every new R session

# If we can get the loci `start` and `end` coordinates from txdb we will not 
# need the `annotations` table above.
```

# GET PROMOTERS

``` r
# get_promoter_sequences <- function(
#     .distances,
#     .folder     = paste0(
#       "../../../",
#       "rcs-gedo2-team_coldstorage/LAB_Share/RNAseq/barley/genome_GP/"
#     ),
#     .txdb       = NULL,
#     .pb         = NULL,
#     .pb_format  = pb_format,
#     .locus_var  = "locus_tag",
#     .chr_var    = "chr",
#     .dist_var   = "promoter_size",
#     .downstream = 0
# ) {
#   # This is not working as expected, need to check whether the txdb exists/was
#   # assigned; tryCatch() maybe
#   # if(!quote(.txdb) %>% as.character() %>% exists()) stop(".txdb is not set")
#   if(is.null(.txdb)) stop(".txdb is NULL")
# 
#   .sequences <- .distances %>%
#     dplyr::group_by(across(all_of(.chr_var))) %>%
#     dplyr::group_map(.f = function(.distances_curr_chr, .keys) {
#       # Initiate (or not) the progress bar - can also receive external object
#       if(is.null(.pb)) {
#         .pb <- progress::progress_bar$new(
#           format = .pb_format,
#           total  = .distances_curr_chr %>%
#             dplyr::select(all_of(.locus_var)) %>%
#             unique() %>%
#             nrow()
#         )
#       }
# 
#       # Start the progress bar
#       .pb$tick(0)
# 
#       # Assign the current chromosome name to an object ("character")
#       .curr_chr <- .keys %>% dplyr::pull(.chr_var)
# 
#       # Should trimming of promoter sizes be moved here???
# 
#       # Open FASTA file for each chromosome b4 extracting promoters in each
#       # - This is highly problematic! Currently assumes the input files will
#       #   have a given filename structure. There must be better ways to do that!
#       .fasta_file <- paste0(
#         .folder,
#         paste0("Hordeum_vulgare.refseq[", .curr_chr, "].fasta")
#       ) %>%
#         Rsamtools::FaFile() %>%
#         open()
# 
#       # Get promoters
#       .sequences_curr_chr <- .distances_curr_chr %>%
#         dplyr::group_by(across(all_of(.locus_var))) %>%
#         dplyr::group_map(.f = function(.distance_curr_locus, .each_keys) {
#           # Assign the current locus locus name to an object
#           .curr_locus    <- .each_keys %>% dplyr::pull(.locus_var)
# 
#           # Assign the upstream and downstream length of the promoter sequence
#           # to get (should be numbers).
#           .len_downstream <- .downstream
#           .len_upstream   <- .distance_curr_locus %>%
#             dplyr::pull(.dist_var) %>%
#             unique()
# 
#           # get promoter sequence for each locus
#           .each_promoters <- .txdb %>%
#             GenomicFeatures::transcripts(
#               filter = list("tx_name" = .curr_locus)
#             ) %>%
#             GenomicFeatures::getPromoterSeq(
#               subject    = .fasta_file,
#               upstream   = .len_upstream,
#               downstream = .len_downstream
#             ) %>%
#             magrittr::set_names(.curr_locus)
# 
#           # progress bar oogie boogie
#           .pb$tick(
#             tokens = list(what = paste(.curr_chr, .curr_locus, sep = " - "))
#           )
# 
#           .each_promoters
#         })
# 
#       # Close FASTA file
#       close(.fasta_file)
# 
#       # End progress bar
#       rm(.pb)
# 
#       .sequences_curr_chr %>%
#         do.call(what = c, arg = .)
#     }) %>%
#     do.call(what = c, arg = .)
# 
#   .sequences
# }
```

# Explicit pipeline

``` r
annotations %>%
  # Choose which loci to use
  CSCpromoters::filter_locus(
    .keep = c("chr1Hg0000021", "chr1Hg0000031", "chr1Hg0000041")
  ) %>%
  # Calculate distances to closest upstream locus
  CSCpromoters::get_promoter_distances() %>%
  # Trim found upstream distances and define promoter lengths
  CSCpromoters::trim_distances(.min_size = 100, .max_size = 2000) %>%
  # Get promoter sequences
  CSCpromoters::get_promoter_sequences(
    .txdb   = txdb,
    .folder = paste0(
      "../../../",
      "rcs-gedo2-team_coldstorage/LAB_Share/RNAseq/barley/genome_GP/"
    )
  )
```

    ## DNAStringSet object of length 2:
    ##     width seq                                               names               
    ## [1]  2000 ATTGCGCTGTTTTCACATGAAAA...AGAGGAACAGGTGTTGGAGAGTG chr1Hg0000031
    ## [2]  1585 ACTAACACATGTACTCCTCCATG...GCCCGTAGGATGTGCTAAGCGTA chr1Hg0000041

# Wrap-up function

``` r
annotations %>%
  # Wrap-up function
  CSCpromoters::get_promoters(
    .keep       = c("chr1Hg0000021", "chr1Hg0000031", "chr1Hg0000041"),
    .txdb       = txdb,
    .folder     = paste0(
      "../../../",
      "rcs-gedo2-team_coldstorage/LAB_Share/RNAseq/barley/genome_GP/"
    )
)
```

    ## DNAStringSet object of length 2:
    ##     width seq                                               names               
    ## [1]  2000 ATTGCGCTGTTTTCACATGAAAA...AGAGGAACAGGTGTTGGAGAGTG chr1Hg0000031
    ## [2]  1585 ACTAACACATGTACTCCTCCATG...GCCCGTAGGATGTGCTAAGCGTA chr1Hg0000041

# EXPORT SEQUENCES

``` r
# write_fasta <- function(.sequences, .output_file) {
#   Biostrings::writeXStringSet(.sequences, filepath = .output_file)
# } # This possibly doesn't belongs to the `CSCpromoters` scope
# 
# write_rda <- function(.sequences, .output_file) {
#   save(.sequences, file = .output_file)
# }
# 
# # Not being used currently
# write_sequences <- function(.sequences, .output_file) {
#   switch (
#     tools::file_ext(.output_file),
#     fasta = write_fasta(.sequences, .output_file),
#     fa    = write_fasta(.sequences, .output_file),
#     rda   = write_rda(  .sequences, .output_file)
#   )
# }
# 
# my_promoters %>% write_fasta("data/promoters_Min-Yao.fasta")
# my_promoters %>% write_rda("data/promoters_Min-Yao.rda")
```
