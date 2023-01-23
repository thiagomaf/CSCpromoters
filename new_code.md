CSCpromoters - new code
================

# SETUP

## Load libraries

``` r
library(magrittr)
library(CSCpromoters)
```

## Define constants

``` r
gff_folder <- paste0(
  "../../../",
  "Data/FASTA sequences/Genomes/Hordeum vulgare - Golden Promise/"
)

fasta_list <- (function(
    .chr_list = c(
      "chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H", "chrUn"
    ),
    .folder   = paste0(
      "../../../",
      "rcs-gedo2-team_coldstorage/LAB_Share/RNAseq/barley/genome_GP/"
    )
) {
  .folder %>% 
    paste0(paste0("Hordeum_vulgare.refseq[", .chr_list, "].fasta")) %>% 
    magrittr::set_names(.chr_list)
})()
```

# LOAD DATA

## H. vulgare cv. Golden Promise TxDB from GFF file

``` r
# # paste0(
# #   gff_folder,
# #   "Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.gff3"
# # ) %>% 
# #   make_txdb(
# #     .data_source = "Hv - Golden Promise",
# #     .organism    = "Hordeum vulgare"
# #   ) %>% 
# #   saveDb(
# #     file = paste0(
# #       gff_folder,
# #       "Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.sqlite"
# #     )
# #   )
# 
# # For some reason this cannot be properly loaded from an .RData or .rda file,
# # must be run on every new R session
# 
# 
# txdb <- paste0(
#   gff_folder,
#   "Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.sqlite"
# ) %>% AnnotationDbi::loadDb()

txdb <- paste0(
  gff_folder,
  "Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.gff3"
) %>%
  make_txdb(
    .data_source = "Hv - Golden Promise",
    .organism    = "Hordeum vulgare"
  )
```

    ## Import genomic features from the file as a GRanges object ... OK
    ## Prepare the 'metadata' data frame ... OK
    ## Make the TxDb object ... OK

## H. vulgare cv. Golden Promise annotations

``` r
# annotations <- paste0(
#   gff_folder,
#   "Annotation_Golden_Promise_v1r1_Apollo_300620_mRNA.fasta"
# ) %>% 
#   CSCfasta::load_fasta() %>%
#   CSCfasta::get_fasta_annotation()

annotations <- txdb %>%
  get_txdb_annotation()
```

# GET PROMOTERS

## Explicit pipeline

``` r
annotations %>%
  # Choose which loci to use
  filter_locus(
    # .keep = c("chr1Hg0000021", "chr1Hg0000031", "chr1Hg0000041")
    # .keep = "chr1Hg0040651"
    .keep = c(
      "chr6Hg0657981",
      "chr1Hg0040651",
      "chr3Hg0328761",
      "chr4Hg0417171",
      "chr3Hg0345511",
      "chr4Hg0383761",
      "chr5Hg0547031",
      "chr4Hg0380431",
      "chr4Hg0421931",
      "chr7Hg0679581"
    )
    # .keep = 1:3
  ) %>% 
  #attr(which = "keep") %>%
  # Calculate distances to closest upstream locus
  get_promoter_distances(.debug = F) %>%
  # Trim found upstream distances and define promoter lengths
  trim_distances(.min_size = 100, .max_size = 2000) %>%
  # Get promoter sequences
  get_promoter_sequences(
    .txdb = txdb, .FASTA_list = fasta_list, .parallel = F
  )
```

    ## DNAStringSet object of length 8:
    ##     width seq                                               names               
    ## [1]  2000 GTAACAATAGTAACAAGGTGCAC...CCAGCCGTACAGACGATATTTCA chr1Hg0040651
    ## [2]  2000 AACTAATCTGTGGTTGGATGACT...ATTGTTCCGGCACGGGGCTGGGG chr3Hg0328761
    ## [3]  2000 ACATAGAAAGTATGCACATGACA...CCTCCCTCCCTCCCTCCCCCCAA chr3Hg0345511
    ## [4]  2000 TGCATTTGACACATCAGATTTGG...ATGGCGCGGCTACGTCTGCTACG chr4Hg0380431
    ## [5]  2000 TTTAATGCAATGTATGAATATGA...AAGACTGCATAAAGTTTGGATCA chr4Hg0383761
    ## [6]  2000 GAGAGCCTTTGGTCAGCAAAGAT...GCTCTCGATCGCATGGAAGAAAA chr4Hg0417171
    ## [7]  1698 GACATATATGTGTCTCATAATGA...TACCCGCGCTACTATTACAGGAA chr5Hg0547031
    ## [8]  2000 GACGGACGGATGGATCATGGATG...TGTGAGCCTGAGATGCAGGGGAA chr6Hg0657981

## Wrap-up function

``` r
annotations %>%
  # Wrap-up function
  get_promoters(
    # .keep       = c("chr1Hg0000021", "chr1Hg0000031", "chr1Hg0000041"),
    .keep       = c(
      "chr6Hg0657981",
      "chr1Hg0040651",
      "chr3Hg0328761",
      "chr4Hg0417171",
      "chr3Hg0345511",
      "chr4Hg0383761",
      "chr5Hg0547031",
      "chr4Hg0380431",
      "chr4Hg0421931",
      "chr7Hg0679581"
    ),
    .txdb       = txdb,
    .FASTA_list = fasta_list
  )
```

    ## DNAStringSet object of length 8:
    ##     width seq                                               names               
    ## [1]  2000 GTAACAATAGTAACAAGGTGCAC...CCAGCCGTACAGACGATATTTCA chr1Hg0040651
    ## [2]  2000 AACTAATCTGTGGTTGGATGACT...ATTGTTCCGGCACGGGGCTGGGG chr3Hg0328761
    ## [3]  2000 ACATAGAAAGTATGCACATGACA...CCTCCCTCCCTCCCTCCCCCCAA chr3Hg0345511
    ## [4]  2000 TGCATTTGACACATCAGATTTGG...ATGGCGCGGCTACGTCTGCTACG chr4Hg0380431
    ## [5]  2000 TTTAATGCAATGTATGAATATGA...AAGACTGCATAAAGTTTGGATCA chr4Hg0383761
    ## [6]  2000 GAGAGCCTTTGGTCAGCAAAGAT...GCTCTCGATCGCATGGAAGAAAA chr4Hg0417171
    ## [7]  1698 GACATATATGTGTCTCATAATGA...TACCCGCGCTACTATTACAGGAA chr5Hg0547031
    ## [8]  2000 GACGGACGGATGGATCATGGATG...TGTGAGCCTGAGATGCAGGGGAA chr6Hg0657981

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
