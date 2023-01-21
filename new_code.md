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
  make_txdb(.data_source = "Hv - Golden Promise", .organism = "Hordeum vulgare")
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

# Explicit pipeline

``` r
annotations %>%
  # Choose which loci to use
  filter_locus(
    .keep = c("chr1Hg0000021", "chr1Hg0000031", "chr1Hg0000041")
  ) %>%
  # Calculate distances to closest upstream locus
  get_promoter_distances() %>%
  # Trim found upstream distances and define promoter lengths
  trim_distances(.min_size = 100, .max_size = 2000) %>%
  # Get promoter sequences
  get_promoter_sequences(
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
  get_promoters(
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
