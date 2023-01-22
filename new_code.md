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
# paste0(
#   gff_folder,
#   "Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.gff3"
# ) %>% 
#   make_txdb(
#     .data_source = "Hv - Golden Promise",
#     .organism    = "Hordeum vulgare"
#   ) %>% 
#   saveDb(
#     file = paste0(
#       gff_folder,
#       "Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.sqlite"
#     )
#   )

# For some reason this cannot be properly loaded from an .RData or .rda file,
# must be run on every new R session

txdb <- paste0(
  gff_folder,
  "Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.sqlite"
) %>% AnnotationDbi::loadDb()
```

    ## Carregando pacotes exigidos: GenomicFeatures

    ## Carregando pacotes exigidos: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Carregando pacotes exigidos: S4Vectors

    ## Carregando pacotes exigidos: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Carregando pacotes exigidos: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Carregando pacotes exigidos: GenomeInfoDb

    ## Carregando pacotes exigidos: GenomicRanges

    ## 
    ## Attaching package: 'GenomicRanges'

    ## The following object is masked from 'package:magrittr':
    ## 
    ##     subtract

    ## Carregando pacotes exigidos: AnnotationDbi

    ## Carregando pacotes exigidos: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

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
    .keep = c("chr1Hg0000021", "chr1Hg0000031", "chr1Hg0000041")
  ) %>%
  # Calculate distances to closest upstream locus
  get_promoter_distances() %>%
  # Trim found upstream distances and define promoter lengths
  trim_distances(.min_size = 100, .max_size = 2000) %>%
  # Get promoter sequences
  get_promoter_sequences2(
    .txdb = txdb, .FASTA_list = fasta_list, .parallel = F
  )
```

    ## DNAStringSet object of length 2:
    ##     width seq                                               names               
    ## [1]  2000 ATTGCGCTGTTTTCACATGAAAA...AGAGGAACAGGTGTTGGAGAGTG chr1Hg0000031
    ## [2]  1585 ACTAACACATGTACTCCTCCATG...GCCCGTAGGATGTGCTAAGCGTA chr1Hg0000041

## Wrap-up function

``` r
# annotations %>%
#   # Wrap-up function
#   get_promoters(
#     .keep       = c("chr1Hg0000021", "chr1Hg0000031", "chr1Hg0000041"),
#     .txdb       = txdb,
#     .folder     = folder2
#   )
```

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
