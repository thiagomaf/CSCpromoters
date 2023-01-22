CSCpromoters - Genomic Features
================

# SETUP

## Load libraries

``` r
library(magrittr)
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.2.2

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(CSCpromoters)
```

## Define constants

``` r
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

# TEST

``` r
GenomicFeatures::seqinfo(txdb)
```

    ## Seqinfo object with 8 sequences from an unspecified genome; no seqlengths:
    ##   seqnames seqlengths isCircular genome
    ##   chr1H            NA         NA   <NA>
    ##   chr2H            NA         NA   <NA>
    ##   chr3H            NA         NA   <NA>
    ##   chr4H            NA         NA   <NA>
    ##   chr5H            NA         NA   <NA>
    ##   chr6H            NA         NA   <NA>
    ##   chr7H            NA         NA   <NA>
    ##   chrUn            NA         NA   <NA>

``` r
GenomeInfoDb::seqlevels(txdb)
```

    ## [1] "chr1H" "chr2H" "chr3H" "chr4H" "chr5H" "chr6H" "chr7H" "chrUn"

``` r
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
    dplyr::mutate(tx_name = stringr::str_remove(tx_name, ".*\\|")) %>%
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

annotations_txdb <- txdb %>%
  get_txdb_annotation()
```

``` r
annotations_txdb %>%
  # Choose which loci to use
  filter_locus(
    .keep = c("chr1Hg0000021", "chr1Hg0000031", "chr1Hg0000041")
  ) %>%
  # Calculate distances to closest upstream locus
  get_promoter_distances() %>% 
  # Trim found upstream distances and define promoter lengths
  trim_distances(.min_size = 100, .max_size = 2000) %>%
  # Get promoter sequences
  get_promoter_sequences(.txdb = txdb, .folder = folder2)
```

    ## DNAStringSet object of length 2:
    ##     width seq                                               names               
    ## [1]  2000 ATTGCGCTGTTTTCACATGAAAA...AGAGGAACAGGTGTTGGAGAGTG chr1Hg0000031
    ## [2]  1585 ACTAACACATGTACTCCTCCATG...GCCCGTAGGATGTGCTAAGCGTA chr1Hg0000041

``` r
annotations_txdb %>%
  # Wrap-up function
  get_promoters(
    .keep       = c("chr1Hg0000021", "chr1Hg0000031", "chr1Hg0000041"),
    .txdb       = txdb,
    .folder     = folder2
)
```

    ## DNAStringSet object of length 2:
    ##     width seq                                               names               
    ## [1]  2000 ATTGCGCTGTTTTCACATGAAAA...AGAGGAACAGGTGTTGGAGAGTG chr1Hg0000031
    ## [2]  1585 ACTAACACATGTACTCCTCCATG...GCCCGTAGGATGTGCTAAGCGTA chr1Hg0000041
