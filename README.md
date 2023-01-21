
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CSCpromoters

<!-- badges: start -->
<!-- badges: end -->

## About

The goal of CSCpromoters is to extract the sequences of promoter regions
from the genome fasta files.

## Description

`CSCpromoters` allows users to extract the sequences of promoter regions
from the genome fasta files. `CSCpromoters` can calculate the length of
promoters based on the distance to the closest upstream gene.
`CSCpromoters` also provides the option for users to enter the
user-defined minimum and maximum promoter length.

## Installation

`CSCpromoters` is part of a larger set of CSC libraries and depends on
`CSCfasta`.

You can install the development version of `CSCpromoters` and `CSCfasta`
from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("thiagomaf/CSCfasta")
devtools::install_github("thiagomaf/CSCpromoters")

# These commands will not work while the repositories are private.
# There are work-around, we can talk about it later.
```

## Usage Example

This is a basic example which shows you how to solve a common problem:

### SETUP

#### Load libraries

``` r
library(magrittr)

#library(CSCfasta)
library(CSCpromoters)
```

#### Define constants

``` r
folder1 <- paste0(
  "./"
)

folder2 <- paste0(
  "./"
)
```

### LOAD DATA

#### H. vulgare cv. Golden Promise annotations

``` r
annotations <- paste0(
  folder1,
  "Annotation_Golden_Promise_v1r1_Apollo_300620_mRNA.fasta"
) %>%
  CSCfasta::load_fasta() %>%
  CSCfasta::get_fasta_annotation()
```

#### H. vulgare cv. Golden Promise TxDB from GFF file

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

### GET PROMOTERS

#### Explicit pipeline

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
  get_promoter_sequences(.txdb = txdb, .folder = folder2)
```

    ## DNAStringSet object of length 2:
    ##     width seq                                               names               
    ## [1]  2000 ATTGCGCTGTTTTCACATGAAAA...AGAGGAACAGGTGTTGGAGAGTG chr1Hg0000031
    ## [2]  1585 ACTAACACATGTACTCCTCCATG...GCCCGTAGGATGTGCTAAGCGTA chr1Hg0000041

#### Wrap-up function

``` r
annotations %>%
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

## Dependencies

### Libraries

- magrittr
- tidyr
- plyr
- dplyr
- purrr
- progress
- GenomicFeatures
- CSCprimers

### Files

- \[folder1\]/Annotation_Golden_Promise_v1r1_Apollo_300620_mRNA.fasta
- \[folder1\]/Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.gff3
- \[folder2\]/Hordeum_vulgare.refseq\[chr1H\].fasta
- \[folder2\]/Hordeum_vulgare.refseq\[chr2H\].fasta
- \[folder2\]/Hordeum_vulgare.refseq\[chr3H\].fasta
- \[folder2\]/Hordeum_vulgare.refseq\[chr4H\].fasta
- \[folder2\]/Hordeum_vulgare.refseq\[chr5H\].fasta
- \[folder2\]/Hordeum_vulgare.refseq\[chr6H\].fasta
- \[folder2\]/Hordeum_vulgare.refseq\[chr7H\].fasta
- \[folder2\]/Hordeum_vulgare.refseq\[chrUn\].fasta
