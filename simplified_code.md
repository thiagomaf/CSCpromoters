CSCpromoters - simplified code
================

# SETUP

## Load libraries

``` r
# library(magrittr)
library(CSCpromoter)
library(CSCprimers)

# library(progress)
# library(GenomicFeatures)
# library(Rsamtools)
# library(tidyr)
# library(dplyr)
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
  CSCprimers::load_fasta() %>%
  CSCprimers::get_fasta_annotation() %>% # Add progress bar to this
  dplyr::mutate(begin  = as.numeric(begin)) %>% 
  dplyr::mutate(end    = as.numeric(end)) %>% 
  dplyr::mutate(len    = as.numeric(len)) %>% 
  dplyr::mutate(strand = as.numeric(strand)) %>% 
  data.table::as.data.table()

# consider moving these `mutate()` to inside the get_fasta_annotation() OR 
# create a test function to be used on the fly - this is assumed to be true for
# downstream analyses.

# Ideally, selection of loci should be done here and be propagated downstream.
```

## H. vulgare cv. Golden Promise TxDB from GFF file

``` r
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

## Calculate distances to closest upstream locus

``` r
get_promoter_distances <- function(
    .annotations,
    .pb        = NULL,
    .pb_format = pb_format,
    .locus_var = "locus_tag",
    .chr_var   = "chr"
) {
  # Initiate (or not) the progress bar - can also receive external object
  if(is.null(.pb)) {
    .pb <- progress::progress_bar$new(
      format = .pb_format,
      total  = .annotations %>% 
        dplyr::select(all_of(c(.locus_var, .chr_var))) %>% 
        unique() %>% 
        nrow()
    )
  }
  
  # Start the progress bar
  .pb$tick(0)
  
  # THIS MUST BE PARALLELIZED
  .distances <- .annotations %>%
    dplyr::group_by(across(all_of(c(.locus_var, .chr_var)))) %>%
    dplyr::group_modify(.f = function(.each_annotation, .keys) {
      # Assign the current locus and chromosome names to objects ("character")
      .curr_locus <- .keys %>% dplyr::pull(.locus_var)
      .curr_chr   <- .keys %>% dplyr::pull(.chr_var)
      
      # Subset relevant loci - a.k.a. "loci of interest"
      .loci_OI <- .annotations %>%
        tidyr::pivot_longer(cols = c("begin", "end")) %>%
        # get all loci in the current chromosome
        #subset(chr == .keys$chr) %>%
        subset(get(.chr_var) == .curr_chr) %>%
        # get all loci but the current locus
        #subset(locus_tag != .keys$locus_tag) %>%
        subset(get(.locus_var) != .curr_locus) %>%
        # get all loci upstream the current locus
        subset(
          dplyr::case_when(
            .each_annotation$strand ==  1 ~ value < .each_annotation$begin,
            .each_annotation$strand == -1 ~ value > .each_annotation$end,
            TRUE               ~ NA
          )
        )
      
      # Initiate the output object
      .each_output <- data.frame(
        closest_locus = NA_character_, 
        dist          = NA
      )
      
      # If there are "loci of interest", updates the output object with the 
      # upstream distance to the closest gene (still, in a given chromosome). 
      # The if() below handles e.g. the first locus in each chromosome which 
      # doesn't have other loci upstream.
      if(nrow(.loci_OI) != 0) {
        .each_output <- .loci_OI %>%
          # get the reference coordinate to the upstream locus closest to the 
          # current locus
          subset(dplyr::case_when(
            .each_annotation$strand ==  1 ~ value == max(value, na.rm = TRUE),
            .each_annotation$strand == -1 ~ value == min(value, na.rm = TRUE)
          )) %>%
          # Calculate the distance to the closest upstream locus we got above
          dplyr::mutate(dist = dplyr::case_when(
            .each_annotation$strand ==  1 ~ .each_annotation$begin - value,
            .each_annotation$strand == -1 ~ value - .each_annotation$end
          )) %>%
          dplyr::select(locus_tag, dist) %>% # locus_tag here must be on the fly
          dplyr::rename(closest_locus = locus_tag) # here too, only locus_tag
      }
      
      # progress bar oogie boogie
      .pb$tick(
        tokens = list(what = paste(.curr_chr, .curr_locus, sep = " - "))
      )
      
      .each_output
    })
  
  # End progress bar
  rm(.pb)
  
  .distances
}

minyao_promoters <- annotations %>%
  #filter_locus(.keep = c(1,2)) %>%
  #filter_locus(.keep = 1:10) %>%
  filter_locus(.keep = c("chr1Hg0000021", "chr1Hg0000031", "chr1Hg0000041")) %>%
  get_promoter_distances()

minyao_promoters
```

    ## # A tibble: 3 × 4
    ## # Groups:   locus_tag, chr [3]
    ##   locus_tag     chr   closest_locus  dist
    ##   <chr>         <chr> <chr>         <dbl>
    ## 1 chr1Hg0000021 chr1H <NA>             NA
    ## 2 chr1Hg0000031 chr1H chr1Hg0000021 15680
    ## 3 chr1Hg0000041 chr1H chr1Hg0000031  1585

## Trim found upstream distances and define promoter lengths

``` r
trim_distances <- function(.distances, .min_size = 100, .max_size = 2000) {
  .distances <- .distances  %>%
    dplyr::ungroup() %>% # breaks the magig if removed!
    subset(dist >= .min_size) %>%
    dplyr::mutate(
      promoter_size = dplyr::case_when(
        dist > .max_size ~ .max_size,
        TRUE ~ dist
      )
    )
  
  .distances # absolutelly useless return but makes it clearer what is happening
}

minyao_promoters2 <- minyao_promoters %>% 
  trim_distances(.min_size = 100, .max_size = 2000)

minyao_promoters2
```

    ## # A tibble: 2 × 5
    ##   locus_tag     chr   closest_locus  dist promoter_size
    ##   <chr>         <chr> <chr>         <dbl>         <dbl>
    ## 1 chr1Hg0000031 chr1H chr1Hg0000021 15680          2000
    ## 2 chr1Hg0000041 chr1H chr1Hg0000031  1585          1585

## Get promoter sequences

``` r
get_promoter_sequences <- function(
    .distances,
    .folder     = paste0(
      "../../../",
      "rcs-gedo2-team_coldstorage/LAB_Share/RNAseq/barley/genome_GP/"
    ),
    .txdb       = NULL,
    .pb         = NULL,
    .pb_format  = pb_format,
    .locus_var  = "locus_tag",
    .chr_var    = "chr",
    .dist_var   = "promoter_size",
    .downstream = 0
) {
  # This is not working as expected, need to check whether the txdb exists/was 
  # assigned; tryCatch() maybe
  # if(!quote(.txdb) %>% as.character() %>% exists()) stop(".txdb is not set") 
  if(is.null(.txdb)) stop(".txdb is NULL")
  
  .sequences <- .distances %>%
    dplyr::group_by(across(all_of(.chr_var))) %>%
    dplyr::group_map(.f = function(.distances_curr_chr, .keys) {
      # Initiate (or not) the progress bar - can also receive external object
      if(is.null(.pb)) {
        .pb <- progress::progress_bar$new(
          format = .pb_format,
          total  = .distances_curr_chr %>%
            dplyr::select(all_of(.locus_var)) %>% 
            unique() %>% 
            nrow()
        )
      }
      
      # Start the progress bar
      .pb$tick(0)
      
      # Assign the current chromosome name to an object ("character")
      .curr_chr <- .keys %>% dplyr::pull(.chr_var)
      
      # Should trimming of promoter sizes be moved here???
      
      # Open FASTA file for each chromosome b4 extracting promoters in each
      # - This is highly problematic! Currently assumes the input files will 
      #   have a given filename structure. There must be better ways to do that!
      .fasta_file <- paste0(
        .folder,
        paste0("Hordeum_vulgare.refseq[", .curr_chr, "].fasta")
      ) %>%
        Rsamtools::FaFile() %>% 
        open()
      
      # Get promoters
      .sequences_curr_chr <- .distances_curr_chr %>%
        dplyr::group_by(across(all_of(.locus_var))) %>%
        dplyr::group_map(.f = function(.distance_curr_locus, .each_keys) {
          # Assign the current locus locus name to an object
          .curr_locus    <- .each_keys %>% dplyr::pull(.locus_var)
          
          # Assign the upstream and downstream length of the promoter sequence 
          # to get (should be numbers).
          .len_downstream <- .downstream
          .len_upstream   <- .distance_curr_locus %>%
            dplyr::pull(.dist_var) %>%
            unique()
          
          # get promoter sequence for each locus
          .each_promoters <- .txdb %>% 
            GenomicFeatures::transcripts(
              filter = list("tx_name" = .curr_locus)
            ) %>%
            GenomicFeatures::getPromoterSeq(
              subject    = .fasta_file,
              upstream   = .len_upstream,
              downstream = .len_downstream
            ) %>%
            magrittr::set_names(.curr_locus)
          
          # progress bar oogie boogie
          .pb$tick(
            tokens = list(what = paste(.curr_chr, .curr_locus, sep = " - "))
          )
          
          .each_promoters
        })
      
      # Close FASTA file
      close(.fasta_file)
      
      # End progress bar
      rm(.pb)
      
      .sequences_curr_chr %>% 
        do.call(what = c, arg = .)
    }) %>% 
    do.call(what = c, arg = .)
  
  .sequences
}

my_promoters <- minyao_promoters2 %>% 
  get_promoter_sequences(.txdb = txdb)

my_promoters
```

    ## DNAStringSet object of length 2:
    ##     width seq                                               names               
    ## [1]  2000 ATTGCGCTGTTTTCACATGAAAA...AGAGGAACAGGTGTTGGAGAGTG chr1Hg0000031
    ## [2]  1585 ACTAACACATGTACTCCTCCATG...GCCCGTAGGATGTGCTAAGCGTA chr1Hg0000041

# EXPORT SEQUENCES

``` r
write_fasta <- function(.sequences, .output_file) {
  Biostrings::writeXStringSet(.sequences, filepath = .output_file)
} # This possibly doesn't belongs to the `CSCpromoters` scope

write_rda <- function(.sequences, .output_file) {
  save(.sequences, file = .output_file)
}

# Not being used currently
write_sequences <- function(.sequences, .output_file) {
  switch (
    tools::file_ext(.output_file),
    fasta = write_fasta(.sequences, .output_file),
    fa    = write_fasta(.sequences, .output_file),
    rda   = write_rda(  .sequences, .output_file)
  )
}

my_promoters %>% write_fasta("data/promoters_Min-Yao.fasta")
my_promoters %>% write_rda("data/promoters_Min-Yao.rda")
```

# TO DELETE LATER (TAM)

``` r
my_promoters5 <- annotations %>%
  get_promoters(
    .keep       = c("chr1Hg0000021", "chr1Hg0000031", "chr1Hg0000041"),
    .txdb       = txdb,
    .folder     = paste0(
      "../../../",
      "rcs-gedo2-team_coldstorage/LAB_Share/RNAseq/barley/genome_GP/"
    )
)
```
