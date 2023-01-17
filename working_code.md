CSCpromoters - working code
================

# SETUP

``` r
library(magrittr)
# library(progress)
library(CSCprimers)
# library(GenomicFeatures)
# library(Rsamtools)
# library(tidyr)
# library(dplyr)

# There is an issue with the loading of libraries. Some conflict with dplyr (or tidyverse as a whole) and either the GenomicFeatures or Rsamtools
# Omit messages!???
```

# LOAD DATA

## H. vulgare cv. Golden Promise annotations

``` r
annotations <- load_fasta(
  .filename = paste0(
    "../../../Data/FASTA sequences/Genomes/Hordeum vulgare - Golden Promise/",
    "Annotation_Golden_Promise_v1r1_Apollo_300620_mRNA.fasta"
  )
) %>%
  get_fasta_annotation() %>% # Add progress bar
  dplyr::mutate(begin  = as.numeric(begin)) %>% 
  dplyr::mutate(end    = as.numeric(end)) %>% 
  dplyr::mutate(len    = as.numeric(len)) %>% 
  dplyr::mutate(strand = as.numeric(strand))
# consider moving these `mutate()` to inside the get_fasta_annotation() OR 
# create a test function to be used on the fly - this is assumed to be true for
# downstream analyses.
```

## H. vulgare cv. Golden Promise TxDB from GFF file

``` r
# gffFile <- paste0(
#   # "E:/OneDrive - University of Cambridge/CSC/",
#   "C:/Users/ta507/OneDrive - University of Cambridge/CSC/",
#   "Data/FASTA sequences/Genomes/Hordeum vulgare - Golden Promise/",
#   "Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.gff3"
# )
# 
# txdb <- makeTxDbFromGFF(
#   file       = gffFile,
#   dataSource = "Hv - Golden Promise",
#   organism   = "Hordeum vulgare"
# )

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


# For some reason this cannot be properly loaded from an .RData or .rda file,
# must be run on every new R session
txdb <- paste0(
  "../../../",
  "Data/FASTA sequences/Genomes/Hordeum vulgare - Golden Promise/",
  "Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.gff3"
) %>% 
  make_txdb(.data_source = "Hv - Golden Promise", .organism = "Hordeum vulgare")
```

    ## Import genomic features from the file as a GRanges object ... OK
    ## Prepare the 'metadata' data frame ... OK
    ## Make the TxDb object ... OK

# GET PROMOTERS

``` r
# # pb <- progress_bar$new(
# #   total = annotations %>% 
# #     dplyr::select("locus_tag", "chr") %>% 
# #     unique() %>% 
# #     nrow()
# # )
# 
# # This should become a function
# # THIS MUST BE PARALLELIZED
# minyao_promoters <- annotations %>%
#   group_by(locus_tag, chr) %>%
#   (function(.data) {
#     pb$tick(0)
#     .data
#   }) %>% 
#   group_modify(.f = function(.each_locus, .keys) {
#     .results <- annotations %>%
#       pivot_longer(cols = c("begin", "end")) %>% # PIVOT TO LONGER - BEGIN AND END VALUES
#       subset(chr == .keys$chr) %>%
#       subset(locus_tag != .keys$locus_tag) %>%
#       subset(case_when(
#         # .each_locus$strand ==  1 ~ value < max(.each_locus$begin, .each_locus$end, na.rm = TRUE),
#         # .each_locus$strand == -1 ~ value > min(.each_locus$begin, .each_locus$end, na.rm = TRUE),
#         .each_locus$strand ==  1 ~ value < .each_locus$begin,
#         .each_locus$strand == -1 ~ value > .each_locus$end,
#         TRUE               ~ NA
#       ))
#     
#     .output <- if(nrow(.results) != 0) {
#       .results <- .results %>%
#         subset(case_when(
#           .each_locus$strand ==  1 ~ value == max(value, na.rm = TRUE),
#           .each_locus$strand == -1 ~ value == min(value, na.rm = TRUE)
#         )) %>%
#         mutate(dist = case_when(
#           .each_locus$strand ==  1 ~ .each_locus$begin - value,
#           .each_locus$strand == -1 ~ value - .each_locus$end
#         )) %>% 
#         select(locus_tag, dist)
#       
#       data.frame(
#         closest_locus = .results$locus_tag, 
#         dist          = .results$dist
#       )
#     } else {
#       data.frame(
#         closest_locus = NA_character_, 
#         dist          = NA
#       )
#     }
#     
#     pb$tick()
#     
#     .output
#   })
# 
# rm(pb)

get_promoter_distances <- function(.annotations, .pb = NULL, .locus_var = "locus_tag", .chr_var = "chr") {
  # Initiate (or not) the progress bar - can also receive external object
  if(is.null(.pb)) {
    .pb <- progress::progress_bar$new(
      total = .annotations %>% 
        dplyr::select(all_of(c(.locus_var, .chr_var))) %>% 
        unique() %>% 
        nrow()
    )
  }
  
  # THIS MUST BE PARALLELIZED
  .distances <- .annotations %>%
    #group_by(locus_tag, chr) %>% # This should probably be set on the fly
    dplyr::group_by(across(all_of(c(.locus_var, .chr_var)))) %>%
    (function(.data) {
      .pb$tick(0)
      .data
    }) %>% 
    dplyr::group_modify(.f = function(.each_locus, .keys) {
      .results <- .annotations %>%
        tidyr::pivot_longer(cols = c("begin", "end")) %>% # PIVOT TO LONGER - BEGIN AND END VALUES
        subset(chr == .keys$chr) %>%             # This should be done on the fly, merging probably
        subset(locus_tag != .keys$locus_tag) %>% # This should be done on the fly, merging probably
        subset(dplyr::case_when(                 # This should be done on the fly, merging probably
          # .each_locus$strand ==  1 ~ value < max(.each_locus$begin, .each_locus$end, na.rm = TRUE),
          # .each_locus$strand == -1 ~ value > min(.each_locus$begin, .each_locus$end, na.rm = TRUE),
          .each_locus$strand ==  1 ~ value < .each_locus$begin,
          .each_locus$strand == -1 ~ value > .each_locus$end,
          TRUE               ~ NA
        ))
      
      .each_output <- if(nrow(.results) != 0) {
        .results <- .results %>%
          subset(dplyr::case_when(
            .each_locus$strand ==  1 ~ value == max(value, na.rm = TRUE),
            .each_locus$strand == -1 ~ value == min(value, na.rm = TRUE)
          )) %>%
          dplyr::mutate(dist = dplyr::case_when(
            .each_locus$strand ==  1 ~ .each_locus$begin - value,
            .each_locus$strand == -1 ~ value - .each_locus$end
          )) %>% 
          dplyr::select(locus_tag, dist)
        
        data.frame(
          closest_locus = .results$locus_tag, 
          dist          = .results$dist
        )
      } else {
        data.frame(
          closest_locus = NA_character_, 
          dist          = NA
        )
      }
      
      .pb$tick()
      
      .each_output
    })
  
  rm(.pb)
  
  .distances
}

minyao_promoters <- annotations %>% 
  dplyr::slice(1:10) %>%
  get_promoter_distances()
  
minyao_promoters
```

    ## # A tibble: 10 × 4
    ## # Groups:   locus_tag, chr [10]
    ##    locus_tag     chr   closest_locus   dist
    ##    <chr>         <chr> <chr>          <dbl>
    ##  1 chr1Hg0000001 chr1H <NA>              NA
    ##  2 chr1Hg0000011 chr1H chr1Hg0000021     68
    ##  3 chr1Hg0000021 chr1H chr1Hg0000011      1
    ##  4 chr1Hg0000031 chr1H chr1Hg0000021  15680
    ##  5 chr1Hg0000041 chr1H chr1Hg0000031   1585
    ##  6 chr1Hg0000051 chr1H chr1Hg0000041   1749
    ##  7 chr1Hg0000061 chr1H chr1Hg0000071   2010
    ##  8 chr1Hg0000071 chr1H chr1Hg0000101 222620
    ##  9 chr1Hg0000101 chr1H chr1Hg0000071 222620
    ## 10 chr1Hg0000111 chr1H chr1Hg0000101   3030

``` r
# This should become a function
# minyao_promoters2 <- minyao_promoters %>% 
#   ungroup() %>%
#   (function(.data, .min_size = 100, .max_size = 2000) {
#     .data %>%
#       subset(dist >= .min_size) %>% 
#       mutate(promoter_size = case_when(
#         dist > .max_size ~ .max_size,
#         TRUE ~ dist
#       ))
#   })

trim_distances <- function(.distances, .min_size = 100, .max_size = 2000) {
  .distances <- .distances  %>% 
    dplyr::ungroup() %>%
    (function(.data) {
      .data %>%
        subset(dist >= .min_size) %>% 
        dplyr::mutate(promoter_size = dplyr::case_when(
          dist > .max_size ~ .max_size,
          TRUE ~ dist
        ))
    })
  
  .distances # absolutelly useless return, but make it clearer what is happening
}

minyao_promoters2 <- minyao_promoters %>% 
  trim_distances(.min_size = 100, .max_size = 2000)

minyao_promoters2
```

    ## # A tibble: 7 × 5
    ##   locus_tag     chr   closest_locus   dist promoter_size
    ##   <chr>         <chr> <chr>          <dbl>         <dbl>
    ## 1 chr1Hg0000031 chr1H chr1Hg0000021  15680          2000
    ## 2 chr1Hg0000041 chr1H chr1Hg0000031   1585          1585
    ## 3 chr1Hg0000051 chr1H chr1Hg0000041   1749          1749
    ## 4 chr1Hg0000061 chr1H chr1Hg0000071   2010          2000
    ## 5 chr1Hg0000071 chr1H chr1Hg0000101 222620          2000
    ## 6 chr1Hg0000101 chr1H chr1Hg0000071 222620          2000
    ## 7 chr1Hg0000111 chr1H chr1Hg0000101   3030          2000

# GET PROMOTERS

``` r
# pb2 <- progress_bar$new(
#   total = minyao_promoters2 %>%
#     dplyr::select(locus_tag, chr) %>% 
#     unique() %>% 
#     nrow()
# )
# 
# my_promoters <- minyao_promoters2 %>%
#   dplyr::slice(1:100) %>% # Taking too long to do the whole table
#   # (function(.data) {
#   #   pb2$tick(0)
#   #   .data
#   # }) %>% 
#   #dplyr::group_by(chr) %>% # 'chr' should probably be set on the fly
#   dplyr::group_by(across(all_of("chr"))) %>%
#   dplyr::group_map(.f = function(
#     .data, .keys, .txdb = txdb,
#     #.min_size = 100, .max_size = 2000, # not being used yet
#     .debug = TRUE
#   ) {
#     # This is not working as expected, need to check whether the txdb exists/was assigned; tryCatch() maybe
#     # if(!quote(.txdb) %>% as.character() %>% exists()) stop(".txdb is not set") 
#     if(is.null(.txdb)) stop(".txdb is NULL")
#     
#     .each_chr <- .keys %>% dplyr::pull("chr") %>% unique()
#     # if(.debug) print(.each_chr) # became native to the function operation, below
#     
#     # Talk to me!
#     print(.each_chr)
#     
#     # Trim promoter sizes
#     #... to-do
#     
#     # Open FASTA file for each chromossome prior to extracting promoters in each
#     .fasta_file <- paste0(
#       # "E:/OneDrive - University of Cambridge/CSC/",
#       "C:/Users/ta507/OneDrive - University of Cambridge/CSC/",
#       "rcs-gedo2-team_coldstorage/LAB_Share/RNAseq/barley/genome_GP/",
#       paste0("Hordeum_vulgare.refseq[", .each_chr, "].fasta")
#     ) %>%
#       Rsamtools::FaFile() %>% 
#       open()
#     
#     # Set progress bar
#     .pb <- progress::progress_bar$new(
#       total = .data %>%
#         dplyr::select(locus_tag) %>% 
#         unique() %>% 
#         nrow()
#     )
#     
#     # Initiate progress bar
#     # Get promoters for each loci
#     .promoters <- .data %>%
#       (function(.data) {
#         .pb$tick(0)
#         .data
#       }) %>% 
#       #dplyr::group_by("locus_tag") %>% # 'locus_tag' must be assign on the fly
#       dplyr::group_by(across(all_of("locus_tag"))) %>%
#       dplyr::group_map(.f = function(.each_data, .each_keys) {
#         # 'locus_tag' and 'promoter_size' must be assign on the fly
#         .each_loci       <- .each_keys %>% dplyr::pull(locus_tag) %>% unique()
#         .each_upstream   <- .each_data %>% dplyr::pull(promoter_size) %>% unique()
#         .each_downstream <- 0 # implement properly lazy ass! Either function argument or data column
#         
#         .each_promoters <- .txdb %>% 
#           GenomicFeatures::transcripts(filter = list("tx_name" = .each_loci)) %>%
#           GenomicFeatures::getPromoterSeq(
#             subject    = .fasta_file,
#             upstream   = .each_upstream,
#             downstream = .each_downstream
#           ) %>% 
#           set_names(.each_loci)
#         
#         .pb$tick()
#         
#         .each_promoters
#       })
#     
#     # Close FASTA file
#     close(.fasta_file)
#     
#     # End progress bar
#     # .pb$tick()
#     # rm(.pb)
#     
#     .promoters %>% 
#       do.call(what = c, arg = .)
#   }) %>% 
#   do.call(what = c, arg = .)

get_promoters <- function(.distances, .txdb = NULL, .locus_var = "locus_tag", .chr_var = "chr", .dist_var = "promoter_size") {
  # This is not working as expected, need to check whether the txdb exists/was assigned; tryCatch() maybe
  # if(!quote(.txdb) %>% as.character() %>% exists()) stop(".txdb is not set") 
  if(is.null(.txdb)) stop(".txdb is NULL")
  
  .sequences <- .distances %>%
    # dplyr::slice(1:100) %>% # Taking too long to do the whole table
    dplyr::group_by(across(all_of(.chr_var))) %>%
    dplyr::group_map(.f = function(
    .data, .keys,
    #.min_size = 100, .max_size = 2000, # not being used yet, might never be
    .debug = TRUE
    ) {
      .each_chr <- .keys %>% dplyr::pull(.chr_var) %>% unique()

      # Talk to me!
      print(.each_chr)
      
      # Should trimming of promoter sizes be moved here???
      
      # Open FASTA file for each chromossome b4 extracting promoters in each
      # - This is highly problematic! Currently assumes the input files will 
      #   have a given filename structure. There must be better ways to do that!
      .fasta_file <- paste0(
        "../../../",
        "rcs-gedo2-team_coldstorage/LAB_Share/RNAseq/barley/genome_GP/",
        paste0("Hordeum_vulgare.refseq[", .each_chr, "].fasta")
      ) %>%
        Rsamtools::FaFile() %>% 
        open()
      
      # Set progress bar
      .pb <- progress::progress_bar$new(
        total = .data %>%
          dplyr::select(all_of(.locus_var)) %>% 
          unique() %>% 
          nrow()
      )
      
      # Initiate progress bar
      # Get promoters for each loci
      .promoters <- .data %>%
        (function(.data) {
          .pb$tick(0)
          .data
        }) %>% 
        #dplyr::group_by("locus_tag") %>% # 'locus_tag' must be assign on the fly
        dplyr::group_by(across(all_of(.locus_var))) %>%
        dplyr::group_map(.f = function(.each_data, .each_keys) {
          # 'locus_tag' and 'promoter_size' must be assign on the fly
          # .each_loci       <- .each_keys %>% dplyr::pull(locus_tag) %>% unique()
          # .each_upstream   <- .each_data %>% dplyr::pull(promoter_size) %>% unique()
          .each_loci       <- .each_keys %>% dplyr::pull(.locus_var) %>% unique()
          .each_upstream   <- .each_data %>% dplyr::pull(.dist_var) %>% unique()
          .each_downstream <- 0 # implement properly lazy ass! Either function argument or data column
          
          .each_promoters <- .txdb %>% 
            GenomicFeatures::transcripts(filter = list("tx_name" = .each_loci)) %>%
            GenomicFeatures::getPromoterSeq(
              subject    = .fasta_file,
              upstream   = .each_upstream,
              downstream = .each_downstream
            ) %>% 
            set_names(.each_loci)
          
          .pb$tick()
          
          .each_promoters
        })
      
      # Close FASTA file
      close(.fasta_file)
      
      # End progress bar - NOT WORKING
      # rm(.pb)
      
      .promoters %>% 
        do.call(what = c, arg = .)
    }) %>% 
    do.call(what = c, arg = .)
  
  .sequences
}

my_promoters <- minyao_promoters2 %>% 
  get_promoters(.txdb = txdb)
```

    ## [1] "chr1H"

``` r
my_promoters
```

    ## DNAStringSet object of length 7:
    ##     width seq                                               names               
    ## [1]  2000 ATTGCGCTGTTTTCACATGAAAA...AGAGGAACAGGTGTTGGAGAGTG chr1Hg0000031
    ## [2]  1585 ACTAACACATGTACTCCTCCATG...GCCCGTAGGATGTGCTAAGCGTA chr1Hg0000041
    ## [3]  1749 AGACTGCATCTAATATAAATTAG...GAGTTGGACAGGTTAGATTGTAT chr1Hg0000051
    ## [4]  2000 GCCACATGGGACACATGGAAGTT...ATCAAGTTAACCTGGAACTCTGC chr1Hg0000061
    ## [5]  2000 AATAAACCCCAAAACCACAAAAC...TTCTTCTACAAAGTTAAGTTAAG chr1Hg0000071
    ## [6]  2000 GTTGTTTGCAATGATAAAATCAC...GCTAGCGGTGGTAGCCCCCCTCG chr1Hg0000101
    ## [7]  2000 TTATATTTAAAGATGTAAATGTT...GCACCCAATGCAATGAGGCACTG chr1Hg0000111

# EXPORT SEQUENCES

``` r
# output_file <- "data/promoters_Min-Yao.fasta"
# Biostrings::writeXStringSet(my_promoters, filepath = output_file)
#
# save(my_promoters, file = "data/my_promoter.rda")

# This possibly doesn't belongs to the `CSCpromoters` scope
write_fasta <- function(.sequences, .output_file) {
  Biostrings::writeXStringSet(.sequences, filepath = .output_file)
}

write_rda <- function(.sequences, .output_file) {
  save(.sequences, file = .output_file)
}

# Not being used currently
write_sequences <- function(.sequences, .output_file) {
  switch (
    tools::file_ext(.output_file),
    fasta = write_fasta(.sequences, .output_file),
    fa    = write_fasta(.sequences, .output_file),
    rda   = write_rda(.sequences, .output_file)
  )
}

my_promoters %>% write_fasta("data/promoters_Min-Yao.fasta")
my_promoters %>% write_rda("data/promoters_Min-Yao.rda")
```
