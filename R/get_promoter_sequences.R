#' @title Extract promoter sequences from the final promoter sizes
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' @inheritParams get_promoters
#' @inheritDotParams get_promoters
#' 
#' @return DNAStringSet
#' @export
#' 
get_promoter_sequences <- function(
    .distances,
    .FASTA_list = NULL,
    .txdb       = NULL,
    .pb         = NULL,
    .pb_format  = ":what - [:bar] :percent (:spin)",
    .locus_var  = "locus_tag",
    .chr_var    = "chr",
    .dist_var   = "promoter_size",
    .downstream = 0,
    .parallel   = FALSE,
    ...
) {
  .fun <- ifelse(
    .parallel,
    no  = {get_promoter_sequences.serial},
    yes = {get_promoter_sequences.parallel}
  )
  
  .fun(
    .distances,
    .FASTA_list,
    .txdb,
    .pb,
    .pb_format,
    .locus_var,
    .chr_var,
    .dist_var,
    .downstream,
    ...
  )
}


#' Extract promoter sequences from the final promoter sizes (NEW - SERIAL)
#'
#' @inheritParams get_promoters
#' @inheritDotParams get_promoters
#'
#' @return DNAStringSet
#' @export
#'
get_promoter_sequences.serial <- function(
    .distances,
    .FASTA_list = NULL,
    .txdb       = NULL,
    .pb         = NULL,
    .pb_format  = ":what - [:bar] :percent (:spin)",
    .locus_var  = "locus_tag",
    .chr_var    = "chr",
    .dist_var   = "promoter_size",
    .downstream = 0,
    ...
) {
  # This is not working as expected, need to check whether the txdb exists/was 
  # assigned; tryCatch() maybe
  # if(!quote(.txdb) %>% as.character() %>% exists()) stop(".txdb is not set") 
  if(is.null(.txdb)) stop(".txdb is NULL")
  
  .sequences <- .distances %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(.chr_var))) %>%
    dplyr::group_map(.f = function(.distances_curr_chr, .keys) {
      # Initiate (or not) the progress bar - can also receive external object
      if(is.null(.pb)) {
        .pb <- progress::progress_bar$new(
          format = .pb_format,
          total  = .distances_curr_chr %>%
            dplyr::select(dplyr::all_of(.locus_var)) %>% 
            unique() %>% 
            nrow()
        )
      }
      
      # Start the progress bar
      .pb$tick(0)
      
      # Assign the current chromosome name to an object ("character")
      #.curr_chr <- .keys %>% dplyr::pull(.chr_var) # Issue #25
      .curr_chr <- .keys %>%
        dplyr::pull(.chr_var) %>% 
        as.character()

      # Should trimming of promoter sizes be moved here???
      
      # Open FASTA file for each chromosome b4 extracting promoters in each
      # - This is highly problematic! Currently assumes the input files will 
      #   have a given filename structure. There must be better ways to do that!
      #.fasta_file <- .FASTA_list[[.curr_chr]] %>% # Issue #25
      .fasta_file <- .FASTA_list %>%
        purrr::pluck(.curr_chr) %>%
        Rsamtools::FaFile() %>% 
        open()
      
      # Get promoters
      .sequences_curr_chr <- .distances_curr_chr %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(.locus_var))) %>%
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
            # GenomicFeatures::transcripts(
            #   filter = list("tx_name" = .curr_locus)
            # ) %>%
            GenomicFeatures::transcripts() %>% 
            subset(stringr::str_detect(tx_name, .curr_locus)) %>% 
            # print() %>%
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
        do.call(what = c, args = .)
    }) %>% 
    do.call(what = c, args = .)
  
  .sequences
}


#' Extract promoter sequences from the final promoter sizes (NEW - PARALLEL)
#'
#' @inheritParams get_promoters
#' @inheritDotParams get_promoters
#'
#' @return DNAStringSet
#' @export
#'
get_promoter_sequences.parallel <- function(
    .distances,
    .FASTA_list = NULL,
    .txdb       = NULL,
    .pb         = NULL,
    .pb_format  = ":what - [:bar] :percent (:spin)",
    .locus_var  = "locus_tag",
    .chr_var    = "chr",
    .dist_var   = "promoter_size",
    .downstream = 0,
    ...
) {
  # This is not working as expected, need to check whether the txdb exists/was 
  # assigned; tryCatch() maybe
  # if(!quote(.txdb) %>% as.character() %>% exists()) stop(".txdb is not set") 
  if(is.null(.txdb)) stop(".txdb is NULL")
  
  # (parallel::detectCores()-1) %>% 
  #   print()
  # stop()
  
  .cluster <- multidplyr::new_cluster(n = 4)
  .cluster %>% multidplyr::cluster_library("dplyr")
  .cluster %>% multidplyr::cluster_library("GenomicFeatures")
  
  .sequences <- .distances %>%
    # plyr::dlply(.variables = .chr_var, .fun = function(.distances_curr_chr) {
    #   # Get keys table
    #   .keys <- .distances_curr_chr %>%
    #     dplyr::select(dplyr::all_of(.chr_var)) %>%
    #     unique()
    # }) #%>% #do.call(what = c, args = .)
    dplyr::group_by(dplyr::across(dplyr::all_of(.chr_var))) %>%
    dplyr::group_map(.f = function(.distances_curr_chr, .keys) {
      # Initiate (or not) the progress bar - can also receive external object
      if(is.null(.pb)) {
        .pb <- progress::progress_bar$new(
          format = .pb_format,
          total  = .distances_curr_chr %>%
            dplyr::select(dplyr::all_of(.locus_var)) %>% 
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
      # .fasta_file <- paste0(
      #   .folder,
      #   paste0("Hordeum_vulgare.refseq[", .curr_chr, "].fasta")
      # ) %>%
      .fasta_file <- .FASTA_list[[.curr_chr]] %>%
        Rsamtools::FaFile() %>% 
        open()
      
      # Get promoters
      .sequences_curr_chr <- .distances_curr_chr %>%
        # plyr::dlply(
        #   .variables = .locus_var,
        #   .fun       = function(.distance_curr_locus) {
        #     # Get keys table
        #     .each_keys <- .distance_curr_locus %>%
        #       dplyr::select(dplyr::all_of(.locus_var)) %>%
        #       unique()
        #   }
        # )
        dplyr::group_by(dplyr::across(dplyr::all_of(.locus_var))) %>%
        multidplyr::partition(.cluster) %>% 
        # dplyr::group_map(.f = function(.distance_curr_locus, .each_keys) {
        dplyr::do((function(.distance_curr_locus, .txdb) {
          # Get keys table
          .each_keys <- .distance_curr_locus %>%
            dplyr::select(dplyr::all_of(.locus_var)) %>%
            unique()
          
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
        })(., .txdb)) %>% 
        dplyr::collect()
      
      # Close FASTA file 
      close(.fasta_file)
      
      # End progress bar
      rm(.pb)
      
      .sequences_curr_chr %>%
        do.call(what = c, args = .)
    }) %>%
    # dplyr::collect() %>% 
    do.call(what = c, args = .)
  
  .sequences
}

# get_promoter_sequences <- function(
#     .distances,
#     .folder     = NULL,
#     .txdb       = NULL,
#     .pb         = NULL,
#     .pb_format  = ":what - [:bar] :percent (:spin)",
#     .locus_var  = "locus_tag",
#     .chr_var    = "chr",
#     .dist_var   = "promoter_size",
#     .downstream = 0,
#     ...
# ) {
#   # This is not working as expected, need to check whether the txdb exists/was 
#   # assigned; tryCatch() maybe
#   # if(!quote(.txdb) %>% as.character() %>% exists()) stop(".txdb is not set") 
#   if(is.null(.txdb)) stop(".txdb is NULL")
#   
#   .sequences <- .distances %>%
#     dplyr::group_by(dplyr::across(dplyr::all_of(.chr_var))) %>%
#     dplyr::group_map(.f = function(.distances_curr_chr, .keys) {
#       # Initiate (or not) the progress bar - can also receive external object
#       if(is.null(.pb)) {
#         .pb <- progress::progress_bar$new(
#           format = .pb_format,
#           total  = .distances_curr_chr %>%
#             dplyr::select(dplyr::all_of(.locus_var)) %>% 
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
#         dplyr::group_by(dplyr::across(dplyr::all_of(.locus_var))) %>%
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
#         do.call(what = c, args = .)
#     }) %>% 
#     do.call(what = c, args = .)
#   
#   .sequences
# }