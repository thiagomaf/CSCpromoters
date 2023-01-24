#' Wrap-up function to query promoter sequences
#'
#' @param .annotations Preferably a `data.table` (but possibly a `data.frame`) 
#'                     of annotations as described {HERE}. `data.table` 
#'                     computations are normally faster.
#'                     
#' @param .distances   {TODO}
#'                     
#' @param .keep        Either a `character` vector of loci names (present in the
#'                     `.annotations` table) or a `numeric` vector of indexes 
#'                     (of the `.annotations` table).
#'                     
#' @param .txdb        {TODO}
#' 
# @param .folder      A `character` indicating the folder in which the FASTA 
#                     files containing the genome sequences of each individual
#                     chromosome. Currently, the library assumes that the
#                     filename will have the following format: 
#                     `paste0("Hordeum_vulgare.refseq[", .curr_chr, "].fasta")`
#                     .
#                     
#' @param .FASTA_list  A named `character` list. Each element is a full path and
#'                     filename to each chromosome. Element names indicate 
#'                     chromosome names matching those found in the annotation 
#'                     table.
#'                     
#' @param .min_size    An integer indicating the minimum promoter size. Any 
#'                     _loci_ identified to have sequence length shorter than 
#'                     `.min_size` will be ignored.
#'                     
#' @param .max_size    An integer indicating the maximum promoter size. Any 
#'                     _loci_ identified to have sequence length longer than 
#'                     `.max_size` will be set (trimmed) to equals `.max_size`.
#'                     
#' @param .downstream  An integer indicating the sequence length downstream the 
#'                     start site to be included in the output promoter 
#'                     sequence. Presumably, will be reduced from the upstream 
#'                     sequence length but have never tested. Defaults to `0`.
#'                     
#' @param .pb          An optional progress bar object created via 
#'                     `progress::progress_bar()`. Defaults to `NULL`, case in
#'                     which a progress bar will be generated internally.
#' 
#' @param .pb_format   A `character` indicating the format of the progress bar.
#'                     A number of tokens can be used here, see them below. It 
#'                     defaults to `":what - [:bar] :percent (:spin)"`, which 
#'                     means that information on the locus being handled in 
#'                     shown on the left (`:what`), followed by the progress 
#'                     bar is within brackets (`[:bar]`), and ended by the 
#'                     percentage and a spinning character on the right 
#'                     (`:percent (:spin)`).
#' 
#' @param .locus_var   A `character` indicating the name of the column in the 
#'                     `.annotation` table containing loci names. This can be 
#'                     altered to fit input annotation tables with different 
#'                     column names (not recommended). Defaults to `"locus_tag"`
#'                     .
#'                     
#' @param .closest_var A `character` indicating the name of the column in the 
#'                     `.annotation` table containing closest loci names. This 
#'                     can be altered to fit input annotation tables with 
#'                     different column names (not recommended). Defaults to 
#'                     `"locus_tag"`.
#'                      
#' @param .chr_var     A `character` indicating the name of the column in the 
#'                     `.annotation` table containing chromosome names. This 
#'                     can be altered to fit input annotation tables with 
#'                     different column names (not recommended). Defaults to 
#'                     `"chr"`.
#' @param .dist_var    A `character` indicating the name of the output column 
#'                     table containing the final promoter lengths. There is 
#'                     absolutely no reason you would need to change that. 
#'                     Defaults to `"promoter_size"`.
#'                     
#' @param .strand_var  A `character` indicating the name of the column in the 
#'                     `.annotation` table containing each loci strand 
#'                     direction. This can be altered to fit input annotation 
#'                     tables with different column names (not recommended). 
#'                     Defaults to `"strand"`.
#'                     
#' @param .start_var   A `character` indicating the name of the column in the 
#'                     `.annotation` table containing each loci start 
#'                     coordinate. This can be altered to fit input annotation 
#'                     tables with different column names (not recommended). 
#'                     Defaults to `"begin"`.
#'                     
#' @param .end_var     A `character` indicating the name of the column in the 
#'                     `.annotation` table containing each loci end coordinate.
#'                     This can be altered to fit input annotation tables with 
#'                     different column names (not recommended). Defaults to 
#'                     `"end"`.
#'                     
#' @param .parallel    A `logical`. If TRUE, apply computations in parallel. 
#'                     Defauts to `FALSE`.
#'                     
#' @param ...          Future-proof stuff.
#'
#' @inheritSection progress::progress_bar Tokens
#'
#' @return DNAStringSet
#' @export
#'
get_promoters <- function(
    .annotations,
    .keep        = NULL,
    .txdb        = NULL,
    # .folder     = NULL,
    .FASTA_list  = NULL,
    .min_size    = 100,
    .max_size    = 2000,
    .downstream  = 0,
    .pb          = NULL,
    .pb_format   = ":what - [:bar] :percent (:spin)",
    .locus_var   = "locus_tag",
    .closest_var = "closest_locus",
    .chr_var     = "chr",
    .strand_var  = "strand",
    .start_var   = "begin",
    .end_var     = "end",
    .dist_var    = "promoter_size",
    .distances   = NULL,
    .parallel    = FALSE,
    ...
) {
  .annotations %>%
    filter_locus(.keep) %>%
    get_promoter_distances() %>% 
    set_promoter_sizes(.min_size, .max_size) %>% 
    get_promoter_sequences(.FASTA_list, .txdb)
}