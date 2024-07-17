gff_folder <- paste0(
  Sys.getenv("OneDrive"),
  "/CSC/Data/FASTA sequences/Genomes/Hordeum vulgare - Golden Promise/"
)

fasta_list <- (function(
    .chr_list = c(
      "chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H", "chrUn"
    ),
    .folder   = paste0(
      Sys.getenv("OneDrive"),
      "/CSC/Data/FASTA sequences/Genomes/Hordeum vulgare - Golden Promise/"
    )
) {
  .folder |>
    paste0(paste0("Hordeum_vulgare.refseq[", .chr_list, "].fasta")) |>
    magrittr::set_names(.chr_list)
})()


#-------------------------------------------------------------------------------
test_that("Loading GFF", {
  txdb <- paste0(
    gff_folder,
    "Horvul_GP_v1r1_Apollo_30_06_20_named_product_GO.gff3"
  ) |>
    make_txdb(
      .data_source = "Hv - Golden Promise",
      .organism    = "Hordeum vulgare"
    )
})


#-------------------------------------------------------------------------------
test_that("Load annotations", {
  annotations <- txdb |>
    get_txdb_annotation()
})


#-------------------------------------------------------------------------------
test_that("Wrap-up function", {
  annotations |>
    get_promoters(
      .txdb       = txdb,
      .FASTA_list = fasta_list,
      .keep       = c("chr1Hg0040651")
    )
})


#-------------------------------------------------------------------------------
test_that("Explicit pipeline", {
  promoter_sizes <- annotations |>
    filter_locus(
      .keep = c("chr1Hg0040651")
    ) |> 
    get_promoter_distances() |>
    set_promoter_sizes(.min_size = 100, .max_size = 2000)
})


#-------------------------------------------------------------------------------
test_that("Get plots", {
  promoter_sizes |>
    plot_promoter_maps(annotations)
})



#-------------------------------------------------------------------------------
test_that("Get sequences", {
  promoter_sizes |> 
    get_promoter_sequences(.txdb = txdb, .FASTA_list = fasta_list)
})

