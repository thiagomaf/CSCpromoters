% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_promoter_sequences.R
\name{get_promoter_sequences.parallel}
\alias{get_promoter_sequences.parallel}
\title{Extract promoter sequences from the final promoter sizes (NEW - PARALLEL)}
\usage{
get_promoter_sequences.parallel(
  .distances,
  .FASTA_list = NULL,
  .txdb = NULL,
  .pb = NULL,
  .pb_format = ":what - [:bar] :percent (:spin)",
  .locus_var = "locus_tag",
  .chr_var = "chr",
  .dist_var = "promoter_size",
  .downstream = 0,
  ...
)
}
\arguments{
\item{.distances}{{TODO}}

\item{.FASTA_list}{A named \code{character} list. Each element is a full path and
filename to each chromosome. Element names indicate
chromosome names matching those found in the annotation
table.}

\item{.txdb}{{TODO}}

\item{.pb}{An optional progress bar object created via
\code{progress::progress_bar()}. Defaults to \code{NULL}, case in
which a progress bar will be generated internally.}

\item{.pb_format}{A \code{character} indicating the format of the progress bar.
A number of tokens can be used here, see them below. It
defaults to \code{":what - [:bar] :percent (:spin)"}, which
means that information on the locus being handled in
shown on the left (\verb{:what}), followed by the progress
bar is within brackets (\verb{[:bar]}), and ended by the
percentage and a spinning character on the right
(\verb{:percent (:spin)}).}

\item{.locus_var}{A \code{character} indicating the name of the column in the
\code{.annotation} table containing loci names. This can be
altered to fit input annotation tables with different
column names (not recommended). Defaults to \code{"locus_tag"}
.}

\item{.chr_var}{A \code{character} indicating the name of the column in the
\code{.annotation} table containing chromosome names. This
can be altered to fit input annotation tables with
different column names (not recommended). Defaults to
\code{"chr"}.}

\item{.dist_var}{A \code{character} indicating the name of column used to
output the distance from each locus to its closest
upstream locus. I would bet you don't have a reason to
change that. Defaults to \code{"dist"}.}

\item{.downstream}{An integer indicating the sequence length downstream the
start site to be included in the output promoter
sequence. Presumably, will be reduced from the upstream
sequence length but have never tested. Defaults to \code{0}.}

\item{...}{
  Arguments passed on to \code{\link[=get_promoters]{get_promoters}}
  \describe{
    \item{\code{.annotations}}{Preferably a \code{data.table} (but possibly a \code{data.frame})
of annotations as described (HERE). \code{data.table}
computations are normally faster.}
    \item{\code{.keep}}{Either a \code{character} vector of loci names (present in the
\code{.annotations} table) or a \code{numeric} vector of indexes
(of the \code{.annotations} table).}
    \item{\code{.min_size}}{An integer indicating the minimum promoter size. Any
\emph{loci} identified to have sequence length shorter than
\code{.min_size} will be ignored.}
    \item{\code{.max_size}}{An integer indicating the maximum promoter size. Any
\emph{loci} identified to have sequence length longer than
\code{.max_size} will be set (trimmed) to equals \code{.max_size}.}
    \item{\code{.closest_var}}{A \code{character} indicating the name of the column in the
\code{.annotation} table containing closest loci names. This
can be altered to fit input annotation tables with
different column names (not recommended). Defaults to
\code{"locus_tag"}.}
    \item{\code{.size_var}}{A \code{character} indicating the name of the output column
table containing the final promoter lengths. I still
doubt you need to change that. Defaults to
\code{"promoter_size"}.}
    \item{\code{.strand_var}}{A \code{character} indicating the name of the column in the
\code{.annotation} table containing each loci strand
direction. This can be altered to fit input annotation
tables with different column names (not recommended).
Defaults to \code{"strand"}.}
    \item{\code{.start_var}}{A \code{character} indicating the name of the column in the
\code{.annotation} table containing each loci start
coordinate. This can be altered to fit input annotation
tables with different column names (not recommended).
Defaults to \code{"begin"}.}
    \item{\code{.end_var}}{A \code{character} indicating the name of the column in the
\code{.annotation} table containing each loci end coordinate.
This can be altered to fit input annotation tables with
different column names (not recommended). Defaults to
\code{"end"}.}
    \item{\code{.parallel}}{A \code{logical}. If TRUE, apply computations in parallel.
Defauts to \code{FALSE}.}
  }}
}
\value{
DNAStringSet
}
\description{
Extract promoter sequences from the final promoter sizes (NEW - PARALLEL)
}
