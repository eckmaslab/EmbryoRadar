#' Calculate TPM (Transcripts Per Million) from raw read counts
#'
#' This function converts a raw read count matrix into TPM (Transcripts Per Million)
#' expression values using supplied gene lengths. TPM normalization accounts for
#' both sequencing depth and gene length, enabling direct comparison of expression
#' values across genes and samples.
#'
#' @param raw_counts A numeric matrix of raw read counts with genes in rows and
#'   samples in columns. Rownames and column names should correspond to gene
#'   identifiers and sample IDs, respectively.
#' @param gene_lengths A numeric vector of gene lengths in base pairs corresponding
#'   to the rows of \code{raw_counts}.
#'
#' @details
#' The TPM for each gene \(i\) in sample \(j\) is computed as:
#' \deqn{
#'   TPM_{ij} = \frac{RPK_{ij}}{\sum_i RPK_{ij}} \times 10^6
#' }
#' where \eqn{RPK_{ij} = \frac{counts_{ij}}{gene\_length_i / 1000}}.
#'
#' Genes with missing or zero lengths are excluded from scaling but retained
#' as \code{NA} in the output matrix.
#'
#' @return A numeric matrix of TPM-normalized expression values with the same
#' dimensions and dimnames as \code{raw_counts}.
#'
#' @examples
#' \dontrun{
#' # Example: Convert raw counts to TPM
#' counts <- matrix(sample(1:1000, 20, replace = TRUE), nrow = 5)
#' lengths <- c(1500, 2000, 800, 1000, 2200)
#' tpm <- calculate_tpm(counts, lengths)
#' head(tpm)
#' }
#'
#' @import data.table
#' @seealso \code{\link{compute_TE_singscores}}, \code{\link{load_TE_gene_sets}}
#' @export
#'
calculate_tpm <- function(raw_counts, gene_lengths) {
  library(data.table)
  if (!is.matrix(raw_counts)) raw_counts <- as.matrix(raw_counts)
  gene_lengths_kb <- gene_lengths / 1000
  valid <- which(!is.na(gene_lengths_kb) & gene_lengths_kb > 0)
  counts_valid <- raw_counts[valid, , drop = FALSE]
  lengths_valid <- gene_lengths_kb[valid]
  rpk <- counts_valid / lengths_valid
  scaling <- colSums(rpk, na.rm=TRUE) / 1e6
  tpm_valid <- sweep(rpk, 2, scaling, "/")
  tpm <- matrix(NA_real_, nrow = nrow(raw_counts), ncol = ncol(raw_counts),
                dimnames = dimnames(raw_counts))
  tpm[valid, ] <- tpm_valid
  tpm
}
