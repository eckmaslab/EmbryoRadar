#' Compute transcriptional enrichment (TE) singscores for bulk or single-cell data
#'
#' This function computes singscores across predefined or user-supplied gene sets
#' associated with transcriptional programs, developmental stages, and trophoblast
#' lineage states. It supports both normalized TPM matrices and raw read counts
#' (which are internally converted to TPM). The resulting singscores quantify
#' per-sample activation of each gene set, enabling downstream visualization
#' or integration with PIE classifier outputs.
#'
#' @param tpm A numeric matrix of transcript-per-million (TPM) expression values
#'   with genes in rows and samples in columns. Either \code{tpm} or
#'   \code{raw_counts} must be provided.
#' @param gene_names Optional character vector of gene symbols corresponding to
#'   the rownames of \code{tpm} or \code{raw_counts}. Used for gene-length matching
#'   and mapping to the packaged gene annotations.
#' @param ensembl Optional character vector of Ensembl gene IDs corresponding to
#'   the rownames of \code{tpm} or \code{raw_counts}. If provided, gene sets will
#'   be automatically mapped to Ensembl identifiers.
#' @param raw_counts Optional raw count matrix (genes × samples). If supplied,
#'   TPM values will be computed internally using gene lengths from the packaged
#'   annotation data.
#' @param additional_gene_sets Optional named list of user-defined gene sets.
#'   Each element should be a character vector of gene symbols or Ensembl IDs.
#'   These will be appended to or replace the default packaged gene sets.
#' @param additional_gene_sets_only Logical; if \code{TRUE}, only the
#'   \code{additional_gene_sets} will be used (no default gene sets are loaded).
#'   Default is \code{FALSE}.
#'
#' @details
#' Internally, this function:
#' \enumerate{
#'   \item Loads the packaged developmental and TE-related gene sets via
#'         \code{\link{load_TE_gene_sets}} unless \code{additional_gene_sets_only = TRUE}.
#'   \item Converts raw counts to TPM if requested.
#'   \item Ranks genes using \code{singscore::rankGenes()}.
#'   \item Computes per-sample singscores using
#'         \code{singscore::simpleScore()} for each up/down gene set pair.
#' }
#'
#' The output can be directly merged with PIE classifier predictions or visualized
#' using standard plotting tools.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{singscores}}{A \code{data.table} containing per-sample scores
#'         for each gene set (columns = gene sets, rows = samples).}
#'   \item{\code{geneSets}}{The list of gene sets actually used in scoring.}
#' }
#'
#' @import data.table
#' @import dplyr
#' @import magrittr
#' @import singscore
#' @import readxl
#'
#' @examples
#' \dontrun{
#' # Example: compute singscores from TPM matrix
#' tpm <- matrix(rexp(5000), nrow = 1000, ncol = 5,
#'               dimnames = list(paste0("Gene", 1:1000),
#'                               paste0("Sample", 1:5)))
#' res <- compute_TE_singscores(tpm = tpm, gene_names = rownames(tpm))
#' head(res$singscores)
#' }
#'
#' @seealso \code{\link{load_TE_gene_sets}}, \code{\link{singscore::simpleScore}}
#'
#' @export

compute_TE_singscores <- function(tpm = NULL, gene_names = NULL, ensembl = NULL, raw_counts = NULL,
                               additional_gene_sets = NULL, additional_gene_sets_only = FALSE) {
  library(data.table); library(dplyr); library(magrittr); library(singscore); library(readxl)

  data_dir <- system.file("extdata/classifier2.0", package = "EmbryoRadar")

  geneAnnPIE <- readRDS(file.path(data_dir, "geneAnnPIE.rds"))
  geneAnnPIE$ensemble <- gsub("\\..*", "", geneAnnPIE$gene_id)

  # test <- readRDS("inst/extdata/classifier2.0/geneAnnPIE.rds")

  gsls <- if (isTRUE(additional_gene_sets_only)) {
    c(additional_gene_sets)
  } else {
    c(additional_gene_sets, load_TE_gene_sets())
  }

  if (!is.null(raw_counts)) {
    message("calculating tpm")
    if (!is.null(gene_names)) {
      gene_lengths <- geneAnnPIE$width[match(gene_names, geneAnnPIE$gene_name)]
    } else if (!is.null(ensembl)) {
      gene_lengths <- geneAnnPIE$width[match(ensembl, geneAnnPIE$ensemble)]
    } else {
      stop("Provide either gene_names or ensembl with raw_counts.")
    }
    tpm <- calculate_tpm(raw_counts, gene_lengths)
  }

  rankData <- singscore::rankGenes(tpm)

  if (is.null(gene_names) & is.null(ensembl)) {
    rownames(rankData) <- rownames(tpm)
  } else if (!is.null(gene_names)) {
    rownames(rankData) <- gene_names
  } else if (!is.null(ensembl)) {
    rownames(rankData) <- ensembl
  }

  # Map gene sets to Ensembl if user supplied ensembl ids
  if (is.null(gene_names) & !is.null(ensembl)) {
    gene_map <- setNames(nm = geneAnnPIE$gene_name, object = geneAnnPIE$ensemble) # gene_name -> ensembl
    gsls <- lapply(gsls, function(genes) {
      mapped <- gene_map[genes]; unname(mapped[!is.na(mapped)])
    })
  }
  gsid <- unlist(lapply(gsls, function(x) {length(x) > 1}))
  if (any(!gsid)) {
    message("geneset: ",paste0(paste0(names(gsls[!gsid]), collapse = "|")," lost during gene mapping"))
  }
  gsls <- gsls[gsid]

  pb <- txtProgressBar(min = 0, max = length(gsls), style = 3)
  scc <- suppressWarnings(lapply(seq_along(gsls), function(i) {
    gs <- names(gsls)[i]; setTxtProgressBar(pb, i)
    base <- sub("_up$|_dn$", "", gs)
    if (grepl("_dn$", gs)) return(NULL)
    if ((any(grepl(paste0(base,"_dn$"), names(gsls)))) & (grepl("_up$", gs))) {
      s <- singscore::simpleScore(rankData,
                                  upSet = unique(gsls[[gs]]),
                                  downSet = unique(gsls[[paste0(base, "_dn")]]))
      data.table::as.data.table(dplyr::select(s, TotalScore)) %>% magrittr::set_colnames(base)
    } else if (all(!(grepl(paste0(base,"_dn$"), names(gsls))))) {
      s <- singscore::simpleScore(rankData, upSet = unique(gsls[[gs]]))
      data.table::as.data.table(dplyr::select(s, TotalScore)) %>% magrittr::set_colnames(gs)
    }
  }))
  close(pb)
  scc <- scc[!(unlist(lapply(scc, nrow)) != ncol(rankData))]
  scc <- dplyr::bind_cols(sampleID = colnames(rankData), scc)

  list(singscores = data.table::data.table(scc), geneSets = gsls)
}


