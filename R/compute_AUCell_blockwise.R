#' Compute AUCell Gene-Set Scores in Sample Blocks
#'
#' @description
#' Computes AUCell enrichment scores from a gene-by-sample expression matrix
#' while processing samples in blocks to reduce peak memory usage. Gene sets
#' can be supplied as individual upregulated sets or as paired upregulated and
#' downregulated sets. For paired sets, the final composite score is calculated
#' as the upregulated AUCell score minus the downregulated AUCell score.
#'
#' @details
#' The expression input must contain genes in rows and samples or cells in
#' columns. Every expression row must have a unique, non-missing gene
#' identifier.
#'
#' Gene sets must be supplied as a named list. Set names ending in
#' `up_suffix` and `dn_suffix` are paired by removing those suffixes. A gene set
#' without either suffix is treated as an upregulated component.
#'
#' For example, gene sets named `"TE_up"` and `"TE_dn"` produce the following
#' output columns:
#'
#' * `TE_up`: AUCell score for the upregulated component.
#' * `TE_dn`: AUCell score for the downregulated component.
#' * `TE`: composite score calculated as `TE_up - TE_dn`.
#'
#' If only an upregulated component is available, its AUCell score is used as
#' the composite score. If only a downregulated component is available, its
#' negative AUCell score is used as the composite score.
#'
#' When `gene_names` is supplied, it is assigned directly to the expression
#' row names. In this mode, the identifiers in `gene_sets` must use the same
#' identifier system as `gene_names`.
#'
#' When `ensembl` is supplied, Ensembl version suffixes are removed before the
#' identifiers are assigned as expression row names. The supplied gene sets
#' are then assumed to contain gene symbols and are mapped internally to
#' Ensembl identifiers using the EmbryoRadar gene annotation.
#'
#' Supply only one of `gene_names` and `ensembl`. If neither is supplied,
#' `tpm` must already have valid row names.
#'
#' When `raw_counts` is supplied, counts are converted to TPM using gene
#' lengths from the EmbryoRadar annotation. Either `gene_names` or `ensembl`
#' must also be supplied so that gene lengths can be matched.
#'
#' Gene sets are restricted to genes present in the expression matrix. Empty
#' sets and set pairs without sufficient overlapping genes are removed before
#' scoring.
#'
#' @param tpm A numeric matrix, data frame, or `data.table` containing
#'   TPM-like expression values. Genes must be in rows and samples or cells in
#'   columns. Row identifiers must be unique.
#'
#' @param gene_names An optional character vector of gene identifiers with
#'   length equal to `nrow(tpm)` or `nrow(raw_counts)`. These identifiers are
#'   assigned directly as expression row names. Gene sets must use the same
#'   identifier system.
#'
#' @param ensembl An optional character vector of Ensembl gene identifiers
#'   with length equal to `nrow(tpm)` or `nrow(raw_counts)`. Version suffixes,
#'   such as `".12"`, are removed. When this argument is used, gene sets are
#'   assumed to contain gene symbols and are mapped internally to Ensembl
#'   identifiers.
#'
#' @param raw_counts An optional numeric gene-by-sample matrix of raw counts.
#'   When provided, counts are converted to TPM before AUCell scoring. Either
#'   `gene_names` or `ensembl` must also be supplied.
#'
#' @param gene_sets A named list of character vectors containing gene-set
#'   members. Names may use `up_suffix` and `dn_suffix` to define paired
#'   components. Gene identifiers must be compatible with the identifiers used
#'   for the expression matrix.
#'
#' @param chunk_size A positive integer specifying the maximum number of
#'   samples or cells processed in each block. If `NULL`, non-positive, or
#'   greater than or equal to the number of samples, all samples are processed
#'   in one block. Defaults to `4000`.
#'
#' @param verbose Logical; if `TRUE`, progress messages are printed. Defaults
#'   to `TRUE`.
#'
#' @param up_suffix Character suffix identifying upregulated gene-set
#'   components. Defaults to `"_up"`.
#'
#' @param dn_suffix Character suffix identifying downregulated gene-set
#'   components. Defaults to `"_dn"`.
#'
#' @param auc_max_rank_frac Numeric value between 0 and 1 specifying the
#'   fraction of ranked genes used as `aucMaxRank` by
#'   [AUCell::AUCell_calcAUC()]. Defaults to `0.05`.
#'
#' @return
#' A `data.table` with one row per sample or cell. The first column,
#' `sampleID`, contains the original expression column names. Remaining columns
#' contain composite scores followed by the available upregulated and
#' downregulated component scores.
#'
#' @seealso
#' [AUCell::AUCell_buildRankings()], [AUCell::AUCell_calcAUC()],
#' [AUCell::getAUC()]
#'
#' @examples
#' \dontrun{
#' expression_matrix <- matrix(
#'   runif(1000),
#'   nrow = 100,
#'   ncol = 10,
#'   dimnames = list(
#'     paste0("Gene", seq_len(100)),
#'     paste0("Sample", seq_len(10))
#'   )
#' )
#'
#' example_sets <- list(
#'   TE_up = paste0("Gene", 1:20),
#'   TE_dn = paste0("Gene", 21:40),
#'   EPI_up = paste0("Gene", 41:60)
#' )
#'
#' scores <- compute_AUCell_blockwise(
#'   tpm = expression_matrix,
#'   gene_sets = example_sets,
#'   chunk_size = 5
#' )
#' }
#'
#' @importFrom data.table data.table
#' @export

compute_AUCell_blockwise <- function (tpm = NULL, gene_names = NULL, ensembl = NULL, raw_counts = NULL,
                                      gene_sets = NULL, chunk_size = 4000, verbose = TRUE, up_suffix = "_up",
                                      dn_suffix = "_dn", auc_max_rank_frac = 0.05)
{
  suppressPackageStartupMessages({
    library(data.table)
    library(AUCell)
  })
  if (is.null(gene_sets)) {
    stop("Please provide `gene_sets` as a named list of gene vectors.")
  }
  if (!is.list(gene_sets)) {
    stop("`gene_sets` must be a named list.")
  }
  if (is.null(names(gene_sets)) || any(names(gene_sets) ==
                                       "")) {
    stop("`gene_sets` must be a named list.")
  }
  gsls <- gene_sets
  data_dir <- system.file("extdata/classifier2.0", package = "EmbryoRadar")
  geneAnnPIE <- readRDS(file.path(data_dir, "geneAnnPIE.rds"))
  geneAnnPIE$ensemble <- gsub("\\..*", "", geneAnnPIE$gene_id)
  if (!is.null(raw_counts)) {
    if (verbose) {
      message("[PIE] Converting raw counts to TPM ...")
    }
    if (!is.null(gene_names)) {
      gene_lengths <- geneAnnPIE$width[match(gene_names,
                                             geneAnnPIE$gene_name)]
    }
    else if (!is.null(ensembl)) {
      ensembl_clean <- gsub("\\..*", "", ensembl)
      gene_lengths <- geneAnnPIE$width[match(ensembl_clean,
                                             geneAnnPIE$ensemble)]
    }
    else {
      stop("Provide either `gene_names` or `ensembl` with `raw_counts`.")
    }
    tpm <- calculate_tpm(raw_counts, gene_lengths)
  }
  if (is.null(tpm)) {
    stop("Provide either `tpm` or `raw_counts`.")
  }
  if (!is.null(gene_names)) {
    rownames(tpm) <- gene_names
  }
  else if (!is.null(ensembl)) {
    rownames(tpm) <- gsub("\\..*", "", ensembl)
  }
  else if (is.null(rownames(tpm))) {
    stop("If `gene_names` and `ensembl` are both NULL, `tpm` must already have rownames.")
  }
  if (!is.null(ensembl)) {
    gene_symbol_to_ensembl <- setNames(geneAnnPIE$ensemble,
                                       geneAnnPIE$gene_name)
    gsls <- lapply(gsls, function(genes) {
      genes <- unique(as.character(genes))
      mapped <- gene_symbol_to_ensembl[genes]
      unname(mapped[!is.na(mapped)])
    })
  }
  gsid_ok <- vapply(gsls, function(x) length(x) > 0, logical(1))
  if (any(!gsid_ok) && verbose) {
    message("[PIE] Dropping empty gene sets: ", paste0(names(gsls)[!gsid_ok],
                                                       collapse = ", "))
  }
  gsls <- gsls[gsid_ok]
  if (length(gsls) == 0) {
    stop("No usable gene sets remain after removing empty sets.")
  }
  all_names <- names(gsls)
  base_from <- function(nm) {
    sub(paste0("(", up_suffix, "|", dn_suffix, ")$"), "",
        nm)
  }
  bases <- unique(base_from(all_names))
  pair_list <- lapply(bases, function(base) {
    nm_up <- paste0(base, up_suffix)
    nm_dn <- paste0(base, dn_suffix)
    nm_raw <- base
    up <- NULL
    dn <- NULL
    if (nm_up %in% all_names) {
      up <- unique(gsls[[nm_up]])
    }
    if (nm_raw %in% all_names) {
      up <- unique(c(up, gsls[[nm_raw]]))
    }
    if (nm_dn %in% all_names) {
      dn <- unique(gsls[[nm_dn]])
    }
    if (is.null(up) && is.null(dn)) {
      return(NULL)
    }
    list(base = base, up = up, dn = dn)
  })
  names(pair_list) <- bases
  pair_list <- pair_list[!vapply(pair_list, is.null, logical(1))]
  if (!length(pair_list)) {
    stop("No usable gene-set bases found.")
  }
  present_genes <- rownames(tpm)
  pair_list <- lapply(pair_list, function(p) {
    p$up <- if (!is.null(p$up)) {
      intersect(unique(p$up), present_genes)
    }
    else {
      NULL
    }
    p$dn <- if (!is.null(p$dn)) {
      intersect(unique(p$dn), present_genes)
    }
    else {
      NULL
    }
    has_up <- !is.null(p$up) && length(p$up) > 1
    has_dn <- !is.null(p$dn) && length(p$dn) > 1
    if (has_up || has_dn) {
      return(p)
    }
    NULL
  })
  pair_list <- pair_list[!vapply(pair_list, is.null, logical(1))]
  if (!length(pair_list)) {
    stop("After trimming to present genes, no gene sets remain.")
  }
  n_cells <- ncol(tpm)
  sample_id <- colnames(tpm)
  if (is.null(sample_id)) {
    sample_id <- paste0("Sample", seq_len(n_cells))
  }
  if (is.null(chunk_size) || chunk_size <= 0L || chunk_size >=
      n_cells) {
    starts <- 1L
    ends <- n_cells
  }
  else {
    starts <- seq.int(1L, n_cells, by = as.integer(chunk_size))
    ends <- pmin(starts + as.integer(chunk_size) - 1L, n_cells)
  }
  if (verbose) {
    message("[PIE] AUCell scoring ", length(pair_list),
            " base set(s) across ", n_cells, " cells/samples in ",
            length(starts), " block(s).")
  }
  base_names <- vapply(pair_list, `[[`, character(1), "base")
  up_cols <- paste0(base_names, up_suffix)
  dn_cols <- paste0(base_names, dn_suffix)
  comp_cols <- base_names
  block_tbls <- vector("list", length(starts))
  for (bi in seq_along(starts)) {
    start_i <- starts[bi]
    end_i <- ends[bi]
    idx <- start_i:end_i
    if (verbose) {
      message(sprintf("[PIE] Block %d/%d: cells/samples %d-%d",
                      bi, length(starts), start_i, end_i))
    }
    if (inherits(tpm, "data.table")) {
      expr_block <- as.matrix(tpm[, ..idx, drop = FALSE])
    }
    else {
      expr_block <- as.matrix(tpm[, idx, drop = FALSE])
    }
    colnames(expr_block) <- sample_id[idx]
    rownames(expr_block) <- rownames(tpm)
    rankings <- AUCell::AUCell_buildRankings(expr_block,
                                             plotStats = FALSE, verbose = FALSE)
    aucMaxRank <- ceiling(nrow(expr_block) * auc_max_rank_frac)
    aucMaxRank <- max(1L, min(aucMaxRank, nrow(expr_block)))
    auc_geneSets <- list()
    for (p in pair_list) {
      b <- p$base
      if (!is.null(p$up) && length(p$up) > 0) {
        auc_geneSets[[paste0(b, up_suffix)]] <- unique(p$up)
      }
      if (!is.null(p$dn) && length(p$dn) > 0) {
        auc_geneSets[[paste0(b, dn_suffix)]] <- unique(p$dn)
      }
    }
    auc_obj <- AUCell::AUCell_calcAUC(auc_geneSets, rankings,
                                      aucMaxRank = aucMaxRank, verbose = FALSE)
    auc_mat <- AUCell::getAUC(auc_obj)
    dt <- data.table(sampleID = colnames(expr_block))
    for (nm in rownames(auc_mat)) {
      dt[[nm]] <- as.numeric(auc_mat[nm, ])
    }
    for (b in base_names) {
      nm_up <- paste0(b, up_suffix)
      nm_dn <- paste0(b, dn_suffix)
      has_up <- nm_up %in% names(dt)
      has_dn <- nm_dn %in% names(dt)
      if (has_up && has_dn) {
        dt[[b]] <- dt[[nm_up]] - dt[[nm_dn]]
      }
      else if (has_up) {
        dt[[b]] <- dt[[nm_up]]
      }
      else if (has_dn) {
        dt[[b]] <- -dt[[nm_dn]]
      }
      else {
        dt[[b]] <- NA_real_
      }
    }
    keep_up <- up_cols[up_cols %in% names(dt)]
    keep_dn <- dn_cols[dn_cols %in% names(dt)]
    data.table::setcolorder(dt, c("sampleID", comp_cols,
                                  keep_up, keep_dn))
    block_tbls[[bi]] <- dt
  }
  scores_all <- data.table::rbindlist(block_tbls, use.names = TRUE,
                                      fill = TRUE)
  return(scores_all)
}
