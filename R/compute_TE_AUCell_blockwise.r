compute_TE_AUCell_blockwise <- function (tpm = NULL, gene_names = NULL, ensembl = NULL, raw_counts = NULL,
                                         additional_gene_sets = NULL, additional_gene_sets_only = FALSE,
                                         chunk_size = 4000, verbose = TRUE, up_suffix = "_up", dn_suffix = "_dn",
                                         auc_max_rank_frac = 0.05)
{
  library(data.table)
  library(dplyr)
  library(magrittr)
  library(AUCell)
  library(readxl)
  data_dir <- system.file("extdata/classifier2.0", package = "EmbryoRadar")
  geneAnnPIE <- readRDS(file.path(data_dir, "geneAnnPIE.rds"))
  geneAnnPIE$ensemble <- gsub("\\..*", "", geneAnnPIE$gene_id)
  if (additional_gene_sets_only) {
    gsls <- list()
  } else {
    gsls <- load_TE_gene_sets()
  }
  if (!is.null(raw_counts)) {
    if (verbose)
      message("[PIE] Converting raw counts to TPM ...")
    if (!is.null(gene_names)) {
      gene_lengths <- geneAnnPIE$width[match(gene_names,
                                             geneAnnPIE$gene_name)]
    }
    else if (!is.null(ensembl)) {
      gene_lengths <- geneAnnPIE$width[match(ensembl,
                                             geneAnnPIE$ensemble)]
    }
    else {
      stop("Provide either gene_names or ensembl with raw_counts.")
    }
    tpm <- calculate_tpm(raw_counts, gene_lengths)
  }
  if (is.null(tpm))
    stop("Provide either `tpm` or `raw_counts`.")
  if (is.null(gene_names) && is.null(ensembl)) {
    rownames(tpm) <- rownames(tpm)
  } else if (!is.null(gene_names)) {
    rownames(tpm) <- gene_names
  } else if (!is.null(ensembl)) {
    rownames(tpm) <- ensembl
    gene_map <- setNames(nm = geneAnnPIE$gene_name, geneAnnPIE$ensemble)
    gsls <- lapply(gsls, function(genes) {
      mapped <- gene_map[unlist(genes)]
      unname(mapped[!is.na(mapped)])
    })
  }

  gsls <- if (isTRUE(additional_gene_sets_only)) {
    c(additional_gene_sets)
  } else {
    c(additional_gene_sets, load_TE_gene_sets())
  }

  gsid_ok <- vapply(gsls, function(x) length(x) > 0, logical(1))
  if (any(!gsid_ok) && verbose) {
    message("[PIE] Dropping tiny gene sets: ", paste0(names(gsls)[!gsid_ok],
                                                      collapse = ", "))
  }
  gsls <- gsls[gsid_ok]
  all_names <- names(gsls)
  base_from <- function(nm) sub("_(up|dn)$", "", nm)
  bases <- unique(base_from(all_names))
  pair_list <- lapply(bases, function(base) {
    nm_up <- paste0(base, up_suffix)
    nm_dn <- paste0(base, dn_suffix)
    nm_raw <- base
    up <- NULL
    dn <- NULL
    if (nm_up %in% all_names)
      up <- unique(gsls[[nm_up]])
    if (nm_raw %in% all_names)
      up <- unique(c(up, gsls[[nm_raw]]))
    if (nm_dn %in% all_names)
      dn <- unique(gsls[[nm_dn]])
    if (is.null(up) && is.null(dn))
      return(NULL)
    list(base = base, up = up, dn = dn)
  })
  names(pair_list) <- bases
  pair_list <- pair_list[!vapply(pair_list, is.null, logical(1))]
  if (!length(pair_list))
    stop("No usable gene-set bases found.")
  present_genes <- rownames(tpm)
  pair_list <- lapply(pair_list, function(p) {
    if (!is.null(ensembl) & !is.null(p$up)) {
      p$up <- geneAnnPIE$ensemble[match(p$up, geneAnnPIE$gene_name)]
    }
    if (!is.null(ensembl) & !is.null(p$dn)) {
      p$dn <- geneAnnPIE$ensemble[match(p$dn, geneAnnPIE$gene_name)]
    }
    p$up <- if (!is.null(p$up))
      intersect(p$up, present_genes)
    else NULL
    p$dn <- if (!is.null(p$dn))
      intersect(p$dn, present_genes)
    else NULL
    if (is.null(p$up) && is.null(p$dn))
      return(NULL)
    if (!is.null(p$up) && length(p$up) > 1)
      return(p)
    if (!is.null(p$dn) && length(p$dn) > 1)
      return(p)
    NULL
  })
  pair_list <- pair_list[!vapply(pair_list, is.null, logical(1))]
  if (!length(pair_list))
    stop("After trimming to present genes, no sets remain.")
  n_cells <- ncol(tpm)
  sample_id <- colnames(tpm)
  if (is.null(sample_id))
    sample_id <- paste0("Sample", seq_len(n_cells))
  if (is.null(chunk_size) || chunk_size <= 0L || chunk_size >=
      n_cells) {
    starts <- 1L
    ends <- n_cells
  } else {
    starts <- seq.int(1L, n_cells, by = as.integer(chunk_size))
    ends <- pmin(starts + as.integer(chunk_size) - 1L, n_cells)
  }
  if (verbose) {
    message("[PIE] AUCell scoring ", length(pair_list),
            " base set(s) across ", n_cells, " cells in ", length(starts),
            " block(s).")
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
    if (verbose)
      message(sprintf("[PIE] Block %d/%d: cells %d-%d",
                      bi, length(starts), start_i, end_i))
    if (class(tpm)[1] == "data.table") {
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
      if (!is.null(p$up) && length(p$up) > 0)
        auc_geneSets[[paste0(b, up_suffix)]] <- unique(p$up)
      if (!is.null(p$dn) && length(p$dn) > 0)
        auc_geneSets[[paste0(b, dn_suffix)]] <- unique(p$dn)
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
    setcolorder(dt, c("sampleID", comp_cols, keep_up, keep_dn))
    block_tbls[[bi]] <- dt
  }
  scores_all <- data.table::rbindlist(block_tbls, use.names = TRUE,
                                      fill = TRUE)
}
