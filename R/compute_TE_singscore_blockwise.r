compute_TE_singscores_blockwise <- function(
    tpm = NULL, gene_names = NULL, ensembl = NULL, raw_counts = NULL,
    additional_gene_sets = NULL, additional_gene_sets_only = FALSE,
    chunk_size = 4000, verbose = TRUE,
    up_suffix = "_up", dn_suffix = "_dn"   # configurable suffixes
) {
  # deps
  library(data.table)
  library(dplyr)
  library(magrittr)
  library(singscore)
  library(readxl)

  # --------------------------
  # 0) Load resources / gene sets
  # --------------------------
  data_dir <- system.file("extdata/classifier2.0", package = "EmbryoRadar")
  geneAnnPIE <- readRDS(file.path(data_dir, "geneAnnPIE.rds"))
  geneAnnPIE$ensemble <- gsub("\\..*", "", geneAnnPIE$gene_id)

  gsls <- if (isTRUE(additional_gene_sets_only)) {
    c(additional_gene_sets)
  } else {
    c(additional_gene_sets, load_TE_gene_sets())
  }

  # --------------------------
  # 1) TPM derivation if raw counts
  # --------------------------
  if (!is.null(raw_counts)) {
    if (verbose) message("[PIE] Converting raw counts to TPM ...")
    if (!is.null(gene_names)) {
      gene_lengths <- geneAnnPIE$width[match(gene_names, geneAnnPIE$gene_name)]
    } else if (!is.null(ensembl)) {
      gene_lengths <- geneAnnPIE$width[match(ensembl, geneAnnPIE$ensemble)]
    } else {
      stop("Provide either gene_names or ensembl with raw_counts.")
    }
    tpm <- calculate_tpm(raw_counts, gene_lengths)
  }
  if (is.null(tpm)) stop("Provide either `tpm` or `raw_counts`.")

  # --------------------------
  # 2) Row names (gene ids used in scoring)
  # --------------------------
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

  # Drop too-small sets (<=1 gene)
  gsid_ok <- vapply(gsls, function(x) length(x) > 0, logical(1))
  if (any(!gsid_ok) && verbose) {
    message("[PIE] Dropping tiny gene sets: ", paste0(names(gsls)[!gsid_ok], collapse = ", "))
  }
  gsls <- gsls[gsid_ok]

  # --------------------------
  # 3) Build base list: handle pairs and UNSUFFIXED (up-only) gracefully
  # --------------------------
  all_names  <- names(gsls)
  base_from  <- function(nm) sub("_(up|dn)$", "", nm)
  bases      <- unique(base_from(all_names))

  pair_list <- lapply(bases, function(base) {
    nm_up  <- paste0(base, up_suffix)
    nm_dn  <- paste0(base, dn_suffix)
    nm_raw <- base                       # UNSUFFIXED

    up <- NULL; dn <- NULL
    if (nm_up %in% all_names)  up <- unique(gsls[[nm_up]])
    if (nm_raw %in% all_names) up <- unique(c(up, gsls[[nm_raw]]))  # merge unsuffixed into UP
    if (nm_dn %in% all_names)  dn <- unique(gsls[[nm_dn]])

    if (is.null(up) && is.null(dn)) return(NULL)
    list(base = base, up = up, dn = dn)
  })
  names(pair_list) <- bases
  pair_list <- pair_list[!vapply(pair_list, is.null, logical(1))]
  if (!length(pair_list)) stop("No usable gene-set bases found.")

  # --------------------------
  # 4) Trim sets to present genes; drop empty/tiny after trim
  # --------------------------
  present_genes <- rownames(tpm)
  pair_list <- lapply(pair_list, function(p) {
    p$up <- if (!is.null(p$up)) intersect(p$up, present_genes) else NULL
    p$dn <- if (!is.null(p$dn)) intersect(p$dn, present_genes) else NULL
    if (is.null(p$up) && is.null(p$dn)) return(NULL)
    # At least one side must have >1 gene
    if (!is.null(p$up) && length(p$up) > 1) return(p)
    if (!is.null(p$dn) && length(p$dn) > 1) return(p)
    NULL
  })
  pair_list <- pair_list[!vapply(pair_list, is.null, logical(1))]
  if (!length(pair_list)) stop("After trimming to present genes, no sets remain.")

  # --------------------------
  # 5) Block-wise ranking + scoring
  # --------------------------
  n_cells   <- ncol(tpm)
  sample_id <- colnames(tpm)
  if (is.null(sample_id)) sample_id <- paste0("Sample", seq_len(n_cells))

  if (is.null(chunk_size) || chunk_size <= 0L || chunk_size >= n_cells) {
    starts <- 1L; ends <- n_cells
  } else {
    starts <- seq.int(1L, n_cells, by = as.integer(chunk_size))
    ends   <- pmin(starts + as.integer(chunk_size) - 1L, n_cells)
  }

  if (verbose) {
    message("[PIE] Scoring ", length(pair_list), " base set(s) across ",
            n_cells, " cells in ", length(starts), " block(s).")
  }

  block_tbls <- vector("list", length(starts))

  for (bi in seq_along(starts)) {
    start_i <- starts[bi]; end_i <- ends[bi]
    idx <- start_i:end_i
    if (verbose) message(sprintf("[PIE] Block %d/%d: cells %d-%d", bi, length(starts), start_i, end_i))

    tpm_block <- as.matrix(tpm[, idx, drop = FALSE])
    colnames(tpm_block) <- sample_id[idx]

    rank_block <- singscore::rankGenes(tpm_block)
    rownames(rank_block) <- rownames(tpm_block)

    block_scores <- lapply(pair_list, function(p) {
      if (!is.null(p$up) && !is.null(p$dn)) {
        ss <- singscore::simpleScore(rank_block, upSet = p$up, downSet = p$dn)
      } else if (!is.null(p$up)) {
        ss <- singscore::simpleScore(rank_block, upSet = p$up)       # UNSUFFIXED handled here
      } else {
        ss <- singscore::simpleScore(rank_block, downSet = p$dn)     # dn-only (rare)
        ss$TotalScore <- -ss$TotalScore                              # flip so higher = "more base"
      }
      unname(ss$TotalScore)
    })

    block_dt <- data.table::as.data.table(block_scores)
    block_dt$sampleID <- colnames(tpm_block)
    setcolorder(block_dt, c("sampleID", vapply(pair_list, `[[`, character(1), "base")))

    block_tbls[[bi]] <- block_dt
  }

  scc_all <- data.table::rbindlist(block_tbls, use.names = TRUE, fill = TRUE)

  list(
    singscores = scc_all,   # data.table: sampleID + one column per base name
    geneSets   = gsls       # the (trimmed) raw gene-set list used
  )
}
