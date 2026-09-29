#' Load packaged gene sets used for transcriptional and developmental scoring
#'
#' @param id_type Character scalar. Which identifier to return in the gene sets.
#'   One of: "gene_name" (HGNC symbols) or "ensembl" (Ensembl gene IDs without version).
#' @param keep_unmapped Logical. If FALSE (default), drop genes that cannot be mapped
#'   to the requested id_type. If TRUE, keep as-is (not recommended for scoring).
#'
#' @return A named list of gene sets in the requested identifier type.
#'
#' @export
load_Maternal_gene_sets <- function(
    id_type = c("gene_name", "ensembl"),
    keep_unmapped = F
) {
  library(data.table); library(dplyr); library(readxl); library(magrittr)

  id_type <- match.arg(id_type)

  data_dir <- system.file("extdata/classifier2.0", package = "EmbryoRadar")
  geneAnnPIE <- readRDS(file.path(data_dir, "geneAnnPIE.rds"))
  geneAnnPIE$ensemble <- gsub("\\..*", "", geneAnnPIE$gene_id)

  # ---- mapping helpers ----
  ens_from_ens <- function(x) gsub("\\..*$", "", x)
  ens_from_sym <- function(x) {
    # map gene_name -> ensembl (drop versions)
    idx <- match(x, geneAnnPIE$gene_name)
    out <- geneAnnPIE$ensemble[idx]
    ens_from_ens(out)
  }
  sym_from_ens <- function(x) {
    x0 <- ens_from_ens(x)
    idx <- match(x0, geneAnnPIE$ensemble)
    geneAnnPIE$gene_name[idx]
  }
  coerce_ids <- function(genes, from = c("gene_name", "ensembl"), to = c("gene_name", "ensembl")) {
    from <- match.arg(from); to <- match.arg(to)
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & genes != ""]
    if (from == to) return(genes)

    if (from == "gene_name" && to == "ensembl") {
      mapped <- ens_from_sym(genes)
    } else if (from == "ensembl" && to == "gene_name") {
      mapped <- sym_from_ens(genes)
    } else {
      stop("Unsupported mapping.")
    }

    if (isTRUE(keep_unmapped)) {
      # keep original where mapping failed
      mapped2 <- mapped
      mapped2[is.na(mapped2)] <- genes[is.na(mapped2)]
      mapped <- mapped2
    } else {
      mapped <- mapped[!is.na(mapped)]
    }
    unique(mapped)
  }

  # ---- load inputs (mostly gene symbols unless explicitly ensembl) ----
  Mat_sc <- readRDS(file.path(data_dir, "top_table_deglist_Maternal_KIRC_as_num_trimmed_200_rules.rds"))
  tempclus <- readRDS(file.path(data_dir, "markers_up_dn_separated.rds"))
  gsls <- list()
  gsls$Mat_sc_up <- Mat_sc %>% filter(FDR < 5e-11, logFC > 10) %>% pull(gene_name)
  gsls$Mat_sc_dn <- Mat_sc %>% filter(FDR < 5e-11, logFC < -10) %>% pull(gene_name)

  # developmental sets: tempclus stores ensembl; currently mapped to gene_name
  gsls$mat_dev_up <- sym_from_ens(tempclus$Maternal$up)
  gsls$mat_dev_dn <- sym_from_ens(tempclus$Maternal$dn)
  gsls$zga_dev_up <- sym_from_ens(tempclus$ZGA$up)
  gsls$zga_dev_dn <- sym_from_ens(tempclus$ZGA$dn)
  gsls$mor_dev_up <- sym_from_ens(tempclus$Morula$up)
  gsls$mor_dev_dn <- sym_from_ens(tempclus$Morula$dn)
  gsls$icm_dev_up <- sym_from_ens(tempclus$ICM$up)
  gsls$icm_dev_dn <- sym_from_ens(tempclus$ICM$dn)
  gsls$TE_dev_up  <- sym_from_ens(tempclus$TE$up)
  gsls$TE_dev_dn  <- sym_from_ens(tempclus$TE$dn)

  # ---- infer the CURRENT id type in each set, then convert if requested ----
  # Known sources:
  # - tempclus/rambow are converted to gene_name above
  # - most others are gene symbols already
  # If id_type == "ensembl", map gene_name -> ensembl for all sets.
  if (id_type == "ensembl") {
    gsls <- lapply(gsls, function(v) coerce_ids(v, from = "gene_name", to = "ensembl"))
  } else {
    gsls <- lapply(gsls, function(v) coerce_ids(v, from = "gene_name", to = "gene_name"))
  }

  # Drop empties after mapping unless keep_unmapped=TRUE
  if (!isTRUE(keep_unmapped)) {
    gsls <- gsls[vapply(gsls, function(x) length(x) > 0, logical(1))]
  }

  return(gsls)
}
