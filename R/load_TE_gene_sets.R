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
load_TE_gene_sets <- function(
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
  cor_res_expr <- fread(file.path(data_dir, "cor_TE_genes_tcga_skcm.csv"))
  TE_sc <- readRDS(file.path(data_dir, "top_table_deglist_te_skcm_as_num_trimmed_200_rules.rds"))
  TE_num_skcm_200 <- readRDS(file.path(data_dir, "top_table_deglist_skcm_te_num_200_rules.rds"))
  TE_num_uvm <- readRDS(file.path(data_dir, "TT_TCGA_UVM.rds"))
  TE_num_uvm_up <- TE_num_uvm %>% filter(FDR < .00005, logFC > 0) %>% pull(gene_name)
  TE_num_uvm_dn <- TE_num_uvm %>% filter(FDR < .00005, logFC < 0) %>% pull(gene_name)
  TE_num_tcga <- readRDS(file.path(data_dir, "top_table_deglist_te_skcm_all_contrast_tcga.rds"))
  TE_num_ccle <- readRDS(file.path(data_dir, "top_table_deglist_te_skcm_all_contrast_ccle.rds"))

  tempclus <- readRDS(file.path(data_dir, "markers_up_dn_separated.rds"))
  bajpai_pigmentation_regulator <- readxl::read_xlsx(file.path(data_dir, "science.ade6289_table_s5.xlsx"))
  bajpai_pigmentation_regulator_up <- bajpai_pigmentation_regulator %>%
    filter(q_value < 0.05, Spearman_correlation_coefficient < 0) %>% pull(Screen_hit)
  bajpai_pigmentation_regulator_dn <- bajpai_pigmentation_regulator %>%
    filter(q_value < 0.05, Spearman_correlation_coefficient > 0) %>% pull(Screen_hit)

  tretinoin <- readRDS(file.path(data_dir, "tretinoin_gs.rds"))
  reactome_keratinization <- readRDS(file.path(data_dir, "reactome_keratinization.rds"))

  Zhu_gs <- readxl::read_xlsx(file.path(data_dir, "1-s2.0-S1535610818301223-mmc4.xlsx"), skip = 1)
  Zhu_gs <- split(Zhu_gs$Gene, Zhu_gs$Signature)
  names(Zhu_gs) <- paste0("Zhu_", gsub(" |-","_", names(Zhu_gs)))

  pozniak_genesets <- readxl::read_xlsx(file.path(data_dir, "mmc3-2.xlsx"), sheet = 3)
  pozniak_genesets <- split(pozniak_genesets$gene, pozniak_genesets$cluster)
  names(pozniak_genesets) <- paste0("Pozniak_", gsub(" ","_", gsub("\\(|\\)","", names(pozniak_genesets))))

  gsls <- list()
  gsls <- c(gsls, pozniak_genesets, Zhu_gs)
  gsls$reactome_keratinization <- reactome_keratinization$reactome_keratinization
  gsls$martens_tretinoin_response_up <- tretinoin$MARTENS_TRETINOIN_RESPONSE_UP
  gsls$martens_tretinoin_response_dn <- tretinoin$MARTENS_TRETINOIN_RESPONSE_DN
  gsls$bajpai_pigmentation_regulator_up <- bajpai_pigmentation_regulator_up
  gsls$bajpai_pigmentation_regulator_dn <- bajpai_pigmentation_regulator_dn

  gsls$TE_sc_up <- TE_sc %>% filter(FDR < 5e-30, logFC > 0) %>% pull(gene_name)
  gsls$TE_sc_dn <- TE_sc %>% filter(FDR < 5e-30, logFC < 0) %>% pull(gene_name)

  gsls$TE_cor_tcga_up <- cor_res_expr %>% filter(p < .000001, r > 0) %>% .$gene %>% unique
  gsls$TE_cor_tcga_dn <- cor_res_expr %>% filter(p < .000001, r < 0) %>% .$gene %>% unique

  gsls$TE_kmcut_up <- TE_num_tcga %>% filter(comparison == "highvsmidlow", FDR < 1e-10, logFC > 0) %>% .$gene_name %>% unique
  gsls$TE_kmcut_dn <- TE_num_tcga %>% filter(comparison == "highvsmidlow", FDR < 1e-10, logFC < 0) %>% .$gene_name %>% unique

  gsls$TE_num_tcga_up <- TE_num_tcga %>% filter(comparison == "TE_num", FDR < 1e-10, logFC > 10) %>% .$gene_name %>% unique
  gsls$TE_num_tcga_dn <- TE_num_tcga %>% filter(comparison == "TE_num", FDR < 1e-10, logFC < -10) %>% .$gene_name %>% unique

  gsls$TE_num_ccle_up <- TE_num_ccle %>% filter(comparison == "TE_num", FDR < 1e-10, logFC > 10) %>% .$gene_name %>% unique
  gsls$TE_num_ccle_dn <- TE_num_ccle %>% filter(comparison == "TE_num", FDR < 1e-10, logFC < -10) %>% .$gene_name %>% unique

  gsls$TE_num_tcgaxccle_up <- intersect(
    TE_num_tcga %>% filter(comparison == "TE_num", FDR < .05, logFC > 2) %>% .$gene_name %>% unique,
    TE_num_ccle %>% filter(comparison == "TE_num", FDR < .05, logFC > 2) %>% .$gene_name %>% unique
  )
  gsls$TE_num_tcgaxccle_dn <- intersect(
    TE_num_tcga %>% filter(comparison == "TE_num", FDR < .05, logFC < -2) %>% .$gene_name %>% unique,
    TE_num_ccle %>% filter(comparison == "TE_num", FDR < .05, logFC < -2) %>% .$gene_name %>% unique
  )

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

  rambow_gs <- readxl::read_xlsx(file.path(data_dir, "1-s2.0-S0092867418307931-mmc1.xlsx"), sheet = 4, skip = 1)
  colnames(rambow_gs) <- paste0("Rambow_", gsub("-","_", gsub(" \\(.*","", colnames(rambow_gs))))
  gsls <- c(gsls, setNames(lapply(seq_len(ncol(rambow_gs)), function(j) {
    unique(sym_from_ens(na.omit(rambow_gs[[j]])))
  }), nm = gsub("\\(|\\)","", colnames(rambow_gs))))

  ICR_res <- suppressMessages(readxl::read_xlsx(file.path(data_dir, "1-s2.0-S0092867418311784-mmc4.xlsx"),
                                                sheet = 1, skip = 2) %>%
                                set_colnames(rep("gene", ncol(.))))
  ICR_res <- ICR_res[,1] %>% mutate(geneset = "ICI_resistance_up") %>%
    rbind(ICR_res[,9] %>% mutate(geneset = "ICI_resistance_dn"))
  ICR_res <- split(ICR_res$gene, ICR_res$geneset)
  gsls$ICR_resistance_up <- na.omit(ICR_res$ICI_resistance_up)
  gsls$ICR_resistance_dn <- na.omit(ICR_res$ICI_resistance_dn)

  e2f <- readRDS(file.path(data_dir, "E2F_target_gs.rds"))
  gsls$E2F1_targets_up <- unlist(e2f$e2f1$E2F1_UP.V1_UP)
  gsls$E2F1_targets_dn <- unlist(e2f$e2f1$E2F1_UP.V1_DN)
  gsls$E2F2_targets    <- unlist(e2f$e2f2$E2F2_TARGET_GENES)
  gsls$E2F3_targets_up <- unlist(e2f$e2f3$E2F3_UP.V1_UP)
  gsls$E2F3_targets_dn <- unlist(e2f$e2f3$E2F3_UP.V1_DN)

  rb <- readRDS(file.path(data_dir, "rb_target_gs.rds"))
  gsls$rb_targets_up <- unlist(rb$rb$VERNELL_RETINOBLASTOMA_PATHWAY_UP)
  gsls$rb_targets_dn <- unlist(rb$rb$VERNELL_RETINOBLASTOMA_PATHWAY_DN)

  gsls$TE_num_uvm_up <- TE_num_uvm_up
  gsls$TE_num_uvm_dn <- TE_num_uvm_dn

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
