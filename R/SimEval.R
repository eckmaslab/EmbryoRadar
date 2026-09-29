#' Evaluate EmbryoRadar robustness under target-dataset technical conditions
#'
#' Simulates paired observations from the same latent embryo-class profile:
#' an ideal full-length-like library and a target-like 3'/UMI library. Target-like
#' depth and ambient contamination are sampled jointly from real target cells.
#' The primary outputs quantify paired prediction-score, margin, class-call and
#' k-TSP rule-state stability. This is a robustness analysis around idealised
#' class centroids; it is not validation of biological labels in the target data.
#'
#' @param target_counts Genes x target cells filtered count matrix. If NULL,
#'   cells are identified from raw_counts using emptyDrops.
#' @param raw_counts Optional genes x droplets raw matrix. When supplied with
#'   target_counts, target barcodes are removed before using remaining droplets
#'   as DecontX background.
#' @param auto_filter_from_raw Run emptyDrops when target_counts is NULL.
#' @param target_name Label used in output file names.
#' @param n_per_class Number of simulated cells per class and replicate.
#' @param n_reps Number of paired simulation replicates.
#' @param seed Random seed.
#' @param cal_subsample_n Maximum target cells used for calibration.
#' @param block_size EmbryoRadar prediction block size.
#' @param input_gene_annotations One of gene_name, ensembl, or original.
#' @param out_dir Output directory.
#' @param decontx_z Optional DecontX cluster labels, aligned to target_counts.
#' @param decontx_max_cells Maximum calibration cells passed to DecontX.
#' @param ambient_cap Maximum contamination fraction used in simulation.
#' @param emptydrops_lower,emptydrops_fdr emptyDrops parameters.
#' @param filtered_matrix Informational flag retained for compatibility.
#' @param pseudobulk If TRUE, skip DecontX and set contamination to zero.
#' @param save_all Save the complete result object and per-cell predictions.
#' @param ideal_depth Scalar ideal-library depth.
#' @param model EmbryoRadar model passed to load_models/prediction.
#' @param classes EmbryoRadar classes to evaluate.
#' @return Invisibly, a list containing calibration, paired predictions,
#'   per-cell rule diagnostics, and stability summaries.
SimEval <- function(
    target_counts = NULL,
    raw_counts = NULL,
    auto_filter_from_raw = FALSE,
    target_name = "target",
    n_per_class = 200,
    n_reps = 3,
    seed = 1,
    cal_subsample_n = 10000,
    block_size = 1000,
    input_gene_annotations = c("gene_name", "ensembl", "original"),
    out_dir = "results_exact/",
    decontx_z = NULL,
    decontx_max_cells = 20000,
    ambient_cap = 0.3,
    emptydrops_lower = 500,
    emptydrops_fdr = 0.001,
    filtered_matrix = NA,
    pseudobulk = FALSE,
    save_all = FALSE,
    ideal_depth = 10000,
    model = "trimmed",
    classes = c("Maternal", "ZGA", "Morula", "EPI", "TE", "Adult")) {

  required <- c("Matrix", "dplyr", "ggplot2", "tidyr", "EmbryoRadar")
  if (!isTRUE(pseudobulk)) {
    required <- c(required, "SingleCellExperiment", "SummarizedExperiment", "celda")
  }
  if (isTRUE(auto_filter_from_raw) && is.null(target_counts)) {
    required <- c(required, "DropletUtils")
  }
  required <- unique(required)
  missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
  }

  input_gene_annotations <- match.arg(input_gene_annotations)
  stopifnot(length(n_per_class) == 1L, n_per_class > 0,
            length(n_reps) == 1L, n_reps > 0,
            length(ideal_depth) == 1L, ideal_depth > 0)
  n_per_class <- as.integer(n_per_class)
  n_reps <- as.integer(n_reps)
  ideal_depth <- as.integer(ideal_depth)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  as_dgC <- function(x) {
    if (inherits(x, "dgCMatrix")) return(x)
    if (inherits(x, "sparseMatrix")) return(methods::as(x, "CsparseMatrix"))
    methods::as(Matrix::Matrix(x, sparse = TRUE), "CsparseMatrix")
  }

  check_count_matrix <- function(x, label) {
    if (is.null(rownames(x))) stop(label, " must have gene row names.")
    if (is.null(colnames(x))) stop(label, " must have cell/barcode column names.")
    if (anyNA(x@x) || any(!is.finite(x@x)) || any(x@x < 0)) {
      stop(label, " must contain finite, non-negative counts.")
    }
    x
  }

  # Sum duplicate features without densifying a count matrix.
  collapse_rows_sum <- function(x, ids) {
    keep <- !is.na(ids) & nzchar(ids)
    x <- x[keep, , drop = FALSE]
    ids <- ids[keep]
    lev <- unique(ids)
    if (length(lev) == length(ids)) {
      rownames(x) <- ids
      return(x)
    }
    map <- Matrix::sparseMatrix(
      i = match(ids, lev), j = seq_along(ids), x = 1,
      dims = c(length(lev), length(ids)),
      dimnames = list(lev, NULL)
    )
    methods::as(map %*% x, "CsparseMatrix")
  }

  standardise_count_ids <- function(x) {
    ids <- rownames(x)
    if (input_gene_annotations == "ensembl") ids <- sub("\\..*$", "", ids)
    collapse_rows_sum(x, ids)
  }

  standardise_probability_ids <- function(x, gene_annotation) {
    ids <- rownames(x)
    if (input_gene_annotations == "gene_name") {
      ids <- gene_annotation$gene_name[match(ids, gene_annotation$gene_id)]
    } else if (input_gene_annotations == "ensembl") {
      ids <- sub("\\..*$", "", ids)
    }
    x <- as_dgC(x)
    x <- collapse_rows_sum(x, ids)
    cs <- Matrix::colSums(x)
    if (any(!is.finite(cs) | cs <= 0)) stop("Invalid embryo class probability column.")
    sweep(x, 2, cs, "/")
  }

  log1p_cpm_sparse <- function(counts) {
    lib <- Matrix::colSums(counts)
    lib[!is.finite(lib) | lib <= 0] <- 1
    log1p(sweep(counts, 2, lib / 1e6, "/"))
  }

  safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
      return(NA_real_)
    }
    suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
  }

  data_dir <- system.file("extdata/classifier2.0", package = "EmbryoRadar")
  if (!nzchar(data_dir)) stop("Could not find EmbryoRadar extdata/classifier2.0.")
  gene_annotation <- readRDS(file.path(data_dir, "geneAnnPIE.rds"))
  embryo_prob <- readRDS(file.path(data_dir, "embryo_prob_by_class.rds"))
  embryo_prob <- standardise_probability_ids(embryo_prob, gene_annotation)
  missing_classes <- setdiff(classes, colnames(embryo_prob))
  if (length(missing_classes)) {
    stop("Classes absent from embryo_prob_by_class: ", paste(missing_classes, collapse = ", "))
  }

  model_object <- EmbryoRadar::load_models(
    model = model, input_gene_annotations = input_gene_annotations
  )$model
  missing_model_classes <- setdiff(classes, names(model_object$classifiers))
  if (length(missing_model_classes)) {
    stop("Classes absent from EmbryoRadar model: ", paste(missing_model_classes, collapse = ", "))
  }
  rule_table <- dplyr::bind_rows(lapply(classes, function(cl) {
    z <- model_object$classifiers[[cl]]$TSPs
    data.frame(
      class = cl,
      rule_index = seq_len(nrow(z)),
      gene1 = z[, 1],
      gene2 = z[, 2],
      rule_id = paste(cl, seq_len(nrow(z)), z[, 1], z[, 2], sep = "::"),
      stringsAsFactors = FALSE
    )
  }))

  raw_counts_all <- NULL
  raw_counts_background <- NULL
  target_source <- NULL
  if (!is.null(raw_counts)) {
    raw_counts_all <- check_count_matrix(as_dgC(raw_counts), "raw_counts")
    raw_counts_all <- standardise_count_ids(raw_counts_all)
  }

  if (is.null(target_counts)) {
    if (!isTRUE(auto_filter_from_raw) || is.null(raw_counts_all)) {
      stop("Provide target_counts, or set auto_filter_from_raw=TRUE with raw_counts.")
    }
    message("Running emptyDrops on raw_counts to derive target cells...")
    set.seed(seed)
    ed <- DropletUtils::emptyDrops(raw_counts_all, lower = emptydrops_lower)
    keep <- !is.na(ed$FDR) & ed$FDR <= emptydrops_fdr
    if (!any(keep)) stop("emptyDrops found no cells at the chosen FDR threshold.")
    target_counts <- raw_counts_all[, keep, drop = FALSE]
    raw_counts_background <- raw_counts_all[, !keep, drop = FALSE]
    target_source <- "derived_from_raw_counts_via_emptyDrops"
    filtered_matrix <- FALSE
  } else {
    target_counts <- check_count_matrix(as_dgC(target_counts), "target_counts")
    target_counts <- standardise_count_ids(target_counts)
    target_source <- "user_supplied_target_counts"
    if (is.na(filtered_matrix)) filtered_matrix <- TRUE
    if (!is.null(raw_counts_all)) {
      if (is.null(colnames(target_counts)) || is.null(colnames(raw_counts_all))) {
        stop("Barcode column names are required when raw_counts supplies DecontX background.")
      }
      background_barcodes <- setdiff(colnames(raw_counts_all), colnames(target_counts))
      if (!length(background_barcodes)) {
        stop("No empty-droplet barcodes remain after removing target cell barcodes.")
      }
      raw_counts_background <- raw_counts_all[, background_barcodes, drop = FALSE]
    }
  }
  if (is.na(filtered_matrix)) filtered_matrix <- FALSE

  common_genes <- intersect(rownames(target_counts), rownames(embryo_prob))
  if (length(common_genes) < 5000L) {
    warning("Low target-embryo gene overlap: ", length(common_genes))
  }
  if (length(intersect(unique(c(rule_table$gene1, rule_table$gene2)), common_genes)) < 200L) {
    stop("Too few k-TSP genes overlap the target and embryo probability matrices.")
  }
  message("Gene overlap (target and embryo): ", length(common_genes))

  # Calibration is sampled in target-cell space. Depth and contamination retain
  # the same cell identifiers so they can be sampled as paired technical states.
  set.seed(seed)
  cal_idx <- if (ncol(target_counts) > cal_subsample_n) {
    sample(seq_len(ncol(target_counts)), cal_subsample_n)
  } else seq_len(ncol(target_counts))
  target_cal_full <- target_counts[, cal_idx, drop = FALSE]
  target_cal_common <- target_cal_full[common_genes, , drop = FALSE]
  total_depth <- Matrix::colSums(target_cal_full)
  compatible_depth <- Matrix::colSums(target_cal_common)
  compatible_fraction <- compatible_depth / pmax(total_depth, 1)

  if (!is.null(decontx_z)) {
    if (length(decontx_z) != ncol(target_counts)) {
      stop("Length of decontx_z must equal ncol(target_counts).")
    }
    z_cal <- decontx_z[cal_idx]
  } else z_cal <- NULL

  dx_idx <- seq_len(ncol(target_cal_full))
  if (length(dx_idx) > decontx_max_cells) {
    set.seed(seed + 999L)
    dx_idx <- sample(dx_idx, decontx_max_cells)
  }
  target_dx <- target_cal_full[, dx_idx, drop = FALSE]
  z_dx <- if (is.null(z_cal)) NULL else z_cal[dx_idx]
  contamination <- setNames(rep(NA_real_, ncol(target_cal_full)), colnames(target_cal_full))
  soup_full <- setNames(rep(0, nrow(target_cal_full)), rownames(target_cal_full))

  if (isTRUE(pseudobulk)) {
    message("Skipping DecontX for pseudobulk input; contamination set to zero.")
    contamination[] <- 0
  } else if (!is.null(raw_counts_background)) {
    shared <- intersect(rownames(target_dx), rownames(raw_counts_background))
    if (length(shared) < 1000L) stop("Too few genes shared with DecontX background.")
    x_sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = target_dx[shared, , drop = FALSE])
    )
    bg <- raw_counts_background[shared, , drop = FALSE]
    bg_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = bg))
    message("Running DecontX with empty-droplet background...")
    fit <- celda::decontX(x = x_sce, background = bg_sce, z = z_dx, verbose = FALSE)
    c_dx <- SummarizedExperiment::colData(fit)$decontX_contamination
    contamination[colnames(target_dx)] <- as.numeric(c_dx)
    soup <- Matrix::rowSums(bg)
    if (sum(soup) > 0) soup_full[names(soup)] <- soup / sum(soup)
  } else {
    message("Running DecontX without external background...")
    x_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = target_dx))
    fit <- celda::decontX(x = x_sce, z = z_dx, verbose = FALSE)
    c_dx <- SummarizedExperiment::colData(fit)$decontX_contamination
    contamination[colnames(target_dx)] <- as.numeric(c_dx)
    native <- tryCatch(as_dgC(celda::decontXcounts(fit)), error = function(e) NULL)
    observed <- tryCatch(as_dgC(SummarizedExperiment::assay(fit, "counts")), error = function(e) NULL)
    if (!is.null(native) && !is.null(observed)) {
      contam_counts <- observed - native
      contam_counts@x[contam_counts@x < 0] <- 0
      soup <- Matrix::rowSums(Matrix::drop0(contam_counts))
    } else soup <- Matrix::rowSums(target_dx)
    if (sum(soup) <= 0) soup <- Matrix::rowSums(target_dx)
    soup_full[names(soup)] <- soup / sum(soup)
  }

  contamination <- pmin(pmax(contamination, 0), ambient_cap)
  contam_median <- stats::median(contamination, na.rm = TRUE)
  if (!is.finite(contam_median)) contam_median <- 0

  soup_profile <- soup_full[common_genes]
  soup_profile[!is.finite(soup_profile) | soup_profile < 0] <- 0
  if (sum(soup_profile) <= 0) {
    soup_profile <- Matrix::rowSums(target_cal_common)
  }
  if (sum(soup_profile) <= 0) soup_profile[] <- 1
  soup_profile <- as.numeric(soup_profile / sum(soup_profile))
  names(soup_profile) <- common_genes

  technical_states <- data.frame(
    target_cell = colnames(target_cal_full),
    total_depth = as.numeric(total_depth),
    compatible_depth = as.integer(round(compatible_depth)),
    compatible_fraction = as.numeric(compatible_fraction),
    ambient_fraction = as.numeric(contamination[colnames(target_cal_full)]),
    stringsAsFactors = FALSE
  )
  technical_states <- technical_states[
    is.finite(technical_states$compatible_depth) & technical_states$compatible_depth > 0 &
      is.finite(technical_states$ambient_fraction), , drop = FALSE
  ]
  if (!nrow(technical_states)) stop("No usable target technical states remain.")

  embryo_prob_common <- embryo_prob[common_genes, classes, drop = FALSE]
  embryo_prob_common <- sweep(embryo_prob_common, 2, Matrix::colSums(embryo_prob_common), "/")

  simulate_paired <- function(rep_seed) {
    set.seed(rep_seed)
    n_total <- n_per_class * length(classes)
    truth <- rep(classes, each = n_per_class)
    sample_id <- paste0("rep", rep_seed, "_", truth, "_", ave(seq_len(n_total), truth, FUN = seq_along))
    tech_index <- sample(seq_len(nrow(technical_states)), n_total, replace = TRUE)
    tech <- technical_states[tech_index, , drop = FALSE]
    ideal <- matrix(0L, nrow = length(common_genes), ncol = n_total,
                    dimnames = list(common_genes, sample_id))
    distorted <- ideal
    realised_ambient <- integer(n_total)

    for (j in seq_len(n_total)) {
      # Both observations are drawn independently from the same latent profile;
      # the realised ideal counts are never used as target probabilities.
      latent <- as.numeric(embryo_prob_common[, truth[j]])
      latent <- latent / sum(latent)
      ideal[, j] <- as.integer(stats::rmultinom(1L, ideal_depth, latent)[, 1])

      final_depth <- max(1L, as.integer(tech$compatible_depth[j]))
      ambient_n <- stats::rbinom(1L, size = final_depth, prob = tech$ambient_fraction[j])
      native_n <- final_depth - ambient_n
      native_counts <- if (native_n > 0L) {
        as.integer(stats::rmultinom(1L, native_n, latent)[, 1])
      } else integer(length(latent))
      ambient_counts <- if (ambient_n > 0L) {
        as.integer(stats::rmultinom(1L, ambient_n, soup_profile)[, 1])
      } else integer(length(latent))
      distorted[, j] <- native_counts + ambient_counts
      realised_ambient[j] <- ambient_n
    }
    tech$sample_id <- sample_id
    tech$truth <- truth
    tech$simulated_depth <- Matrix::colSums(distorted)
    tech$ambient_umis <- realised_ambient
    tech$realised_ambient_fraction <- realised_ambient / pmax(tech$simulated_depth, 1)
    list(
      ideal = methods::as(Matrix::Matrix(ideal, sparse = TRUE), "CsparseMatrix"),
      distorted = methods::as(Matrix::Matrix(distorted, sparse = TRUE), "CsparseMatrix"),
      metadata = tech
    )
  }

  ktsp_predict <- function(counts) {
    logmat <- log1p_cpm_sparse(counts)
    suppressMessages(EmbryoRadar::EmbryoRadar_pred_blockwise(
      ExpressionSet = logmat,
      model = model,
      block_size = block_size,
      input_gene_annotations = input_gene_annotations,
      classes = classes
    )$prediction)
  }

  add_prediction_summary <- function(pred) {
    score <- as.matrix(pred[, classes, drop = FALSE])
    storage.mode(score) <- "double"
    ord <- t(apply(score, 1, order, decreasing = TRUE, na.last = TRUE))
    top_i <- ord[, 1]
    second_i <- if (ncol(score) > 1L) ord[, 2] else ord[, 1]
    data.frame(
      sample_id = rownames(pred),
      predicted_top = classes[top_i],
      top_score = score[cbind(seq_len(nrow(score)), top_i)],
      second_score = score[cbind(seq_len(nrow(score)), second_i)],
      prediction_margin = score[cbind(seq_len(nrow(score)), top_i)] -
        score[cbind(seq_len(nrow(score)), second_i)],
      stringsAsFactors = FALSE
    )
  }

  rule_diagnostics_one <- function(counts, rules, representation, rep_id) {
    groups <- c(split(rules, rules$class), list(ALL = rules))
    dplyr::bind_rows(lapply(names(groups), function(cl) {
      rr <- groups[[cl]]
      present <- rr$gene1 %in% rownames(counts) & rr$gene2 %in% rownames(counts)
      rr2 <- rr[present, , drop = FALSE]
      if (!nrow(rr2)) return(NULL)
      a <- as.matrix(counts[rr2$gene1, , drop = FALSE])
      b <- as.matrix(counts[rr2$gene2, , drop = FALSE])
      both_detected <- a > 0 & b > 0
      both_zero <- a == 0 & b == 0
      exactly_one <- xor(a > 0, b > 0)
      tied_positive <- a == b & a > 0
      evaluable <- a != b
      denom <- nrow(rr2)
      data.frame(
        sample_id = colnames(counts), class = cl, representation = representation,
        rep = rep_id, n_rules_total = nrow(rr), n_rules_present = denom,
        fraction_rules_missing = 1 - denom / nrow(rr),
        n_both_detected = colSums(both_detected),
        fraction_both_detected = colSums(both_detected) / denom,
        n_both_zero = colSums(both_zero),
        fraction_both_zero = colSums(both_zero) / denom,
        n_exactly_one_detected = colSums(exactly_one),
        fraction_exactly_one_detected = colSums(exactly_one) / denom,
        n_tied_positive = colSums(tied_positive),
        fraction_tied_positive = colSums(tied_positive) / denom,
        n_evaluable_rules = colSums(evaluable),
        fraction_evaluable_rules = colSums(evaluable) / denom,
        stringsAsFactors = FALSE
      )
    }))
  }

  paired_rule_stability <- function(ideal, distorted, rules, rep_id) {
    groups <- c(split(rules, rules$class), list(ALL = rules))
    dplyr::bind_rows(lapply(names(groups), function(cl) {
      rr <- groups[[cl]]
      present <- rr$gene1 %in% rownames(ideal) & rr$gene2 %in% rownames(ideal) &
        rr$gene1 %in% rownames(distorted) & rr$gene2 %in% rownames(distorted)
      rr <- rr[present, , drop = FALSE]
      if (!nrow(rr)) return(NULL)
      ai <- as.matrix(ideal[rr$gene1, , drop = FALSE])
      bi <- as.matrix(ideal[rr$gene2, , drop = FALSE])
      ad <- as.matrix(distorted[rr$gene1, , drop = FALSE])
      bd <- as.matrix(distorted[rr$gene2, , drop = FALSE])
      si <- sign(ai - bi)
      sd <- sign(ad - bd)
      both_directional <- si != 0 & sd != 0
      directional_den <- colSums(both_directional)
      data.frame(
        sample_id = colnames(ideal), class = cl, rep = rep_id,
        n_rules_compared = nrow(rr),
        rule_state_concordance = colMeans(si == sd),
        n_directional_in_both = directional_den,
        directional_concordance = ifelse(
          directional_den > 0,
          colSums(si == sd & both_directional) / directional_den,
          NA_real_
        ),
        directional_flip_rate = ifelse(
          directional_den > 0,
          colSums(si == -sd & both_directional) / directional_den,
          NA_real_
        ),
        tie_introduced_rate = colMeans(si != 0 & sd == 0),
        tie_resolved_rate = colMeans(si == 0 & sd != 0),
        stringsAsFactors = FALSE
      )
    }))
  }

  prediction_list <- vector("list", n_reps)
  raw_prediction_list <- if (isTRUE(save_all)) vector("list", n_reps) else NULL
  technical_list <- vector("list", n_reps)
  rule_qc_list <- vector("list", n_reps)
  rule_stability_list <- vector("list", n_reps)

  for (r in seq_len(n_reps)) {
    current_seed <- seed + (r - 1L) * 1000L
    message("Running paired replicate ", r, "/", n_reps, " (seed=", current_seed, ")")
    sim <- simulate_paired(current_seed + 777L)
    pred_ideal <- ktsp_predict(sim$ideal)
    pred_distorted <- ktsp_predict(sim$distorted)
    if (!identical(rownames(pred_ideal), rownames(pred_distorted))) {
      stop("Ideal and distorted predictions are not aligned by sample ID.")
    }

    s_ideal <- add_prediction_summary(pred_ideal)
    s_dist <- add_prediction_summary(pred_distorted)
    paired <- data.frame(
      sample_id = rownames(pred_ideal), rep = r,
      truth = sim$metadata$truth[match(rownames(pred_ideal), sim$metadata$sample_id)],
      ideal_top = s_ideal$predicted_top,
      distorted_top = s_dist$predicted_top,
      class_agreement = s_ideal$predicted_top == s_dist$predicted_top,
      ideal_threshold_call = if ("prediction" %in% colnames(pred_ideal)) pred_ideal$prediction else NA_character_,
      distorted_threshold_call = if ("prediction" %in% colnames(pred_distorted)) pred_distorted$prediction else NA_character_,
      ideal_margin = s_ideal$prediction_margin,
      distorted_margin = s_dist$prediction_margin,
      margin_delta = s_dist$prediction_margin - s_ideal$prediction_margin,
      stringsAsFactors = FALSE
    )
    paired$threshold_call_agreement <-
      paired$ideal_threshold_call == paired$distorted_threshold_call
    tech_match <- match(paired$sample_id, sim$metadata$sample_id)
    paired$total_depth <- sim$metadata$total_depth[tech_match]
    paired$compatible_depth <- sim$metadata$compatible_depth[tech_match]
    paired$ambient_fraction <- sim$metadata$ambient_fraction[tech_match]
    paired$realised_ambient_fraction <- sim$metadata$realised_ambient_fraction[tech_match]
    for (cl in classes) {
      paired[[paste0(cl, "_ideal")]] <- pred_ideal[[cl]]
      paired[[paste0(cl, "_distorted")]] <- pred_distorted[[cl]]
      paired[[paste0(cl, "_delta")]] <- pred_distorted[[cl]] - pred_ideal[[cl]]
    }

    technical_list[[r]] <- dplyr::mutate(sim$metadata, rep = r)
    prediction_list[[r]] <- paired
    rule_qc_list[[r]] <- dplyr::bind_rows(
      rule_diagnostics_one(sim$ideal, rule_table, "ideal", r),
      rule_diagnostics_one(sim$distorted, rule_table, "distorted", r)
    )
    rule_stability_list[[r]] <- paired_rule_stability(sim$ideal, sim$distorted, rule_table, r)

    if (isTRUE(save_all)) {
      raw_prediction_list[[r]] <- list(ideal = pred_ideal, distorted = pred_distorted)
    }
  }

  paired_predictions <- dplyr::bind_rows(prediction_list)
  simulated_technical_states <- dplyr::bind_rows(technical_list)
  rule_diagnostics <- dplyr::bind_rows(rule_qc_list)
  rule_stability <- dplyr::bind_rows(rule_stability_list)

  score_stability <- dplyr::bind_rows(lapply(classes, function(cl) {
    x <- paired_predictions[[paste0(cl, "_ideal")]]
    y <- paired_predictions[[paste0(cl, "_distorted")]]
    data.frame(
      class = cl,
      spearman_rho = safe_cor(x, y),
      mean_absolute_score_delta = mean(abs(y - x), na.rm = TRUE),
      median_score_delta = stats::median(y - x, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  prediction_stability <- data.frame(
    target_name = target_name,
    n_simulated_pairs = nrow(paired_predictions),
    class_agreement = mean(paired_predictions$class_agreement, na.rm = TRUE),
    threshold_call_agreement = mean(
      paired_predictions$threshold_call_agreement, na.rm = TRUE
    ),
    mean_absolute_margin_delta = mean(abs(paired_predictions$margin_delta), na.rm = TRUE),
    median_margin_delta = stats::median(paired_predictions$margin_delta, na.rm = TRUE),
    margin_spearman_rho = safe_cor(
      paired_predictions$ideal_margin, paired_predictions$distorted_margin
    ),
    stringsAsFactors = FALSE
  )
  all_rule_stability <- rule_stability |>
    dplyr::filter(class == "ALL") |>
    dplyr::select(sample_id, rep, rule_state_concordance,
                  directional_concordance, tie_introduced_rate)
  technical_outcomes <- paired_predictions |>
    dplyr::left_join(all_rule_stability, by = c("sample_id", "rep")) |>
    dplyr::mutate(
      absolute_margin_delta = abs(margin_delta),
      class_switch = as.numeric(!class_agreement),
      threshold_call_switch = as.numeric(!threshold_call_agreement),
      log10_compatible_depth = log10(pmax(compatible_depth, 1))
    )
  technical_associations <- dplyr::bind_rows(lapply(
    c("absolute_margin_delta", "class_switch", "threshold_call_switch", "rule_state_concordance",
      "directional_concordance", "tie_introduced_rate"),
    function(outcome) {
      dplyr::bind_rows(lapply(
        c("log10_compatible_depth", "ambient_fraction"),
        function(covariate) data.frame(
          outcome = outcome,
          technical_covariate = covariate,
          spearman_rho = safe_cor(
            technical_outcomes[[outcome]], technical_outcomes[[covariate]]
          ),
          stringsAsFactors = FALSE
        )
      ))
    }
  ))
  rule_stability_summary <- rule_stability |>
    dplyr::group_by(class) |>
    dplyr::summarise(
      mean_rule_state_concordance = mean(rule_state_concordance, na.rm = TRUE),
      mean_directional_concordance = mean(directional_concordance, na.rm = TRUE),
      mean_directional_flip_rate = mean(directional_flip_rate, na.rm = TRUE),
      mean_tie_introduced_rate = mean(tie_introduced_rate, na.rm = TRUE),
      mean_tie_resolved_rate = mean(tie_resolved_rate, na.rm = TRUE),
      .groups = "drop"
    )
  rule_diagnostics_summary <- rule_diagnostics |>
    dplyr::group_by(representation, class) |>
    dplyr::summarise(
      mean_fraction_both_detected = mean(fraction_both_detected, na.rm = TRUE),
      mean_fraction_both_zero = mean(fraction_both_zero, na.rm = TRUE),
      mean_fraction_exactly_one_detected = mean(fraction_exactly_one_detected, na.rm = TRUE),
      mean_fraction_evaluable = mean(fraction_evaluable_rules, na.rm = TRUE),
      median_n_evaluable = stats::median(n_evaluable_rules, na.rm = TRUE),
      .groups = "drop"
    )

  stability_plot <- ggplot2::ggplot(
    score_stability,
    ggplot2::aes(x = stats::reorder(class, spearman_rho), y = spearman_rho, fill = class)
  ) +
    ggplot2::geom_col(colour = "grey30", linewidth = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::labs(
      title = paste0("Paired EmbryoRadar score stability - ", target_name),
      subtitle = paste0("Same latent profiles; ideal versus target-like | n = ",
                        nrow(paired_predictions), " paired cells"),
      x = NULL, y = "Spearman correlation"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank())

  rule_plot_data <- rule_diagnostics_summary |>
    dplyr::filter(class == "ALL") |>
    tidyr::pivot_longer(
      cols = c(mean_fraction_both_detected, mean_fraction_both_zero,
               mean_fraction_exactly_one_detected, mean_fraction_evaluable),
      names_to = "metric", values_to = "value"
    )
  rule_plot <- ggplot2::ggplot(
    rule_plot_data, ggplot2::aes(representation, value, fill = representation)
  ) +
    ggplot2::geom_col(colour = "grey30", linewidth = 0.25) +
    ggplot2::facet_wrap(~metric, nrow = 1) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = paste0("k-TSP rule evaluability - ", target_name),
      x = NULL, y = "Mean fraction across cells"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

  file_prefix <- gsub("[^A-Za-z0-9_.-]+", "_", target_name)
  paths <- list(
    paired_predictions = file.path(out_dir, paste0("paired_predictions_", file_prefix, ".csv")),
    score_stability = file.path(out_dir, paste0("score_stability_", file_prefix, ".csv")),
    prediction_stability = file.path(out_dir, paste0("prediction_stability_", file_prefix, ".csv")),
    technical_associations = file.path(out_dir, paste0("technical_associations_", file_prefix, ".csv")),
    rule_diagnostics = file.path(out_dir, paste0("rule_diagnostics_per_cell_", file_prefix, ".csv")),
    rule_stability = file.path(out_dir, paste0("rule_stability_per_cell_", file_prefix, ".csv")),
    rule_stability_summary = file.path(out_dir, paste0("rule_stability_summary_", file_prefix, ".csv")),
    rule_diagnostics_summary = file.path(out_dir, paste0("rule_diagnostics_summary_", file_prefix, ".csv")),
    technical_states = file.path(out_dir, paste0("simulated_technical_states_", file_prefix, ".csv")),
    stability_plot = file.path(out_dir, paste0("paired_score_stability_", file_prefix, ".pdf")),
    rule_plot = file.path(out_dir, paste0("rule_evaluability_", file_prefix, ".pdf")),
    result_rds = file.path(out_dir, paste0("sim_eval_results_", file_prefix, ".rds"))
  )
  utils::write.csv(paired_predictions, paths$paired_predictions, row.names = FALSE)
  utils::write.csv(score_stability, paths$score_stability, row.names = FALSE)
  utils::write.csv(prediction_stability, paths$prediction_stability, row.names = FALSE)
  utils::write.csv(technical_associations, paths$technical_associations, row.names = FALSE)
  utils::write.csv(rule_diagnostics, paths$rule_diagnostics, row.names = FALSE)
  utils::write.csv(rule_stability, paths$rule_stability, row.names = FALSE)
  utils::write.csv(rule_stability_summary, paths$rule_stability_summary, row.names = FALSE)
  utils::write.csv(rule_diagnostics_summary, paths$rule_diagnostics_summary, row.names = FALSE)
  utils::write.csv(simulated_technical_states, paths$technical_states, row.names = FALSE)
  ggplot2::ggsave(paths$stability_plot, stability_plot, width = 5.5, height = 4.2)
  ggplot2::ggsave(paths$rule_plot, rule_plot, width = 9, height = 3.8)

  result <- list(
    target_name = target_name,
    target_source = target_source,
    filtered_matrix = filtered_matrix,
    pseudobulk = pseudobulk,
    common_genes = common_genes,
    calibration = list(
      technical_states = technical_states,
      ambient_cap = ambient_cap,
      ambient_median = contam_median,
      soup_profile = soup_profile
    ),
    paired_predictions = paired_predictions,
    raw_predictions = raw_prediction_list,
    score_stability = score_stability,
    prediction_stability = prediction_stability,
    technical_associations = technical_associations,
    rule_diagnostics = rule_diagnostics,
    rule_diagnostics_summary = rule_diagnostics_summary,
    rule_stability = rule_stability,
    rule_stability_summary = rule_stability_summary,
    simulated_technical_states = simulated_technical_states,
    plots = list(score_stability = stability_plot, rule_evaluability = rule_plot),
    paths = paths
  )
  if (isTRUE(save_all)) saveRDS(result, paths$result_rds)

  message("Saved paired stability outputs to: ", normalizePath(out_dir, mustWork = FALSE))
  invisible(result)
}
