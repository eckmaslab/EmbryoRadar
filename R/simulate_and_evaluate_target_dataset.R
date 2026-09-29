simulate_and_evaluate_target_dataset <- function(
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
    plot_width = 4.5,
    plot_height = 3.5,
    decontx_z = NULL,
    decontx_max_cells = 20000,
    ambient_cap = 0.30,
    emptydrops_lower = 500,
    emptydrops_fdr = 0.001,
    filtered_matrix = NA,
    pseudobulk = FALSE,
    save_all = FALSE,
    ideal_depth = 10000
) {
  suppressPackageStartupMessages({
    library(Matrix)
    library(dplyr)
    library(ggplot2)
    library(EmbryoRadar)
    library(tidyr)
    library(pROC)
    library(ComplexHeatmap)
    library(viridis)
    library(grid)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(celda)
    library(DropletUtils)
  })

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  input_gene_annotations <- match.arg(input_gene_annotations)

  ## ============================================================
  ## Helpers
  ## ============================================================
  confusionHeatmap_local <- function(tb) {
    hm <- prop.table(tb, margin = 2)
    Heatmap(
      hm,
      name = "proportion",
      col = viridis::viridis(100),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_side = "left",
      row_title = "Prediction",
      column_title = "Reference",
      column_title_side = "bottom",
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.0f", tb[i, j]), x, y,
                  gp = gpar(fontsize = 10, col = "white"))
      }
    )
  }

  as_dgC <- function(x) {
    if (inherits(x, "dgCMatrix")) return(x)
    Matrix(as.matrix(x), sparse = TRUE)
  }

  log1p_cpm_sparse <- function(counts) {
    lib <- Matrix::colSums(counts)
    lib[lib == 0] <- 1
    log1p(1e6 * sweep(counts, 2, lib, "/"))
  }

  rowVars_sparse <- function(m) {
    mu <- Matrix::rowMeans(m)
    mu2 <- Matrix::rowMeans(m^2)
    as.numeric(mu2 - mu^2)
  }

  fbeta <- function(prec, rec, beta = 1) {
    b2 <- beta^2
    if (is.na(prec) || is.na(rec) || (b2 * prec + rec) == 0) return(NA_real_)
    (1 + b2) * prec * rec / (b2 * prec + rec)
  }

  compute_prediction_spearman <- function(pred_ref, pred, class_levels) {
    common_classes <- intersect(class_levels, intersect(colnames(pred_ref), colnames(pred)))

    out <- lapply(common_classes, function(cl) {
      ok <- is.finite(pred_ref[[cl]]) & is.finite(pred[[cl]])

      if (sum(ok) < 3) {
        return(data.frame(
          class = cl,
          rho = NA_real_,
          p.value = NA_real_,
          stringsAsFactors = FALSE
        ))
      }

      ct <- suppressWarnings(
        cor.test(
          x = pred_ref[[cl]][ok],
          y = pred[[cl]][ok],
          method = "spearman",
          exact = FALSE
        )
      )

      data.frame(
        class = cl,
        rho = unname(ct$estimate),
        p.value = ct$p.value,
        stringsAsFactors = FALSE
      )
    })

    out <- dplyr::bind_rows(out)
    out$p.adj <- p.adjust(out$p.value, method = "BH")
    out
  }

  eval_metrics <- function(pred, truth, class_levels) {
    truth <- factor(truth, levels = class_levels)
    pred_lab <- factor(pred$max_score, levels = class_levels)

    bind_rows(lapply(class_levels, function(cl) {
      pos <- truth == cl
      pred_pos <- pred_lab == cl

      TP <- sum(pos & pred_pos, na.rm = TRUE)
      FP <- sum(!pos & pred_pos, na.rm = TRUE)
      FN <- sum(pos & !pred_pos, na.rm = TRUE)
      TN <- sum(!pos & !pred_pos, na.rm = TRUE)

      sens <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
      spec <- if ((TN + FP) == 0) NA_real_ else TN / (TN + FP)
      prec <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
      f1 <- fbeta(prec, sens, beta = 1)
      f2 <- fbeta(prec, sens, beta = 2)

      auc <- NA_real_
      if (cl %in% colnames(pred) && sum(pos) >= 5 && sum(!pos) >= 5) {
        rr <- tryCatch(
          pROC::roc(
            response = as.logical(pos),
            predictor = pred[[cl]],
            quiet = TRUE,
            direction = "<"
          ),
          error = function(e) NULL
        )
        if (!is.null(rr)) auc <- as.numeric(pROC::auc(rr))
      }

      data.frame(
        class = cl,
        sensitivity = sens,
        specificity = spec,
        precision = prec,
        F1 = f1,
        F2 = f2,
        AUC = auc,
        TP = TP, FP = FP, FN = FN, TN = TN,
        stringsAsFactors = FALSE
      )
    }))
  }

  simulate_pseudocells_per_class_ideal <- function(prob_by_class,
                                                   n_per_class = 1000,
                                                   depth = 10000,
                                                   seed = 1) {
    set.seed(seed)
    genes <- rownames(prob_by_class)
    classes <- colnames(prob_by_class)

    mats <- vector("list", length(classes))
    names(mats) <- classes

    for (cl in classes) {
      p <- as.numeric(prob_by_class[, cl])
      p <- p / sum(p + 1e-12)

      M <- matrix(0L, nrow = length(genes), ncol = n_per_class)
      for (j in seq_len(n_per_class)) {
        M[, j] <- as.integer(rmultinom(1, size = depth, prob = p)[, 1])
      }

      M <- Matrix::Matrix(M, sparse = TRUE)
      colnames(M) <- paste0(cl, "_", seq_len(n_per_class))
      mats[[cl]] <- M
    }

    out <- mats[[1]]
    if (length(mats) > 1) {
      for (k in 2:length(mats)) out <- cbind(out, mats[[k]])
    }
    rownames(out) <- genes
    out
  }

  ## ============================================================
  ## Load EmbryoRadar references
  ## ============================================================
  data.dir <- file.path(system.file("extdata/classifier2.0", package = "EmbryoRadar"))
  if (data.dir == "") stop("Could not find EmbryoRadar extdata/classifier2.0")

  geneAnnPIE <- readRDS(file.path(data.dir, "geneAnnPIE.rds"))
  embryo_prob_by_class <- readRDS(file.path(data.dir, "embryo_prob_by_class.rds"))
  feats <- readRDS(file.path(data.dir, "feats.rds"))
  class_levels <- c("Maternal", "ZGA", "Morula", "EPI", "TE", "Adult")

  if (input_gene_annotations == "gene_name") {
    rownames(embryo_prob_by_class) <- geneAnnPIE$gene_name[match(rownames(embryo_prob_by_class), geneAnnPIE$gene_id)]
    feats <- geneAnnPIE$gene_name[match(feats, geneAnnPIE$gene_id)]
  } else if (input_gene_annotations == "ensembl") {
    rownames(embryo_prob_by_class) <- gsub("\\..*$", "", rownames(embryo_prob_by_class))
    feats <- gsub("\\..*$", "", feats)
  }

  feats <- unique(feats[!is.na(feats)])
  rownames(embryo_prob_by_class) <- make.unique(rownames(embryo_prob_by_class))

  ## ============================================================
  ## Input handling
  ## ============================================================
  raw_counts_background <- NULL
  raw_counts_all <- NULL
  target_source <- NULL

  if (!is.null(raw_counts)) {
    raw_counts_all <- as_dgC(raw_counts)
    rownames(raw_counts_all) <- make.unique(rownames(raw_counts_all))
  }

  if (is.null(target_counts)) {
    if (!isTRUE(auto_filter_from_raw)) {
      stop("target_counts is NULL. Set auto_filter_from_raw = TRUE with raw_counts, or provide target_counts directly.")
    }
    if (is.null(raw_counts_all)) {
      stop("auto_filter_from_raw = TRUE requires raw_counts.")
    }

    message("Running emptyDrops on raw_counts to derive target cells...")
    set.seed(seed)
    ed <- DropletUtils::emptyDrops(raw_counts_all, lower = emptydrops_lower)
    keep <- !is.na(ed$FDR) & ed$FDR <= emptydrops_fdr

    if (sum(keep) == 0) {
      stop("emptyDrops found no cells at the chosen FDR threshold.")
    }

    target_counts <- raw_counts_all[, keep, drop = FALSE]
    raw_counts_background <- raw_counts_all[, !keep, drop = FALSE]
    target_source <- "derived_from_raw_counts_via_emptyDrops"

    message("emptyDrops retained ", sum(keep), " droplets as target cells.")
    message("Background droplets retained for DecontX: ", ncol(raw_counts_background))

    filtered_matrix <- FALSE
  } else {
    target_counts <- as_dgC(target_counts)
    rownames(target_counts) <- make.unique(rownames(target_counts))
    target_source <- "user_supplied_target_counts"

    if (!is.null(raw_counts_all)) {
      raw_counts_background <- raw_counts_all
    }

    if (is.na(filtered_matrix)) filtered_matrix <- TRUE
  }

  if (is.na(filtered_matrix)) filtered_matrix <- FALSE

  ## ============================================================
  ## Harmonise target to embryo reference space
  ## ============================================================
  common_genes_full <- intersect(rownames(target_counts), rownames(embryo_prob_by_class))
  if (length(common_genes_full) < 5000) {
    warning("Low target-embryo gene overlap: ", length(common_genes_full))
  }
  message("Gene overlap (target and embryo): ", length(common_genes_full))

  target_counts_full <- target_counts[common_genes_full, , drop = FALSE]
  embryo_prob_by_class_full <- embryo_prob_by_class[common_genes_full, , drop = FALSE]

  ## ============================================================
  ## Calibration subset
  ## ============================================================
  set.seed(seed)
  if (ncol(target_counts_full) > cal_subsample_n) {
    idx <- sample(seq_len(ncol(target_counts_full)), cal_subsample_n)
    target_cal <- target_counts_full[, idx, drop = FALSE]
  } else {
    idx <- seq_len(ncol(target_counts_full))
    target_cal <- target_counts_full
  }

  if (!is.null(decontx_z)) {
    if (length(decontx_z) != ncol(target_counts)) {
      stop("Length of decontx_z must equal ncol(target_counts) after target cell definition.")
    }
    decontx_z_cal <- decontx_z[idx]
  } else {
    decontx_z_cal <- NULL
  }

  ## ============================================================
  ## Depth
  ## ============================================================
  target_lib <- Matrix::colSums(target_cal)
  depth_med <- as.numeric(stats::median(target_lib))

  ## ============================================================
  ## Ambient estimation
  ## - if pseudobulk = TRUE, skip DecontX entirely
  ## - otherwise:
  ##     * with raw background: DecontX(x, background = raw_bg)
  ##     * without raw background: DecontX(x) on filtered cells only
  ## ============================================================
  soup_profile <- setNames(rep(0, nrow(target_cal)), rownames(target_cal))
  ambient_med <- 0
  ambient_med_raw <- 0
  decontx_contam <- NULL
  decontx_contam_use <- NULL

  target_dx <- target_cal
  z_dx <- decontx_z_cal

  if (ncol(target_dx) > decontx_max_cells) {
    set.seed(seed + 999)
    dx_idx_local <- sample(seq_len(ncol(target_dx)), decontx_max_cells)
    target_dx <- target_dx[, dx_idx_local, drop = FALSE]
    z_dx <- if (!is.null(decontx_z_cal)) decontx_z_cal[dx_idx_local] else NULL
  }

  if (isTRUE(pseudobulk)) {
    message("Skipping DecontX because target matrix looks pseudobulked or too small (n = ", ncol(target_dx), ").")
    message("Setting ambient contamination to 0 and using zero soup profile.")

    soup_profile <- setNames(rep(0, nrow(target_cal)), rownames(target_cal))
    ambient_med_raw <- 0
    ambient_med <- 0
    decontx_contam <- rep(0, ncol(target_dx))
    names(decontx_contam) <- colnames(target_dx)
    decontx_contam_use <- decontx_contam

  } else if (!is.null(raw_counts_background)) {
    raw_common <- intersect(rownames(target_dx), rownames(raw_counts_background))
    if (length(raw_common) < 1000) {
      stop("Too few shared genes between target_dx and raw_counts background for DecontX.")
    }

    target_dx2 <- target_dx[raw_common, , drop = FALSE]
    raw_bg <- raw_counts_background[raw_common, , drop = FALSE]

    if (!is.null(colnames(raw_bg)) && !is.null(colnames(target_dx2))) {
      raw_bg <- raw_bg[, setdiff(colnames(raw_bg), colnames(target_dx2)), drop = FALSE]
    }

    if (ncol(raw_bg) == 0) {
      stop("raw_counts background became empty after removing overlapping barcodes.")
    }

    sce_x <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = target_dx2)
    )
    sce_bg <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = raw_bg)
    )

    message("Running DecontX with raw background for ambient estimation...")
    dx_fit <- celda::decontX(
      x = sce_x,
      background = sce_bg,
      z = z_dx,
      verbose = FALSE
    )

    decontx_contam <- SummarizedExperiment::colData(dx_fit)$decontX_contamination
    names(decontx_contam) <- colnames(dx_fit)

    decontx_contam_use <- pmin(pmax(as.numeric(decontx_contam), 0), ambient_cap)
    decontx_contam_use <- decontx_contam_use[is.finite(decontx_contam_use)]

    ambient_med_raw <- as.numeric(stats::median(decontx_contam, na.rm = TRUE))
    ambient_med <- as.numeric(stats::median(decontx_contam_use, na.rm = TRUE))

    soup_vec <- Matrix::rowSums(raw_bg)
    soup_profile_sub <- as.numeric(soup_vec / sum(soup_vec + 1e-12))
    names(soup_profile_sub) <- rownames(raw_bg)
    soup_profile[names(soup_profile_sub)] <- soup_profile_sub

    message("DecontX median contamination (raw): ", round(ambient_med_raw, 4))
    message("ambient median: ", round(ambient_med, 4))

  } else {
    message("Running DecontX on filtered matrix only (no raw background provided)...")

    sce_x <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = target_dx)
    )

    dx_fit <- celda::decontX(
      x = sce_x,
      z = z_dx,
      verbose = FALSE
    )

    decontx_contam <- SummarizedExperiment::colData(dx_fit)$decontX_contamination
    names(decontx_contam) <- colnames(dx_fit)

    decontx_contam_use <- pmin(pmax(as.numeric(decontx_contam), 0), ambient_cap)
    decontx_contam_use <- decontx_contam_use[is.finite(decontx_contam_use)]

    ambient_med_raw <- as.numeric(stats::median(decontx_contam, na.rm = TRUE))
    ambient_med <- as.numeric(stats::median(decontx_contam_use, na.rm = TRUE))

    native_counts <- tryCatch(
      as(celda::decontXcounts(dx_fit), "dgCMatrix"),
      error = function(e) NULL
    )

    obs_counts <- tryCatch(
      as(SummarizedExperiment::assay(dx_fit, "counts"), "dgCMatrix"),
      error = function(e) NULL
    )

    if (!is.null(native_counts) && !is.null(obs_counts)) {
      contam_counts <- obs_counts - native_counts
      contam_counts@x[contam_counts@x < 0] <- 0
      contam_counts <- Matrix::drop0(contam_counts)
      soup_vec <- Matrix::rowSums(contam_counts)
    } else {
      soup_vec <- rep(0, nrow(target_dx))
      names(soup_vec) <- rownames(target_dx)
    }

    if (sum(soup_vec) <= 0) {
      soup_vec <- Matrix::rowSums(target_dx)
    }

    soup_profile_sub <- as.numeric(soup_vec / sum(soup_vec + 1e-12))
    names(soup_profile_sub) <- names(soup_vec)
    soup_profile[names(soup_profile_sub)] <- soup_profile_sub

    message("DecontX median contamination (raw): ", round(ambient_med_raw, 4))
    message("ambient median: ", round(ambient_med, 4))
  }

  ## ============================================================
  ## Technical parameter estimation
  ## ============================================================
  b_g_raw <- Matrix::rowMeans(target_cal > 0)
  b_g_raw <- pmin(pmax(as.numeric(b_g_raw), 1e-6), 1)
  names(b_g_raw) <- rownames(target_cal)

  target_mean <- as.numeric(Matrix::rowMeans(target_cal))
  target_var  <- rowVars_sparse(target_cal)
  size_g_raw  <- (target_mean^2) / pmax(target_var - target_mean, 1e-8)
  size_g_raw  <- pmin(pmax(size_g_raw, 1e-2), 1e6)
  names(size_g_raw) <- rownames(target_cal)

  target_det <- as.numeric(Matrix::rowMeans(target_cal > 0))
  obs_zero   <- 1 - target_det
  pois_zero  <- exp(-target_mean)
  pi_g_raw   <- pmin(pmax(pmax(obs_zero - pois_zero, 0), 0), 0.95)
  names(pi_g_raw) <- rownames(target_cal)

  stopifnot(all(names(b_g_raw) %in% common_genes_full))
  stopifnot(all(names(size_g_raw) %in% common_genes_full))
  stopifnot(all(names(pi_g_raw) %in% common_genes_full))

  baseline_target_median <- list(
    depth        = depth_med,
    bias_alpha   = 1.0,
    disp_scale   = 1.0,
    pi_scale     = 1.0,
    ambient_frac = ambient_med
  )

  ## ============================================================
  ## Apply technical distortion to an existing ideal matrix
  ## ============================================================
  apply_technical_distortion <- function(
    counts_mat,
    depth_target,
    bias_alpha = 1,
    disp_scale = 1,
    pi_scale = 1,
    ambient_frac = 0,
    seed = 1
  ) {
    set.seed(seed)

    genes <- rownames(counts_mat)
    G <- nrow(counts_mat)
    N <- ncol(counts_mat)

    out <- matrix(0L, nrow = G, ncol = N)
    rownames(out) <- genes
    colnames(out) <- colnames(counts_mat)

    b <- (b_g_raw[genes] ^ bias_alpha)
    b[!is.finite(b)] <- 1

    size_use <- size_g_raw[genes] / disp_scale
    size_use <- pmin(pmax(size_use, 1e-2), 1e6)

    pi_use <- if (pi_scale > 0) pmin(pmax(pi_scale * pi_g_raw[genes], 0), 0.99) else rep(0, G)

    if (length(ambient_frac) == 1) {
      ambient_frac <- rep(ambient_frac, N)
    } else if (length(ambient_frac) != N) {
      stop("ambient_frac must have length 1 or ncol(counts_mat)")
    }

    ambient_frac <- pmin(pmax(as.numeric(ambient_frac), 0), ambient_cap)

    for (j in seq_len(N)) {
      v <- counts_mat[, j]

      if (sum(v) == 0) next

      p <- as.numeric(v) * b
      if (sum(p) <= 0) {
        p <- rep(1 / G, G)
      } else {
        p <- p / sum(p)
      }

      mu <- depth_target * p
      v2 <- rnbinom(n = G, mu = mu, size = size_use)

      if (pi_scale > 0) {
        zmask <- runif(G) < pi_use
        v2[zmask] <- 0L
      }

      ambient_umis <- as.integer(round(depth_target * ambient_frac[j]))
      if (ambient_umis > 0) {
        soup_lambda <- ambient_umis * soup_profile[genes]
        v2 <- v2 + rpois(G, lambda = soup_lambda)
      }

      out[, j] <- v2
    }

    Matrix::Matrix(out, sparse = TRUE)
  }

  ## ============================================================
  ## Prediction wrapper
  ## ============================================================
  ktsp_predict <- function(counts_mat_full, model = "trimmed", block_size = 1000) {
    feats_use <- intersect(feats, rownames(counts_mat_full))
    if (length(feats_use) < 200) stop("Too few k-TSP genes present in counts.")
    counts_sub <- counts_mat_full[feats_use, , drop = FALSE]
    logmat <- log1p_cpm_sparse(counts_sub)
    suppressMessages(
      EmbryoRadar_pred_blockwise(
        ExpressionSet = logmat,
        model = model,
        block_size = block_size,
        input_gene_annotations = input_gene_annotations,
        classes = class_levels
      )$prediction
    )
  }

  ## ============================================================
  ## Main replicate loop
  ## ============================================================
  pred_list <- vector("list", n_reps)
  truth_list <- vector("list", n_reps)
  tb_list <- vector("list", n_reps)
  reftb_list <- vector("list", n_reps)
  metrics_list <- vector("list", n_reps)
  pred_spearman_list <- vector("list", n_reps)

  for (r in seq_len(n_reps)) {
    current_seed <- seed + (r - 1) * 1000
    message("Running replicate ", r, "/", n_reps, " seed=", current_seed)

    sim_counts_ideal <- simulate_pseudocells_per_class_ideal(
      prob_by_class = embryo_prob_by_class_full,
      n_per_class   = n_per_class,
      depth         = ideal_depth,
      seed          = current_seed + 777
    )

    if (!is.null(decontx_contam_use) && length(decontx_contam_use) > 0) {
      ambient_frac_sim <- sample(
        decontx_contam_use,
        size = ncol(sim_counts_ideal),
        replace = TRUE
      )
    } else {
      ambient_frac_sim <- baseline_target_median$ambient_frac
    }

    sim_counts <- apply_technical_distortion(
      counts_mat    = sim_counts_ideal,
      depth_target  = baseline_target_median$depth,
      bias_alpha    = baseline_target_median$bias_alpha,
      disp_scale    = baseline_target_median$disp_scale,
      pi_scale      = baseline_target_median$pi_scale,
      ambient_frac  = ambient_frac_sim,
      seed          = current_seed
    )

    pred_ref <- ktsp_predict(sim_counts_ideal, model = "trimmed", block_size = block_size)
    pred <- ktsp_predict(sim_counts, model = "trimmed", block_size = block_size)

    pred_spearman <- compute_prediction_spearman(
      pred_ref = pred_ref,
      pred = pred,
      class_levels = class_levels
    )
    pred_spearman$rep <- r
    pred_spearman_list[[r]] <- pred_spearman

    truth <- gsub("_.*", "", colnames(sim_counts))

    tb <- table(
      prediction = factor(pred$max_score, levels = class_levels),
      reference  = factor(truth, levels = class_levels)
    )

    ref_tb <- table(
      prediction = factor(pred_ref$max_score, levels = class_levels),
      reference  = factor(truth, levels = class_levels)
    )

    met <- eval_metrics(pred, truth, class_levels)
    met$rep <- r

    pred_list[[r]] <- pred
    truth_list[[r]] <- truth
    tb_list[[r]] <- tb
    reftb_list[[r]] <- ref_tb
    metrics_list[[r]] <- met
  }

  pred_spearman_long <- bind_rows(pred_spearman_list)
  metrics_long <- bind_rows(metrics_list)
  metrics_long <- left_join(metrics_long, pred_spearman_long, by = c("class", "rep"))

  ## ============================================================
  ## QC summary
  ## ============================================================
  qc_metrics_long <- metrics_long %>%
    dplyr::select(class, rep, sensitivity, specificity, rho) %>%
    tidyr::pivot_longer(
      cols = c(sensitivity, specificity, rho),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      class = factor(class, levels = class_levels),
      metric = factor(metric, levels = c("sensitivity", "specificity", "rho"))
    )

  qc_metrics_summary <- qc_metrics_long %>%
    dplyr::group_by(class, metric) %>%
    dplyr::summarise(
      mean_value = mean(value, na.rm = TRUE),
      min_value  = min(value, na.rm = TRUE),
      max_value  = max(value, na.rm = TRUE),
      .groups = "drop"
    )

  qc_plot <- ggplot(qc_metrics_summary, aes(x = class, y = mean_value, fill = class)) +
    geom_col(width = 0.75, colour = "grey30", linewidth = 0.3) +
    geom_errorbar(aes(ymin = min_value, ymax = max_value), width = 0.2, linewidth = 0.5) +
    facet_wrap(~ metric, nrow = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = c(
      "Maternal" = "pink",
      "ZGA"      = "navy",
      "Morula"   = "cyan",
      "EPI"      = "lightgreen",
      "TE"       = "red",
      "Adult"    = "grey60"
    )) +
    labs(
      title = paste0("QC metrics: ", target_name),
      subtitle = paste0("Bar = mean | Error bars = range | n_reps = ", n_reps),
      x = NULL, y = "Score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#F0F4F8"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

  ## ============================================================
  ## Save outputs
  ## ============================================================
  reftb_sum <- Reduce("+", reftb_list)
  tb_sum <- Reduce("+", tb_list)
  ref_cf_hm <- confusionHeatmap_local(reftb_sum)
  cf_hm <- confusionHeatmap_local(tb_sum)

  heatmap_path <- file.path(out_dir, paste0("confusionHM_simulated_at_", target_name, "_level.pdf"))
  pdf(heatmap_path, height = 3, width = 6)
  print(ref_cf_hm + cf_hm)
  dev.off()

  qc_plot_path <- file.path(out_dir, paste0("QC_metrics_", target_name, ".pdf"))
  suppressWarnings(ggsave(qc_plot_path, qc_plot, width = 5.5, height = 4.5, limitsize = FALSE))

  if (isTRUE(save_all)) {
    utils::write.csv(
      metrics_long,
      file.path(out_dir, paste0("prediction_metrics_full_", target_name, ".csv")),
      row.names = FALSE
    )
    utils::write.csv(
      qc_metrics_long,
      file.path(out_dir, paste0("QC_metrics_long_", target_name, ".csv")),
      row.names = FALSE
    )
    saveRDS(
      list(
        target_name = target_name,
        target_source = target_source,
        filtered_matrix = filtered_matrix,
        pseudobulk = pseudobulk,
        auto_filter_from_raw = auto_filter_from_raw,
        common_genes_full = common_genes_full,
        calibration = list(
          depth_med = depth_med,
          ambient_med = ambient_med,
          ambient_med_raw = ambient_med_raw,
          soup_profile = soup_profile,
          decontx_contam = decontx_contam,
          decontx_contam_use = decontx_contam_use,
          b_g_raw = b_g_raw,
          size_g_raw = size_g_raw,
          pi_g_raw = pi_g_raw
        ),
        baseline_target_median = baseline_target_median,
        prediction = pred_list,
        truth = truth_list,
        ref_confusion_table = reftb_list,
        confusion_table = tb_list,
        metrics = metrics_long,
        qc_metrics_long = qc_metrics_long,
        qc_metrics_summary = qc_metrics_summary,
        heatmap_path = heatmap_path,
        qc_plot_path = qc_plot_path
      ),
      file = file.path(out_dir, paste0("sim_eval_results_", target_name, ".rds"))
    )
  }

  message("\nSaved confusion heatmap to: ", heatmap_path)
  message("Saved QC metrics plot to: ", qc_plot_path)

  invisible(list(
    target_name = target_name,
    target_source = target_source,
    filtered_matrix = filtered_matrix,
    auto_filter_from_raw = auto_filter_from_raw,
    common_genes_full = common_genes_full,
    calibration = list(
      depth_med = depth_med,
      ambient_med = ambient_med,
      ambient_med_raw = ambient_med_raw,
      soup_profile = soup_profile,
      decontx_contam = decontx_contam,
      decontx_contam_use = decontx_contam_use,
      b_g_raw = b_g_raw,
      size_g_raw = size_g_raw,
      pi_g_raw = pi_g_raw
    ),
    baseline_target_median = baseline_target_median,
    prediction = pred_list,
    truth = truth_list,
    ref_confusion_table = reftb_list,
    confusion_table = tb_list,
    metrics = metrics_long,
    qc_metrics_long = qc_metrics_long,
    qc_metrics_summary = qc_metrics_summary,
    qc_plot = qc_plot,
    ref_confusionHM = ref_cf_hm,
    confusionHM = cf_hm,
    heatmap_path = heatmap_path,
    qc_plot_path = qc_plot_path
  ))
}
