## ============================
## Drop all Biobase handling.
## - ExpressionSet can be a matrix/data.frame/etc.
## - We coerce ONCE to a matrix (expr_mat) and then block over columns.
## - Each block_mat passed to multiclassPairs is a MATRIX.
## ============================

EmbryoRadar_pred_blockwise <- function(ExpressionSet = NULL, classes = NULL,
                             model = c("trimmed","200","30","10"),
                             default_threshold = .5,
                             max_score = FALSE, singscore = FALSE,
                             weighted = FALSE, sumddCT = FALSE, ref = NULL,
                             input_gene_annotations = c("original", "ensembl", "gene_name"),
                             threshold = c("default","simulation_2.0","refined","stringent"),
                             tissue_origin = NULL,
                             block_size = 1000) {

  require(data.table)
  require(multiclassPairs)

  model <- model[1]
  threshold <- threshold[1]
  input_gene_annotations <- input_gene_annotations[1]

  data_dir <- system.file("extdata/classifier2.0", package = "EmbryoRadar")

  thr_file <- switch(threshold,
                     "refined" = file.path(data_dir, "refined_threshold.rds"),
                     "stringent" = file.path(data_dir, "threshold_ver1.rds"),
                     "default" = file.path(data_dir, "embryo_threshold_roc_30.csv"),
                     "simulation_2.0" = file.path(data_dir, "simulation_2.0_threshold.csv"))
  thr <- if (grepl("\\.rds$", thr_file)) readRDS(thr_file) else data.table::fread(thr_file)

  lr <- load_models(model = model, input_gene_annotations = input_gene_annotations)$model
  features <- unique(unlist(lapply(lr$classifiers, function(x) unique(c(x$TSPs[,1], x$TSPs[,2])))))
  if (is.null(classes)) classes <- names(lr$classifiers)

  message("Predicting ...")

  ## ---------------------------
  ## sumddCT and singscore branches unchanged (not shown here)
  ## ---------------------------
  ## ... keep your existing sumddCT and singscore code above this point ...

  ## ---------------------------
  ## Blockwise multiclassPairs prediction (MATRIX input)
  ## ---------------------------
  if (is.null(ExpressionSet)) stop("ExpressionSet is NULL")

  expr_mat <- as.matrix(ExpressionSet)  ## <- the only coercion we do
  if (is.null(rownames(expr_mat))) stop("ExpressionSet must have rownames (genes).")
  if (is.null(colnames(expr_mat))) colnames(expr_mat) <- paste0("sample_", seq_len(ncol(expr_mat)))

  run_block <- function(block_mat, verbose = FALSE) {
    multiclassPairs::predict_one_vs_rest_TSP(
      classifier = lr,
      Data = block_mat,                   ## <- guaranteed MATRIX
      tolerate_missed_genes = TRUE,
      weighted_votes = TRUE,
      classes = classes,
      verbose = verbose
    )
  }

  if (is.null(block_size) || block_size >= ncol(expr_mat)) {

    lr_res_pred <- run_block(expr_mat, verbose = TRUE)

  } else {

    stopifnot(is.numeric(block_size), length(block_size) == 1, block_size >= 1)
    block_size <- as.integer(block_size)

    ns <- ncol(expr_mat)
    idx_list <- split(seq_len(ns), ceiling(seq_len(ns) / block_size))

    message("Blockwise prediction: ", ns, " samples in ", length(idx_list),
            " blocks (block_size=", block_size, ")")

    pb <- txtProgressBar(min = 0, max = length(idx_list), style = 3)
    pred_blocks <- vector("list", length(idx_list))

    for (b in seq_along(idx_list)) {
      setTxtProgressBar(pb, b)
      idx <- idx_list[[b]]
      block_mat <- expr_mat[, idx, drop = FALSE]   ## <- MATRIX
      pred_blocks[[b]] <- run_block(block_mat, verbose = FALSE)
    }
    close(pb)

    lr_res_pred <- do.call(rbind, pred_blocks)
    lr_res_pred <- lr_res_pred[colnames(expr_mat), , drop = FALSE]
  }

  ## ---------------------------
  ## Thresholding unchanged
  ## ---------------------------
  mat <- lr_res_pred$Maternal > thr$threshold[thr$class=="Maternal"]
  zga <- lr_res_pred$ZGA > thr$threshold[thr$class=="ZGA"]
  mor <- lr_res_pred$Morula > thr$threshold[thr$class=="Morula"]
  EPI <- lr_res_pred$EPI > thr$threshold[thr$class=="ICM"]
  te  <- lr_res_pred$TE > thr$threshold[thr$class=="TE"]
  nor <- lr_res_pred$Adult > thr$threshold[thr$class=="Normal"]

  lr_res_pred$prediction <- vapply(seq_along(mat), function(i) {
    if (all(!mat[i],!zga[i],!mor[i],!EPI[i],!te[i],!nor[i])) "unclassified" else
      paste0(c("Maternal","ZGA","Morula","EPI","TE","Adult")[c(mat[i],zga[i],mor[i],EPI[i],te[i],nor[i])], collapse="&")
  }, character(1))

  class(lr_res_pred) <- "data.frame"
  list(model = lr, prediction = lr_res_pred, features = features)
}

## Example::
## res <- EmbryoRadar_pred(ExpressionSet = expr_mat, model="200", block_size = 1000)
