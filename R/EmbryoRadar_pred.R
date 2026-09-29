#' Predict developmental lineage or transcriptional class using the PIE classifier
#'
#' This is the main prediction function of the \pkg{EmbryoRadar} package.
#' It applies the Pre-Implantation Embryo (PIE) classifier to expression data
#' (bulk or single-cell) to estimate similarity with developmental or
#' transcriptional lineages such as Maternal, ZGA, Morula, EPI, TE, or Adult.
#'
#' @param ExpressionSet A \code{Biobase::ExpressionSet} object containing
#'   log-transformed or Adultized expression data. Row names should be gene
#'   identifiers and column names should correspond to samples.
#' @param classes Optional character vector specifying which classes to predict.
#'   Defaults to all classes available in the selected model.
#' @param model Character string specifying the trained model to use.
#'   One of \code{"200"}, \code{"30"}, or \code{"10"}.
#'   \itemize{
#'     \item \code{"200"} – Classifier trained with 200 k-TSP rules (full model)
#'     \item \code{"30"} – Classifier trained with 30 k-TSP rules (bulk/well-based behaviour; supports tissue thresholds)
#'     \item \code{"10"} – Classifier trained with 10 k-TSP rules (mini model)
#'   }
#' @param default_threshold Numeric; default threshold for classification when
#'   tissue-specific thresholds are unavailable. Default is \code{0.5}.
#' @param max_score Logical; if \code{TRUE}, assigns each sample to the class
#'   with the highest score passing threshold criteria. Default is \code{FALSE}.
#' @param singscore Logical; if \code{TRUE}, applies \code{\link{singscore}}
#'   to rank-based scoring of TSP rules instead of pairwise comparison.
#' @param weighted Logical; if \code{TRUE}, uses classifier rule weights for
#'   scoring aggregation (available when \code{singscore = TRUE}).
#' @param sumddCT Logical; if \code{TRUE}, computes the mean ΔΔCT score for
#'   each TSP rule (used in qPCR-based experiments).
#' @param ref Optional reference sample or group for ΔΔCT Adultization.
#' @param input_gene_annotations Character string specifying how input genes
#'   are annotated: one of \code{"original"}, \code{"ensembl"}, or
#'   \code{"gene_name"}.
#' @param threshold Character string specifying which threshold set to use for
#'   classification. One of \code{"default"}, \code{"simulation_2.0"},
#'   \code{"refined"}, or \code{"stringent"}.
#' @param tissue_origin Optional character string specifying the tissue type
#'   for applying tissue-specific thresholds (relevant for the \code{"30"} model).
#'
#' @details
#' This function loads the trained multiclass pairwise comparison models via
#' \code{\link{load_models}} and applies them to an input expression dataset.
#' Depending on the parameters, it supports three main scoring modes:
#' \enumerate{
#'   \item **TSP prediction (default):** Standard multiclassPairs rule-based
#'         scoring using log-expression values.
#'   \item **Singscore mode:** Rank-based scoring via
#'         \code{singscore::simpleScore()} for robustness in single-cell data.
#'   \item **ΔΔCT mode:** Mean ΔΔCT calculation for rule pairs, suited for qPCR.
#' }
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{model}}{The trained classifier object loaded via
#'         \code{\link{load_models}}.}
#'   \item{\code{prediction}}{A \code{data.frame} containing class probabilities
#'         or scores for each sample.}
#'   \item{\code{features}}{A vector of all genes used by the model.}
#' }
#'
#' @examples
#' \dontrun{
#' library(Biobase)
#' data <- readRDS(system.file("extdata/classifier2.0/GSE232692_exprset.rds",
#'                             package = "EmbryoRadar"))
#' res <- predPIE2.0(ExpressionSet = data,
#'                   model = "200",
#'                   input_gene_annotations = "ensembl",
#'                   threshold = "simulation_2.0")
#' head(res$prediction)
#' }
#'
#' @import data.table
#' @import multiclassPairs
#' @seealso
#' \code{\link{load_models}}, \code{\link{compute_TE_singscores}},
#' \code{\link{singscore::simpleScore}}
#' @export
EmbryoRadar_pred <- function(ExpressionSet = NULL, classes = NULL,
                       model = c("trimmed","200","30","10"),
                       default_threshold = .5,
                       max_score = FALSE, singscore = FALSE,
                       weighted = FALSE, sumddCT = FALSE, ref = NULL,
                       input_gene_annotations = c("original", "ensembl", "gene_name"),
                       threshold = c("default","simulation_2.0","refined","stringent"),
                       tissue_origin = NULL) {
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

  if (sumddCT) {
    message("Calculating sumddCT for each rule")
    out <- setNames(lapply(classes, function(cl) {
      message(paste0("Scoring ", cl))
      pb <- txtProgressBar(min = 0, max = nrow(lr$classifiers[[cl]]$TSPs),
                           style = 3)
      score_dat <- lapply(seq_len(nrow(lr$classifiers[[cl]]$TSPs)),
                          function(i) {
                            print(i)
                            hi <- lr$classifiers[[cl]]$TSPs[i, 1]
                            lo <- lr$classifiers[[cl]]$TSPs[i, 2]
                            cn <- paste0(hi,"-",lo)
                            setTxtProgressBar(pb, i)
                            if (all(hi %in% rownames(ExpressionSet), lo %in% rownames(ExpressionSet))) {
                              if ("ExpressionSet" %in% class(ExpressionSet)) {
                                df <- data.frame(ExpressionSet[hi, ]@assayData$exprs -
                                                   ExpressionSet[lo, ]@assayData$exprs)
                              } else {
                                df <- data.frame(ExpressionSet[hi, ] - ExpressionSet[lo, ])
                                colnames(df) <- cn
                              }
                            } else {
                              df <- data.frame(rep(NA, ncol(ExpressionSet)))
                              colnames(df) <- cn
                            }
                            return(df)
                          })
      close(pb)
      score_dat <- suppressMessages(dplyr::bind_cols(score_dat))
      score_dat <- data.frame(sample = colnames(ExpressionSet),
                              mean_sumddCT = rowMeans(score_dat, na.rm = TRUE),
                              score_dat)
      list(model = lr, prediction = score_dat, features = c(lr$classifiers[[cl]]$TSPs[,
                                                                                      1], lr$classifiers[[cl]]$TSPs[, 2]))
    }), nm = classes)
    return(out)
  }

  if (singscore) {
    require(singscore)
    message("Applying rules with singscore")
    rankData <- singscore::rankGenes(ExpressionSet)
    out <- setNames(lapply(classes, function(cl){
      message(paste0("Scoring ",cl))
      pb <- txtProgressBar(min = 0, max = nrow(lr$classifiers[[cl]]$TSPs), style = 3)
      score_dat <- lapply(seq_len(nrow(lr$classifiers[[cl]]$TSPs)), function(i) {
        hi <- lr$classifiers[[cl]]$TSPs[i,1]; lo <- lr$classifiers[[cl]]$TSPs[i,2]
        setTxtProgressBar(pb, i)
        if (hi %in% rownames(rankData) & lo %in% rownames(rankData)) {
          s <- singscore::simpleScore(rankData, upSet = hi, downSet = lo)
        } else if (hi %in% rownames(rankData)) {
          s <- singscore::simpleScore(rankData, upSet = hi)
        } else {
          s <- singscore::simpleScore(rankData, upSet = NA, downSet = NA)
        }
        s$Weighted_score <- s$TotalScore * lr$classifiers[[cl]]$score[i]
        data.table::data.table(score = if (weighted) s$Weighted_score else s$TotalScore)
      })
      close(pb)
      score_dat <- suppressMessages(setNames(dplyr::bind_cols(score_dat),
                                             paste0(lr$classifiers[[cl]]$TSPs[,1], "-", lr$classifiers[[cl]]$TSPs[,2])))
      score_dat <- cbind(sample = colnames(rankData),
                         mean_singscore = rowMeans(score_dat, na.rm = TRUE), score_dat)
      list(model = lr, prediction = score_dat,
           features = c(lr$classifiers[[cl]]$TSPs[,1], lr$classifiers[[cl]]$TSPs[,2]))
    }), nm = classes)
    return(out)
  }

  lr_res_pred <- multiclassPairs::predict_one_vs_rest_TSP(
    classifier = lr, Data = ExpressionSet,
    tolerate_missed_genes = TRUE, weighted_votes = TRUE,
    classes = classes, verbose = TRUE
  )

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
  # colnames(lr_res_pred)[which(colnames(lr_res_pred) == "EPI")] <- "EPI"
  # colnames(lr_res_pred)[which(colnames(lr_res_pred) == "Adult")] <- "Adult"
  list(model = lr, prediction = lr_res_pred, features = features)
}


