#' Load trained PIE classifier models and gene annotations
#'
#' This function loads a pre-trained \pkg{multiclassPairs} classifier and its
#' corresponding gene annotations packaged within \pkg{EmbryoRadar}. It supports
#' multiple model configurations and optionally converts model gene identifiers
#' to Ensembl IDs or gene symbols to match user input.
#'
#' @param model Character string specifying which PIE classifier model to load.
#'   One of:
#'   \itemize{
#'     \item \code{"200"} – Classifier trained with 200 k-TSP rules
#'     \item \code{"30"} – Classifier trained with 30 k-TSP rules
#'     \item \code{"10"} – Classifier trained with 10 k-TSP rules
#'   }
#'   Default is \code{"200"}.
#'
#' @param input_gene_annotations Character string specifying how the loaded
#'   model should represent gene identifiers. One of:
#'   \itemize{
#'     \item \code{"original"} – Uses original identifiers as trained
#'     \item \code{"ensembl"} – Converts gene IDs to Ensembl IDs (stripped of version)
#'     \item \code{"gene_name"} – Converts gene IDs to gene symbols
#'   }
#'   Default is \code{"original"}.
#'
#' @details
#' The function loads serialized \pkg{multiclassPairs} models and their
#' corresponding gene annotations from the package's
#' \code{inst/extdata/classifier2.0/} directory.
#' Each model contains a list of classifiers with Top Scoring Pairs (TSPs) used
#' to distinguish developmental lineages or transcriptional states.
#'
#' When an alternative \code{input_gene_annotations} format is specified,
#' the function maps each gene pair (\code{TSP}) in the classifier to either
#' Ensembl IDs or official gene symbols using the \code{geneAnnPIE.rds}
#' annotation table.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{model}}{The loaded \pkg{multiclassPairs} classifier object.}
#'   \item{\code{features}}{A character vector of all unique genes used in the model.}
#' }
#'
#' @examples
#' \dontrun{
#' # Load the default model (200 rules)
#' mod <- load_models(model = "200", input_gene_annotations = "ensembl")
#' names(mod$model$classifiers)
#' length(mod$features)
#' }
#'
#' @import data.table
#' @import multiclassPairs
#' @seealso
#' \code{\link{predPIE2.0}}, \code{\link{compute_TE_singscores}},
#' \code{\link{load_TE_gene_sets}}
#' @export
load_models <- function(model = c("trimmed","200","30","10"),
                        input_gene_annotations = c("original","ensembl","gene_name")) {
  require(data.table)
  require(multiclassPairs)

  data_dir <- system.file("extdata/classifier2.0", package = "EmbryoRadar")

  model_file <- switch(model[1],
                       "trimmed" = "trimmed_EmbryoRadar_with_sc_model_hopefully_final_260220.rds",
                       "200" = "all_gene_ground_truth_classifier2.0_200_with_Normal.rds",
                       "30"  = "all_gene_ground_truth_classifier2.0_30_with_Normal.rds",
                       "10"  = "all_gene_ground_truth_classifier2.0_10_with_Normal.rds"
  )

  lr <- readRDS(file.path(data_dir, model_file))
  geneAnnPIE <- readRDS(file.path(data_dir, "geneAnnPIE.rds"))

  if (input_gene_annotations[1] == "ensembl") {
    geneAnnPIE$ensembl <- gsub("\\..*$","",geneAnnPIE$gene_id)
    for (x in seq_along(lr$classifiers)) {
      ind1 <- match(lr$classifiers[[x]]$TSPs[,1], geneAnnPIE$gene_id)
      ind2 <- match(lr$classifiers[[x]]$TSPs[,2], geneAnnPIE$gene_id)
      lr$classifiers[[x]]$TSPs[,1] <- geneAnnPIE$ensembl[ind1]
      lr$classifiers[[x]]$TSPs[,2] <- geneAnnPIE$ensembl[ind2]
      names(lr$classifiers[[x]]$tieVote) <- names(lr$classifiers[[x]]$score) <-
        rownames(lr$classifiers[[x]]$TSPs) <- paste0(lr$classifiers[[x]]$TSPs[,1],",",lr$classifiers[[x]]$TSPs[,2])
    }
  } else if (input_gene_annotations[1] == "gene_name") {
    for (x in seq_along(lr$classifiers)) {
      ind1 <- match(lr$classifiers[[x]]$TSPs[,1], geneAnnPIE$gene_id)
      ind2 <- match(lr$classifiers[[x]]$TSPs[,2], geneAnnPIE$gene_id)
      lr$classifiers[[x]]$TSPs[,1] <- geneAnnPIE$gene_name[ind1]
      lr$classifiers[[x]]$TSPs[,2] <- geneAnnPIE$gene_name[ind2]
      names(lr$classifiers[[x]]$tieVote) <- names(lr$classifiers[[x]]$score) <-
        rownames(lr$classifiers[[x]]$TSPs) <- paste0(lr$classifiers[[x]]$TSPs[,1],",",lr$classifiers[[x]]$TSPs[,2])
    }
  }

  features <- unique(unlist(lapply(lr$classifiers, function(x) unique(c(x$TSPs[,1], x$TSPs[,2])))))
  list(model = lr, features = features)
}
