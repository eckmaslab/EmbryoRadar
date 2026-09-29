audit_kTSP_technical <- function(
  object,
  technical_cols = c(
    "nCount_RNA", "nFeature_RNA", "nCount_SCT",
    "nFeature_SCT", "percent.mito"
  ),
  adjust_cols = character(),
  group_col = NULL,
  assay = "SCT",
  layer = "data",
  score_col = "TE",
  model_name = "trimmed",
  n_shuffle = 100,
  seed = 1,
  min_cells_per_group = 30,
  out_dir = NULL) {

  required <- c("SeuratObject", "EmbryoRadar", "ggplot2")
  missing_packages <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages)) {
    stop(
      "Missing packages: ",
      paste(missing_packages, collapse = ", ")
    )
  }

  ## Use score_col as a safe filename prefix.
  score_prefix <- gsub("[^A-Za-z0-9._-]+", "_", score_col)
  score_prefix <- gsub("^_+|_+$", "", score_prefix)

  if (!nzchar(score_prefix)) {
    stop("score_col does not produce a valid filename prefix.")
  }

  ## Derive the output directory from score_col unless explicitly supplied.
  if (is.null(out_dir)) {
    out_dir <- paste0(score_prefix, "_technical_audit")
  }

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  safe_spearman <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)

    if (
      sum(ok) < 3L ||
        length(unique(x[ok])) < 2L ||
        length(unique(y[ok])) < 2L
    ) {
      return(NA_real_)
    }

    suppressWarnings(
      stats::cor(
        x[ok],
        y[ok],
        method = "spearman"
      )
    )
  }

  safe_cor_test <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)

    if (
      sum(ok) < 3L ||
        length(unique(x[ok])) < 2L ||
        length(unique(y[ok])) < 2L
    ) {
      return(
        c(
          rho = NA_real_,
          p_value = NA_real_
        )
      )
    }

    z <- suppressWarnings(
      stats::cor.test(
        x[ok],
        y[ok],
        method = "spearman",
        exact = FALSE
      )
    )

    c(
      rho = unname(z$estimate),
      p_value = z$p.value
    )
  }

  ## ============================================================
  ## 1. Extract expression matrix and metadata
  ## ============================================================

  expression_matrix <- tryCatch(
    SeuratObject::LayerData(
      object,
      assay = assay,
      layer = layer
    ),
    error = function(e) {
      SeuratObject::GetAssayData(
        object,
        assay = assay,
        slot = layer
      )
    }
  )

  metadata <- object[[]]

  missing_cells <- setdiff(
    colnames(expression_matrix),
    rownames(metadata)
  )

  if (length(missing_cells)) {
    stop(
      "Some expression-matrix cells are absent from the metadata: ",
      paste(head(missing_cells), collapse = ", ")
    )
  }

  metadata <- metadata[
    colnames(expression_matrix),
    ,
    drop = FALSE
  ]

  technical_cols <- intersect(
    technical_cols,
    colnames(metadata)
  )

  if (!length(technical_cols)) {
    stop("None of technical_cols are present in the metadata.")
  }

  missing_adjust <- setdiff(
    adjust_cols,
    colnames(metadata)
  )

  if (length(missing_adjust)) {
    stop(
      "Missing adjustment columns: ",
      paste(missing_adjust, collapse = ", ")
    )
  }

  if (
    !is.null(group_col) &&
      !group_col %in% colnames(metadata)
  ) {
    stop(
      "group_col is not present in the metadata: ",
      group_col
    )
  }

  if (!score_col %in% colnames(metadata)) {
    stop(
      "Score column is not present in metadata: ",
      score_col
    )
  }

  score_metadata <- suppressWarnings(
    as.numeric(metadata[[score_col]])
  )

  if (all(!is.finite(score_metadata))) {
    stop(
      "The selected score column contains no finite numeric values: ",
      score_col
    )
  }

  ## ============================================================
  ## 2. Load the requested one-vs-rest k-TSP classifier
  ## ============================================================

  loaded_model <- EmbryoRadar::load_models(
    model = model_name,
    input_gene_annotations = "gene_name"
  )$model

  classifier_names <- names(
    loaded_model$classifiers
  )

  if (!score_col %in% classifier_names) {
    stop(
      "No classifier named '",
      score_col,
      "' was found in model '",
      model_name,
      "'. Available classifiers: ",
      paste(classifier_names, collapse = ", ")
    )
  }

  score_classifier <- loaded_model$classifiers[[score_col]]

  if (
    is.null(score_classifier$TSPs) ||
      nrow(score_classifier$TSPs) == 0L
  ) {
    stop(
      "The classifier contains no k-TSP rules: ",
      score_col
    )
  }

  if (
    is.null(score_classifier$score) ||
      length(score_classifier$score) != nrow(score_classifier$TSPs)
  ) {
    stop(
      "Classifier weights are missing or do not match the number of rules."
    )
  }

  rules <- data.frame(
    rule_index = seq_len(
      nrow(score_classifier$TSPs)
    ),
    gene1 = score_classifier$TSPs[, 1],
    gene2 = score_classifier$TSPs[, 2],
    weight = as.numeric(score_classifier$score),
    stringsAsFactors = FALSE
  )

  rules$rule <- paste0(
    rules$gene1,
    " > ",
    rules$gene2
  )

  rules$present <- (
    rules$gene1 %in% rownames(expression_matrix) &
      rules$gene2 %in% rownames(expression_matrix)
  )

  rules_used <- rules[
    rules$present,
    ,
    drop = FALSE
  ]

  if (nrow(rules_used) < 2L) {
    stop(
      "Fewer than two ",
      score_col,
      " rules are evaluable in the selected assay/layer."
    )
  }

  if (
    any(!is.finite(rules_used$weight)) ||
      sum(rules_used$weight) == 0
  ) {
    stop(
      "The evaluable classifier weights are non-finite or sum to zero."
    )
  }

  ## ============================================================
  ## 3. Recompute the selected k-TSP score
  ## ============================================================

  score_rules <- function(gene1, gene2, weights) {
    a <- as.matrix(
      expression_matrix[
        gene1,
        ,
        drop = FALSE
      ]
    )

    b <- as.matrix(
      expression_matrix[
        gene2,
        ,
        drop = FALSE
      ]
    )

    ## Follow the EmbryoRadar/switchBox convention:
    ## gene1 > gene2 gives a vote of 1.
    ## gene1 <= gene2, including ties, gives a vote of 0.
    votes <- a > b

    as.numeric(
      crossprod(
        weights / sum(weights),
        votes
      )
    )
  }

  score_recomputed <- score_rules(
    gene1 = rules_used$gene1,
    gene2 = rules_used$gene2,
    weights = rules_used$weight
  )

  validation <- data.frame(
    score = score_col,
    model = model_name,
    assay = assay,
    layer = layer,
    rules_in_model = nrow(rules),
    rules_evaluated = nrow(rules_used),
    rules_missing = sum(!rules$present),
    max_score_difference = max(
      abs(score_recomputed - score_metadata),
      na.rm = TRUE
    ),
    score_correlation = safe_spearman(
      score_recomputed,
      score_metadata
    ),
    stringsAsFactors = FALSE
  )

  if (
    is.finite(validation$max_score_difference) &&
      validation$max_score_difference > 1e-8
  ) {
    warning(
      "Recomputed ",
      score_col,
      " score does not exactly match the metadata score. ",
      "Check assay, layer, gene annotation and model."
    )
  }

  ## ============================================================
  ## 4. Observed correlations with technical variables
  ## ============================================================

  observed_correlations <- do.call(
    rbind,
    lapply(
      technical_cols,
      function(v) {
        technical_value <- suppressWarnings(
          as.numeric(metadata[[v]])
        )

        z <- safe_cor_test(
          score_recomputed,
          technical_value
        )

        data.frame(
          score = score_col,
          technical_variable = v,
          rho = unname(z["rho"]),
          p_value = unname(z["p_value"]),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  observed_correlations$p_adj_BH <- stats::p.adjust(
    observed_correlations$p_value,
    method = "BH"
  )

  ## ============================================================
  ## 5. Shuffled-rule null test
  ## ============================================================

  random_pairs <- function(
    gene1,
    gene2,
    max_tries = 10000L) {

    pool <- c(gene1, gene2)
    n_rules <- length(gene1)

    for (i in seq_len(max_tries)) {
      shuffled_pool <- sample(
        pool,
        length(pool),
        replace = FALSE
      )

      left <- shuffled_pool[
        seq_len(n_rules)
      ]

      right <- shuffled_pool[
        n_rules + seq_len(n_rules)
      ]

      if (all(left != right)) {
        return(
          list(
            gene1 = left,
            gene2 = right
          )
        )
      }
    }

    stop(
      "Unable to create shuffled rules without self-pairs."
    )
  }

  set.seed(seed)

  shuffled_correlations <- vector(
    "list",
    n_shuffle
  )

  for (i in seq_len(n_shuffle)) {
    shuffled <- random_pairs(
      rules_used$gene1,
      rules_used$gene2
    )

    shuffled_score <- score_rules(
      gene1 = shuffled$gene1,
      gene2 = shuffled$gene2,
      weights = rules_used$weight
    )

    shuffled_correlations[[i]] <- do.call(
      rbind,
      lapply(
        technical_cols,
        function(v) {
          data.frame(
            score = score_col,
            shuffle = i,
            technical_variable = v,
            rho = safe_spearman(
              shuffled_score,
              suppressWarnings(
                as.numeric(metadata[[v]])
              )
            ),
            stringsAsFactors = FALSE
          )
        }
      )
    )
  }

  shuffled_correlations <- do.call(
    rbind,
    shuffled_correlations
  )

  null_summary <- do.call(
    rbind,
    lapply(
      technical_cols,
      function(v) {
        null_r <- shuffled_correlations$rho[
          shuffled_correlations$technical_variable == v
        ]

        observed_r <- observed_correlations$rho[
          observed_correlations$technical_variable == v
        ]

        data.frame(
          score = score_col,
          technical_variable = v,
          observed_rho = observed_r,
          shuffled_median = stats::median(
            null_r,
            na.rm = TRUE
          ),
          shuffled_q025 = unname(
            stats::quantile(
              null_r,
              0.025,
              na.rm = TRUE
            )
          ),
          shuffled_q975 = unname(
            stats::quantile(
              null_r,
              0.975,
              na.rm = TRUE
            )
          ),
          empirical_p = (
            1 +
              sum(
                abs(null_r) >= abs(observed_r),
                na.rm = TRUE
              )
          ) / (
            1 +
              sum(is.finite(null_r))
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  null_summary$empirical_p_adj_BH <- stats::p.adjust(
    null_summary$empirical_p,
    method = "BH"
  )

  ## ============================================================
  ## 6. Partial Spearman correlations
  ## ============================================================

  partial_correlations <- do.call(
    rbind,
    lapply(
      technical_cols,
      function(v) {
        d <- data.frame(
          rank_score = rank(
            score_recomputed,
            ties.method = "average",
            na.last = "keep"
          ),
          rank_technical = rank(
            suppressWarnings(
              as.numeric(metadata[[v]])
            ),
            ties.method = "average",
            na.last = "keep"
          )
        )

        for (j in seq_along(adjust_cols)) {
          value <- metadata[[adjust_cols[j]]]

          d[[paste0("adjust", j)]] <- if (
            is.character(value) ||
              is.factor(value)
          ) {
            factor(value)
          } else {
            value
          }
        }

        d <- d[
          stats::complete.cases(d),
          ,
          drop = FALSE
        ]

        if (nrow(d) < 3L) {
          partial_r <- NA_real_
        } else {
          terms <- paste0(
            "adjust",
            seq_along(adjust_cols)
          )

          rhs <- if (length(terms)) {
            paste(terms, collapse = " + ")
          } else {
            "1"
          }

          m_score <- tryCatch(
            stats::lm(
              stats::as.formula(
                paste("rank_score ~", rhs)
              ),
              data = d
            ),
            error = function(e) NULL
          )

          m_technical <- tryCatch(
            stats::lm(
              stats::as.formula(
                paste("rank_technical ~", rhs)
              ),
              data = d
            ),
            error = function(e) NULL
          )

          if (
            is.null(m_score) ||
              is.null(m_technical)
          ) {
            partial_r <- NA_real_
          } else {
            partial_r <- suppressWarnings(
              stats::cor(
                stats::residuals(m_score),
                stats::residuals(m_technical),
                method = "pearson"
              )
            )
          }
        }

        data.frame(
          score = score_col,
          technical_variable = v,
          raw_rho = safe_spearman(
            score_recomputed,
            suppressWarnings(
              as.numeric(metadata[[v]])
            )
          ),
          partial_rho = partial_r,
          adjustment_columns = paste(
            adjust_cols,
            collapse = ";"
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  ## ============================================================
  ## 7. Rule-level technical audit
  ## ============================================================

  a <- as.matrix(
    expression_matrix[
      rules_used$gene1,
      ,
      drop = FALSE
    ]
  )

  b <- as.matrix(
    expression_matrix[
      rules_used$gene2,
      ,
      drop = FALSE
    ]
  )

  rule_votes <- a > b

  rule_correlations <- do.call(
    rbind,
    lapply(
      technical_cols,
      function(v) {
        technical_value <- suppressWarnings(
          as.numeric(metadata[[v]])
        )

        rule_rho <- apply(
          rule_votes,
          1,
          safe_spearman,
          y = technical_value
        )

        data.frame(
          score = score_col,
          rules_used,
          technical_variable = v,
          rule_rho = rule_rho,
          fraction_both_zero = rowMeans(
            a == 0 & b == 0,
            na.rm = TRUE
          ),
          fraction_exactly_one_detected = rowMeans(
            xor(a > 0, b > 0),
            na.rm = TRUE
          ),
          fraction_both_detected = rowMeans(
            a > 0 & b > 0,
            na.rm = TRUE
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  ## ============================================================
  ## 8. Optional within-group correlations
  ## ============================================================

  stratified_correlations <- NULL

  if (!is.null(group_col)) {
    group_value <- as.character(
      metadata[[group_col]]
    )

    groups <- split(
      seq_len(nrow(metadata)),
      group_value
    )

    stratified_correlations <- do.call(
      rbind,
      lapply(
        names(groups),
        function(g) {
          idx <- groups[[g]]

          do.call(
            rbind,
            lapply(
              technical_cols,
              function(v) {
                z <- safe_cor_test(
                  score_recomputed[idx],
                  suppressWarnings(
                    as.numeric(metadata[[v]][idx])
                  )
                )

                data.frame(
                  score = score_col,
                  group = g,
                  n_cells = length(idx),
                  eligible = length(idx) >= min_cells_per_group,
                  technical_variable = v,
                  rho = unname(z["rho"]),
                  p_value = unname(z["p_value"]),
                  stringsAsFactors = FALSE
                )
              }
            )
          )
        }
      )
    )

    stratified_correlations$p_adj_BH <- stats::p.adjust(
      stratified_correlations$p_value,
      method = "BH"
    )
  }

  ## ============================================================
  ## 9. Plot shuffled-rule null distributions
  ## ============================================================

  null_plot_data <- merge(
    shuffled_correlations,
    observed_correlations[
      ,
      c(
        "technical_variable",
        "rho"
      )
    ],
    by = "technical_variable",
    suffixes = c(
      "_shuffled",
      "_observed"
    )
  )

  null_plot <- ggplot2::ggplot(
    null_plot_data,
    ggplot2::aes(
      x = .data$rho_shuffled
    )
  ) +
    ggplot2::geom_histogram(
      bins = 30,
      fill = "grey65",
      colour = "white"
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(
        xintercept = .data$rho_observed
      ),
      colour = "red",
      linewidth = 0.7
    ) +
    ggplot2::facet_wrap(
      ~technical_variable,
      scales = "free_y"
    ) +
    ggplot2::labs(
      title = paste0(
        score_col,
        " shuffled-rule null test"
      ),
      subtitle = paste0(
        "Red line: observed ",
        score_col,
        " score correlation"
      ),
      x = "Shuffled-score Spearman rho",
      y = "Shuffle count"
    ) +
    ggplot2::theme_bw(
      base_size = 11
    )

  ## ============================================================
  ## 10. Output paths
  ## ============================================================

  paths <- list(
    validation = file.path(
      out_dir,
      paste0(
        score_prefix,
        "_score_validation.csv"
      )
    ),
    observed = file.path(
      out_dir,
      paste0(
        score_prefix,
        "_technical_correlations.csv"
      )
    ),
    shuffled = file.path(
      out_dir,
      paste0(
        score_prefix,
        "_shuffled_correlations.csv"
      )
    ),
    null_summary = file.path(
      out_dir,
      paste0(
        score_prefix,
        "_null_summary.csv"
      )
    ),
    partial = file.path(
      out_dir,
      paste0(
        score_prefix,
        "_partial_correlations.csv"
      )
    ),
    rules = file.path(
      out_dir,
      paste0(
        score_prefix,
        "_rule_correlations.csv"
      )
    ),
    stratified = file.path(
      out_dir,
      paste0(
        score_prefix,
        "_stratified_correlations.csv"
      )
    ),
    null_plot = file.path(
      out_dir,
      paste0(
        score_prefix,
        "_shuffled_null_plot.pdf"
      )
    ),
    result = file.path(
      out_dir,
      paste0(
        score_prefix,
        "_technical_audit.rds"
      )
    )
  )

  utils::write.csv(
    validation,
    paths$validation,
    row.names = FALSE
  )

  utils::write.csv(
    observed_correlations,
    paths$observed,
    row.names = FALSE
  )

  utils::write.csv(
    shuffled_correlations,
    paths$shuffled,
    row.names = FALSE
  )

  utils::write.csv(
    null_summary,
    paths$null_summary,
    row.names = FALSE
  )

  utils::write.csv(
    partial_correlations,
    paths$partial,
    row.names = FALSE
  )

  utils::write.csv(
    rule_correlations,
    paths$rules,
    row.names = FALSE
  )

  if (!is.null(stratified_correlations)) {
    utils::write.csv(
      stratified_correlations,
      paths$stratified,
      row.names = FALSE
    )
  }

  ggplot2::ggsave(
    filename = paths$null_plot,
    plot = null_plot,
    width = 8,
    height = 6
  )

  result <- list(
    score_col = score_col,
    model_name = model_name,
    validation = validation,
    observed_correlations = observed_correlations,
    null_summary = null_summary,
    shuffled_correlations = shuffled_correlations,
    partial_correlations = partial_correlations,
    rule_correlations = rule_correlations,
    stratified_correlations = stratified_correlations,
    plot = null_plot,
    paths = paths
  )

  # saveRDS(
  #   result,
  #   paths$result
  # )

  invisible(result)
}
