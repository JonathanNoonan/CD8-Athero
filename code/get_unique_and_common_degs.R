# Unique vs Common Genes
# Required for FindMarkers(test.use='MAST')
library(MAST)
library(S4Vectors)
library(SummarizedExperiment)
library(SingleCellExperiment)

get_unique_and_common_degs <- function(
    seu,
    groups,
    ident_col      = "Tissue",
    unique_groups  = groups,
    common_sets    = list(),
    min.pct        = 0.25,
    logfc.threshold= 0.2,
    pval_cutoff    = 0.05,
    test.use       = "MAST",
    only.pos       = TRUE
){
  # 1) set identities
  Idents(seu) <- ident_col

  # 2) pairwise comparison helper
  compare_pairs <- function(ident1, ident2) {
    pair <- FindMarkers(
      seu,
      ident.1        = ident1,
      ident.2        = ident2,
      min.pct        = min.pct,
      logfc.threshold= logfc.threshold,
      only.pos       = only.pos,
      test.use       = test.use,
      verbose        = FALSE
    )
    if (nrow(pair) == 0) return(NULL)

    pair$gene <- rownames(pair)
    rownames(pair) <- NULL
    pair <- pair[pair$p_val_adj < pval_cutoff, , drop = FALSE]
    if (nrow(pair) == 0) return(NULL)

    pair$contrast <- paste(ident1, "vs", ident2)
    pair$ident1   <- ident1
    pair$ident2   <- ident2
    pair
  }

  # 3) run all pairwise comparisons for 'groups'
  pairwise_results <- list()
  for (i in seq_along(groups)) {
    for (j in seq_along(groups)) {
      if (i == j) next
      nm <- paste(groups[i], "vs", groups[j])
      pairwise_results[[nm]] <- compare_pairs(groups[i], groups[j])
    }
  }

  # 4) collapse to "up in group vs all others"
  group_up <- list()
  for (g in groups) {
    other <- setdiff(groups, g)
    res_g <- NULL
    for (o in other) {
      nm <- paste(g, "vs", o)
      df <- pairwise_results[[nm]]
      if (is.null(df) || nrow(df) == 0) next

      if (is.null(res_g)) {
        res_g <- df
      } else {
        new_genes <- !df$gene %in% res_g$gene
        if (any(new_genes)) {
          res_g <- rbind(res_g, df[new_genes, , drop = FALSE])
        }
      }
    }
    group_up[[g]] <- res_g
  }

  # 5) unique genes per selected group
  unique_list <- list()
  for (g in unique_groups) {
    df_g <- group_up[[g]]
    if (is.null(df_g) || nrow(df_g) == 0) next

    others <- setdiff(unique_groups, g)
    genes_others <- unique(unlist(lapply(group_up[others], function(x) if (is.null(x)) character(0) else x$gene)))

    unique_g <- df_g[!df_g$gene %in% genes_others, , drop = FALSE]
    if (nrow(unique_g) == 0) next

    unique_g$cluster <- g
    unique_list[[g]] <- unique_g
  }
  unique_combined <- if (length(unique_list) > 0) do.call(rbind, unique_list) else NULL

  # 6) common genes for specified combos
  common_list <- list()
  if (length(common_sets) > 0) {
    for (i in seq_along(common_sets)) {
      gs <- common_sets[[i]]
      label <- paste(gs, collapse = "-")

      gene_lists <- lapply(group_up[gs], function(df) if (is.null(df)) character(0) else df$gene)
      if (any(lengths(gene_lists) == 0)) next  # skip if any group has no genes

      common_genes <- Reduce(intersect, gene_lists)
      if (length(common_genes) == 0) next

      base_df <- group_up[[gs[1]]]
      df_common <- base_df[base_df$gene %in% common_genes, , drop = FALSE]
      df_common$cluster <- label
      common_list[[label]] <- df_common
    }
  }
  common_combined <- if (length(common_list) > 0) do.call(rbind, common_list) else NULL

  # 7) final combined rbind
  if (!is.null(unique_combined) && !is.null(common_combined)) {
    out <- rbind(unique_combined, common_combined)
  } else if (!is.null(unique_combined)) {
    out <- unique_combined
  } else if (!is.null(common_combined)) {
    out <- common_combined
  } else {
    out <- NULL
  }

  # clean rownames
  if (!is.null(out)) rownames(out) <- NULL

  out
}
