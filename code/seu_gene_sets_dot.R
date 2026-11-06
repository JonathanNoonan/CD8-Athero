seu_gene_sets_dot <- function(
    seu,
    gene_sets, # named list: list("Group A"=c("G1","G2(+)"), ...)
    strip_gene = FALSE,
    group_by = "cluster",              # meta.data column (Seurat mode)
    assay    = NULL,                   # NULL => DefaultAssay(seu)
    use_scaled_for_fill = TRUE,        # Seurat mode: fill = avg.exp.scaled if TRUE else avg.exp

    # ---------- FindAllMarkers / DF support ----------
    markers_df = NULL,                 # pass a FindAllMarkers-style data.frame to use "markers" mode
    markers_gene_col   = "gene",
    markers_group_col  = "cluster",
    markers_lfc_col    = "avg_log2FC",
    markers_padj_col   = "p_val_adj",  # will fall back to "p_adj" if not present
    markers_pct_col    = NULL,         # if present (e.g., "pct.1"), used for size instead of -log10(p)

    # aesthetics
    size_range   = c(1, 20),
    size_breaks  = c(5, 25, 50, 100),  # when size is % expressing
    size_limits  = c(0, 100),          # when size is % expressing; also used to clamp -log10(p)
    fill_limits  = c(-1, 1),           # colour limits for fill_val in BOTH modes
    lfc_limits   = c(-2, 2),           # kept for backwards compat (not used directly)
    fill_palette = c("blue", "lightgrey", "darkred"),
    x_label      = "gene",
    base_size    = 14,
    x_angle      = 45,
    x_text_size  = 12,
    y_text_size  = 10,
    x_title_size = 14,
    bracket_text_size = 5,
    bracket_tick = 0.25,

    # show an X for gene/cluster combos that are *absent* from the DF
    show_missing    = FALSE,
    missing_label   = "×",
    missing_size    = 3,
    missing_color   = "black",
    missing_alpha   = 0.8,

    # layout/headroom
    y_expand_add_bottom = 0.60,
    y_expand_add_top    = 0.25,
    y_bracket_offset    = 0.40,
    y_label_offset      = 0.18,
    plot_margin_cm      = c(1.2, 0.6, 0.6, 0.6),

    # clustering
    cluster_groups = FALSE,
    cluster_genes  = c("none","within_group","global"),
    dist_method    = "euclidean",
    hclust_method  = "complete"
){
  # ---- helpers ----
  clamp_vec <- function(x, lim) {
    if (is.null(lim) || length(lim) != 2 || any(!is.finite(lim))) return(as.numeric(x))
    pmin(pmax(as.numeric(x), lim[1]), lim[2])
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  strip_paren <- function(x) {
    if (strip_gene) {
      gsub("\\s*\\([^)]*\\)", "", x)
    } else {
      x
    }
  }

  contiguous_runs <- function(idx) {
    if (!length(idx)) return(list())
    idx <- sort(idx)
    cuts <- c(0, which(diff(idx) > 1), length(idx))
    runs <- Map(function(a,b) idx[(a+1):b], head(cuts, -1), tail(cuts, -1))
    lapply(runs, function(run) c(run[1], run[length(run)]))
  }

  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Please install 'reshape2' (install.packages('reshape2')).")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install 'dplyr' (install.packages('dplyr')).")
  }

  # ---- inputs / guards ----
  cluster_genes <- match.arg(cluster_genes)
  stopifnot(is.list(gene_sets), length(gene_sets) > 0)

  if (is.null(markers_df)) {
    stopifnot(inherits(seu, "Seurat"))
    if (!group_by %in% colnames(seu@meta.data)) {
      stop(sprintf("group_by '%s' not found in seu@meta.data.", group_by))
    }
  } else {
    stopifnot(is.data.frame(markers_df))
    cols_needed <- c(markers_gene_col, markers_group_col, markers_lfc_col)
    if (!all(cols_needed %in% colnames(markers_df))) {
      stop("markers_df must contain columns: ", paste(cols_needed, collapse = ", "))
    }
    if (!markers_padj_col %in% colnames(markers_df)) {
      if ("p_adj" %in% colnames(markers_df)) {
        markers_padj_col <- "p_adj"
      } else {
        stop("Adjusted p-value column not found in markers_df (tried '",
             markers_padj_col, "' and 'p_adj').")
      }
    }
  }

  if (is.null(assay) && inherits(seu, "Seurat")) assay <- Seurat::DefaultAssay(seu)

  # 1) Order-preserving, de-duplicated display list
  display_groups <- lapply(gene_sets, function(v) v[!duplicated(v)])
  display_all <- unlist(display_groups, use.names = FALSE)
  match_all   <- strip_paren(display_all)

  # 2) Keep genes present in the current data source
  if (is.null(markers_df)) {
    present_mask <- match_all %in% rownames(seu)
  } else {
    present_mask <- match_all %in% markers_df[[markers_gene_col]]
  }
  if (!any(present_mask)) stop("None of the requested genes were found.")
  display_all <- display_all[present_mask]
  match_all   <- match_all[present_mask]

  # prune per-group blocks & remove empties
  if (is.null(markers_df)) {
    group_blocks <- lapply(display_groups, function(v) v[strip_paren(v) %in% rownames(seu)])
  } else {
    group_blocks <- lapply(display_groups, function(v) v[strip_paren(v) %in% markers_df[[markers_gene_col]]])
  }
  nonempty <- vapply(group_blocks, length, integer(1)) > 0
  group_blocks <- group_blocks[nonempty]
  group_names  <- names(gene_sets)[nonempty]

  # ensure duplicates only once while keeping order
  first_keep <- !duplicated(strip_paren(display_all))
  features_ordered <- strip_paren(display_all[first_keep])   # matching / x mapping
  display_ordered  <- display_all[first_keep]                # x-axis labels
  lab_map <- setNames(display_ordered, nm = features_ordered)

  # ---------- Build plotting data 'df' ----------
  if (is.null(markers_df)) {
    # -------- Seurat-expression mode --------
    dp <- Seurat::DotPlot(
      object   = seu,
      features = features_ordered,
      group.by = group_by,
      assay    = assay,
      scale    = use_scaled_for_fill
    )
    df <- dp$data  # features.plot, id, pct.exp, avg.exp, avg.exp.scaled

    # factor levels for y → keep all groups present in metadata
    groups <- seu@meta.data[[group_by]]
    if (!is.factor(groups)) groups <- factor(groups, levels = unique(groups))
    df$id <- factor(as.character(df$id), levels = levels(groups))

    # x order + labels
    df$features.plot <- factor(df$features.plot, levels = features_ordered)

    # fill column (create a clamped version used for plotting)
    if (use_scaled_for_fill) {
      df$fill_val <- clamp_vec(df$avg.exp.scaled, fill_limits)
      fill_scale  <- ggplot2::scale_fill_gradientn(
        colours = fill_palette,
        limits  = fill_limits,
        oob     = scales::squish,
        name    = "avg expr (scaled)"
      )
    } else {
      df$fill_val <- as.numeric(df$avg.exp)
      fill_scale  <- ggplot2::scale_fill_gradientn(
        colours = fill_palette,
        limits  = fill_limits,
        oob     = scales::squish,
        name    = "avg expr"
      )
    }

    aes_size <- "pct.exp"

  } else {
    # -------- FindAllMarkers / DF mode --------
    suppressPackageStartupMessages({ library(dplyr) })

    df0 <- markers_df %>%
      dplyr::rename(
        gene    = !!markers_gene_col,
        cluster = !!markers_group_col,
        lfc     = !!markers_lfc_col,
        padj    = !!markers_padj_col
      ) %>%
      dplyr::filter(gene %in% features_ordered)

    if (!nrow(df0)) {
      stop("No rows in markers_df for the selected genes.")
    }

    # size variable: user-supplied pct column OR -log10(padj) (with protection against 0)
    if (!is.null(markers_pct_col) && markers_pct_col %in% colnames(markers_df)) {
      df0 <- df0 %>%
        dplyr::mutate(pct.exp = as.numeric(.data[[markers_pct_col]]))
      size_var <- "pct.exp"
    } else {
      df0 <- df0 %>%
        dplyr::mutate(
          neglog10p_raw = -log10(pmax(padj, .Machine$double.xmin)),  # avoid Inf
          neglog10p     = clamp_vec(neglog10p_raw, size_limits)      # keep within limits
        )
      size_var <- "neglog10p"
    }

    df <- df0 %>%
      dplyr::transmute(
        features.plot = gene,
        id            = cluster,
        avg_log2FC    = lfc,
        !!size_var    := .data[[size_var]]
      )

    # keep ALL groups from markers_df, even if they have no selected genes
    grp_col <- markers_df[[markers_group_col]]
    if (is.factor(grp_col)) {
      all_groups <- levels(grp_col)     # respect your existing factor levels
    } else {
      all_groups <- unique(grp_col)     # fallback: order of appearance
    }

    df$id <- factor(as.character(df$id), levels = all_groups)
    df$features.plot <- factor(df$features.plot, levels = features_ordered)

    # fill (clamped for colour scale; using fill_limits so behaviour matches Seurat mode)
    df$fill_val <- clamp_vec(df$avg_log2FC, fill_limits)
    fill_scale  <- ggplot2::scale_fill_gradientn(
      colours = fill_palette,
      limits  = fill_limits,
      oob     = scales::squish,
      name    = "avg log2FC"
    )

    aes_size <- size_var
  }

  # ---------- clustering matrix ----------
  mat <- reshape2::acast(
    df,
    formula = features.plot ~ id,
    value.var = "fill_val",
    fun.aggregate = mean,
    na.rm = TRUE
  )
  mat[is.na(mat)] <- 0

  # ---------- optional clustering ----------
  if (cluster_groups && ncol(mat) > 1) {
    dend_g <- stats::hclust(stats::dist(t(mat), method = dist_method), method = hclust_method)
    new_group_order <- colnames(mat)[dend_g$order]
    df$id <- factor(as.character(df$id), levels = new_group_order)
    mat <- mat[, new_group_order, drop = FALSE]
  }

  apply_gene_order <- function(new_order) {
    new_order <- intersect(new_order, rownames(mat))
    df$features.plot <<- factor(as.character(df$features.plot), levels = new_order)
  }

  if (cluster_genes == "within_group") {
    new_gene_order <- character(0)
    for (i in seq_along(group_blocks)) {
      block <- group_blocks[[i]]
      if (!length(block)) next
      syms <- strip_paren(block)
      syms <- syms[syms %in% rownames(mat)]
      if (length(syms) <= 1) {
        new_gene_order <- c(new_gene_order, syms)
      } else {
        dend_b <- stats::hclust(stats::dist(mat[syms, , drop = FALSE], method = dist_method),
                                method = hclust_method)
        new_gene_order <- c(new_gene_order, syms[dend_b$order])
      }
    }
    apply_gene_order(new_gene_order)
  } else if (cluster_genes == "global" && nrow(mat) > 1) {
    dend_all <- stats::hclust(stats::dist(mat, method = dist_method), method = hclust_method)
    new_gene_order <- rownames(mat)[dend_all$order]
    apply_gene_order(new_gene_order)
  }
  # if "none": keep original order

  # ---------- build grid of ALL gene × cluster combos (for X overlay) ----------
  x_levels <- levels(df$features.plot)
  y_levels <- levels(df$id)
  grid_all <- expand.grid(
    features.plot = factor(x_levels, levels = x_levels),
    id            = factor(y_levels, levels = y_levels)
  )

  df_full <- dplyr::left_join(
    grid_all,
    df,
    by = c("features.plot", "id")
  )

  # ---------- plot ----------
  suppressPackageStartupMessages({ library(ggplot2) })
  base_theme <- if (requireNamespace("ggpubr", quietly = TRUE)) {
    ggpubr::theme_pubr(base_size = base_size)
  } else {
    ggplot2::theme_minimal(base_size = base_size)
  }

  # clamp % expressing, if applicable
  if (aes_size == "pct.exp" && "pct.exp" %in% names(df)) {
    df$pct.exp <- clamp_vec(df$pct.exp, size_limits)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = features.plot, y = id)) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data[[aes_size]], fill = .data[["fill_val"]]),
      shape = 21, color = "black", stroke = 0.6
    ) +
    ggplot2::scale_x_discrete(labels = lab_map[levels(df$features.plot)]) +
    {
      # size scale
      if (is.null(markers_df)) {
        ggplot2::scale_size(
          range  = size_range,
          limits = size_limits,
          breaks = size_breaks,
          name   = "% expressing"
        )
      } else if (aes_size == "pct.exp") {
        ggplot2::scale_size(
          range  = size_range,
          limits = size_limits,
          breaks = size_breaks,
          name   = "% expressing"
        )
      } else {
        ggplot2::scale_size_continuous(
          range  = size_range,
          limits = size_limits %||% range(df[[aes_size]], na.rm = TRUE),
          name   = expression(-log[10]("adj p"))
        )
      }
    } +
    fill_scale +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(order = 2),
      size = ggplot2::guide_legend(order = 1)
    ) +
    ggplot2::labs(x = x_label, y = NULL) +
    base_theme +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(add = c(y_expand_add_bottom, y_expand_add_top))
    ) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = x_angle, hjust = 1, vjust = 1, size = x_text_size),
      axis.text.y  = ggplot2::element_text(size = y_text_size),
      axis.title.x = ggplot2::element_text(size = x_title_size, vjust = 1),
      panel.grid   = ggplot2::element_blank(),
      legend.box   = "horizontal",
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 10),
      legend.text  = ggplot2::element_text(size = 9),
      plot.margin  = grid::unit(plot_margin_cm, "cm")
    )

  # ---------- Optional: add X for gene/cluster combos ABSENT from df ----------
  if (show_missing && nrow(df_full) > 0) {
    missing_df <- df_full[is.na(df_full[[aes_size]]), c("features.plot", "id"), drop = FALSE]

    if (nrow(missing_df) > 0) {
      p <- p +
        ggplot2::geom_text(
          data = missing_df,
          ggplot2::aes(x = features.plot, y = id),
          label = missing_label,
          size  = missing_size,
          colour = missing_color,
          alpha  = missing_alpha,
          inherit.aes = FALSE
        )
    }
  }

  # ---------- Brackets + labels ----------
  y_top    <- length(levels(df$id)) + y_bracket_offset
  text_of  <- y_label_offset

  seg_list <- list()
  lab_list <- list()
  for (i in seq_along(group_blocks)) {
    block <- group_blocks[[i]]
    if (!length(block)) next
    syms <- strip_paren(block)
    pos  <- match(syms, x_levels)
    pos  <- pos[!is.na(pos)]
    if (!length(pos)) next

    runs <- contiguous_runs(pos)
    for (r in runs) {
      x1 <- r[1]; x2 <- r[2]; xm <- (x1 + x2)/2
      seg_list[[length(seg_list)+1]] <- data.frame(x1=x1, x2=x2)
      lab_list[[length(lab_list)+1]] <- data.frame(x=xm, label=names(group_blocks)[i], stringsAsFactors = FALSE)
    }
  }
  seg_df   <- if (length(seg_list)) do.call(rbind, seg_list) else data.frame(x1=numeric(0), x2=numeric(0))
  label_df <- if (length(lab_list)) do.call(rbind, lab_list) else data.frame(x=numeric(0), label=character(0))

  if (nrow(seg_df)) {
    p <- p +
      ggplot2::annotate("segment",
                        x = seg_df$x1, xend = seg_df$x2,
                        y = y_top, yend = y_top, size = 0.5) +
      ggplot2::annotate("segment",
                        x = seg_df$x1, xend = seg_df$x1,
                        y = y_top, yend = y_top - bracket_tick, size = 0.5) +
      ggplot2::annotate("segment",
                        x = seg_df$x2, xend = seg_df$x2,
                        y = y_top, yend = y_top - bracket_tick, size = 0.5)
  }
  if (nrow(label_df)) {
    p <- p +
      ggplot2::geom_text(
        data = transform(label_df, y = y_top + text_of),
        ggplot2::aes(x = x, y = y, label = label),
        vjust = -0.2, hjust = 0.5, size = bracket_text_size,
        inherit.aes = FALSE
      )
  }

  return(p)
}
