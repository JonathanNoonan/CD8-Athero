count_and_plot_deg_lollipop <- function(
    deg_list,
    contrast_levels = c("Plaque", "Blood"),   # 1st → +x, 2nd → -x
    pval_col = "p_val_adj",
    lfc_col  = "avg_log2FC",
    alpha = 0.05,
    lfc_thresh = 0.25,
    gene_col = "gene",
    aspect = 1,
    palette = c(Plaque = "#C00000", Blood = "#1F4E79"),
    x_limits = NULL,
    # NEW: font controls
    font_family = "Times",
    txt_pt = 10,
    pt_size = 4
){
  stopifnot(length(contrast_levels) == 2)

  # pt -> ggplot text size (mm)
  txt_gg <- txt_pt / 2.95

  # 1. Count up-regulated genes per direction
  counts_df <- purrr::map2_dfr(deg_list, names(deg_list), function(df, nm){
    if(!all(c(pval_col, lfc_col, gene_col) %in% names(df)))
      stop(paste("Missing columns in", nm))

    df <- df %>%
      dplyr::filter(!is.na(.data[[pval_col]]),
                    !is.na(.data[[lfc_col]])) %>%
      dplyr::filter(.data[[pval_col]] < alpha,
                    abs(.data[[lfc_col]]) >= lfc_thresh)

    tibble::tibble(
      Item  = nm,
      Group = contrast_levels[c(1,2)],
      Count = c(
        sum(df[[lfc_col]] > 0),
        sum(df[[lfc_col]] < 0)
      )
    )
  }) %>%
    dplyr::mutate(
      CountSigned = ifelse(Group == contrast_levels[1], Count, -Count),
      Group = factor(Group, levels = contrast_levels)
    )

  # 2. Make symmetrical axis limits
  if (is.null(x_limits)) {
    m <- max(abs(counts_df$CountSigned), na.rm = TRUE)
    x_limits <- c(-m * 1.5, m * 1.5)
  }

  # 3. Plot
  counts_df <- counts_df %>%
    dplyr::mutate(Item = factor(Item, levels = rev(names(deg_list))))

  ggplot2::ggplot(counts_df, ggplot2::aes(x = CountSigned, y = Item, fill = Group)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = CountSigned, yend = Item),
      color = "black", linewidth = 1.2
    ) +
    ggplot2::geom_point(shape = 21, size = pt_size, stroke = 1, color = "black") +
    ggplot2::geom_text(
      ggplot2::aes(
        label = Count,
        hjust = ifelse(CountSigned > 0, -0.4, 1.4)
      ),
      size = txt_gg,                 # <-- uses your pt standard
      family = font_family,          # <-- forces Times on labels
      color = "black",
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(
      limits = x_limits,
      breaks = pretty(x_limits, n = 5),
      labels = function(x) abs(x)
    ) +
    ggplot2::labs(x = "# up-regulated genes", y = NULL) +
    ggpubr::theme_pubr(base_size = txt_pt, base_family = font_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = font_family),
      aspect.ratio = aspect,
      axis.text.y = ggplot2::element_text(size = txt_pt, color = "black"),
      axis.text.x = ggplot2::element_text(size = txt_pt, color = "black"),
      axis.title.x = ggplot2::element_text(size = txt_pt + 2, color = "black"),
      legend.position = "none",
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )
}
