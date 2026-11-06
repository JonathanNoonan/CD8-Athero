grabMeta <- function(sc) {
  if (inherits(x=sc, what ="Seurat")) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[length(meta)] <- "cluster.active.ident"
    } else {
      colnames(meta)[length(meta)] <- "cluster"
    }
  }
  else if (inherits(x=sc, what ="SummarizedExperiment")){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[clu] <- "cluster.active.idents"
    } else {
      colnames(meta)[clu] <- "cluster"
    }
  }
  return(meta)
}


theCall <- function(x) {
  if (x == "gene") {
    x <- "CTgene"
  } else if(x == "nt") {
    x <- "CTnt"
  } else if (x == "aa") {
    x <- "CTaa"
  } else if (x == "gene+nt") {
    x <- "CTstrict"
  }
  return(x)
}


recalculateFrequency <- function(sc, cloneCall ="gene+nt",
                                 groupBy = NULL, proportion = TRUE, filterNA = FALSE,
                                 cloneTypes=c(Rare = 1e-4, Small = 0.001,
                                              Medium = 0.01, Large = 0.1, Hyperexpanded = 1)) {
  Con.df <- NULL
  master.index <- NULL
  cloneCall <- theCall(cloneCall)
  options( dplyr.summarise.inform = FALSE )
  cloneTypes <- c(None = 0, cloneTypes)
  meta <- grabMeta(sc)
  location <- which(colnames(meta) == "Frequency")
  meta <- meta[,-location]
  meta$barcode <- rownames(meta)
  cell.names <- rownames(meta)
  dat2 <- na.omit(unique(meta[,c("barcode", cloneCall, groupBy)]))
  dat2 <- as.data.frame(dat2 %>% group_by(dat2[,cloneCall],
                                          dat2[,groupBy]) %>% summarise(Frequency = n()))
  colnames(dat2)[c(1,2)] <- c(cloneCall, groupBy)
  x <- unique(meta[,groupBy])
  for (i in seq_along(x)) {
    sub1 <- subset(meta, meta[,groupBy] == x[i])
    sub2 <- subset(dat2, dat2[,groupBy] == x[i])
    sub1 <- sub1[,-which(colnames(sub1) == groupBy)]
    merge <- merge(sub1, sub2, by=cloneCall)
    if (proportion == TRUE) {
      index <- 1/length(merge$Frequency)
      merge$Frequency <- merge$Frequency/length(merge$Frequency)
      merge$Frequency <- ifelse(merge$Frequency == index, "Single", merge$Frequency)
    }
    Con.df <- rbind.data.frame(Con.df, merge)
  }

  Con.df$cloneType <- NA
  for (x in seq_along(cloneTypes)) { names(cloneTypes)[x] <-
    paste0(names(cloneTypes[x]), ' (', cloneTypes[x-1],
           ' < X <= ', cloneTypes[x], ')') }
  for (i in 2:length(cloneTypes)) { Con.df$cloneType <-
    ifelse(Con.df$Frequency > cloneTypes[i-1] & Con.df$Frequency
           <= cloneTypes[i], names(cloneTypes[i]), Con.df$cloneType) }
  Con.df$cloneType <- ifelse(Con.df$Frequency == "Single", "Single Clone", Con.df$cloneType)
  PreMeta <- Con.df[,c("barcode", "CTgene", "CTnt",
                       "CTaa", "CTstrict", "Frequency", "cloneType")]
  location <- which(colnames(PreMeta) == "Frequency")
  colnames(PreMeta)[location] <- paste0("recalFrequency.", groupBy)
  dup <- PreMeta$barcode[which(duplicated(PreMeta$barcode))]
  `%!in%` = Negate(`%in%`)
  PreMeta <- PreMeta[PreMeta$barcode %!in% dup,]
  rownames(PreMeta) <- PreMeta$barcode
  if (inherits(x=sc, what ="Seurat")) {
    col.name <- names(PreMeta) %||% colnames(PreMeta)
    sc[[col.name]] <- PreMeta
  } else {
    rownames <- rownames(colData(sc))
    colData(sc) <- cbind(colData(sc), PreMeta[rownames,])[, union(colnames(colData(sc)),  colnames(PreMeta))]
    rownames(colData(sc)) <- rownames
  }
  if (filterNA == TRUE) { sc <- filteringNA(sc) }
  return(sc)
}


prepare_circlize_data2 <- function(
  seu_obj,
  scale =0.8
  ){
  celltypes<-levels(seu_obj)
  cell_colors <- scales::hue_pal()(length(celltypes))
  data_plot <- get_metadata2(seu_obj, color = cell_colors, coord_scale = scale)
  data_plot <- cell_order(data_plot)
  data_plot$x_polar2 <- log10(data_plot$x_polar)
  data_plot
}


get_metadata2 <- function(
  seu_obj,
  reductions = "eff",
  coord_scale = 0.8,
  color
  ){
  metadata<-seu_obj@meta.data
  metadata$Cluster<-seu_obj@active.ident
  metadata$dim1<-as.numeric(seu_obj[[reductions]]@cell.embeddings[,1])
  metadata$dim2<-as.numeric(seu_obj[[reductions]]@cell.embeddings[,2])
  metadata$x<-transform_coordinates(metadata$dim1, zoom = coord_scale)
  metadata$y<-transform_coordinates(metadata$dim2, zoom = coord_scale)
  color_df<-data.frame(Cluster=levels(seu_obj), Colors=color)
  cellnames<-rownames(metadata)
  metadata$cells<-rownames(metadata)
  metadata<-merge(metadata, color_df, by='Cluster')
  rownames(metadata)<-metadata$cells
  metadata<-metadata[cellnames,]
  metadata
}


extract_counts <- function(df) {
  df$count <- as.numeric(sub("/.*", "", df$Overlap))
  return(df)
}


insert_line_breaks <- function(text, n = 3) {
  if (!is.character(text)) {
    text <- as.character(text)
  }

  words <- strsplit(text, " ")[[1]]

  if (length(words) > n) {
    split_words <- split(words, ceiling(seq_along(words) / n))
    broken_text <- sapply(split_words, paste, collapse = " ")
    return(paste(broken_text, collapse = "\n"))
  } else {
    return(text)
  }
}


library(data.table)
library(data.table)
library(svglite)
library(igraph)

## ---------- helper to compute network stats per tissue ----------

compute_network_stats <- function(df,
                                  tissue_label,
                                  chain     = "TRA",
                                  sequence  = "aa",
                                  threshold = 0.75) {
  # only PIDs actually present in this tissue subset
  ids <- unique(df$PID)

  out_list <- vector("list", length(ids))
  names(out_list) <- ids

  for (i in seq_along(ids)) {
    pid   <- ids[i]
    sub_i <- df[df$PID == pid, , drop = FALSE]
    if (nrow(sub_i) == 0) next

    # be defensive around clonalCluster
    Clusters <- try(
      clonalCluster(
        sub_i,
        chain    = chain,
        sequence = sequence,
        threshold = threshold
      ),
      silent = TRUE
    )
    if (inherits(Clusters, "try-error")) next

    # pick the right cluster column for this chain (TRA_cluster or TRB_cluster)
    cluster_col <- paste0(chain, "_cluster")
    if (!cluster_col %in% colnames(Clusters)) {
      warning("Column ", cluster_col, " not found for PID = ", pid,
              " (chain = ", chain, "); skipping.")
      next
    }

    x <- as.data.frame(table(Clusters[[cluster_col]]))
    y <- as.data.frame(table(Clusters$cdr3_aa1))

    if (nrow(x) == 0 || nrow(y) == 0) next

    out_list[[i]] <- data.frame(
      Median.Cluster    = median(x$Freq),
      Mean.Cluster      = mean(x$Freq),
      Median.Clone.Size = median(y$Freq),
      Mean.Clone.Size   = mean(y$Freq),
      Tissue            = tissue_label,
      PID               = pid,
      chain             = chain,
      stringsAsFactors  = FALSE
    )
  }

  out_list <- Filter(Negate(is.null), out_list)
  if (!length(out_list)) {
    return(data.frame())
  }

  data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

## ---------- helper to make SVG igraph plots ----------
plot_clone_svg <- function(df, tissue_label, outfile, color, chain,
                           sequence = "aa", threshold = 0.75) {

  # allow vector tissue_label (e.g., c("Stable","Unstable"))
  Tissue_df <- subset(df, Tissue %in% tissue_label)
  if (!nrow(Tissue_df)) {
    warning("No rows for Tissue %in% ", paste(tissue_label, collapse = ", "),
            " – skipping plot for ", outfile)
    return(invisible(NULL))
  }

  clust <- clonalCluster(
    Tissue_df,
    chain       = chain,
    sequence    = sequence,
    threshold   = threshold,
    exportGraph = TRUE
  )

  svglite(outfile, width = 7.5, height = 7.5)
  plot(
    clust,
    vertex.size     = sqrt(igraph::V(clust)$size) * 1.5,
    vertex.label    = NA,
    edge.arrow.size = 0.2,
    vertex.color    = color,
    margin          = 0.5
  )
  invisible(dev.off())
}
