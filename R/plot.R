#' Get the palette for the given number of categories
#'
#' @keywords internal
#'
get_palette <- function(n) {
  RColorBrewer::brewer.pal(n = n, name = "Accent")
}


#' Plot a similarity matrix heatmap
#'
#' @keywords internal
#'
similarity.heatmap <- function(
                               S,
                               samples = NULL,
                               label = NULL,
                               cluster = NULL,
                               color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
                               label_colors = NULL,
                               label_palette = "Set3",
                               cluster_colors = NULL,
                               cluster_palette = "Paired",
                               ...) {
  #
  # Check matrix is square
  stopifnot(nrow(S) == ncol(S))
  num <- nrow(S)
  #
  # Use integers as sample names if not supplied
  if (is.null(samples)) {
    samples <- 1:num
  }
  stopifnot(num == length(samples))
  #
  # Configuration and meta data
  meta <- data.frame(sample = samples)
  ann_colors <- list()
  rownames(meta) <- colnames(S) <- rownames(S) <- samples
  #
  # If we have labels configure the plot settings
  if (!is.null(label)) {
    meta$label = factor(label)
    if (is.null(label_colors)) {
      label_colors <- RColorBrewer::brewer.pal(n = nlevels(meta$label), name = label_palette)
    }
    if (length(label_colors) < nlevels(meta$label)) {
      stop(sprintf(
        "Not enough colours in palette (%d) for the number of labels (%d)",
        length(label_colors), nlevels(meta$label)
      ))
    }
    ann_colors$label <- label_colors
    names(ann_colors[["label"]]) <- levels(meta$label)
  }
  #
  # If we have clusters configure the plot settings
  if (!is.null(cluster)) {
    meta$cluster = factor(cluster)
    if (is.null(cluster_colors)) {
      cluster_colors <- RColorBrewer::brewer.pal(n = nlevels(meta$cluster), name = cluster_palette)
    }
    if (length(cluster_colors) < nlevels(meta$cluster)) {
      stop(sprintf(
        "Not enough colours in palette (%d) for the number of clusters (%d)",
        length(cluster_colors), nlevels(meta$cluster)
      ))
    }
    ann_colors$cluster <- cluster_colors
    names(ann_colors$cluster) <- levels(meta$cluster)
  }
  #
  # Make the heatmap
  pheatmap::pheatmap(S,
    color = color,
    # annotation_row = meta,
    treeheight_col = 0, # Don't show dendrogram for columns
    annotation_col = meta %>% dplyr::select(-sample),
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_colors = ann_colors,
    ...
  )
}
