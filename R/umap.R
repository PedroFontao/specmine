############################################################################
################################ UMAP ######################################
############################################################################


# perform umap analysis
# n_components - number of dimensions in the embedding (default: 2)
# n_neighbors  - number of nearest neighbours (controls local vs global structure)
# min_dist     - minimum distance between points in the embedding
# metric       - distance metric: "euclidean", "cosine", etc.
# scale        - if TRUE, scales data before UMAP (recommended when variables have very different ranges)
# ret_model    - if TRUE, returns the UMAP model object for later use with new data
# seed         - random seed for reproducibility
# write.file   - if TRUE, saves the embedding to a CSV file
# file.out     - base name for output file

umap_analysis_dataset = function(dataset, n_components = 2, n_neighbors = 15, min_dist = 0.1, metric = "euclidean", scale = FALSE, ret_model = FALSE, seed = 42, write.file = FALSE, file.out = "umap", ...) {
  mat_check = as.matrix(dataset$data)
  if (any(is.na(mat_check)) || any(is.nan(mat_check)) || any(is.infinite(mat_check))) {
    stop("dataset$data contains NA, NaN or Inf values. Please clean your data first.")
  }
  if (!requireNamespace("uwot", quietly = TRUE)) stop("Package 'uwot' is required. Install it with: install.packages('uwot')")
  mat = t(mat_check)
  if (scale) mat = scale(mat)
  set.seed(seed)
  umap_model = uwot::umap(mat, n_components = n_components, n_neighbors = n_neighbors, min_dist = min_dist, metric = metric, ret_model = ret_model, ...)
  if (ret_model) { embedding = umap_model$embedding } else { embedding = umap_model }
  rownames(embedding) = colnames(dataset$data)
  colnames(embedding) = paste0("UMAP", seq_len(n_components))
  if (write.file) write.csv(embedding, file = paste0(file.out, "_embedding.csv"))
  result = list(embedding = embedding, params = list(n_components = n_components, n_neighbors = n_neighbors, min_dist = min_dist, metric = metric, scale = scale, seed = seed))
  if (ret_model) result$model = umap_model
  return(result)
}


########################## UMAP PLOTS ##################################


# 2d scores plot
# column.class - name of metadata column to colour points by
# dims         - which UMAP dimensions to plot (default: c(1,2))
# labels       - if TRUE, shows sample names next to points
# ellipses     - if TRUE, draws confidence ellipses per group
# bw           - if TRUE, uses black and white with shapes instead of colours
# pallette     - colour palette index for scale_colour_brewer
# leg.pos      - legend position ("right", "left", "top", "bottom")
# xlim         - optional x axis limits c(min, max)
# ylim         - optional y axis limits c(min, max)

umap_scoresplot2D = function(dataset, umap.result, column.class = NULL, dims = c(1,2), labels = FALSE, ellipses = FALSE, bw = FALSE, pallette = 2, leg.pos = "right", xlim = NULL, ylim = NULL) {
  has.legend = FALSE
  emb = umap.result$embedding
  if (ncol(emb) < 2) stop("Embedding must have at least 2 components for a 2D plot.")
  umap.points = data.frame(emb[, dims])
  names(umap.points) = c("x", "y")
  if (is.null(column.class)) {
    group.values = factor(rep(4, ncol(dataset$data)))
  } else {
    group.values = dataset$metadata[, column.class]
    has.legend = TRUE
  }
  umap.points$group = group.values
  umap.points$label = colnames(dataset$data)
  if (bw) shape.values = 1:length(levels(group.values))
  if (bw) umap.plot = ggplot2::ggplot(data = umap.points, ggplot2::aes_string(x = 'x', y = 'y', shape = 'group'))
  else umap.plot = ggplot2::ggplot(data = umap.points, ggplot2::aes_string(x = 'x', y = 'y', colour = 'group'))
  umap.plot = umap.plot + ggplot2::geom_point(size = 3, alpha = 1)
  if (bw) umap.plot = umap.plot + ggplot2::scale_shape_manual(values = shape.values)
  else umap.plot = umap.plot + ggplot2::scale_colour_brewer(type = "qual", palette = pallette)
  umap.plot = umap.plot + ggplot2::xlab(paste0("UMAP", dims[1])) + ggplot2::ylab(paste0("UMAP", dims[2])) + ggplot2::ggtitle("UMAP 2D Scores Plot")
  if (has.legend) {
    if (bw) umap.plot = umap.plot + ggplot2::theme_bw()
    else umap.plot = umap.plot + ggplot2::theme(legend.position = leg.pos)
  }
  if (!is.null(xlim)) umap.plot = umap.plot + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) umap.plot = umap.plot + ggplot2::ylim(ylim[1], ylim[2])
  if (labels) umap.plot = umap.plot + ggplot2::geom_text(data = umap.points, ggplot2::aes_string(x = 'x', y = 'y', label = 'label'), hjust = -0.1, vjust = 0)
  if (!bw & ellipses) {
    df.ellipses = calculate_ellipses(umap.points)
    umap.plot = umap.plot + ggplot2::geom_path(data = df.ellipses, ggplot2::aes_string(x = 'x', y = 'y', colour = 'group'), size = 1, linetype = 2)
  }
  umap.plot
}


# 3d scores plot (uses plotly)
# column.class - name of metadata column to colour points by
# dims         - which UMAP dimensions to plot (default: c(1,2,3))
# title        - plot title

umap_scoresplot3D = function(dataset, umap.result, column.class = NULL, dims = c(1,2,3), title = "UMAP 3D Scores Plot") {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required. Install it with: install.packages('plotly')")
  emb = umap.result$embedding
  if (ncol(emb) < 3) stop("Embedding must have at least 3 components for a 3D plot.")
  df = as.data.frame(emb[, dims])
  colnames(df) = c("UMAP1", "UMAP2", "UMAP3")
  if (!is.null(column.class) && column.class %in% colnames(dataset$metadata)) {
    df$Group = as.factor(dataset$metadata[, column.class])
  } else {
    df$Group = factor(rep("Samples", nrow(df)))
  }
  plotly::plot_ly(df, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, color = ~Group, type = "scatter3d", mode = "markers") %>% plotly::layout(title = title)
}


# pairs plot for multiple UMAP components
# column.class - name of metadata column to colour points by
# dims         - which UMAP dimensions to include (default: first 5)

umap_pairs_plot = function(dataset, umap.result, column.class = NULL, dims = 1:min(5, ncol(umap.result$embedding)), ...) {
  if (!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required. Install it with: install.packages('GGally')")
  emb = umap.result$embedding
  if (is.null(column.class)) {
    group.values = rep(4, ncol(dataset$data))
  } else {
    group.values = dataset$metadata[, column.class]
  }
  pairs.df = data.frame(emb[, dims])
  pairs.df$group = group.values
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color = group), ...)
}


# kmeans clustering with 2 UMAP dimensions
# num.clusters  - number of k-means clusters
# dims          - which UMAP dimensions to plot (default: c(1,2))
# kmeans.result - optional pre-computed clustering result (from clustering())
# labels        - if TRUE, shows sample names next to points
# bw            - if TRUE, uses black and white with shapes instead of colours
# ellipses      - if TRUE, draws confidence ellipses per group
# leg.pos       - legend position
# xlim          - optional x axis limits c(min, max)
# ylim          - optional y axis limits c(min, max)

umap_kmeans_plot2D = function(dataset, umap.result, num.clusters = 3, dims = c(1,2), kmeans.result = NULL, labels = FALSE, bw = FALSE, ellipses = FALSE, leg.pos = "right", xlim = NULL, ylim = NULL) {
  emb = umap.result$embedding
  if (is.null(kmeans.result)) kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  umap.points = data.frame(emb[, dims])
  names(umap.points) = c("x", "y")
  umap.points$group = factor(kmeans.result$cluster)
  umap.points$label = colnames(dataset$data)
  if (bw) shape.values = 1:num.clusters
  if (bw) umap.plot = ggplot2::ggplot(data = umap.points, ggplot2::aes_string(x = 'x', y = 'y', shape = 'group'))
  else umap.plot = ggplot2::ggplot(data = umap.points, ggplot2::aes_string(x = 'x', y = 'y', colour = 'group'))
  umap.plot = umap.plot + ggplot2::geom_point(size = 3, alpha = .6)
  if (bw) umap.plot = umap.plot + ggplot2::scale_shape_manual(values = shape.values)
  else umap.plot = umap.plot + ggplot2::scale_colour_brewer(palette = "Set1")
  umap.plot = umap.plot + ggplot2::xlab(paste0("UMAP", dims[1])) + ggplot2::ylab(paste0("UMAP", dims[2])) + ggplot2::ggtitle("UMAP 2D K-means Plot")
  if (bw) umap.plot = umap.plot + ggplot2::theme_bw()
  else umap.plot = umap.plot + ggplot2::theme(legend.position = leg.pos)
  if (!is.null(xlim)) umap.plot = umap.plot + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) umap.plot = umap.plot + ggplot2::ylim(ylim[1], ylim[2])
  if (labels) umap.plot = umap.plot + ggplot2::geom_text(data = umap.points, ggplot2::aes_string(x = 'x', y = 'y', label = 'label'), hjust = -0.1, vjust = 0, size = 3)
  if (!bw & ellipses) {
    df.ellipses = calculate_ellipses(umap.points)
    umap.plot = umap.plot + ggplot2::geom_path(data = df.ellipses, ggplot2::aes_string(x = 'x', y = 'y', colour = 'group'), size = 1, linetype = 2)
  }
  umap.plot
}


# kmeans clustering with 3 UMAP dimensions (uses plotly)
# num.clusters  - number of k-means clusters
# dims          - which UMAP dimensions to plot (default: c(1,2,3))
# kmeans.result - optional pre-computed clustering result (from clustering())
# title         - plot title

umap_kmeans_plot3D = function(dataset, umap.result, num.clusters = 3, dims = c(1,2,3), kmeans.result = NULL, title = "UMAP 3D K-means Plot") {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required. Install it with: install.packages('plotly')")
  emb = umap.result$embedding
  if (ncol(emb) < 3) stop("Embedding must have at least 3 components for a 3D plot.")
  if (is.null(kmeans.result)) kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  df = as.data.frame(emb[, dims])
  colnames(df) = c("UMAP1", "UMAP2", "UMAP3")
  df$Group = factor(kmeans.result$cluster)
  plotly::plot_ly(df, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, color = ~Group, type = "scatter3d", mode = "markers") %>% plotly::layout(title = title)
}


# pairs plot with kmeans clusters
# num.clusters  - number of k-means clusters
# dims          - which UMAP dimensions to include (default: first 5)
# kmeans.result - optional pre-computed clustering result (from clustering())

umap_pairs_kmeans_plot = function(dataset, umap.result, num.clusters = 3, kmeans.result = NULL, dims = 1:min(5, ncol(umap.result$embedding))) {
  if (!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required. Install it with: install.packages('GGally')")
  if (is.null(kmeans.result)) kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  emb = umap.result$embedding
  pairs.df = data.frame(emb[, dims])
  pairs.df$group = factor(kmeans.result$cluster)
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color = group))
}
