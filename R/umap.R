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

#' UMAP analysis
#'
#' Performs Uniform Manifold Approximation and Projection (UMAP)
#' dimensionality reduction on a specmine dataset.
#'
#' @param dataset Dataset to analyse.
#' @param n_components Number of dimensions in the embedding.
#' @param n_neighbors Number of nearest neighbours.
#' @param min_dist Minimum distance between points in the embedding.
#' @param metric Distance metric passed to \code{uwot::umap()}.
#' @param scale Logical indicating whether the data should be scaled.
#' @param ret_model Logical indicating whether the UMAP model should be returned.
#' @param seed Random seed used for reproducibility.
#' @param write.file Logical indicating whether the embedding should be written
#'   to a CSV file.
#' @param file.out Output file prefix used when \code{write.file = TRUE}.
#' @param ... Additional arguments passed to \code{uwot::umap()}.
#'
#' @return A list containing the UMAP embedding and the analysis parameters.
#'   When \code{ret_model = TRUE}, the fitted UMAP model is also returned.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("uwot", quietly = TRUE)) {
#'   datamat <- matrix(
#'     rnorm(40),
#'     nrow = 5,
#'     dimnames = list(paste0("x", 1:5), paste0("s", 1:8))
#'   )
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(class = factor(rep(c("A", "B"), length.out = 8)))
#'   )
#'   result <- umap_analysis_dataset(
#'     dataset,
#'     n_components = 2,
#'     n_neighbors = 3,
#'     seed = 42
#'   )
#' }
#' }
#'
#' @export
umap_analysis_dataset = function(dataset, n_components = 2, n_neighbors = 15, min_dist = 0.1, metric = "euclidean", scale = FALSE, ret_model = FALSE, seed = 42, write.file = FALSE, file.out = NULL, ...) {
  mat_check = as.matrix(dataset$data)
  if (any(is.na(mat_check)) || any(is.nan(mat_check)) || any(is.infinite(mat_check))) {
    stop("dataset$data contains NA, NaN or Inf values. Please clean your data first.")
  }
  if (!requireNamespace("uwot", quietly = TRUE)) stop("Package 'uwot' is required. Install it with: install.packages('uwot')")
  mat = t(mat_check)
  if (scale) mat = scale(mat)
  set.seed(seed)
  umap_model = uwot::umap(mat, n_components = n_components, n_neighbors = n_neighbors, min_dist = min_dist, metric = metric, ret_model = ret_model, ...)
  if (ret_model) {
    embedding = umap_model$embedding
  } else {
    embedding = umap_model
  }
  rownames(embedding) = colnames(dataset$data)
  colnames(embedding) = paste0("UMAP", seq_len(n_components))
  if (isTRUE(write.file)) {
    if (is.null(file.out) || !nzchar(file.out)) {
      stop("Please provide 'file.out' when write.file = TRUE.")
    }
    utils::write.csv(embedding, file = paste0(file.out, "_embedding.csv"))
  }
  result = list(
    embedding = embedding,
    params = list(
      n_components = n_components,
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      metric = metric,
      scale = scale,
      seed = seed
    )
  )
  if (ret_model) result$model = umap_model
  return(result)
}

########################## UMAP PLOTS ##################################

#' UMAP 2D scores plot
#'
#' Creates a two-dimensional scores plot from a UMAP result.
#'
#' @param dataset Dataset used to calculate the embedding.
#' @param umap.result Result returned by \code{umap_analysis_dataset()}.
#' @param column.class Optional metadata column used to colour or group samples.
#' @param dims Two embedding dimensions to plot.
#' @param labels Logical indicating whether sample labels should be displayed.
#' @param ellipses Logical indicating whether group ellipses should be displayed.
#' @param bw Logical indicating whether a black-and-white plot should be used.
#' @param pallette RColorBrewer palette number.
#' @param leg.pos Legend position.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#' @return A \code{ggplot} object.
#' @examples
#' \donttest{
#' if (requireNamespace("uwot", quietly = TRUE)) {
#'   datamat <- matrix(rnorm(40), nrow = 5)
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), length.out = 8))
#'     )
#'   )
#'   result <- umap_analysis_dataset(dataset, n_neighbors = 3)
#'   umap_scoresplot2D(dataset, result, column.class = "class")
#' }
#' }
#' @export
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
  if (bw) {
    umap.plot = ggplot2::ggplot(data = umap.points, ggplot2::aes_string(x = "x", y = "y", shape = "group"))
  } else {
    umap.plot = ggplot2::ggplot(data = umap.points, ggplot2::aes_string(x = "x", y = "y", colour = "group"))
  }
  umap.plot = umap.plot + ggplot2::geom_point(size = 3, alpha = 1)
  if (bw) {
    umap.plot = umap.plot + ggplot2::scale_shape_manual(values = shape.values)
  } else {
    umap.plot = umap.plot + ggplot2::scale_colour_brewer(type = "qual", palette = pallette)
  }
  umap.plot = umap.plot + ggplot2::xlab(paste0("UMAP", dims[1])) + ggplot2::ylab(paste0("UMAP", dims[2])) + ggplot2::ggtitle("UMAP 2D Scores Plot")
  if (has.legend) {
    if (bw) umap.plot = umap.plot + ggplot2::theme_bw()
    else umap.plot = umap.plot + ggplot2::theme(legend.position = leg.pos)
  }
  if (!is.null(xlim)) umap.plot = umap.plot + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) umap.plot = umap.plot + ggplot2::ylim(ylim[1], ylim[2])
  if (labels) umap.plot = umap.plot + ggplot2::geom_text(data = umap.points, ggplot2::aes_string(x = "x", y = "y", label = "label"), hjust = -0.1, vjust = 0)
  if (!bw & ellipses) {
    df.ellipses = calculate_ellipses(umap.points)
    umap.plot = umap.plot + ggplot2::geom_path(data = df.ellipses, ggplot2::aes_string(x = "x", y = "y", colour = "group"), size = 1, linetype = 2)
  }
  umap.plot
}

#' UMAP 3D scores plot
#'
#' Creates an interactive three-dimensional scores plot from a UMAP result.
#'
#' @param dataset Dataset used to calculate the embedding.
#' @param umap.result Result returned by \code{umap_analysis_dataset()}.
#' @param column.class Optional metadata column used to colour samples.
#' @param dims Three embedding dimensions to plot.
#' @param title Plot title.
#' @return A \code{plotly} htmlwidget.
#' @examples
#' \donttest{
#' if (requireNamespace("uwot", quietly = TRUE) &&
#'     requireNamespace("plotly", quietly = TRUE)) {
#'   datamat <- matrix(rnorm(40), nrow = 5)
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), length.out = 8))
#'     )
#'   )
#'   result <- umap_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     n_neighbors = 3
#'   )
#'   umap_scoresplot3D(dataset, result, column.class = "class")
#' }
#' }
#' @export
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
  plotly::layout(
    plotly::plot_ly(df, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, color = ~Group, type = "scatter3d", mode = "markers"),
    title = title
  )
}

#' UMAP pairs plot
#'
#' Creates a pairs plot for selected UMAP dimensions.
#'
#' @param dataset Dataset used to calculate the embedding.
#' @param umap.result Result returned by \code{umap_analysis_dataset()}.
#' @param column.class Optional metadata column used to colour samples.
#' @param dims Embedding dimensions to include in the pairs plot.
#' @param ... Additional arguments passed to \code{GGally::ggpairs()}.
#' @return A \code{GGally} pairs plot object.
#' @examples
#' \donttest{
#' if (requireNamespace("uwot", quietly = TRUE) &&
#'     requireNamespace("GGally", quietly = TRUE)) {
#'   datamat <- matrix(rnorm(40), nrow = 5)
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), length.out = 8))
#'     )
#'   )
#'   result <- umap_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     n_neighbors = 3
#'   )
#'   umap_pairs_plot(
#'     dataset,
#'     result,
#'     column.class = "class",
#'     dims = 1:3
#'   )
#' }
#' }
#' @export
umap_pairs_plot = function(dataset, umap.result, column.class = NULL, dims = 1:min(5, ncol(umap.result$embedding)), ...) {
  if (!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required. Install it with: install.packages('GGally')")
  emb = umap.result$embedding
  if (is.null(column.class)) group.values = rep(4, ncol(dataset$data))
  else group.values = dataset$metadata[, column.class]
  pairs.df = data.frame(emb[, dims])
  pairs.df$group = group.values
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color = group), ...)
}

#' UMAP 2D k-means plot
#'
#' Creates a two-dimensional UMAP plot coloured by k-means clusters.
#'
#' @param dataset Dataset to cluster.
#' @param umap.result Result returned by \code{umap_analysis_dataset()}.
#' @param num.clusters Number of k-means clusters.
#' @param dims Two embedding dimensions to plot.
#' @param kmeans.result Optional result returned by \code{clustering()}.
#' @param labels Logical indicating whether sample labels should be displayed.
#' @param bw Logical indicating whether a black-and-white plot should be used.
#' @param ellipses Logical indicating whether group ellipses should be displayed.
#' @param leg.pos Legend position.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' if (requireNamespace("uwot", quietly = TRUE)) {
#'   datamat <- matrix(
#'     rnorm(40),
#'     nrow = 5,
#'     dimnames = list(
#'       paste0("x", 1:5),
#'       paste0("s", 1:8)
#'     )
#'   )
#'
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), length.out = 8))
#'     )
#'   )
#'
#'   umap.result <- umap_analysis_dataset(
#'     dataset,
#'     n_components = 2,
#'     n_neighbors = 3,
#'     seed = 42
#'   )
#'
#'   umap_kmeans_plot2D(
#'     dataset,
#'     umap.result,
#'     num.clusters = 2,
#'     dims = c(1, 2)
#'   )
#' }
#'
#' @export
umap_kmeans_plot2D = function(dataset, umap.result, num.clusters = 3, dims = c(1,2), kmeans.result = NULL, labels = FALSE, bw = FALSE, ellipses = FALSE, leg.pos = "right", xlim = NULL, ylim = NULL) {
  emb = umap.result$embedding
  if (is.null(kmeans.result)) kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  umap.points = data.frame(emb[, dims])
  names(umap.points) = c("x", "y")
  umap.points$group = factor(kmeans.result$cluster)
  umap.points$label = colnames(dataset$data)
  if (bw) shape.values = 1:num.clusters
  if (bw) umap.plot = ggplot2::ggplot(data = umap.points, ggplot2::aes_string(x = "x", y = "y", shape = "group"))
  else umap.plot = ggplot2::ggplot(data = umap.points, ggplot2::aes_string(x = "x", y = "y", colour = "group"))
  umap.plot = umap.plot + ggplot2::geom_point(size = 3, alpha = .6)
  if (bw) umap.plot = umap.plot + ggplot2::scale_shape_manual(values = shape.values)
  else umap.plot = umap.plot + ggplot2::scale_colour_brewer(palette = "Set1")
  umap.plot = umap.plot + ggplot2::xlab(paste0("UMAP", dims[1])) + ggplot2::ylab(paste0("UMAP", dims[2])) + ggplot2::ggtitle("UMAP 2D K-means Plot")
  if (bw) umap.plot = umap.plot + ggplot2::theme_bw()
  else umap.plot = umap.plot + ggplot2::theme(legend.position = leg.pos)
  if (!is.null(xlim)) umap.plot = umap.plot + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) umap.plot = umap.plot + ggplot2::ylim(ylim[1], ylim[2])
  if (labels) umap.plot = umap.plot + ggplot2::geom_text(data = umap.points, ggplot2::aes_string(x = "x", y = "y", label = "label"), hjust = -0.1, vjust = 0, size = 3)
  if (!bw & ellipses) {
    df.ellipses = calculate_ellipses(umap.points)
    umap.plot = umap.plot + ggplot2::geom_path(data = df.ellipses, ggplot2::aes_string(x = "x", y = "y", colour = "group"), size = 1, linetype = 2)
  }
  umap.plot
}

#' UMAP 3D k-means plot
#'
#' Creates an interactive three-dimensional UMAP plot coloured by k-means
#' clusters.
#'
#' @param dataset Dataset to cluster.
#' @param umap.result Result returned by \code{umap_analysis_dataset()}.
#' @param num.clusters Number of k-means clusters.
#' @param dims Three embedding dimensions to plot.
#' @param kmeans.result Optional result returned by \code{clustering()}.
#' @param title Plot title.
#'
#' @return A \code{plotly} htmlwidget.
#'
#' @examples
#' if (requireNamespace("uwot", quietly = TRUE) &&
#'     requireNamespace("plotly", quietly = TRUE)) {
#'   datamat <- matrix(
#'     rnorm(40),
#'     nrow = 5,
#'     dimnames = list(
#'       paste0("x", 1:5),
#'       paste0("s", 1:8)
#'     )
#'   )
#'
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), length.out = 8))
#'     )
#'   )
#'
#'   umap.result <- umap_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     n_neighbors = 3,
#'     seed = 42
#'   )
#'
#'   umap_kmeans_plot3D(
#'     dataset,
#'     umap.result,
#'     num.clusters = 2,
#'     dims = c(1, 2, 3)
#'   )
#' }
#'
#' @export
umap_kmeans_plot3D = function(dataset, umap.result, num.clusters = 3, dims = c(1,2,3), kmeans.result = NULL, title = "UMAP 3D K-means Plot") {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required. Install it with: install.packages('plotly')")
  emb = umap.result$embedding
  if (ncol(emb) < 3) stop("Embedding must have at least 3 components for a 3D plot.")
  if (is.null(kmeans.result)) kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  df = as.data.frame(emb[, dims])
  colnames(df) = c("UMAP1", "UMAP2", "UMAP3")
  df$Group = factor(kmeans.result$cluster)
  plotly::layout(
    plotly::plot_ly(df, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, color = ~Group, type = "scatter3d", mode = "markers"),
    title = title
  )
}

#' UMAP pairs k-means plot
#'
#' Creates a pairs plot for UMAP dimensions coloured by k-means clusters.
#'
#' @param dataset Dataset to cluster.
#' @param umap.result Result returned by \code{umap_analysis_dataset()}.
#' @param num.clusters Number of k-means clusters.
#' @param kmeans.result Optional result returned by \code{clustering()}.
#' @param dims Embedding dimensions to include in the pairs plot.
#'
#' @return A \code{GGally} pairs plot object.
#'
#' @examples
#' if (requireNamespace("uwot", quietly = TRUE) &&
#'     requireNamespace("GGally", quietly = TRUE)) {
#'   datamat <- matrix(
#'     rnorm(40),
#'     nrow = 5,
#'     dimnames = list(
#'       paste0("x", 1:5),
#'       paste0("s", 1:8)
#'     )
#'   )
#'
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), length.out = 8))
#'     )
#'   )
#'
#'   umap.result <- umap_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     n_neighbors = 3,
#'     seed = 42
#'   )
#'
#'   umap_pairs_kmeans_plot(
#'     dataset,
#'     umap.result,
#'     num.clusters = 2,
#'     dims = 1:3
#'   )
#' }
#'
#' @export
umap_pairs_kmeans_plot = function(dataset, umap.result, num.clusters = 3, kmeans.result = NULL, dims = 1:min(5, ncol(umap.result$embedding))) {
  if (!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required. Install it with: install.packages('GGally')")
  if (is.null(kmeans.result)) kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  emb = umap.result$embedding
  pairs.df = data.frame(emb[, dims])
  pairs.df$group = factor(kmeans.result$cluster)
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color = group))
}
