############################################################################
################################ t-SNE #####################################
############################################################################
# perform t-SNE analysis 
# n_components - number of dimensions in the embedding (default: 2) 
# perplexity - balances local vs global structure (typical range: 5-50) 
# max_iter - maximum number of iterations (default: 1000) 
# theta - speed/accuracy trade-off (0 = exact t-SNE, default: 0.5) 
# eta - learning rate (default: "auto") 
# scale - if TRUE, scales data before t-SNE 
# seed - random seed for reproducibility 
# write.file - if TRUE, saves the embedding to a CSV file 
# file.out - base name for output file


#' t-SNE analysis
#'
#' Performs t-distributed Stochastic Neighbor Embedding (t-SNE)
#' dimensionality reduction on a specmine dataset.
#'
#' @param dataset Dataset to analyse.
#' @param n_components Number of dimensions in the embedding.
#' @param perplexity Perplexity parameter controlling the balance between
#'   local and global structure. It must be sufficiently small relative to
#'   the number of samples.
#' @param max_iter Maximum number of optimization iterations.
#' @param theta Speed-accuracy trade-off parameter. A value of zero performs
#'   exact t-SNE, whereas larger values are faster but approximate.
#' @param eta Learning rate. The default value \code{"auto"} is converted to
#'   a numeric value based on the number of samples.
#' @param scale Logical indicating whether the data should be scaled before
#'   t-SNE.
#' @param seed Random seed used for reproducibility.
#' @param write.file Logical indicating whether the embedding should be
#'   written to a CSV file.
#' @param file.out Output file prefix used when
#'   \code{write.file = TRUE}.
#' @param ... Additional arguments passed to \code{Rtsne::Rtsne()}.
#'
#' @return A list containing the t-SNE embedding and the analysis parameters.
#'
#' @examples
#' if (requireNamespace("Rtsne", quietly = TRUE)) {
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
#'   tsne.result <- tsne_analysis_dataset(
#'     dataset,
#'     n_components = 2,
#'     perplexity = 2,
#'     max_iter = 250,
#'     seed = 42
#'   )
#' }
#'
#' @export

tsne_analysis_dataset = function(dataset, n_components = 2, perplexity = 30, max_iter = 1000, theta = 0.5, eta = "auto", scale = FALSE, seed = 42, write.file = FALSE, file.out = NULL, ...) {
  if (!requireNamespace("Rtsne", quietly = TRUE)) stop("Package 'Rtsne' is required. Install it with: install.packages('Rtsne')")
  mat_check = as.matrix(dataset$data)
  if (any(is.na(mat_check)) || any(is.nan(mat_check)) || any(is.infinite(mat_check))) {
    stop("dataset$data contains NA, NaN or Inf values. Please clean your data first.")
  }
  mat = t(mat_check)
  if (scale) mat = scale(mat)
  if (any(duplicated(mat))) {
    warning("Duplicate rows detected. Adding small jitter to avoid t-SNE errors.")
    set.seed(seed)
    mat = mat + matrix(rnorm(nrow(mat) * ncol(mat), 0, 1e-6), nrow = nrow(mat))
  }
  max_perplexity = floor((nrow(mat) - 1) / 3)
  if (perplexity > max_perplexity) {
    warning(paste("perplexity", perplexity, "is too large for", nrow(mat), "samples. Setting perplexity to", max_perplexity))
    perplexity = max_perplexity
  }
  if (is.character(eta) && eta == "auto") eta = max(200, nrow(mat) / 12)
  set.seed(seed)
  tsne_model = Rtsne::Rtsne(mat, dims = n_components, perplexity = perplexity, max_iter = max_iter, theta = theta, eta = eta, check_duplicates = FALSE, ...)
  embedding = tsne_model$Y
  rownames(embedding) = colnames(dataset$data)
  colnames(embedding) = paste0("tSNE", seq_len(n_components))
  if (isTRUE(write.file)) {
    if (is.null(file.out) || !nzchar(file.out)) {
      stop("Please provide 'file.out' when write.file = TRUE.")
    }
    utils::write.csv(embedding, file = paste0(file.out, "_embedding.csv"))
  }
  result = list(embedding = embedding, params = list(n_components = n_components, perplexity = perplexity, max_iter = max_iter, theta = theta, eta = eta, scale = scale, seed = seed))
  return(result)
}

########################## t-SNE PLOTS ##################################

#' t-SNE 2D scores plot
#'
#' Creates a two-dimensional scores plot from a t-SNE result.
#'
#' @param dataset Dataset used to calculate the embedding.
#' @param tsne.result Result returned by \code{tsne_analysis_dataset()}.
#' @param column.class Optional metadata column used to colour or group
#'   samples.
#' @param dims Two embedding dimensions to plot.
#' @param labels Logical indicating whether sample labels should be displayed.
#' @param ellipses Logical indicating whether group ellipses should be
#'   displayed.
#' @param bw Logical indicating whether a black-and-white plot should be used.
#' @param pallette RColorBrewer palette number.
#' @param leg.pos Legend position.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' if (requireNamespace("Rtsne", quietly = TRUE)) {
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
#'   tsne.result <- tsne_analysis_dataset(
#'     dataset,
#'     perplexity = 2,
#'     max_iter = 250,
#'     seed = 42
#'   )
#'
#'   tsne_scoresplot2D(
#'     dataset,
#'     tsne.result,
#'     column.class = "class"
#'   )
#' }
#'
#' @export
tsne_scoresplot2D = function(dataset, tsne.result, column.class = NULL, dims = c(1,2), labels = FALSE, ellipses = FALSE, bw = FALSE, pallette = 2, leg.pos = "right", xlim = NULL, ylim = NULL) {
  has.legend = FALSE
  emb = tsne.result$embedding
  if (ncol(emb) < 2) stop("Embedding must have at least 2 components for a 2D plot.")
  tsne.points = data.frame(emb[, dims])
  names(tsne.points) = c("x", "y")
  if (is.null(column.class)) {
    group.values = factor(rep(4, ncol(dataset$data)))
  } else {
    group.values = dataset$metadata[, column.class]
    has.legend = TRUE
  }
  tsne.points$group = group.values
  tsne.points$label = colnames(dataset$data)
  if (bw) shape.values = 1:length(levels(group.values))
  if (bw) tsne.plot = ggplot2::ggplot(data = tsne.points, ggplot2::aes_string(x = 'x', y = 'y', shape = 'group'))
  else tsne.plot = ggplot2::ggplot(data = tsne.points, ggplot2::aes_string(x = 'x', y = 'y', colour = 'group'))
  tsne.plot = tsne.plot + ggplot2::geom_point(size = 3, alpha = 1)
  if (bw) tsne.plot = tsne.plot + ggplot2::scale_shape_manual(values = shape.values)
  else tsne.plot = tsne.plot + ggplot2::scale_colour_brewer(type = "qual", palette = pallette)
  tsne.plot = tsne.plot + ggplot2::xlab(paste0("tSNE", dims[1])) + ggplot2::ylab(paste0("tSNE", dims[2])) + ggplot2::ggtitle("t-SNE 2D Scores Plot")
  if (has.legend) {
    if (bw) tsne.plot = tsne.plot + ggplot2::theme_bw()
    else tsne.plot = tsne.plot + ggplot2::theme(legend.position = leg.pos)
  }
  if (!is.null(xlim)) tsne.plot = tsne.plot + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) tsne.plot = tsne.plot + ggplot2::ylim(ylim[1], ylim[2])
  if (labels) tsne.plot = tsne.plot + ggplot2::geom_text(data = tsne.points, ggplot2::aes_string(x = 'x', y = 'y', label = 'label'), hjust = -0.1, vjust = 0)
  if (!bw & ellipses) {
    df.ellipses = calculate_ellipses(tsne.points)
    tsne.plot = tsne.plot + ggplot2::geom_path(data = df.ellipses, ggplot2::aes_string(x = 'x', y = 'y', colour = 'group'), size = 1, linetype = 2)
  }
  tsne.plot
}

#' t-SNE 3D scores plot
#'
#' Creates an interactive three-dimensional scores plot from a t-SNE result.
#'
#' @param dataset Dataset used to calculate the embedding.
#' @param tsne.result Result returned by \code{tsne_analysis_dataset()}.
#' @param column.class Optional metadata column used to colour samples.
#' @param dims Three embedding dimensions to plot.
#' @param title Plot title.
#'
#' @return A \code{plotly} htmlwidget.
#'
#' @examples
#' if (requireNamespace("Rtsne", quietly = TRUE) &&
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
#'   tsne.result <- tsne_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     perplexity = 2,
#'     max_iter = 250,
#'     seed = 42
#'   )
#'
#'   tsne_scoresplot3D(
#'     dataset,
#'     tsne.result,
#'     column.class = "class"
#'   )
#' }
#'
#' @export
tsne_scoresplot3D = function(dataset, tsne.result, column.class = NULL, dims = c(1,2,3), title = "t-SNE 3D Scores Plot") {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required. Install it with: install.packages('plotly')")
  emb = tsne.result$embedding
  if (ncol(emb) < 3) stop("Embedding must have at least 3 components for a 3D plot.")
  df = as.data.frame(emb[, dims])
  colnames(df) = c("tSNE1", "tSNE2", "tSNE3")
  if (!is.null(column.class) && column.class %in% colnames(dataset$metadata)) {
    df$Group = as.factor(dataset$metadata[, column.class])
  } else {
    df$Group = factor(rep("Samples", nrow(df)))
  }
  plotly::layout(
    plotly::plot_ly(
      df,
      x = ~tSNE1,
      y = ~tSNE2,
      z = ~tSNE3,
      color = ~Group,
      type = "scatter3d",
      mode = "markers"
    ),
    title = title
  )
}

#' t-SNE pairs plot
#'
#' Creates a pairs plot for selected t-SNE dimensions.
#'
#' @param dataset Dataset used to calculate the embedding.
#' @param tsne.result Result returned by \code{tsne_analysis_dataset()}.
#' @param column.class Optional metadata column used to colour samples.
#' @param dims Embedding dimensions to include in the pairs plot.
#' @param ... Additional arguments passed to \code{GGally::ggpairs()}.
#'
#' @return A \code{GGally} pairs plot object.
#'
#' @examples
#' if (requireNamespace("Rtsne", quietly = TRUE) &&
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
#'   tsne.result <- tsne_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     perplexity = 2,
#'     max_iter = 250,
#'     seed = 42
#'   )
#'
#'   tsne_pairs_plot(
#'     dataset,
#'     tsne.result,
#'     column.class = "class",
#'     dims = 1:3
#'   )
#' }
#'
#' @export
tsne_pairs_plot = function(dataset, tsne.result, column.class = NULL, dims = 1:min(5, ncol(tsne.result$embedding)), ...) {
  if (!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required. Install it with: install.packages('GGally')")
  emb = tsne.result$embedding
  if (is.null(column.class)) {
    group.values = rep(4, ncol(dataset$data))
  } else {
    group.values = dataset$metadata[, column.class]
  }
  pairs.df = data.frame(emb[, dims])
  pairs.df$group = group.values
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color = group), ...)
}

#' t-SNE 2D k-means plot
#'
#' Creates a two-dimensional t-SNE plot coloured by k-means clusters.
#'
#' @param dataset Dataset to cluster.
#' @param tsne.result Result returned by \code{tsne_analysis_dataset()}.
#' @param num.clusters Number of k-means clusters.
#' @param dims Two embedding dimensions to plot.
#' @param kmeans.result Optional result returned by \code{clustering()}.
#' @param labels Logical indicating whether sample labels should be displayed.
#' @param bw Logical indicating whether a black-and-white plot should be used.
#' @param ellipses Logical indicating whether group ellipses should be
#'   displayed.
#' @param leg.pos Legend position.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' if (requireNamespace("Rtsne", quietly = TRUE)) {
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
#'   tsne.result <- tsne_analysis_dataset(
#'     dataset,
#'     perplexity = 2,
#'     max_iter = 250,
#'     seed = 42
#'   )
#'
#'   kmeans.result <- list(
#'     cluster = rep(1:2, each = 4)
#'   )
#'
#'   tsne_kmeans_plot2D(
#'     dataset,
#'     tsne.result,
#'     num.clusters = 2,
#'     kmeans.result = kmeans.result
#'   )
#' }
#'
#' @export
tsne_kmeans_plot2D = function(dataset, tsne.result, num.clusters = 3, dims = c(1,2), kmeans.result = NULL, labels = FALSE, bw = FALSE, ellipses = FALSE, leg.pos = "right", xlim = NULL, ylim = NULL) {
  emb = tsne.result$embedding
  if (is.null(kmeans.result)) kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  tsne.points = data.frame(emb[, dims])
  names(tsne.points) = c("x", "y")
  tsne.points$group = factor(kmeans.result$cluster)
  tsne.points$label = colnames(dataset$data)
  if (bw) shape.values = 1:num.clusters
  if (bw) tsne.plot = ggplot2::ggplot(data = tsne.points, ggplot2::aes_string(x = 'x', y = 'y', shape = 'group'))
  else tsne.plot = ggplot2::ggplot(data = tsne.points, ggplot2::aes_string(x = 'x', y = 'y', colour = 'group'))
  tsne.plot = tsne.plot + ggplot2::geom_point(size = 3, alpha = .6)
  if (bw) tsne.plot = tsne.plot + ggplot2::scale_shape_manual(values = shape.values)
  else tsne.plot = tsne.plot + ggplot2::scale_colour_brewer(palette = "Set1")
  tsne.plot = tsne.plot + ggplot2::xlab(paste0("tSNE", dims[1])) + ggplot2::ylab(paste0("tSNE", dims[2])) + ggplot2::ggtitle("t-SNE 2D K-means Plot")
  if (bw) tsne.plot = tsne.plot + ggplot2::theme_bw()
  else tsne.plot = tsne.plot + ggplot2::theme(legend.position = leg.pos)
  if (!is.null(xlim)) tsne.plot = tsne.plot + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) tsne.plot = tsne.plot + ggplot2::ylim(ylim[1], ylim[2])
  if (labels) tsne.plot = tsne.plot + ggplot2::geom_text(data = tsne.points, ggplot2::aes_string(x = 'x', y = 'y', label = 'label'), hjust = -0.1, vjust = 0, size = 3)
  if (!bw & ellipses) {
    df.ellipses = calculate_ellipses(tsne.points)
    tsne.plot = tsne.plot + ggplot2::geom_path(data = df.ellipses, ggplot2::aes_string(x = 'x', y = 'y', colour = 'group'), size = 1, linetype = 2)
  }
  tsne.plot
}

#' t-SNE 3D k-means plot
#'
#' Creates an interactive three-dimensional t-SNE plot coloured by k-means
#' clusters.
#'
#' @param dataset Dataset to cluster.
#' @param tsne.result Result returned by \code{tsne_analysis_dataset()}.
#' @param num.clusters Number of k-means clusters.
#' @param dims Three embedding dimensions to plot.
#' @param kmeans.result Optional result returned by \code{clustering()}.
#' @param title Plot title.
#'
#' @return A \code{plotly} htmlwidget.
#'
#' @examples
#' if (requireNamespace("Rtsne", quietly = TRUE) &&
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
#'   tsne.result <- tsne_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     perplexity = 2,
#'     max_iter = 250,
#'     seed = 42
#'   )
#'
#'   kmeans.result <- list(
#'     cluster = rep(1:2, each = 4)
#'   )
#'
#'   tsne_kmeans_plot3D(
#'     dataset,
#'     tsne.result,
#'     num.clusters = 2,
#'     kmeans.result = kmeans.result
#'   )
#' }
#'
#' @export
tsne_kmeans_plot3D = function(dataset, tsne.result, num.clusters = 3, dims = c(1,2,3), kmeans.result = NULL, title = "t-SNE 3D K-means Plot") {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required. Install it with: install.packages('plotly')")
  emb = tsne.result$embedding
  if (ncol(emb) < 3) stop("Embedding must have at least 3 components for a 3D plot.")
  if (is.null(kmeans.result)) kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  df = as.data.frame(emb[, dims])
  colnames(df) = c("tSNE1", "tSNE2", "tSNE3")
  df$Group = factor(kmeans.result$cluster)
  plotly::layout(
    plotly::plot_ly(
      df,
      x = ~tSNE1,
      y = ~tSNE2,
      z = ~tSNE3,
      color = ~Group,
      type = "scatter3d",
      mode = "markers"
    ),
    title = title
  )
}

#' t-SNE pairs k-means plot
#'
#' Creates a pairs plot for t-SNE dimensions coloured by k-means clusters.
#'
#' @param dataset Dataset to cluster.
#' @param tsne.result Result returned by \code{tsne_analysis_dataset()}.
#' @param num.clusters Number of k-means clusters.
#' @param kmeans.result Optional result returned by \code{clustering()}.
#' @param dims Embedding dimensions to include in the pairs plot.
#'
#' @return A \code{GGally} pairs plot object.
#'
#' @examples
#' if (requireNamespace("Rtsne", quietly = TRUE) &&
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
#'   tsne.result <- tsne_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     perplexity = 2,
#'     max_iter = 250,
#'     seed = 42
#'   )
#'
#'   kmeans.result <- list(
#'     cluster = rep(1:2, each = 4)
#'   )
#'
#'   tsne_pairs_kmeans_plot(
#'     dataset,
#'     tsne.result,
#'     num.clusters = 2,
#'     kmeans.result = kmeans.result,
#'     dims = 1:3
#'   )
#' }
#'
#' @export
tsne_pairs_kmeans_plot = function(dataset, tsne.result, num.clusters = 3, kmeans.result = NULL, dims = 1:min(5, ncol(tsne.result$embedding))) {
  if (!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required. Install it with: install.packages('GGally')")
  if (is.null(kmeans.result)) kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  emb = tsne.result$embedding
  pairs.df = data.frame(emb[, dims])
  pairs.df$group = factor(kmeans.result$cluster)
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color = group))
}