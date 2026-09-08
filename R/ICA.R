############################################################################
################################ ICA #######################################
############################################################################

# perform ICA analysis
# n_components  - number of independent components to extract (default: 2)
# alg.typ       - algorithm type: "parallel" or "deflation" (default: "parallel")
# fun           - contrast function: "logcosh", "exp" (default: "logcosh")
# maxit         - maximum number of iterations (default: 200)
# tol           - convergence tolerance (default: 1e-4)
# scale         - if TRUE, scales data before ICA (recommended)
# seed          - random seed for reproducibility
# write.file    - if TRUE, saves the components to a CSV file
# file.out      - base name for output file

#' ICA analysis
#'
#' Performs Independent Component Analysis (ICA) on a specmine dataset using
#' the FastICA algorithm.
#'
#' @param dataset Dataset to analyse.
#' @param n_components Number of independent components to extract.
#' @param alg.typ Algorithm used by FastICA. Supported values include
#'   \code{"parallel"} and \code{"deflation"}.
#' @param fun Contrast function used by FastICA. Supported values include
#'   \code{"logcosh"} and \code{"exp"}.
#' @param maxit Maximum number of iterations.
#' @param tol Convergence tolerance.
#' @param scale Logical indicating whether the data should be scaled before
#'   ICA.
#' @param seed Random seed used for reproducibility.
#' @param write.file Logical indicating whether the independent components
#'   should be written to a CSV file.
#' @param file.out Output file prefix used when
#'   \code{write.file = TRUE}.
#' @param ... Additional arguments passed to \code{fastICA::fastICA()}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{embedding}: the independent component scores;
#'   \item \code{S}: the estimated source matrix;
#'   \item \code{A}: the estimated mixing matrix;
#'   \item \code{K}: the whitening matrix;
#'   \item \code{W}: the estimated unmixing matrix;
#'   \item \code{loadings}: the component loadings;
#'   \item \code{params}: the analysis parameters.
#' }
#'
#' @examples
#' if (requireNamespace("fastICA", quietly = TRUE)) {
#'   datamat <- matrix(
#'     rnorm(240),
#'     nrow = 8,
#'     dimnames = list(
#'       paste0("feature", 1:8),
#'       paste0("sample", 1:30)
#'     )
#'   )
#'
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), each = 15))
#'     )
#'   )
#'
#'   ica.result <- ica_analysis_dataset(
#'     dataset,
#'     n_components = 2,
#'     maxit = 100,
#'     seed = 42
#'   )
#'
#'   dim(ica.result$embedding)
#' }
#'
#' @export
ica_analysis_dataset = function(dataset, n_components = 2, alg.typ = "parallel", fun = "logcosh", maxit = 200, tol = 1e-4, scale = FALSE, seed = 42, write.file = FALSE, file.out = NULL, ...) {
  if (!requireNamespace("fastICA", quietly = TRUE)) stop("Package 'fastICA' is required. Install it with: install.packages('fastICA')")
  mat_check = as.matrix(dataset$data)
  if (any(is.na(mat_check)) || any(is.nan(mat_check)) || any(is.infinite(mat_check))) {
    stop("dataset$data contains NA, NaN or Inf values. Please clean your data first.")
  }
  mat = t(mat_check)
  if (scale) mat = scale(mat)
  if (n_components > ncol(mat)) {
    stop(paste("n_components (", n_components, ") cannot exceed the number of variables (", ncol(mat), ")."))
  }
  if (n_components > nrow(mat)) {
    stop(paste("n_components (", n_components, ") cannot exceed the number of samples (", nrow(mat), ")."))
  }
  set.seed(seed)
  ica_model = fastICA::fastICA(mat, n.comp = n_components, alg.typ = alg.typ, fun = fun, maxit = maxit, tol = tol, ...)
  embedding = ica_model$S
  
  rownames(embedding) = colnames(dataset$data)
  colnames(embedding) = paste0("IC", seq_len(n_components))
  
  colnames(ica_model$S) = colnames(embedding)
  colnames(ica_model$A) = rownames(dataset$data)
  
  if (isTRUE(write.file)) {
    if (is.null(file.out) || !nzchar(file.out)) {
      stop("Please provide 'file.out' when write.file = TRUE.")
    }
    
    utils::write.csv(
      embedding,
      file = paste0(file.out, "_components.csv")
    )
  }
  
  result = list(
    embedding = embedding,
    S = ica_model$S,
    A = ica_model$A,
    K = ica_model$K,
    W = ica_model$W,
    loadings = ica_model$A,
    params = list(
      n_components = n_components,
      alg.typ = alg.typ,
      fun = fun,
      maxit = maxit,
      tol = tol,
      scale = scale,
      seed = seed
    )
  )
  return(result)
}

########################## ICA PLOTS ##################################

#' ICA 2D scores plot
#'
#' Creates a two-dimensional scores plot from an ICA result.
#'
#' @param dataset Dataset used to calculate the independent components.
#' @param ica.result Result returned by \code{ica_analysis_dataset()}.
#' @param column.class Optional metadata column used to colour or group
#'   samples.
#' @param dims Two independent components to plot.
#' @param labels Logical indicating whether sample labels should be displayed.
#' @param ellipses Logical indicating whether group ellipses should be
#'   displayed.
#' @param bw Logical indicating whether a black-and-white plot should be used.
#' @param palette RColorBrewer palette number.
#' @param leg.pos Legend position.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' if (requireNamespace("fastICA", quietly = TRUE)) {
#'   datamat <- matrix(rnorm(240), nrow = 8)
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), each = 15))
#'     )
#'   )
#'
#'   ica.result <- ica_analysis_dataset(
#'     dataset,
#'     n_components = 2,
#'     maxit = 100,
#'     seed = 42
#'   )
#'
#'   ica_scoresplot2D(
#'     dataset,
#'     ica.result,
#'     column.class = "class"
#'   )
#' }
#'
#' @export
ica_scoresplot2D = function(dataset, ica.result, column.class = NULL, dims = c(1,2), labels = FALSE, ellipses = FALSE, bw = FALSE, palette = 2, leg.pos = "right", xlim = NULL, ylim = NULL) {
  has.legend = FALSE
  emb = ica.result$embedding
  if (ncol(emb) < 2) stop("Embedding must have at least 2 components for a 2D plot.")
  if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
  ica.points = data.frame(emb[, dims])
  names(ica.points) = c("x", "y")
  if (is.null(column.class)) {
    group.values = factor(rep(4, ncol(dataset$data)))
  } else {
    group.values = as.factor(dataset$metadata[, column.class])
    has.legend = TRUE
  }
  ica.points$group = group.values
  ica.points$label = colnames(dataset$data)
  if (bw) shape.values = 1:length(levels(group.values))
  if (bw) ica.plot = ggplot2::ggplot(data = ica.points, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], shape = .data[["group"]]))
  else ica.plot = ggplot2::ggplot(data = ica.points, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["group"]]))
  ica.plot = ica.plot + ggplot2::geom_point(size = 3, alpha = 1)
  if (bw) ica.plot = ica.plot + ggplot2::scale_shape_manual(values = shape.values)
  else ica.plot = ica.plot + ggplot2::scale_colour_brewer(type = "qual", palette = palette)
  ica.plot = ica.plot + ggplot2::xlab(paste0("IC", dims[1])) + ggplot2::ylab(paste0("IC", dims[2])) + ggplot2::ggtitle("ICA 2D Scores Plot")
  ica.plot = ica.plot + ggplot2::theme_bw()
  if (has.legend) ica.plot = ica.plot + ggplot2::theme(legend.position = leg.pos)
  if (!is.null(xlim)) ica.plot = ica.plot + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) ica.plot = ica.plot + ggplot2::ylim(ylim[1], ylim[2])
  if (labels) ica.plot = ica.plot + ggplot2::geom_text(data = ica.points, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]), hjust = -0.1, vjust = 0)
  if (!bw & ellipses) {
    df.ellipses = calculate_ellipses(ica.points)
    ica.plot = ica.plot + ggplot2::geom_path(data = df.ellipses, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["group"]]), linewidth = 1, linetype = 2)
  }
  ica.plot
}

#' ICA 3D scores plot
#'
#' Creates an interactive three-dimensional scores plot from an ICA result.
#'
#' @param dataset Dataset used to calculate the independent components.
#' @param ica.result Result returned by \code{ica_analysis_dataset()}.
#' @param column.class Optional metadata column used to colour samples.
#' @param dims Three independent components to plot.
#' @param title Plot title.
#'
#' @return A \code{plotly} htmlwidget.
#'
#' @examples
#' if (requireNamespace("fastICA", quietly = TRUE) &&
#'     requireNamespace("plotly", quietly = TRUE)) {
#'   datamat <- matrix(rnorm(240), nrow = 8)
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), each = 15))
#'     )
#'   )
#'
#'   ica.result <- ica_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     maxit = 100,
#'     seed = 42
#'   )
#'
#'   ica_scoresplot3D(
#'     dataset,
#'     ica.result,
#'     column.class = "class"
#'   )
#' }
#'
#' @export
ica_scoresplot3D = function(dataset, ica.result, column.class = NULL, dims = c(1,2,3), title = "ICA 3D Scores Plot") {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required. Install it with: install.packages('plotly')")
  emb = ica.result$embedding
  if (ncol(emb) < 3) stop("Embedding must have at least 3 components for a 3D plot.")
  if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
  df = as.data.frame(emb[, dims])
  colnames(df) = c("IC1", "IC2", "IC3")
  if (!is.null(column.class) && column.class %in% colnames(dataset$metadata)) {
    df$Group = as.factor(dataset$metadata[, column.class])
  } else {
    df$Group = factor(rep("Samples", nrow(df)))
  }
  plotly::layout(
    plotly::plot_ly(
      df,
      x = ~IC1,
      y = ~IC2,
      z = ~IC3,
      color = ~Group,
      type = "scatter3d",
      mode = "markers"
    ),
    title = title
  )
}

#' ICA pairs plot
#'
#' Creates a pairs plot for selected independent components.
#'
#' @param dataset Dataset used to calculate the independent components.
#' @param ica.result Result returned by \code{ica_analysis_dataset()}.
#' @param column.class Optional metadata column used to colour samples.
#' @param dims Independent components to include in the pairs plot.
#' @param ... Additional arguments passed to \code{GGally::ggpairs()}.
#'
#' @return A \code{GGally} pairs plot object.
#'
#' @examples
#' if (requireNamespace("fastICA", quietly = TRUE) &&
#'     requireNamespace("GGally", quietly = TRUE)) {
#'   datamat <- matrix(rnorm(240), nrow = 8)
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), each = 15))
#'     )
#'   )
#'
#'   ica.result <- ica_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     maxit = 100,
#'     seed = 42
#'   )
#'
#'   ica_pairs_plot(
#'     dataset,
#'     ica.result,
#'     column.class = "class",
#'     dims = 1:3
#'   )
#' }
#'
#' @export
ica_pairs_plot = function(dataset, ica.result, column.class = NULL, dims = 1:min(5, ncol(ica.result$embedding)), ...) {
  if (!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required. Install it with: install.packages('GGally')")
  emb = ica.result$embedding
  if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
  if (is.null(column.class)) {
    group.values = rep(4, ncol(dataset$data))
  } else {
    group.values = as.factor(dataset$metadata[, column.class])
  }
  pairs.df = data.frame(emb[, dims])
  pairs.df$group = group.values
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color = group), ...)
}

#' ICA loadings plot
#'
#' Creates a two-dimensional plot of the loadings associated with two
#' independent components.
#'
#' @param ica.result Result returned by \code{ica_analysis_dataset()}.
#' @param dims Two independent components whose loadings should be plotted.
#' @param top.n Optional number of variables with the largest loading
#'   magnitude to display.
#' @param labels Logical indicating whether variable labels should be
#'   displayed.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' if (requireNamespace("fastICA", quietly = TRUE)) {
#'   datamat <- matrix(
#'     rnorm(240),
#'     nrow = 8,
#'     dimnames = list(
#'       paste0("feature", 1:8),
#'       paste0("sample", 1:30)
#'     )
#'   )
#'
#'   dataset <- list(data = datamat)
#'
#'   ica.result <- ica_analysis_dataset(
#'     dataset,
#'     n_components = 2,
#'     maxit = 100,
#'     seed = 42
#'   )
#'
#'   ica_loadingsplot(
#'     ica.result,
#'     dims = c(1, 2),
#'     labels = TRUE
#'   )
#' }
#'
#' @export
ica_loadingsplot = function(
    ica.result,
    dims = c(1, 2),
    top.n = NULL,
    labels = FALSE
) {
  loadings.mat = ica.result$loadings
  
  if (length(dims) != 2) {
    stop("'dims' must contain exactly two component indices.")
  }
  
  if (any(dims < 1) || any(dims > nrow(loadings.mat))) {
    stop(
      paste(
        "dims out of range: only",
        nrow(loadings.mat),
        "components available."
      )
    )
  }
  
  variable.names = colnames(loadings.mat)
  
  if (is.null(variable.names) ||
      length(variable.names) != ncol(loadings.mat)) {
    variable.names = seq_len(ncol(loadings.mat))
  }
  
  df = data.frame(
    variable = variable.names,
    IC_x = as.numeric(loadings.mat[dims[1], ]),
    IC_y = as.numeric(loadings.mat[dims[2], ]),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(top.n)) {
    if (length(top.n) != 1 || top.n < 1) {
      stop("'top.n' must be a positive integer.")
    }
    
    top.n = min(as.integer(top.n), nrow(df))
    importance = sqrt(df$IC_x^2 + df$IC_y^2)
    selected = order(importance, decreasing = TRUE)[seq_len(top.n)]
    df = df[selected, , drop = FALSE]
  }
  
  p = ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["IC_x"]],
      y = .data[["IC_y"]]
    )
  ) +
    ggplot2::geom_point(
      colour = "steelblue",
      alpha = 0.7,
      size = 2
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    ggplot2::xlab(paste0("IC", dims[1], " Loading")) +
    ggplot2::ylab(paste0("IC", dims[2], " Loading")) +
    ggplot2::ggtitle("ICA Loadings Plot") +
    ggplot2::theme_bw()
  
  if (labels) {
    p = p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data[["variable"]]),
        hjust = -0.1,
        vjust = 0,
        size = 3
      )
  }
  
  p
}

#' ICA 2D k-means plot
#'
#' Creates a two-dimensional ICA plot coloured by k-means clusters.
#'
#' @param dataset Dataset to cluster.
#' @param ica.result Result returned by \code{ica_analysis_dataset()}.
#' @param num.clusters Number of k-means clusters.
#' @param dims Two independent components to plot.
#' @param kmeans.result Optional k-means result containing a \code{cluster}
#'   vector.
#' @param use.embedding Logical indicating whether k-means should be applied
#'   to the ICA embedding when \code{kmeans.result} is not supplied.
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
#' if (requireNamespace("fastICA", quietly = TRUE)) {
#'   datamat <- matrix(rnorm(240), nrow = 8)
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), each = 15))
#'     )
#'   )
#'
#'   ica.result <- ica_analysis_dataset(
#'     dataset,
#'     n_components = 2,
#'     maxit = 100,
#'     seed = 42
#'   )
#'
#'   kmeans.result <- list(
#'     cluster = rep(1:2, each = 15)
#'   )
#'
#'   ica_kmeans_plot2D(
#'     dataset,
#'     ica.result,
#'     num.clusters = 2,
#'     kmeans.result = kmeans.result
#'   )
#' }
#'
#' @export
ica_kmeans_plot2D = function(dataset, ica.result, num.clusters = 3, dims = c(1,2), kmeans.result = NULL, use.embedding = TRUE, labels = FALSE, bw = FALSE, ellipses = FALSE, leg.pos = "right", xlim = NULL, ylim = NULL) {
  emb = ica.result$embedding
  if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
  if (is.null(kmeans.result)) {
    if (use.embedding) {
      kmeans.result = kmeans(emb, centers = num.clusters, nstart = 25)
    } else {
      kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
    }
  }
  ica.points = data.frame(emb[, dims])
  names(ica.points) = c("x", "y")
  ica.points$group = factor(kmeans.result$cluster)
  ica.points$label = colnames(dataset$data)
  if (bw) shape.values = 1:num.clusters
  if (bw) ica.plot = ggplot2::ggplot(data = ica.points, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], shape = .data[["group"]]))
  else ica.plot = ggplot2::ggplot(data = ica.points, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["group"]]))
  ica.plot = ica.plot + ggplot2::geom_point(size = 3, alpha = 0.6)
  if (bw) ica.plot = ica.plot + ggplot2::scale_shape_manual(values = shape.values)
  else ica.plot = ica.plot + ggplot2::scale_colour_brewer(palette = "Set1")
  ica.plot = ica.plot + ggplot2::xlab(paste0("IC", dims[1])) + ggplot2::ylab(paste0("IC", dims[2])) + ggplot2::ggtitle("ICA 2D K-means Plot")
  ica.plot = ica.plot + ggplot2::theme_bw()
  ica.plot = ica.plot + ggplot2::theme(legend.position = leg.pos)
  if (!is.null(xlim)) ica.plot = ica.plot + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) ica.plot = ica.plot + ggplot2::ylim(ylim[1], ylim[2])
  if (labels) ica.plot = ica.plot + ggplot2::geom_text(data = ica.points, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]), hjust = -0.1, vjust = 0, size = 3)
  if (!bw & ellipses) {
    df.ellipses = calculate_ellipses(ica.points)
    ica.plot = ica.plot + ggplot2::geom_path(data = df.ellipses, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["group"]]), linewidth = 1, linetype = 2)
  }
  ica.plot
}

#' ICA 3D k-means plot
#'
#' Creates an interactive three-dimensional ICA plot coloured by k-means
#' clusters.
#'
#' @param dataset Dataset to cluster.
#' @param ica.result Result returned by \code{ica_analysis_dataset()}.
#' @param num.clusters Number of k-means clusters.
#' @param dims Three independent components to plot.
#' @param kmeans.result Optional k-means result containing a \code{cluster}
#'   vector.
#' @param use.embedding Logical indicating whether k-means should be applied
#'   to the ICA embedding when \code{kmeans.result} is not supplied.
#' @param title Plot title.
#'
#' @return A \code{plotly} htmlwidget.
#'
#' @examples
#' if (requireNamespace("fastICA", quietly = TRUE) &&
#'     requireNamespace("plotly", quietly = TRUE)) {
#'   datamat <- matrix(rnorm(240), nrow = 8)
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), each = 15))
#'     )
#'   )
#'
#'   ica.result <- ica_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     maxit = 100,
#'     seed = 42
#'   )
#'
#'   kmeans.result <- list(
#'     cluster = rep(1:2, each = 15)
#'   )
#'
#'   ica_kmeans_plot3D(
#'     dataset,
#'     ica.result,
#'     num.clusters = 2,
#'     kmeans.result = kmeans.result
#'   )
#' }
#'
#' @export
ica_kmeans_plot3D = function(dataset, ica.result, num.clusters = 3, dims = c(1,2,3), kmeans.result = NULL, use.embedding = TRUE, title = "ICA 3D K-means Plot") {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required. Install it with: install.packages('plotly')")
  emb = ica.result$embedding
  if (ncol(emb) < 3) stop("Embedding must have at least 3 components for a 3D plot.")
  if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
  if (is.null(kmeans.result)) {
    if (use.embedding) {
      kmeans.result = kmeans(emb, centers = num.clusters, nstart = 25)
    } else {
      kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
    }
  }
  df = as.data.frame(emb[, dims])
  colnames(df) = c("IC1", "IC2", "IC3")
  df$Group = factor(kmeans.result$cluster)
  plotly::layout(
    plotly::plot_ly(
      df,
      x = ~IC1,
      y = ~IC2,
      z = ~IC3,
      color = ~Group,
      type = "scatter3d",
      mode = "markers"
    ),
    title = title
  )
}

#' ICA pairs k-means plot
#'
#' Creates a pairs plot for independent components coloured by k-means
#' clusters.
#'
#' @param dataset Dataset to cluster.
#' @param ica.result Result returned by \code{ica_analysis_dataset()}.
#' @param num.clusters Number of k-means clusters.
#' @param kmeans.result Optional k-means result containing a \code{cluster}
#'   vector.
#' @param use.embedding Logical indicating whether k-means should be applied
#'   to the ICA embedding when \code{kmeans.result} is not supplied.
#' @param dims Independent components to include in the pairs plot.
#'
#' @return A \code{GGally} pairs plot object.
#'
#' @examples
#' if (requireNamespace("fastICA", quietly = TRUE) &&
#'     requireNamespace("GGally", quietly = TRUE)) {
#'   datamat <- matrix(rnorm(240), nrow = 8)
#'   dataset <- list(
#'     data = datamat,
#'     metadata = data.frame(
#'       class = factor(rep(c("A", "B"), each = 15))
#'     )
#'   )
#'
#'   ica.result <- ica_analysis_dataset(
#'     dataset,
#'     n_components = 3,
#'     maxit = 100,
#'     seed = 42
#'   )
#'
#'   kmeans.result <- list(
#'     cluster = rep(1:2, each = 15)
#'   )
#'
#'   ica_pairs_kmeans_plot(
#'     dataset,
#'     ica.result,
#'     num.clusters = 2,
#'     kmeans.result = kmeans.result,
#'     dims = 1:3
#'   )
#' }
#'
#' @export
ica_pairs_kmeans_plot = function(dataset, ica.result, num.clusters = 3, kmeans.result = NULL, use.embedding = TRUE, dims = 1:min(5, ncol(ica.result$embedding))) {
  if (!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required. Install it with: install.packages('GGally')")
  emb = ica.result$embedding
  if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
  if (is.null(kmeans.result)) {
    if (use.embedding) {
      kmeans.result = kmeans(emb, centers = num.clusters, nstart = 25)
    } else {
      kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
    }
  }
  pairs.df = data.frame(emb[, dims])
  pairs.df$group = factor(kmeans.result$cluster)
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color = group))
}