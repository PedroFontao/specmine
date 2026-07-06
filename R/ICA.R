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


ica_analysis_dataset = function(dataset, n_components = 2, alg.typ = "parallel", fun = "logcosh", maxit = 200, tol = 1e-4, scale = FALSE, seed = 42, write.file = FALSE, file.out = "ica", ...) {
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
  if (write.file) write.csv(embedding, file = paste0(file.out, "_components.csv"))
  result = list(
    embedding  = embedding,
    S          = ica_model$S,
    A          = ica_model$A,
    K          = ica_model$K,
    W          = ica_model$W,
    loadings   = t(ica_model$A),
    params     = list(n_components = n_components, alg.typ = alg.typ, fun = fun, maxit = maxit, tol = tol, scale = scale, seed = seed)
  )
  return(result)
}



########################## ICA PLOTS ##################################



# 2d scores plot
# column.class - name of metadata column to colour points by
# dims         - which IC dimensions to plot (default: c(1,2))
# labels       - if TRUE, shows sample names next to points
# ellipses     - if TRUE, draws confidence ellipses per group
# bw           - if TRUE, uses black and white with shapes instead of colours
# palette      - colour palette index for scale_colour_brewer
# leg.pos      - legend position ("right", "left", "top", "bottom")
# xlim         - optional x axis limits c(min, max)
# ylim         - optional y axis limits c(min, max)


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



# 3d scores plot (uses plotly)
# column.class - name of metadata column to colour points by
# dims         - which IC dimensions to plot (default: c(1,2,3))
# title        - plot title


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
  plotly::plot_ly(df, x = ~IC1, y = ~IC2, z = ~IC3, color = ~Group, type = "scatter3d", mode = "markers") %>%
    plotly::layout(title = title)
}



# pairs plot for multiple IC components
# column.class - name of metadata column to colour points by
# dims         - which IC dimensions to include (default: first 5)


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



# loadings plot: variable contributions to each IC
# dims         - which ICs to plot loadings for (default: c(1,2))
# top.n        - if set, shows only the top N variables by absolute loading
# labels       - if TRUE, shows variable names (use only with few variables or top.n)


ica_loadingsplot = function(ica.result, dims = c(1,2), top.n = NULL, labels = FALSE) {
  loadings.mat = ica.result$loadings
  if (any(dims > nrow(loadings.mat))) stop(paste("dims out of range: only", nrow(loadings.mat), "components available."))
  df = data.frame(
    variable = colnames(ica.result$S),
    IC_x     = loadings.mat[dims[1], ],
    IC_y     = loadings.mat[dims[2], ]
  )
  if (is.null(df$variable)) df$variable = seq_len(ncol(loadings.mat))
  if (!is.null(top.n)) {
    importance = sqrt(df$IC_x^2 + df$IC_y^2)
    df = df[order(importance, decreasing = TRUE)[1:min(top.n, nrow(df))], ]
  }
  p = ggplot2::ggplot(df, ggplot2::aes(x = .data[["IC_x"]], y = .data[["IC_y"]])) +
    ggplot2::geom_point(colour = "steelblue", alpha = 0.7, size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::xlab(paste0("IC", dims[1], " Loading")) +
    ggplot2::ylab(paste0("IC", dims[2], " Loading")) +
    ggplot2::ggtitle("ICA Loadings Plot") +
    ggplot2::theme_bw()
  if (labels) p = p + ggplot2::geom_text(ggplot2::aes(label = .data[["variable"]]), hjust = -0.1, vjust = 0, size = 3)
  p
}



# kmeans clustering with 2 IC dimensions
# num.clusters  - number of k-means clusters
# dims          - which IC dimensions to plot (default: c(1,2))
# kmeans.result - optional pre-computed clustering result (from clustering())
# use.embedding - if TRUE, runs k-means on the IC embedding instead of raw data
# labels        - if TRUE, shows sample names next to points
# bw            - if TRUE, uses black and white with shapes instead of colours
# ellipses      - if TRUE, draws confidence ellipses per group
# leg.pos       - legend position
# xlim          - optional x axis limits c(min, max)
# ylim          - optional y axis limits c(min, max)


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



# kmeans clustering with 3 IC dimensions (uses plotly)
# num.clusters  - number of k-means clusters
# dims          - which IC dimensions to plot (default: c(1,2,3))
# kmeans.result - optional pre-computed clustering result (from clustering())
# use.embedding - if TRUE, runs k-means on the IC embedding instead of raw data
# title         - plot title


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
  plotly::plot_ly(df, x = ~IC1, y = ~IC2, z = ~IC3, color = ~Group, type = "scatter3d", mode = "markers") %>%
    plotly::layout(title = title)
}



# pairs plot with kmeans clusters
# num.clusters  - number of k-means clusters
# dims          - which IC dimensions to include (default: first 5)
# kmeans.result - optional pre-computed clustering result (from clustering())
# use.embedding - if TRUE, runs k-means on the IC embedding instead of raw data


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