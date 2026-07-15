############################################################################
########################### DBSCAN / HDBSCAN ###############################
############################################################################



# perform DBSCAN clustering
# eps         - radius of the epsilon neighbourhood
# minPts      - minimum number of points required to form a dense region
# scale       - if TRUE, scales data before clustering
# seed        - random seed for reproducibility
# ret_model   - if TRUE, returns the fitted dbscan model
# write.file  - if TRUE, saves cluster assignments to a CSV file
# file.out    - base name for output file


dbscan_analysis_dataset = function(dataset, eps = 0.5, minPts = 5, scale = FALSE, seed = 42, ret_model = FALSE, write.file = FALSE, file.out = "dbscan", ...) {
  if (!requireNamespace("dbscan", quietly = TRUE)) stop("Package 'dbscan' is required. Install it with: install.packages('dbscan')")
  mat_check = as.matrix(dataset$data)
  if (any(is.na(mat_check)) || any(is.nan(mat_check)) || any(is.infinite(mat_check))) {
    stop("dataset$data contains NA, NaN or Inf values. Please clean your data first.")
  }
  mat = t(mat_check)
  if (scale) mat = base::scale(mat)
  set.seed(seed)
  dbscan_model = dbscan::dbscan(mat, eps = eps, minPts = minPts, ...)
  clusters = dbscan_model$cluster
  names(clusters) = colnames(dataset$data)
  result = list(
    cluster = clusters,
    params = list(eps = eps, minPts = minPts, scale = scale, seed = seed)
  )
  if (!is.null(dbscan_model$isseed)) result$isseed = dbscan_model$isseed
  if (ret_model) result$model = dbscan_model
  if (write.file) {
    out = data.frame(sample = names(clusters), cluster = clusters)
    write.csv(out, file = paste0(file.out, "_clusters.csv"), row.names = FALSE)
  }
  return(result)
}



# perform HDBSCAN clustering
# minPts      - minimum number of points required in the dense neighbourhood
# scale       - if TRUE, scales data before clustering
# seed        - random seed for reproducibility
# ret_model   - if TRUE, returns the fitted hdbscan model
# write.file  - if TRUE, saves cluster assignments to a CSV file
# file.out    - base name for output file


hdbscan_analysis_dataset = function(dataset, minPts = 5, scale = FALSE, seed = 42, ret_model = FALSE, write.file = FALSE, file.out = "hdbscan", ...) {
  if (!requireNamespace("dbscan", quietly = TRUE)) stop("Package 'dbscan' is required. Install it with: install.packages('dbscan')")
  mat_check = as.matrix(dataset$data)
  if (any(is.na(mat_check)) || any(is.nan(mat_check)) || any(is.infinite(mat_check))) {
    stop("dataset$data contains NA, NaN or Inf values. Please clean your data first.")
  }
  mat = t(mat_check)
  if (scale) mat = base::scale(mat)
  set.seed(seed)
  hdbscan_model = dbscan::hdbscan(mat, minPts = minPts, ...)
  clusters = hdbscan_model$cluster
  names(clusters) = colnames(dataset$data)
  result = list(
    cluster = clusters,
    membership_prob = hdbscan_model$membership_prob,
    outlier_scores = hdbscan_model$outlier_scores,
    hc = hdbscan_model$hc,
    params = list(minPts = minPts, scale = scale, seed = seed)
  )
  if (ret_model) result$model = hdbscan_model
  if (write.file) {
    out = data.frame(
      sample = names(clusters),
      cluster = clusters,
      membership_prob = hdbscan_model$membership_prob,
      outlier_scores = hdbscan_model$outlier_scores
    )
    write.csv(out, file = paste0(file.out, "_clusters.csv"), row.names = FALSE)
  }
  return(result)
}



########################## DBSCAN / HDBSCAN PLOTS ##########################



# helper function to choose coordinates for plotting
# method      - one of "pca", "umap", "tsne", "ica"
# result.obj  - result object for methods that require it (umap/tsne/ica)
# dims        - dimensions to plot


.get_density_plot_data = function(dataset, method = "pca", result.obj = NULL, dims = c(1,2), scale = FALSE) {
  method = tolower(method)
  
  if (method == "pca") {
    pca.result = pca(dataset, scale = scale)
    emb = pca.result$scores
    if (any(dims > ncol(emb))) stop(paste("dims out of range: PCA has only", ncol(emb), "components."))
    plot.df = data.frame(emb[, dims])
    names(plot.df) = c("x", "y")
    attr(plot.df, "xlab") = paste0("PC", dims[1])
    attr(plot.df, "ylab") = paste0("PC", dims[2])
  } else if (method %in% c("umap", "tsne", "ica")) {
    if (is.null(result.obj)) stop(paste("result.obj is required when method =", method))
    emb = result.obj$embedding
    if (ncol(emb) < 2) stop("Embedding must have at least 2 components for a 2D plot.")
    if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
    plot.df = data.frame(emb[, dims])
    names(plot.df) = c("x", "y")
    prefix = switch(method, umap = "UMAP", tsne = "tSNE", ica = "IC")
    attr(plot.df, "xlab") = paste0(prefix, dims[1])
    attr(plot.df, "ylab") = paste0(prefix, dims[2])
  } else {
    stop("method must be one of: 'pca', 'umap', 'tsne', 'ica'")
  }
  
  plot.df$label = colnames(dataset$data)
  plot.df
}



# 2d plot for DBSCAN clusters
# method      - coordinates used for plotting: "pca", "umap", "tsne", "ica"
# result.obj  - result object when method is umap/tsne/ica
# dims        - dimensions to plot
# labels      - if TRUE, shows sample names next to points
# bw          - if TRUE, uses black and white with shapes instead of colours
# leg.pos     - legend position
# xlim        - optional x axis limits c(min, max)
# ylim        - optional y axis limits c(min, max)


dbscan_plot2D = function(dataset, dbscan.result, method = "pca", result.obj = NULL, dims = c(1,2), labels = FALSE, bw = FALSE, leg.pos = "right", xlim = NULL, ylim = NULL, scale = FALSE) {
  plot.df = .get_density_plot_data(dataset, method = method, result.obj = result.obj, dims = dims, scale = scale)
  plot.df$group = factor(dbscan.result$cluster)
  levels(plot.df$group)[levels(plot.df$group) == "0"] = "Noise"
  
  if (bw) shape.values = seq_len(length(levels(plot.df$group)))
  
  if (bw) {
    p = ggplot2::ggplot(plot.df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], shape = .data[["group"]]))
  } else {
    p = ggplot2::ggplot(plot.df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["group"]]))
  }
  
  p = p + ggplot2::geom_point(size = 3, alpha = 0.9)
  
  if (bw) p = p + ggplot2::scale_shape_manual(values = shape.values)
  else p = p + ggplot2::scale_colour_brewer(palette = "Set1")
  
  p = p +
    ggplot2::xlab(attr(plot.df, "xlab")) +
    ggplot2::ylab(attr(plot.df, "ylab")) +
    ggplot2::ggtitle(paste0("DBSCAN 2D Plot (", toupper(method), ")")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = leg.pos)
  
  if (!is.null(xlim)) p = p + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) p = p + ggplot2::ylim(ylim[1], ylim[2])
  
  if (labels) {
    p = p + ggplot2::geom_text(ggplot2::aes(label = .data[["label"]]), hjust = -0.1, vjust = 0, size = 3)
  }
  
  p
}



# 2d plot for HDBSCAN clusters
# method      - coordinates used for plotting: "pca", "umap", "tsne", "ica"
# result.obj  - result object when method is umap/tsne/ica
# dims        - dimensions to plot
# labels      - if TRUE, shows sample names next to points
# bw          - if TRUE, uses black and white with shapes instead of colours
# leg.pos     - legend position
# xlim        - optional x axis limits c(min, max)
# ylim        - optional y axis limits c(min, max)


hdbscan_plot2D = function(dataset, hdbscan.result, method = "pca", result.obj = NULL, dims = c(1,2), labels = FALSE, bw = FALSE, leg.pos = "right", xlim = NULL, ylim = NULL, scale = FALSE) {
  plot.df = .get_density_plot_data(dataset, method = method, result.obj = result.obj, dims = dims, scale = scale)
  plot.df$group = factor(hdbscan.result$cluster)
  levels(plot.df$group)[levels(plot.df$group) == "0"] = "Noise"
  
  if (!is.null(hdbscan.result$membership_prob)) {
    plot.df$membership_prob = hdbscan.result$membership_prob
  } else {
    plot.df$membership_prob = 1
  }
  
  if (bw) shape.values = seq_len(length(levels(plot.df$group)))
  
  if (bw) {
    p = ggplot2::ggplot(plot.df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], shape = .data[["group"]]))
  } else {
    p = ggplot2::ggplot(plot.df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["group"]]))
  }
  
  p = p + ggplot2::geom_point(ggplot2::aes(alpha = .data[["membership_prob"]]), size = 3)
  
  if (bw) p = p + ggplot2::scale_shape_manual(values = shape.values)
  else p = p + ggplot2::scale_colour_brewer(palette = "Set1")
  
  p = p +
    ggplot2::scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
    ggplot2::xlab(attr(plot.df, "xlab")) +
    ggplot2::ylab(attr(plot.df, "ylab")) +
    ggplot2::ggtitle(paste0("HDBSCAN 2D Plot (", toupper(method), ")")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = leg.pos)
  
  if (!is.null(xlim)) p = p + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) p = p + ggplot2::ylim(ylim[1], ylim[2])
  
  if (labels) {
    p = p + ggplot2::geom_text(ggplot2::aes(label = .data[["label"]]), hjust = -0.1, vjust = 0, size = 3)
  }
  
  p
}



# 3d plot for DBSCAN clusters (uses plotly)
# method      - coordinates used for plotting: "umap", "tsne", "ica"
# result.obj  - result object with embedding
# dims        - dimensions to plot
# title       - plot title


dbscan_plot3D = function(dataset, dbscan.result, method = "umap", result.obj, dims = c(1,2,3), title = NULL) {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required. Install it with: install.packages('plotly')")
  method = tolower(method)
  if (!(method %in% c("umap", "tsne", "ica"))) stop("method must be one of: 'umap', 'tsne', 'ica'")
  if (missing(result.obj) || is.null(result.obj)) stop("result.obj is required for 3D plots")
  
  emb = result.obj$embedding
  if (ncol(emb) < 3) stop("Embedding must have at least 3 components for a 3D plot.")
  if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
  
  df = as.data.frame(emb[, dims])
  prefix = switch(method, umap = "UMAP", tsne = "tSNE", ica = "IC")
  colnames(df) = paste0(prefix, 1:3)
  df$Group = factor(dbscan.result$cluster)
  levels(df$Group)[levels(df$Group) == "0"] = "Noise"
  
  if (is.null(title)) title = paste0("DBSCAN 3D Plot (", toupper(method), ")")
  
  plotly::plot_ly(
    df,
    x = as.formula(paste0("~", colnames(df)[1])),
    y = as.formula(paste0("~", colnames(df)[2])),
    z = as.formula(paste0("~", colnames(df)[3])),
    color = ~Group,
    type = "scatter3d",
    mode = "markers"
  ) %>% plotly::layout(title = title)
}



# 3d plot for HDBSCAN clusters (uses plotly)
# method      - coordinates used for plotting: "umap", "tsne", "ica"
# result.obj  - result object with embedding
# dims        - dimensions to plot
# title       - plot title


hdbscan_plot3D = function(dataset, hdbscan.result, method = "umap", result.obj, dims = c(1,2,3), title = NULL) {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required. Install it with: install.packages('plotly')")
  method = tolower(method)
  if (!(method %in% c("umap", "tsne", "ica"))) stop("method must be one of: 'umap', 'tsne', 'ica'")
  if (missing(result.obj) || is.null(result.obj)) stop("result.obj is required for 3D plots")
  
  emb = result.obj$embedding
  if (ncol(emb) < 3) stop("Embedding must have at least 3 components for a 3D plot.")
  if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
  
  df = as.data.frame(emb[, dims])
  prefix = switch(method, umap = "UMAP", tsne = "tSNE", ica = "IC")
  colnames(df) = paste0(prefix, 1:3)
  df$Group = factor(hdbscan.result$cluster)
  levels(df$Group)[levels(df$Group) == "0"] = "Noise"
  df$membership_prob = if (!is.null(hdbscan.result$membership_prob)) hdbscan.result$membership_prob else 1
  
  if (is.null(title)) title = paste0("HDBSCAN 3D Plot (", toupper(method), ")")
  
  plotly::plot_ly(
    df,
    x = as.formula(paste0("~", colnames(df)[1])),
    y = as.formula(paste0("~", colnames(df)[2])),
    z = as.formula(paste0("~", colnames(df)[3])),
    color = ~Group,
    size = ~membership_prob,
    type = "scatter3d",
    mode = "markers"
  ) %>% plotly::layout(title = title)
}



# pairs plot for DBSCAN clusters
# method      - coordinates used for plotting: "umap", "tsne", "ica"
# result.obj  - result object with embedding
# dims        - dimensions to include (default: first 5)


dbscan_pairs_plot = function(dataset, dbscan.result, method = "umap", result.obj, dims = 1:min(5, ncol(result.obj$embedding)), ...) {
  if (!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required. Install it with: install.packages('GGally')")
  if (missing(result.obj) || is.null(result.obj)) stop("result.obj is required")
  emb = result.obj$embedding
  if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
  pairs.df = data.frame(emb[, dims])
  pairs.df$group = factor(dbscan.result$cluster)
  levels(pairs.df$group)[levels(pairs.df$group) == "0"] = "Noise"
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color = group), ...)
}



# pairs plot for HDBSCAN clusters
# method      - coordinates used for plotting: "umap", "tsne", "ica"
# result.obj  - result object with embedding
# dims        - dimensions to include (default: first 5)


hdbscan_pairs_plot = function(dataset, hdbscan.result, method = "umap", result.obj, dims = 1:min(5, ncol(result.obj$embedding)), ...) {
  if (!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required. Install it with: install.packages('GGally')")
  if (missing(result.obj) || is.null(result.obj)) stop("result.obj is required")
  emb = result.obj$embedding
  if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
  pairs.df = data.frame(emb[, dims])
  pairs.df$group = factor(hdbscan.result$cluster)
  levels(pairs.df$group)[levels(pairs.df$group) == "0"] = "Noise"
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color = group), ...)
}



# optional prediction for new data using a fitted DBSCAN model
# newdata      - matrix/data.frame with variables in columns and samples in rows
# original.data - matrix/data.frame used to fit the original model


dbscan_predict_newdata = function(dbscan.result, newdata, original.data) {
  if (is.null(dbscan.result$model)) stop("dbscan.result$model is required. Run dbscan_analysis_dataset(..., ret_model = TRUE).")
  stats::predict(dbscan.result$model, newdata = newdata, data = original.data)
}



# optional prediction for new data using a fitted HDBSCAN model
# newdata      - matrix/data.frame with variables in columns and samples in rows
# original.data - matrix/data.frame used to fit the original model


hdbscan_predict_newdata = function(hdbscan.result, newdata, original.data) {
  if (is.null(hdbscan.result$model)) stop("hdbscan.result$model is required. Run hdbscan_analysis_dataset(..., ret_model = TRUE).")
  stats::predict(hdbscan.result$model, newdata = newdata, data = original.data)
}
