############################################################################
################################ GMM #######################################
############################################################################



# perform Gaussian Mixture Model clustering
# num.clusters - number of mixture components; if NULL, selected automatically by BIC
# model.names  - models to be fitted by mclust (default: NULL = all valid models)
# scale        - if TRUE, scales data before clustering
# seed         - random seed for reproducibility
# ret_model    - if TRUE, returns the full mclust model object
# write.file   - if TRUE, saves cluster assignments to a CSV file
# file.out     - base name for output file


gmm_analysis_dataset = function(dataset, num.clusters = NULL, model.names = NULL, scale = FALSE, seed = 42, ret_model = FALSE, write.file = FALSE, file.out = "gmm", ...) {
  if (!requireNamespace("mclust", quietly = TRUE)) stop("Package 'mclust' is required. Install it with: install.packages('mclust')")
  
  mat_check = as.matrix(dataset$data)
  if (any(is.na(mat_check)) || any(is.nan(mat_check)) || any(is.infinite(mat_check))) {
    stop("dataset$data contains NA, NaN or Inf values. Please clean your data first.")
  }
  
  mat = t(mat_check)
  if (scale) mat = base::scale(mat)
  
  set.seed(seed)
  gmm_model = mclust::Mclust(data = mat, G = num.clusters, modelNames = model.names, ...)
  
  clusters = gmm_model$classification
  names(clusters) = colnames(dataset$data)
  
  result = list(
    cluster = clusters,
    z = gmm_model$z,
    uncertainty = gmm_model$uncertainty,
    bic = gmm_model$bic,
    loglik = gmm_model$loglik,
    modelName = gmm_model$modelName,
    G = gmm_model$G,
    parameters = gmm_model$parameters,
    params = list(num.clusters = num.clusters, model.names = model.names, scale = scale, seed = seed)
  )
  
  if (ret_model) result$model = gmm_model
  
  if (write.file) {
    out = data.frame(
      sample = names(clusters),
      cluster = clusters,
      uncertainty = gmm_model$uncertainty
    )
    write.csv(out, file = paste0(file.out, "_clusters.csv"), row.names = FALSE)
  }
  
  return(result)
}



############################## GMM PLOTS ###################################



# helper function to choose coordinates for plotting
# method      - one of "pca", "umap", "tsne", "ica"
# result.obj  - result object for methods that require it (umap/tsne/ica)
# dims        - dimensions to plot


.get_gmm_plot_data = function(dataset, method = "pca", result.obj = NULL, dims = c(1,2), scale = FALSE) {
  method = tolower(method)
  
  if (method == "pca") {
    mat = t(as.matrix(dataset$data))
    if (scale) mat = base::scale(mat)
    pca.result = stats::prcomp(mat, center = TRUE, scale. = FALSE)
    emb = pca.result$x
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



# 2d plot for GMM clusters
# method      - coordinates used for plotting: "pca", "umap", "tsne", "ica"
# result.obj  - result object when method is umap/tsne/ica
# dims        - dimensions to plot
# labels      - if TRUE, shows sample names next to points
# bw          - if TRUE, uses black and white with shapes instead of colours
# leg.pos     - legend position
# xlim        - optional x axis limits c(min, max)
# ylim        - optional y axis limits c(min, max)


gmm_plot2D = function(dataset, gmm.result, method = "pca", result.obj = NULL, dims = c(1,2), labels = FALSE, bw = FALSE, leg.pos = "right", xlim = NULL, ylim = NULL, scale = FALSE) {
  plot.df = .get_gmm_plot_data(dataset, method = method, result.obj = result.obj, dims = dims, scale = scale)
  plot.df$group = factor(gmm.result$cluster)
  plot.df$uncertainty = gmm.result$uncertainty
  
  if (bw) shape.values = seq_len(length(levels(plot.df$group)))
  
  if (bw) {
    p = ggplot2::ggplot(plot.df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], shape = .data[["group"]]))
  } else {
    p = ggplot2::ggplot(plot.df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["group"]]))
  }
  
  p = p + ggplot2::geom_point(ggplot2::aes(alpha = 1 - .data[["uncertainty"]]), size = 3)
  
  if (bw) p = p + ggplot2::scale_shape_manual(values = shape.values)
  else p = p + ggplot2::scale_colour_brewer(palette = "Set1")
  
  p = p +
    ggplot2::scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
    ggplot2::xlab(attr(plot.df, "xlab")) +
    ggplot2::ylab(attr(plot.df, "ylab")) +
    ggplot2::ggtitle(paste0("GMM 2D Plot (", toupper(method), ")")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = leg.pos)
  
  if (!is.null(xlim)) p = p + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) p = p + ggplot2::ylim(ylim[1], ylim[2])
  
  if (labels) {
    p = p + ggplot2::geom_text(ggplot2::aes(label = .data[["label"]]), hjust = -0.1, vjust = 0, size = 3)
  }
  
  p
}



# 3d plot for GMM clusters (uses plotly)
# method      - coordinates used for plotting: "umap", "tsne", "ica"
# result.obj  - result object with embedding
# dims        - dimensions to plot
# title       - plot title


gmm_plot3D = function(dataset, gmm.result, method = "umap", result.obj, dims = c(1,2,3), title = NULL) {
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
  df$Group = factor(gmm.result$cluster)
  df$certainty = 1 - gmm.result$uncertainty
  
  if (is.null(title)) title = paste0("GMM 3D Plot (", toupper(method), ")")
  
  plotly::plot_ly(
    df,
    x = as.formula(paste0("~", colnames(df)[1])),
    y = as.formula(paste0("~", colnames(df)[2])),
    z = as.formula(paste0("~", colnames(df)[3])),
    color = ~Group,
    size = ~certainty,
    type = "scatter3d",
    mode = "markers"
  ) %>% plotly::layout(title = title)
}



# pairs plot for GMM clusters
# result.obj  - result object with embedding
# dims        - dimensions to include (default: first 5)


gmm_pairs_plot = function(dataset, gmm.result, result.obj, dims = 1:min(5, ncol(result.obj$embedding)), ...) {
  if (!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required. Install it with: install.packages('GGally')")
  
  if (missing(result.obj) || is.null(result.obj)) stop("result.obj is required")
  emb = result.obj$embedding
  if (any(dims > ncol(emb))) stop(paste("dims out of range: embedding has only", ncol(emb), "components."))
  
  pairs.df = data.frame(emb[, dims])
  pairs.df$group = factor(gmm.result$cluster)
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color = group), ...)
}



# BIC plot for GMM model selection
# shows the selected BIC surface stored in mclust result


gmm_bic_plot = function(gmm.result) {
  if (is.null(gmm.result$model)) stop("gmm.result$model is required. Run gmm_analysis_dataset(..., ret_model = TRUE).")
  plot(gmm.result$model, what = "BIC")
}



# uncertainty plot for GMM
# method      - coordinates used for plotting: "pca", "umap", "tsne", "ica"
# result.obj  - result object when method is umap/tsne/ica
# dims        - dimensions to plot
# leg.pos     - legend position


gmm_uncertainty_plot2D = function(dataset, gmm.result, method = "pca", result.obj = NULL, dims = c(1,2), leg.pos = "right", xlim = NULL, ylim = NULL, scale = FALSE) {
  plot.df = .get_gmm_plot_data(dataset, method = method, result.obj = result.obj, dims = dims, scale = scale)
  plot.df$uncertainty = gmm.result$uncertainty
  
  p = ggplot2::ggplot(plot.df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["uncertainty"]])) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::scale_colour_gradient(low = "steelblue", high = "red") +
    ggplot2::xlab(attr(plot.df, "xlab")) +
    ggplot2::ylab(attr(plot.df, "ylab")) +
    ggplot2::ggtitle(paste0("GMM Uncertainty Plot (", toupper(method), ")")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = leg.pos)
  
  if (!is.null(xlim)) p = p + ggplot2::xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) p = p + ggplot2::ylim(ylim[1], ylim[2])
  
  p
}
