############################################################################
######################## COMPARE EMBEDDINGS ################################
############################################################################

compare_embeddings = function(dataset, embeddings.list, scale = FALSE, n_neighbors = 5,
                              write.file = FALSE, file.out = "embedding_comparison") {
  if (!is.list(embeddings.list) || length(embeddings.list) == 0) {
    stop("embeddings.list must be a non-empty named list.")
  }
  
  if (is.null(names(embeddings.list)) || any(names(embeddings.list) == "")) {
    stop("embeddings.list must be a named list.")
  }
  
  results = lapply(names(embeddings.list), function(method_name) {
    emb = embeddings.list[[method_name]]
    
    metrics = embedding_quality_metrics(
      dataset = dataset,
      embedding = emb,
      scale = scale,
      n_neighbors = n_neighbors
    )
    
    data.frame(
      method = method_name,
      n_samples = metrics$n_samples,
      n_components = metrics$n_components,
      trustworthiness = metrics$trustworthiness,
      continuity = metrics$continuity,
      stringsAsFactors = FALSE
    )
  })
  
  results_df = do.call(rbind, results)
  rownames(results_df) = NULL
  
  if (write.file) {
    write.csv(results_df, file = paste0(file.out, "_metrics.csv"), row.names = FALSE)
  }
  
  return(results_df)
}


############################################################################
######################## COMPARE CLUSTERINGS ###############################
############################################################################

compare_clusterings = function(dataset, clusterings.list, scale = FALSE, remove.noise = TRUE,
                               write.file = FALSE, file.out = "clustering_comparison") {
  if (!requireNamespace("clusterCrit", quietly = TRUE)) {
    stop("Package 'clusterCrit' is required. Install it with: install.packages('clusterCrit')")
  }
  
  if (!is.list(clusterings.list) || length(clusterings.list) == 0) {
    stop("clusterings.list must be a non-empty named list.")
  }
  
  if (is.null(names(clusterings.list)) || any(names(clusterings.list) == "")) {
    stop("clusterings.list must be a named list.")
  }
  
  results = lapply(names(clusterings.list), function(method_name) {
    cl = clusterings.list[[method_name]]
    
    metrics = cluster_quality_metrics(
      dataset = dataset,
      clusters = cl,
      scale = scale,
      remove.noise = remove.noise
    )
    
    data.frame(
      method = method_name,
      n_samples = metrics$n_samples,
      n_clusters = metrics$n_clusters,
      silhouette = metrics$silhouette,
      calinski_harabasz = metrics$calinski_harabasz,
      davies_bouldin = metrics$davies_bouldin,
      stringsAsFactors = FALSE
    )
  })
  
  results_df = do.call(rbind, results)
  rownames(results_df) = NULL
  
  if (write.file) {
    write.csv(results_df, file = paste0(file.out, "_metrics.csv"), row.names = FALSE)
  }
  
  return(results_df)
}
