############################################################################
######################## QUALITY METRICS ###################################
############################################################################

# clustering metrics:
# - silhouette (higher is better)
# - calinski-harabasz (higher is better)
# - davies-bouldin (lower is better)

cluster_quality_metrics = function(dataset, clusters, scale = FALSE, remove.noise = TRUE,
                                   write.file = FALSE, file.out = "cluster_quality") {
  if (!requireNamespace("clusterCrit", quietly = TRUE)) {
    stop("Package 'clusterCrit' is required. Install it with: install.packages('clusterCrit')")
  }
  
  mat_check = as.matrix(dataset$data)
  if (any(is.na(mat_check)) || any(is.nan(mat_check)) || any(is.infinite(mat_check))) {
    stop("dataset$data contains NA, NaN or Inf values. Please clean your data first.")
  }
  
  x = t(mat_check)
  if (scale) x = base::scale(x)
  
  clusters = as.vector(clusters)
  
  if (length(clusters) != nrow(x)) {
    stop("Length of clusters must match the number of samples.")
  }
  
  if (remove.noise) {
    keep = clusters != 0
    x = x[keep, , drop = FALSE]
    clusters = clusters[keep]
  }
  
  clusters = as.integer(as.factor(clusters))
  n_clusters = length(unique(clusters))
  cluster_sizes = table(clusters)
  
  if (n_clusters < 2) {
    stop("At least 2 clusters are required to compute clustering quality metrics.")
  }
  
  crit = clusterCrit::intCriteria(
    traj = as.matrix(x),
    part = clusters,
    crit = c("Silhouette", "Calinski_Harabasz", "Davies_Bouldin")
  )
  
  result = list(
    n_samples = nrow(x),
    n_clusters = n_clusters,
    cluster_sizes = cluster_sizes,
    silhouette = unname(crit$silhouette),
    calinski_harabasz = unname(crit$calinski_harabasz),
    davies_bouldin = unname(crit$davies_bouldin),
    params = list(scale = scale, remove.noise = remove.noise)
  )
  
  if (write.file) {
    out = data.frame(
      n_samples = nrow(x),
      n_clusters = n_clusters,
      silhouette = unname(crit$silhouette),
      calinski_harabasz = unname(crit$calinski_harabasz),
      davies_bouldin = unname(crit$davies_bouldin)
    )
    write.csv(out, file = paste0(file.out, "_metrics.csv"), row.names = FALSE)
  }
  
  return(result)
}


############################################################################
######################## EMBEDDING METRICS #################################
############################################################################

# helper: get k nearest neighbors (excluding self)
.get_knn_indices = function(x, k) {
  d = as.matrix(stats::dist(x))
  diag(d) = Inf
  t(apply(d, 1, order))[ , 1:k, drop = FALSE]
}

# helper: compute full rank matrix
.get_rank_matrix = function(x) {
  d = as.matrix(stats::dist(x))
  ranks = t(apply(d, 1, rank, ties.method = "average"))
  diag(ranks) = NA
  ranks
}

# trustworthiness
# higher is better, range usually [0,1]
trustworthiness_metric = function(original_data, embedded_data, n_neighbors = 5) {
  original_data = as.matrix(original_data)
  embedded_data = as.matrix(embedded_data)
  
  n = nrow(original_data)
  if (n != nrow(embedded_data)) stop("original_data and embedded_data must have the same number of rows.")
  if (n_neighbors >= n / 2) stop("n_neighbors must be smaller than n_samples / 2.")
  
  orig_ranks = .get_rank_matrix(original_data)
  emb_knn = .get_knn_indices(embedded_data, n_neighbors)
  
  penalty = 0
  for (i in seq_len(n)) {
    for (j in emb_knn[i, ]) {
      r_ij = orig_ranks[i, j]
      if (!is.na(r_ij) && r_ij > n_neighbors) {
        penalty = penalty + (r_ij - n_neighbors)
      }
    }
  }
  
  norm = 2 / (n * n_neighbors * (2 * n - 3 * n_neighbors - 1))
  1 - norm * penalty
}

# continuity
# higher is better, range usually [0,1]
continuity_metric = function(original_data, embedded_data, n_neighbors = 5) {
  original_data = as.matrix(original_data)
  embedded_data = as.matrix(embedded_data)
  
  n = nrow(original_data)
  if (n != nrow(embedded_data)) stop("original_data and embedded_data must have the same number of rows.")
  if (n_neighbors >= n / 2) stop("n_neighbors must be smaller than n_samples / 2.")
  
  emb_ranks = .get_rank_matrix(embedded_data)
  orig_knn = .get_knn_indices(original_data, n_neighbors)
  
  penalty = 0
  for (i in seq_len(n)) {
    for (j in orig_knn[i, ]) {
      r_ij = emb_ranks[i, j]
      if (!is.na(r_ij) && r_ij > n_neighbors) {
        penalty = penalty + (r_ij - n_neighbors)
      }
    }
  }
  
  norm = 2 / (n * n_neighbors * (2 * n - 3 * n_neighbors - 1))
  1 - norm * penalty
}

# global embedding quality function
embedding_quality_metrics = function(dataset, embedding, scale = FALSE, n_neighbors = 5,
                                     write.file = FALSE, file.out = "embedding_quality") {
  mat_check = as.matrix(dataset$data)
  if (any(is.na(mat_check)) || any(is.nan(mat_check)) || any(is.infinite(mat_check))) {
    stop("dataset$data contains NA, NaN or Inf values. Please clean your data first.")
  }
  
  original_data = t(mat_check)
  if (scale) original_data = base::scale(original_data)
  
  embedding = as.matrix(embedding)
  
  if (nrow(original_data) != nrow(embedding)) {
    stop("The number of rows in embedding must match the number of samples.")
  }
  
  trust = trustworthiness_metric(original_data, embedding, n_neighbors = n_neighbors)
  cont  = continuity_metric(original_data, embedding, n_neighbors = n_neighbors)
  
  result = list(
    n_samples = nrow(original_data),
    n_components = ncol(embedding),
    trustworthiness = trust,
    continuity = cont,
    params = list(scale = scale, n_neighbors = n_neighbors)
  )
  
  if (write.file) {
    out = data.frame(
      n_samples = nrow(original_data),
      n_components = ncol(embedding),
      trustworthiness = trust,
      continuity = cont
    )
    write.csv(out, file = paste0(file.out, "_metrics.csv"), row.names = FALSE)
  }
  
  return(result)
}