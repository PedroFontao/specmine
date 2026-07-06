############################################################################
################ TESTE SIMPLES DAS QUALITY METRICS #########################
############################################################################

# instalar se necessário
install.packages(c("specmine.datasets", "clusterCrit", "uwot", "Rtsne"))

library(specmine.datasets)
library(clusterCrit)
library(uwot)
library(Rtsne)
devtools::load_all(".")

data(cachexia)

# dados: amostras em linhas, variáveis em colunas
x = t(as.matrix(cachexia$data))
x = scale(x)

############################################################################
########################### TESTE CLUSTERING ###############################
############################################################################

# clustering simples com k-means
km = kmeans(x, centers = 2, nstart = 25)

# usar a tua função
clust_metrics = cluster_quality_metrics(
  dataset = cachexia,
  clusters = km$cluster,
  scale = TRUE
)

cat("CLUSTER QUALITY METRICS\n")
print(clust_metrics)

cat("\nSilhouette:", clust_metrics$silhouette, "\n")
cat("Calinski-Harabasz:", clust_metrics$calinski_harabasz, "\n")
cat("Davies-Bouldin:", clust_metrics$davies_bouldin, "\n")

############################################################################
########################### TESTE EMBEDDINGS ###############################
############################################################################

# PCA
pca_res = prcomp(x)
pca_emb = pca_res$x[, 1:2]

pca_metrics = embedding_quality_metrics(
  dataset = cachexia,
  embedding = pca_emb,
  scale = TRUE,
  n_neighbors = 5
)

cat("\nPCA EMBEDDING METRICS\n")
print(pca_metrics)

# UMAP
set.seed(42)
umap_emb = uwot::umap(
  x,
  n_components = 2,
  n_neighbors = 15,
  min_dist = 0.1,
  metric = "euclidean"
)

umap_metrics = embedding_quality_metrics(
  dataset = cachexia,
  embedding = umap_emb,
  scale = TRUE,
  n_neighbors = 5
)

cat("\nUMAP EMBEDDING METRICS\n")
print(umap_metrics)

# t-SNE
set.seed(42)
tsne_res = Rtsne::Rtsne(
  x,
  dims = 2,
  perplexity = 20,
  max_iter = 1000,
  theta = 0.5,
  check_duplicates = FALSE
)

tsne_emb = tsne_res$Y

tsne_metrics = embedding_quality_metrics(
  dataset = cachexia,
  embedding = tsne_emb,
  scale = TRUE,
  n_neighbors = 5
)

cat("\nt-SNE EMBEDDING METRICS\n")
print(tsne_metrics)

############################################################################
######################## COMPARAÇÃO RÁPIDA ################################
############################################################################

comparison = data.frame(
  method = c("PCA", "UMAP", "t-SNE"),
  trustworthiness = c(
    pca_metrics$trustworthiness,
    umap_metrics$trustworthiness,
    tsne_metrics$trustworthiness
  ),
  continuity = c(
    pca_metrics$continuity,
    umap_metrics$continuity,
    tsne_metrics$continuity
  )
)

cat("\nEMBEDDING COMPARISON\n")
print(comparison)
