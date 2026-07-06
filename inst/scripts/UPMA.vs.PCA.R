library(specmine)
library(specmine.datasets)
library(rgl)
data("cachexia")

#usar o kmeans para fazer clustering das amostras
km_res <- kmeans_clustering(cachexia,
                            num.clusters = 2,
                            type = "samples")  # ou "variables"


# FAZER O PCA
pca_res <- pca_analysis_dataset(cachexia, scale = TRUE)

pca_kmeans_plot2D(cachexia,
                  pca.result   = pca_res,
                  num.clusters = 4,
                  pcas         = c(1, 2),
                  labels       = FALSE)


pca_kmeans_plot3D(cachexia,
                  pca.result   = pca_res,
                  num.clusters = 2,
                  pcas         = c(1,2,3),
                  kmeans.result = NULL,
                  labels       = FALSE,
                  size         = 1)
rglwidget()


# Correr o UMAP
umap_res <- umap_analysis_dataset(cachexia, n_neighbors = 15, min_dist = 0.1)

umap_scoresplot2D(cachexia, umap_res, color.by = "Muscle.loss")


# Adicionar os clusters ao metadata temporariamente
cachexia_temp <- cachexia
cachexia_temp$metadata$KMeans <- as.character(km_res$cluster)

# Plot UMAP com os clusters do k-means
umap_scoresplot2D(cachexia_temp, umap_res, color.by = "KMeans")


# UMAP 3D
umap_res3D <- umap_analysis_dataset(cachexia, n_components = 3)
#visualização 3d com o muscle loss, ou seja doente vs controlo
umap_scoresplot3D(cachexia, umap_res3D, column.class = "Muscle.loss")
#3d com o kmeans de 2
umap_scoresplot3D(cachexia_temp, umap_res3D, color.class = "KMeans")


# com normalização dos dados scale = TRUE
umap_scaled <- umap_analysis_dataset(cachexia, scale = TRUE)
umap_scoresplot2D(cachexia, umap_scaled, color.class = "Muscle.loss")








