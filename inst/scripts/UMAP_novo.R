library(specmine)
library(specmine.datasets)
library(tidyr)
devtools::load_all(".")
data(cachexia)
propolis = data("propolis")


# analise base 2D e 3D
umap_res   = umap_analysis_dataset(cachexia)
umap_res3D = umap_analysis_dataset(cachexia, n_components = 3)

# scores 2D por grupo clinico
umap_scoresplot2D(cachexia, umap_res, column.class = "Muscle.loss")

# scores 2D com labels e elipses
umap_scoresplot2D(cachexia, umap_res, column.class = "Muscle.loss", labels = TRUE, ellipses = TRUE)

# scores 2D em preto e branco
umap_scoresplot2D(cachexia, umap_res, column.class = "Muscle.loss", bw = TRUE)

# scores 3D por grupo clinico
umap_scoresplot3D(cachexia, umap_res3D, column.class = "Muscle.loss")

# pairs plot com 3 componentes
umap_pairs_plot(cachexia, umap_res3D, column.class = "Muscle.loss")

# kmeans 2D com 2 clusters
umap_kmeans_plot2D(cachexia, umap_res, num.clusters = 2)

# kmeans 2D com 4 clusters e elipses
umap_kmeans_plot2D(cachexia, umap_res, num.clusters = 4, ellipses = TRUE)

# kmeans 3D com 3 clusters
umap_kmeans_plot3D(cachexia, umap_res3D, num.clusters = 4)

# pairs plot com kmeans
umap_pairs_kmeans_plot(cachexia, umap_res3D, num.clusters = 3)

# com scaling
umap_res_scaled = umap_analysis_dataset(cachexia, scale = TRUE)
umap_scoresplot2D(cachexia, umap_res_scaled, column.class = "Muscle.loss")

# guardar embedding em CSV
umap_analysis_dataset(cachexia, write.file = TRUE, file.out = "cachexia_umap")

