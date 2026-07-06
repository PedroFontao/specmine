library(specmine)
library(specmine.datasets)
devtools::load_all(".")
data(cachexia)

# analise base 2D
tsne_res = tsne_analysis_dataset(cachexia)

# analise 3D
tsne_res3D = tsne_analysis_dataset(cachexia, n_components = 3)

# scores 2D por grupo
tsne_scoresplot2D(cachexia, tsne_res, column.class = "Muscle.loss")

# com labels e elipses
tsne_scoresplot2D(cachexia, tsne_res, column.class = "Muscle.loss", labels = TRUE, ellipses = TRUE)

# preto e branco
tsne_scoresplot2D(cachexia, tsne_res, column.class = "Muscle.loss", bw = TRUE)

# scores 3D
tsne_scoresplot3D(cachexia, tsne_res3D, column.class = "Muscle.loss")

# kmeans 2D
tsne_kmeans_plot2D(cachexia, tsne_res, num.clusters = 2)

# kmeans 3D
tsne_kmeans_plot3D(cachexia, tsne_res3D, num.clusters = 3)

# pairs plot
tsne_pairs_plot(cachexia, tsne_res3D, column.class = "Muscle.loss")

# pairs com kmeans
tsne_pairs_kmeans_plot(cachexia, tsne_res3D, num.clusters = 3)
