############################################################################
######################## TESTE SIMPLES GMM - CACHEXIA ######################
############################################################################

# instalar se necessário
install.packages(c("specmine.datasets", "mclust"))

library(specmine.datasets)
library(mclust)
devtools::load_all(".")

data(cachexia)

# amostras em linhas, variáveis em colunas
x = t(as.matrix(cachexia$data))
x = scale(x)

# correr GMM
gmm_res = Mclust(x)

# resultados principais
cat("Modelo escolhido:", gmm_res$modelName, "\n")
cat("Número de clusters escolhido:", gmm_res$G, "\n")

cat("\nClusters:\n")
print(table(gmm_res$classification))

cat("\nResumo da incerteza:\n")
print(summary(gmm_res$uncertainty))

# comparar com a classe real
if ("Muscle.loss" %in% colnames(cachexia$metadata)) {
  col_class = "Muscle.loss"
} else if ("muscle.loss" %in% colnames(cachexia$metadata)) {
  col_class = "muscle.loss"
} else {
  col_class = colnames(cachexia$metadata)[1]
}

real_class = as.factor(cachexia$metadata[, col_class])

cat("\nTabela GMM vs classe real:\n")
print(table(GMM = gmm_res$classification, Class = real_class))

# PCA só para visualização rápida
pca_res = prcomp(x)

par(mfrow = c(1, 2))

plot(
  pca_res$x[,1], pca_res$x[,2],
  col = gmm_res$classification,
  pch = 19,
  xlab = "PC1", ylab = "PC2",
  main = "GMM clusters"
)

plot(
  pca_res$x[,1], pca_res$x[,2],
  col = heat.colors(length(gmm_res$uncertainty))[rank(gmm_res$uncertainty)],
  pch = 19,
  xlab = "PC1", ylab = "PC2",
  main = "GMM uncertainty"
)

# plot nativo do mclust
plot(gmm_res, what = "BIC")
plot(gmm_res, what = "classification")
plot(gmm_res, what = "uncertainty")

shiny::runApp("app.R")

