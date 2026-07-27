############################################################################
################################PCA#########################################
############################################################################

# perform pca analysis - classical
#' 
#' Auto-exported function: pca_analysis_dataset
#'
#' @param dataset Dataset to analyse.
#' @param scale Logical indicating whether variables should be scaled.
#' @param center Logical indicating whether variables should be centered.
#' @param write.file Logical indicating whether scores and loadings should be written to files.
#' @param file.out Output file prefix used when \code{write.file = TRUE}.
#' @param ... Additional arguments passed to \code{prcomp()}.
#'
#' @return An object of class \code{prcomp} containing PCA results. The
#'   \code{x} component stores sample scores and the \code{rotation} component
#'   stores variable loadings.
#'
#' @examples
#' datamat <- matrix(
#'   rnorm(20),
#'   nrow = 4,
#'   dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_res <- pca_analysis_dataset(dataset)
#'
#' @keywords internal
#' @export
pca_analysis_dataset = function(dataset, scale = TRUE, center = TRUE, 
                                write.file = FALSE, file.out = NULL, ...) {
  
  pca.result = prcomp(t(dataset$data), center = center, scale. = scale, ...)
  
  if (isTRUE(write.file)) {
    if (is.null(file.out) || !nzchar(file.out)) {
      stop("Please provide 'file.out' when write.file = TRUE.")
    }
    utils::write.csv(pca.result$x, file = paste0(file.out, "_scores.csv"))
    utils::write.csv(pca.result$rotation, file = paste0(file.out, "_loadings.csv"))
  }
  
  pca.result
}

# returns information about importance of the PC's
# pcs - PCs to get; sd - get std dev; prop - get proportion of variance; cumul - get cumulative
# min.cum - allows to define minimum cumulative % of variance
pca_importance = function(pca.res, pcs = 1:length(pca.res$sdev), sd = TRUE, prop = TRUE, cumul = TRUE, min.cum = NULL)
{
  rows = c()
  if (sd) rows = c(1)
  if (prop) rows = c(rows, 2)
  if (cumul) rows = c(rows, 3)
  
  if (inherits(pca.res, "prcomp")){
    if (!is.null(min.cum)) {
      cum = summary(pca.res)$importance[3,]
      pcs = 1:(min(which(cum > min.cum)))
    }
    
    res = summary(pca.res)$importance
  } else if (inherits(pca.res, "princomp")){
    vars = pca.res$sdev^2
    vars = vars/sum(vars)
    cum = cumsum(vars)
    if (!is.null(min.cum)) {
      pcs = 1:(min(which(cum > min.cum)))
    }
    res = rbind("Standard deviation" = pca.res$sdev, "Proportion of Variance" = vars, "Cumulative Proportion" = cum)
  }
  res[rows, pcs]
}

# robust PCA analysis 

# center - how the data will be centered "mean" or "median" (or NULL if none)
# scale - how the data will be scaled "sd" or "mad" (or NULL if none)
# k - number of PCs to compute

# returns objects of class princomp
#' Robust PCA analysis
#'
#' Performs robust PCA using \code{pcaPP::PCAgrid()}.
#'
#' @param dataset Dataset to analyse.
#' @param center Centering method.
#' @param scale Scaling method.
#' @param k Number of principal components to compute.
#' @param write.file Logical indicating whether scores and loadings should be written to files.
#' @param file.out Output file prefix used when \code{write.file = TRUE}.
#' @param ... Additional arguments passed to \code{pcaPP::PCAgrid()}.
#'
#' @return An object returned by \code{pcaPP::PCAgrid()}.
#'
#' @examples
#' \donttest{
#' datamat <- matrix(
#'   rnorm(30),
#'   nrow = 6,
#'   dimnames = list(paste0("x", 1:6), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_robust(dataset, k = 2)
#' }
#'
#' @export
pca_robust = function(dataset, center = "median", scale = "mad", k = 10,
                      write.file = FALSE, file.out = NULL, ...)
{
  pca.res = pcaPP::PCAgrid(t(dataset$data), k = k, center = center, scale = scale, scores = TRUE, ...)
  
  if (isTRUE(write.file)) {
    if (is.null(file.out) || !nzchar(file.out)) {
      stop("Please provide 'file.out' when write.file = TRUE.")
    }
    utils::write.csv(pca.res$scores, file = paste0(file.out, "_scores.csv"))
    utils::write.csv(pca.res$loadings, file = paste0(file.out, "_loadings.csv"))
  }
  
  pca.res
}


########################## PCA PLOTS ##################################

#scree plot
#' PCA scree plot
#'
#' Draws a scree plot with individual and cumulative explained variance.
#'
#' @param pca.result PCA result object.
#' @param num.pcs Number of principal components to display.
#' @param cex.leg Legend text size.
#' @param leg.pos Legend position.
#' @param lab.text Legend labels.
#' @param fill.col Line colours.
#' @param ylab Y-axis label.
#' @param xlab X-axis label.
#' @param ... Additional graphical parameters.
#'
#' @return No return value, called for its side effect of plotting.
#'
#' @examples
#' datamat <- matrix(
#'   rnorm(20),
#'   nrow = 4,
#'   dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_res <- pca_analysis_dataset(dataset)
#' pca_screeplot(pca_res)
#'
#' @export
pca_screeplot = function(pca.result, num.pcs = NULL, cex.leg = 0.8, leg.pos = "right", 
                         lab.text = c("individual percent","cumulative percent"), 
                         fill.col = c("blue","red"), ylab = "Percentage", xlab = "Principal components",
                         ...){
  importance = pca_importance(pca.result)
  if (is.null(num.pcs)) num.pcs = dim(importance)[2]
  opar=par(mfrow=c(1,1))
  on.exit(par(opar))
  matplot(seq(1, num.pcs), data.frame(t(importance[2:3,1:num.pcs]*100)), type="l", lty=1, col=fill.col, 
          xaxt='n', ylab= ylab,xlab=xlab, ...)
  legend(leg.pos, lab.text, cex=cex.leg, fill=fill.col)
  axis(1, at = 1:dim(importance)[2],labels= colnames(importance))
}

#2d scores plot
#' PCA 2D scores plot
#'
#' Creates a 2D PCA scores plot.
#'
#' @param dataset Dataset used in the PCA.
#' @param pca.result PCA result object.
#' @param column.class Optional metadata column used for colouring groups.
#' @param pcas Principal components to plot.
#' @param labels Logical indicating whether sample labels should be shown.
#' @param ellipses Logical indicating whether group ellipses should be drawn.
#' @param bw Logical indicating whether a black-and-white style should be used.
#' @param pallette Brewer palette identifier.
#' @param leg.pos Legend position.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' datamat <- matrix(
#'   rnorm(20),
#'   nrow = 4,
#'   dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_res <- pca_analysis_dataset(dataset)
#' pca_scoresplot2D(dataset, pca_res, column.class = "class")
#' }
#'
#' @export
pca_scoresplot2D = function(dataset, pca.result, column.class = NULL, pcas = c(1,2), labels = FALSE, 
                            ellipses = FALSE, bw=FALSE, pallette = 2, leg.pos = "right", xlim = NULL, ylim = NULL)
{
  has.legend = FALSE
  if (inherits(pca.result, "prcomp")){
    scores = pca.result$x
  } else if (inherits(pca.result, "princomp")){
    scores = pca.result$scores
  }
  pca.points = data.frame(scores[,pcas])
  names(pca.points) = c("x","y")
  if (is.null(column.class)){
    group.values = factor(rep(4, ncol(dataset$data)))
  } else {
    group.values = dataset$metadata[,column.class]
    has.legend = TRUE
  }
  pca.points$group = group.values
  pca.points$label = colnames(dataset$data)
  if (bw) shape.values = 1:length(levels(group.values))
  if (bw)
    pca.plot = ggplot2::ggplot(data = pca.points, ggplot2::aes_string(x='x', y='y', shape='group'))
  else
    pca.plot = ggplot2::ggplot(data = pca.points, ggplot2::aes_string(x='x', y='y',colour='group')) 
  pca.plot = pca.plot + ggplot2::geom_point(size=3, alpha=1) 
  if (bw) 
    pca.plot = pca.plot + ggplot2::scale_shape_manual(values = shape.values)
  else
    pca.plot = pca.plot + ggplot2::scale_colour_brewer(type = "qual", palette=pallette) 
  pca.plot = pca.plot +
    ggplot2::xlab(paste(paste("PC",pcas[1]," -",sep=""), paste(pca_importance(pca.result, pcas[1], sd=FALSE, prop=TRUE, cumul = FALSE)*100,"%",sep=""))) + 
    ggplot2::ylab(paste(paste("PC",pcas[2]," -",sep=""), paste(pca_importance(pca.result, pcas[2], sd=FALSE, prop=TRUE, cumul = FALSE)*100,"%",sep="")))
  if (has.legend) {
    if (bw) pca.plot = pca.plot + ggplot2::theme_bw() 
    else pca.plot = pca.plot + ggplot2::theme(legend.position = leg.pos)
  }
  if (!is.null(xlim)){
    pca.plot = pca.plot + ggplot2::xlim(xlim[1],xlim[2])
  }
  if (!is.null(ylim)){
    pca.plot = pca.plot + ggplot2::ylim(ylim[1],ylim[2])
  }
  if (labels){
    pca.plot = pca.plot + ggplot2::geom_text(data = pca.points, ggplot2::aes_string(x = 'x',y = 'y',label='label'),hjust=-0.1, vjust=0)
  }
  if (!bw & ellipses){
    df.ellipses = calculate_ellipses(pca.points)
    pca.plot = pca.plot + ggplot2::geom_path(data=df.ellipses, ggplot2::aes_string(x='x', y='y',colour='group'), size=1, linetype=2) 
  }
  pca.plot
}

#3d scores plot
#' PCA 3D scores plot using rgl
#'
#' Creates an interactive 3D PCA scores plot.
#'
#' @param dataset Dataset used in the PCA.
#' @param pca.result PCA result object.
#' @param column.class Optional metadata column used for colouring groups.
#' @param pcas Principal components to plot.
#' @param size Point size.
#' @param labels Logical indicating whether sample labels should be shown.
#'
#' @return No return value, called for its side effect of plotting.
#'
#' @examples
#' \donttest{
#' datamat <- matrix(
#'   rnorm(20),
#'   nrow = 4,
#'   dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_res <- pca_analysis_dataset(dataset)
#' pca_scoresplot3D_rgl(dataset, pca_res, column.class = "class")
#' }
#'
#' @export
pca_scoresplot3D_rgl = function(dataset, pca.result, column.class = NULL, pcas = c(1,2,3), size = 1, 
                                labels = FALSE) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package rgl needed for this function to work. Please install it: install.packages('rgl')",
         call. = FALSE)
  }
  
  if (inherits(pca.result, "prcomp")){
    scores = pca.result$x
  } else if (inherits(pca.result, "princomp")){
    scores = pca.result$scores
  }
  rgl::plot3d(scores[,pcas], type = "s", col = as.integer(dataset$metadata[,column.class]),
              size=size)
  if (labels){
    rgl::text3d(scores[,pcas],texts=colnames(dataset$data), cex=0.6)
  }
}

#' PCA 3D scores plot
#'
#' Creates a static 3D PCA scores plot using \code{scatterplot3d}.
#'
#' @param dataset Dataset used in the PCA.
#' @param pca.result PCA result object.
#' @param column.class Optional metadata column used for colouring groups.
#' @param pcas Principal components to plot.
#'
#' @return No return value, called for its side effect of plotting.
#'
#' @examples
#' \donttest{
#' datamat <- matrix(
#'   rnorm(20),
#'   nrow = 4,
#'   dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_res <- pca_analysis_dataset(dataset)
#' pca_scoresplot3D(dataset, pca_res, column.class = "class")
#' }
#'
#' @export
pca_scoresplot3D = function(dataset, pca.result, column.class = NULL, pcas=c(1,2,3))
{
  if(!requireNamespace("scatterplot3d", quietly = TRUE)){
    stop("Package scatterplot3d needed for this function to work. Please install it: install.packages('scatterplot3d').",
         call. = FALSE)
  }
  
  has.legend = FALSE
  if (inherits(pca.result, "prcomp")){
    scores = pca.result$x
  } else if (inherits(pca.result, "princomp")){
    scores = pca.result$scores
  }
  if (is.null(column.class)){
    group.values = rep(4, ncol(dataset$data))
  } else {
    group.values = as.integer(dataset$metadata[,column.class])
    has.legend = TRUE
  }
  
  scatterplot3d::scatterplot3d(scores[,pcas], color=group.values, pch=17)
  if (has.legend){
    classes = dataset$metadata[,column.class]
    legend(-1.5, 2.5, levels(classes), col = 1:length(classes), cex = 0.7, pt.cex = 1, pch= 17)
  }
}

#biplots
#' PCA biplot
#'
#' Draws a PCA biplot for a \code{prcomp} or \code{princomp} result.
#'
#' @param dataset Dataset used in the PCA.
#' @param pca.result PCA result object.
#' @param cex Text expansion factor.
#' @param legend.cex Legend text size.
#' @param x.colors Colour vector or metadata column name.
#' @param inset Legend inset.
#' @param legend.place Legend position.
#' @param ... Additional arguments passed to the plotting method.
#'
#' @return No return value, called for its side effect of plotting.
#'
#' @examples
#' datamat <- matrix(
#'   rnorm(20),
#'   nrow = 4,
#'   dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_res <- pca_analysis_dataset(dataset)
#' pca_biplot(dataset, pca_res)
#'
#' @export
pca_biplot= function(dataset, pca.result, cex = 0.8, legend.cex = 0.8, x.colors = 1, inset = c(0, 0), legend.place = "topright", ...) {
  x.flag = FALSE
  if (x.colors %in% colnames(dataset$metadata)){
    x.colors.meta = x.colors
    x.colors = as.integer(dataset$metadata[, x.colors])
    x.flag = TRUE
    opar=par(xpd=TRUE, mar=par()$mar+c(0,0,0,6))
    on.exit(par(opar))
  }
  
  if (inherits(pca.result, "prcomp")){
    biplot_prcomp_modified(pca.result, cex = cex, x.colors = x.colors, ...)
  } else if (inherits(pca.result, "princomp")){
    biplot_princomp_modified(pca.result, cex = cex, x.colors = x.colors, ...)
  } else {
    stop("Class not supported")
  } 
  if (x.flag){
    legend(legend.place, inset = inset, levels(dataset$metadata[, x.colors.meta]), cex=legend.cex, bty="n", fill = sort(as.integer(factor(levels(dataset$metadata[, x.colors.meta]))))) 
    opar=par(mar=c(5, 4, 4, 2) + 0.1)
    on.exit(par(opar))
  }
}

biplot_princomp_modified = function (x, x.colors, choices = 1L:2L, scale = 1, pc.biplot = FALSE, ...) 
{
  if (length(choices) != 2L) 
    stop("length of choices must be 2")
  if (!length(scores <- x$scores)) 
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))), 
         domain = NA)
  lam <- x$sdev[choices]
  if (is.null(n <- x$n.obs)) 
    n <- 1
  lam <- lam * sqrt(n)
  if (scale < 0 || scale > 1) 
    warning("'scale' is outside [0, 1]")
  if (scale != 0) 
    lam <- lam^scale
  else lam <- 1
  if (pc.biplot) 
    lam <- lam/sqrt(n)
  biplot_default_modified(t(t(scores[, choices])/lam), t(t(x$loadings[, 
                                                                      choices]) * lam), x.colors = x.colors, ...)
  invisible()
}

biplot_prcomp_modified = function (x, x.colors, choices = 1L:2L, scale = 1, pc.biplot = FALSE, ...) 
{
  if (length(choices) != 2L) 
    stop("length of choices must be 2")
  if (!length(scores <- x$x)) 
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))), 
         domain = NA)
  if (is.complex(scores)) 
    stop("biplots are not defined for complex PCA")
  lam <- x$sdev[choices]
  n <- NROW(scores)
  lam <- lam * sqrt(n)
  if (scale < 0 || scale > 1) 
    warning("'scale' is outside [0, 1]")
  if (scale != 0) 
    lam <- lam^scale
  else lam <- 1
  if (pc.biplot) 
    lam <- lam/sqrt(n)
  biplot_default_modified(t(t(scores[, choices])/lam), t(t(x$rotation[, 
                                                                      choices]) * lam), x.colors = x.colors, ...)
  invisible()
}

biplot_default_modified = function (x, y, var.axes = TRUE, col, x.colors, colors, cex = rep(par("cex"), 2), 
                                    xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL, 
                                    arrow.len = 0.1, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, 
                                    ...) 
{
  n <- nrow(x)
  p <- nrow(y)
  if (missing(xlabs)) {
    xlabs <- dimnames(x)[[1L]]
    if (is.null(xlabs)) 
      xlabs <- 1L:n
  }
  xlabs <- as.character(xlabs)
  dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
  if (missing(ylabs)) {
    ylabs <- dimnames(y)[[1L]]
    if (is.null(ylabs)) 
      ylabs <- paste("Var", 1L:p)
  }
  ylabs <- as.character(ylabs)
  dimnames(y) <- list(ylabs, dimnames(y)[[2L]])
  if (length(cex) == 1L) 
    cex <- c(cex, cex)
  if (missing(col)) {
    col <- par("col")
    if (!is.numeric(col)) 
      col <- match(col, palette(), nomatch = 1L)
    col <- c(col, col + 1L)
  }
  else if (length(col) == 1L) 
    col <- c(col, col)
  unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)), 
                                  abs(max(x, na.rm = TRUE)))
  rangx1 <- unsigned.range(x[, 1L])
  rangx2 <- unsigned.range(x[, 2L])
  rangy1 <- unsigned.range(y[, 1L])
  rangy2 <- unsigned.range(y[, 2L])
  if (missing(xlim) && missing(ylim)) 
    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
  else if (missing(xlim)) 
    xlim <- rangx1
  else if (missing(ylim)) 
    ylim <- rangx2
  ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
  on.exit(par(op))
  op <- par(pty = "s")
  if (!is.null(main)) 
    op <- c(op, par(mar = par("mar") + c(0, 0, 1, 0)))
  plot(x, type = "n", xlim = xlim, ylim = ylim, col = col[1L], 
       xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
  text(x, xlabs, cex = cex[1L], col = x.colors, ...)
  opar=par(new = TRUE)
  on.exit(par(opar))
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim * 
         ratio, xlab = "", ylab = "", col = col[1L], ...)
  axis(3, col = col[2L], ...)
  axis(4, col = col[2L], ...)
  box(col = col[1L])
  text(y, labels = ylabs, cex = cex[2L], col = col[2L], ...)
  if (var.axes) 
    arrows(0, 0, y[, 1L] * 0.8, y[, 2L] * 0.8, col = col[2L], 
           length = arrow.len)
  invisible()
}

#' PCA 3D biplot
#'
#' Creates an interactive 3D PCA biplot using \code{rgl}.
#'
#' @param dataset Dataset used in the PCA.
#' @param pca.result PCA result object.
#' @param column.class Optional metadata column used for colouring groups.
#' @param pcas Principal components to plot.
#'
#' @return No return value, called for its side effect of plotting.
#'
#' @examples
#' \donttest{
#' datamat <- matrix(
#'   rnorm(20),
#'   nrow = 4,
#'   dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_res <- pca_analysis_dataset(dataset)
#' pca_biplot3D(dataset, pca_res, column.class = "class")
#' }
#'
#' @export
pca_biplot3D = function(dataset, pca.result, column.class = NULL, pcas = c(1,2,3)){
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package rgl needed for this function to work. Please install it: install.packages('rgl')",
         call. = FALSE)
  }
  
  if (inherits(pca.result, "prcomp")){
    scores = pca.result$x
    rotation = pca.result$rotation
  } else if (inherits(pca.result, "princomp")){
    scores = pca.result$scores
    rotation = pca.result$loadings
  }  
  pca_scoresplot3D_rgl(dataset, pca.result, column.class, pcas)
  rgl::text3d(scores[,pcas], texts=colnames(dataset$data), cex=0.6)
  rgl::text3d(rotation[,pcas], texts = rownames(rotation), col = "red", cex=0.6)
  coords = NULL
  for (i in 1:nrow(rotation)){
    coords = rbind(coords, rbind(c(0,0,0), rotation[i, pcas]))
  }
  rgl::lines3d(coords, col="red", lwd = 4)
}

#pca pairs plot
#' PCA pairs plot
#'
#' Creates a pairs plot of PCA scores.
#'
#' @param dataset Dataset used in the PCA.
#' @param pca.result PCA result object.
#' @param column.class Optional metadata column used for colouring groups.
#' @param pcas Principal components to include.
#' @param ... Additional arguments passed to \code{GGally::ggpairs()}.
#'
#' @return A \code{ggmatrix} object.
#'
#' @examples
#' \donttest{
#' datamat <- matrix(
#'   rnorm(30),
#'   nrow = 6,
#'   dimnames = list(paste0("x", 1:6), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_res <- pca_analysis_dataset(dataset)
#' pca_pairs_plot(dataset, pca_res, column.class = "class", pcas = 1:3)
#' }
#'
#' @export
pca_pairs_plot = function(dataset, pca.result, column.class = NULL, pcas = c(1,2,3,4,5), ...){
  if (inherits(pca.result, "prcomp")){
    scores = pca.result$x
  } else if (inherits(pca.result, "princomp")){
    scores = pca.result$scores
  }  
  
  if (is.null(column.class)){
    group.values = rep(4, ncol(dataset$data))
  } else {
    group.values = dataset$metadata[,column.class]
  }
  pairs.df = data.frame(scores[,pcas])
  pairs.df$group = group.values
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color=group), ...)
}

#kmeans clustering with 3 PCs
#' PCA 3D k-means plot
#'
#' Plots PCA scores in 3D coloured by k-means cluster assignment.
#'
#' @param dataset Dataset used in the PCA.
#' @param pca.result PCA result object.
#' @param num.clusters Number of clusters.
#' @param pcas Principal components to plot.
#' @param kmeans.result Optional precomputed k-means result.
#' @param labels Logical indicating whether sample labels should be shown.
#' @param size Point size.
#' @param ... Additional arguments passed to \code{rgl::plot3d()}.
#'
#' @return No return value, called for its side effect of plotting.
#'
#' @examples
#' \donttest{
#' datamat <- matrix(
#'   rnorm(20),
#'   nrow = 4,
#'   dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_res <- pca_analysis_dataset(dataset)
#' pca_kmeans_plot3D(dataset, pca_res, num.clusters = 2)
#' }
#'
#' @export
pca_kmeans_plot3D = function(dataset, pca.result, num.clusters = 3, pcas = c(1,2,3), 
                             kmeans.result = NULL, labels = FALSE, size = 1,...) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package rgl needed for this function to work. Please install it: install.packages('rgl')",
         call. = FALSE)
  }
  
  if (inherits(pca.result, "prcomp")){
    scores = pca.result$x
  } else if (inherits(pca.result, "princomp")){
    scores = pca.result$scores
  }
  
  if (is.null(kmeans.result)){
    kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  }
  rgl::plot3d(scores[,pcas], type = "s", col = kmeans.result$cluster, size=size,...)
  if (labels){
    rgl::text3d(scores[,pcas],adj = c(1.2,1.2), texts=colnames(dataset$data), cex=0.6)
  }
}

#kmeans clustering with 2 first PCs
#' PCA 2D k-means plot
#'
#' Plots PCA scores in 2D coloured by k-means cluster assignment.
#'
#' @param dataset Dataset used in the PCA.
#' @param pca.result PCA result object.
#' @param num.clusters Number of clusters.
#' @param pcas Principal components to plot.
#' @param kmeans.result Optional precomputed k-means result.
#' @param labels Logical indicating whether sample labels should be shown.
#' @param bw Logical indicating whether a black-and-white style should be used.
#' @param ellipses Logical indicating whether cluster ellipses should be drawn.
#' @param leg.pos Legend position.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' datamat <- matrix(
#'   rnorm(20),
#'   nrow = 4,
#'   dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_res <- pca_analysis_dataset(dataset)
#' pca_kmeans_plot2D(dataset, pca_res, num.clusters = 2)
#' }
#'
#' @export
pca_kmeans_plot2D = function(dataset, pca.result, num.clusters = 3,
                             pcas = c(1,2), kmeans.result = NULL, labels = FALSE, bw=FALSE,
                             ellipses = FALSE, leg.pos = "right", xlim = NULL, ylim = NULL){
  if (inherits(pca.result, "prcomp")){
    scores = pca.result$x
  } else if (inherits(pca.result, "princomp")){
    scores = pca.result$scores
  }
  
  if (is.null(kmeans.result)){
    kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  }
  pca.points = data.frame(scores[,pcas])
  names(pca.points) = c("x","y")
  pca.points$group = factor(kmeans.result$cluster)
  pca.points$label = colnames(dataset$data)
  if (bw) shape.values = 1:num.clusters
  
  if (bw) pca.plot = ggplot2::ggplot(data = pca.points, ggplot2::aes_string(x='x', y='y',shape='group'))
  else pca.plot = ggplot2::ggplot(data = pca.points, ggplot2::aes_string(x='x', y='y',colour='group')) 
  pca.plot = pca.plot + ggplot2::geom_point(size=3, alpha=.6)
  
  if (bw) pca.plot = pca.plot + ggplot2::scale_shape_manual(values = shape.values)
  else  pca.plot = pca.plot + ggplot2::scale_colour_brewer(palette="Set1") 
  
  pca.plot = pca.plot + ggplot2::xlab(paste("PC",pcas[1],sep="")) + 
    ggplot2::ylab(paste("PC",pcas[2],sep=""))
  
  if (bw) pca.plot = pca.plot + ggplot2::theme_bw()
  else pca.plot = pca.plot + ggplot2::theme(legend.position = leg.pos)
  
  if (!is.null(xlim)){
    pca.plot = pca.plot + ggplot2::xlim(xlim[1],xlim[2])
  }
  if (!is.null(ylim)){
    pca.plot = pca.plot + ggplot2::ylim(ylim[1],ylim[2])
  }
  if (labels){
    pca.plot = pca.plot + ggplot2::geom_text(data = pca.points, ggplot2::aes_string(x='x',y='y',label='label'),hjust=-0.1, vjust=0, size = 3)
  }
  if (!bw & ellipses){
    df.ellipses = calculate_ellipses(pca.points)
    pca.plot = pca.plot + ggplot2::geom_path(data=df.ellipses, ggplot2::aes_string(x='x', y='y',colour='group'), size=1, linetype=2) 
  }
  pca.plot
}

#pca pairs with kmeans clusters plot
#' PCA pairs plot with k-means clusters
#'
#' Creates a pairs plot of PCA scores coloured by k-means cluster assignment.
#'
#' @param dataset Dataset used in the PCA.
#' @param pca.result PCA result object.
#' @param num.clusters Number of clusters.
#' @param kmeans.result Optional precomputed k-means result.
#' @param pcas Principal components to include.
#'
#' @return A \code{ggmatrix} object.
#'
#' @examples
#' \donttest{
#' datamat <- matrix(
#'   rnorm(30),
#'   nrow = 6,
#'   dimnames = list(paste0("x", 1:6), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_res <- pca_analysis_dataset(dataset)
#' pca_pairs_kmeans_plot(dataset, pca_res, num.clusters = 2, pcas = 1:3)
#' }
#'
#' @export
pca_pairs_kmeans_plot = function(dataset, pca.result, num.clusters = 3, kmeans.result = NULL, pcas = c(1,2,3,4,5)){
  if (inherits(pca.result, "prcomp")){
    scores = pca.result$x
  } else if (inherits(pca.result, "princomp")){
    scores = pca.result$scores
  }
  
  if (is.null(kmeans.result)){
    kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  }
  pairs.df = data.frame(scores[,pcas])
  pairs.df$group = factor(kmeans.result$cluster)
  GGally::ggpairs(pairs.df, mapping = ggplot2::aes(color=group))
}

#draw ellipses
calculate_ellipses = function(data){
  df_ell <- data.frame()
  for(g in levels(data$group)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(data[data$group==g,], ellipse::ellipse(cor(x, y), 
                                                                                            scale=c(sd(x),sd(y)), centre=c(mean(x),mean(y))))),group=g))
  }
  df_ell
}