#' Multi-class summary metrics
#'
#' Compute overall and class-averaged performance metrics for multiclass
#' classification models, including ROC AUC and log-loss when class
#' probabilities are available.
#'
#' @param data A data frame containing at least the columns `pred` and `obs`,
#'   plus one probability column per class when ROC or log-loss are needed.
#' @param lev An optional character vector with the class levels.
#' @param model An optional fitted model object passed by `caret`.
#'
#' @returns A named numeric vector with overall and class-averaged
#'   classification statistics.
#'
#' @examples
#' data <- data.frame(
#'   pred = factor(c("A", "B", "A", "B"), levels = c("A", "B")),
#'   obs = factor(c("A", "B", "B", "B"), levels = c("A", "B")),
#'   A = c(0.8, 0.2, 0.7, 0.3),
#'   B = c(0.2, 0.8, 0.3, 0.7)
#' )
#' multiClassSummary(data)
#' @export
multiClassSummary <- function(data, lev = NULL, model = NULL) {
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) {
    stop("levels of observed and predicted data do not match")
  }
  
  prob_stats <- lapply(levels(data[, "pred"]), function(class) {
    pred <- ifelse(data[, "pred"] == class, 1, 0)
    obs <- ifelse(data[, "obs"] == class, 1, 0)
    prob <- data[, class]
    
    cap_prob <- pmin(pmax(prob, 0.000001), 0.999999)
    prob_stats <- c(Metrics::auc(obs, prob), Metrics::logLoss(obs, cap_prob))
    names(prob_stats) <- c("ROC", "logLoss")
    prob_stats
  })
  
  prob_stats <- do.call(rbind, prob_stats)
  rownames(prob_stats) <- paste("Class:", levels(data[, "pred"]))
  
  CM <- caret::confusionMatrix(data[, "pred"], data[, "obs"])
  
  class_stats <- cbind(CM$byClass, prob_stats)
  class_stats <- colMeans(class_stats)
  
  overall_stats <- c(CM$overall)
  
  stats <- c(overall_stats, class_stats)
  stats <- stats[!names(stats) %in% c("AccuracyNull", "Prevalence", "Detection Prevalence")]
  
  names(stats) <- gsub("[[:blank:]]+", "_", names(stats))
  stats
}

# MODEL SELECTION AND PREDICTION

#' Train a classifier and predict new samples
#'
#' Train a classification model from a dataset and use the fitted model to
#' predict the class labels of new samples.
#'
#' @param dataset A dataset object containing data and metadata.
#' @param new.samples A data frame or matrix with new samples to classify.
#' @param column.class The metadata column containing the class labels.
#' @param model A model name accepted by \code{caret::train()}.
#' @param validation Validation method, such as \code{"boot"},
#'   \code{"boot632"}, \code{"cv"}, \code{"repeatedcv"}, \code{"LOOCV"},
#'   \code{"LGOCV"}, or \code{"oob"} where supported.
#' @param num.folds Number of folds used in resampling.
#' @param num.repeats Number of repeats used in repeated resampling.
#' @param tunelength Number of tuning levels evaluated by \code{caret}.
#' @param tunegrid Optional data frame of tuning parameter combinations.
#' @param metric Optional performance metric used for model selection.
#' @param summary.function Summary function passed to
#'   \code{caret::trainControl()}.
#'
#' @returns A list with two elements: \code{train.result}, the fitted training
#'   object, and \code{predictions.result}, a data frame with predicted classes
#'   for \code{new.samples}.
#'
#' @examples
#' \dontrun{
#' datamat <- matrix(
#'   rnorm(24),
#'   nrow = 4,
#'   dimnames = list(paste0("v", 1:4), paste0("s", 1:6))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "A", "B", "B", "B")))
#' dataset <- list(data = datamat, metadata = metadata)
#' new.samples <- datamat[, 1:2, drop = FALSE]
#' train_and_predict(dataset, new.samples, "class", model = "rpart", validation = "cv")
#' }
#' @export
train_and_predict <- function(dataset, new.samples, column.class, model, validation,
                              num.folds = 10, num.repeats = 10, tunelength = 10,
                              tunegrid = NULL, metric = NULL,
                              summary.function = caret::defaultSummary) {
  train.result <- train_classifier(
    dataset, column.class, model, validation, num.folds, num.repeats,
    tunelength, tunegrid, metric, summary.function
  )
  predict.result <- predict_samples(train.result, new.samples)
  result <- list(train.result = train.result, predictions.result = predict.result)
  result
}

#' Train a classifier
#'
#' Train a classifier from a dataset object using metadata or data-derived class
#' labels and a resampling strategy supported by \code{caret}.
#'
#' @param dataset A dataset object.
#' @param column.class The metadata column containing the class labels.
#' @param model A model name accepted by \code{caret::train()}.
#' @param validation Validation method used in training.
#' @param num.folds Number of folds used in resampling.
#' @param num.repeats Number of repeats used in repeated resampling.
#' @param tunelength Number of tuning levels evaluated by \code{caret}.
#' @param tunegrid Optional data frame of tuning parameter combinations.
#' @param metric Optional performance metric used for model selection.
#' @param summary.function Summary function passed to
#'   \code{caret::trainControl()}.
#' @param class.in.metadata Logical; if \code{TRUE}, class labels are taken from
#'   \code{dataset$metadata}, otherwise from \code{dataset$data}.
#'
#' @returns A \code{caret} training object returned by \code{caret::train()}.
#'
#' @examples
#' \dontrun{
#' datamat <- matrix(
#'   rnorm(24),
#'   nrow = 4,
#'   dimnames = list(paste0("v", 1:4), paste0("s", 1:6))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "A", "B", "B", "B")))
#' dataset <- list(data = datamat, metadata = metadata)
#' train_classifier(dataset, "class", model = "rpart", validation = "cv")
#' }
#' @export
train_classifier <- function(dataset, column.class, model, validation,
                             num.folds = 10, num.repeats = 10,
                             tunelength = 10, tunegrid = NULL, metric = NULL,
                             summary.function = caret::defaultSummary,
                             class.in.metadata = TRUE) {
  if (class.in.metadata) {
    train.result <- trainClassifier(
      dataset$data, dataset$metadata[, column.class], model, validation,
      num.folds, num.repeats, tunelength, tunegrid, metric, summary.function
    )
  } else {
    train.result <- trainClassifier(
      dataset$data, dataset$data[column.class, ], model, validation,
      num.folds, num.repeats, tunelength, tunegrid, metric, summary.function
    )
  }
  train.result
}

#' Train a classifier from a data matrix
#'
#' Internal helper used to build a \code{caret} classifier from a numeric data
#' matrix and a vector of sample classes.
#'
#' @param datamat A numeric data matrix with variables in rows and samples in columns.
#' @param sampleclass A vector of class labels for the samples.
#' @param model A model name accepted by \code{caret::train()}.
#' @param validation Validation method used in training.
#' @param num.folds Number of folds used in resampling.
#' @param num.repeats Number of repeats used in repeated resampling.
#' @param tunelength Number of tuning levels evaluated by \code{caret}.
#' @param tunegrid Optional data frame of tuning parameter combinations.
#' @param metric Optional performance metric used for model selection.
#' @param summary.function Summary function passed to \code{caret::trainControl()}.
#' @param class.in.metadata Logical; retained for interface compatibility.
#'
#' @returns A \code{caret} model object containing the fitted classifier,
#'   tuning results, and resampling summaries.
#' @keywords internal
#' @noRd
trainClassifier <- function(datamat, sampleclass, model, validation,
                            num.folds = 10, num.repeats = 10,
                            tunelength = 10, tunegrid = NULL, metric = NULL,
                            summary.function = caret::defaultSummary,
                            class.in.metadata = TRUE) {
  samples.df.ml <- data.frame(t(datamat))
  rnames <- gsub("[-\\ ]", "_", rownames(datamat))
  colnames(samples.df.ml) <- paste("X", rnames, sep = "")
  rownames(samples.df.ml) <- colnames(datamat)
  train.metric <- metric
  
  if (class.in.metadata) {
    samples.df.ml$class <- sampleclass
    if (is.null(metric) && is.factor(samples.df.ml$class)) {
      train.metric <- "Accuracy"
    } else if (is.null(metric) && !is.factor(samples.df.ml$class)) {
      train.metric <- "RMSE"
    }
  } else {
    if (is.null(metric) && is.factor(samples.df.ml$sampleclass)) {
      train.metric <- "Accuracy"
    } else if (is.null(metric) && !is.factor(samples.df.ml$sampleclass)) {
      train.metric <- "RMSE"
    }
  }
  
  class.probs <- isTRUE(train.metric == "ROC")
  
  train.control.args <- list(
    method = validation,
    number = num.folds,
    classProbs = class.probs,
    summaryFunction = summary.function
  )
  
  if (validation %in% c("repeatedcv", "adaptive_cv")) {
    train.control.args$repeats <- num.repeats
  }
  
  train.control <- do.call(caret::trainControl, train.control.args)
  
  if (class.in.metadata) {
    result.train <- caret::train(
      class ~ ., data = samples.df.ml, method = model,
      tuneLength = tunelength, metric = train.metric,
      trControl = train.control, tuneGrid = tunegrid
    )
  } else {
    result.train <- caret::train(
      sampleclass ~ ., data = samples.df.ml, method = model,
      tuneLength = tunelength, metric = train.metric,
      trControl = train.control, tuneGrid = tunegrid
    )
  }
  
  result.train
}

#' Predict sample classes
#'
#' Internal helper that predicts class labels for new samples using a trained
#' classifier.
#'
#' @param train.result A fitted training object.
#' @param new.samples A data frame or matrix containing new samples.
#'
#' @returns A \code{data.frame} with the sample names and predicted class labels.
#' @keywords internal
#' @noRd
predict_samples <- function(train.result, new.samples) {
  new.samples.df <- data.frame(t(new.samples))
  rnames <- gsub("[-\\ ]", "_", rownames(new.samples))
  colnames(new.samples.df) <- paste("X", rnames, sep = "")
  predict.result <- predict(train.result, newdata = new.samples.df)
  result <- data.frame(sample = rownames(new.samples.df), predicted.class = predict.result)
  result$sample <- as.character(result$sample)
  result
}

#' Train multiple models and compare their performance
#'
#' Train a set of models, collect their resampling performance, optionally
#' compute variable importance, and store fitted models and tuning summaries.
#'
#' @param dataset A dataset object.
#' @param models A character vector with model names accepted by \code{caret::train()}.
#' @param column.class The metadata column containing the class labels.
#' @param validation Validation method used in training.
#' @param num.folds Number of folds used in resampling.
#' @param num.repeats Number of repeats used in repeated resampling.
#' @param tunelength Number of tuning levels evaluated by \code{caret}.
#' @param tunegrid Optional list of tuning grids, one per model.
#' @param metric Optional performance metric used for model selection.
#' @param summary.function Summary function, or \code{"default"} to select the
#'   package default.
#' @param class.in.metadata Logical; if \code{TRUE}, class labels are taken from metadata.
#' @param compute.varimp Logical; if \code{TRUE}, variable importance is computed.
#'
#' @returns A list containing model performance, variable importance, full tuning
#'   results, best tuning settings, optional confusion matrices, and final fitted models.
#'
#' @examples
#' \dontrun{
#' datamat <- matrix(
#'   rnorm(24),
#'   nrow = 4,
#'   dimnames = list(paste0("v", 1:4), paste0("s", 1:6))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "A", "B", "B", "B")))
#' dataset <- list(data = datamat, metadata = metadata)
#' train_models_performance(
#'   dataset,
#'   models = c("rpart"),
#'   column.class = "class",
#'   validation = "cv",
#'   compute.varimp = FALSE
#' )
#' }
#' @export
train_models_performance <- function(dataset, models, column.class, validation,
                                     num.folds = 10, num.repeats = 10,
                                     tunelength = 10, tunegrid = NULL,
                                     metric = NULL, summary.function = "default",
                                     class.in.metadata = TRUE,
                                     compute.varimp = TRUE) {
  result.df <- NULL
  classification.flag <- FALSE
  if (compute.varimp) vars.imp <- list()
  final.result <- list()
  full.results <- list()
  
  if (is.factor(dataset$metadata[, column.class])) {
    classification.flag <- TRUE
    confusion.matrices <- list()
  }
  
  if (is.character(summary.function)) {
    if (!is.null(metric) && metric == "ROC" && summary.function == "default") {
      summary.function <- multiClassSummary
    } else if (summary.function == "default") {
      summary.function <- caret::defaultSummary
    }
  }
  
  best.tunes <- list()
  final.models <- list()
  for (i in seq_along(models)) {
    if (!is.null(tunegrid)) {
      tune.grid <- tunegrid[[models[i]]]
    } else {
      tune.grid <- NULL
    }
    
    train.result <- train_classifier(
      dataset, column.class, models[i], validation, num.folds,
      num.repeats, tunelength, tune.grid, metric, summary.function,
      class.in.metadata = class.in.metadata
    )
    
    if (compute.varimp) {
      vips <- var_importance(train.result)
      rownames(vips) <- substring(rownames(vips), 2, nchar(rownames(vips)))
      vips$Mean <- apply(vips, 1, mean)
    }
    
    bestTune <- train.result$bestTune
    result.df <- rbind(
      result.df,
      train.result$result[rownames(bestTune), c("Accuracy", "Kappa", "AccuracySD", "KappaSD")]
    )
    
    if (compute.varimp) {
      vars.imp[[i]] <- vips[order(vips$Mean, decreasing = TRUE), ]
      vips <- NULL
    }
    
    full.results[[i]] <- train.result$results
    if (classification.flag) confusion.matrices[[i]] <- try(caret::confusionMatrix(train.result), TRUE)
    best.tunes[[i]] <- train.result$bestTune
    final.models[[i]] <- train.result$finalModel
  }
  
  rownames(result.df) <- models
  if (compute.varimp) names(vars.imp) <- models
  names(full.results) <- models
  if (classification.flag) names(confusion.matrices) <- models
  names(best.tunes) <- models
  names(final.models) <- models
  
  result.df <- result.df[, colnames(result.df) %in% c(
    "RMSE", "Rsquared", "RMSESD", "RsquaredSD",
    "Accuracy", "AccuracySD", "Kappa", "KappaSD",
    "ROC", "Sensitivity", "Specificity",
    "SensitivitySD", "SpecificitySD", "ROCSD"
  )]
  
  final.result$performance <- result.df
  if (compute.varimp) final.result$vips <- vars.imp
  final.result$full.results <- full.results
  final.result$best.tunes <- best.tunes
  if (classification.flag) final.result$confusion.matrices <- confusion.matrices
  final.result$final.models <- final.models
  final.result
}

# VARIABLE IMPORTANCE

#' Extract variable importance
#'
#' Internal helper that extracts variable-importance values from a trained model.
#'
#' @param train.result A fitted training object.
#'
#' @returns A variable-importance object as returned by \code{caret::varImp()},
#'   typically stored as a data frame of importance scores.
#' @keywords internal
#' @noRd
var_importance <- function(train.result) {
  vip <- caret::varImp(train.result)
  vip$importance
}

#' Summarise variable importance tables
#'
#' Keep the top rows of each variable-importance table produced during model
#' comparison.
#'
#' @param performances A list returned by \code{train_models_performance()}.
#' @param number.rows Number of rows to keep from each variable-importance table.
#'
#' @returns A list of truncated variable-importance tables, one per model.
#'
#' @examples
#' performances <- list(
#'   vips = list(
#'     rf = data.frame(Mean = c(0.9, 0.7, 0.4), row.names = c("v1", "v2", "v3")),
#'     svm = data.frame(Mean = c(0.8, 0.6, 0.5), row.names = c("v1", "v2", "v3"))
#'   )
#' )
#' summary_var_importance(performances, 2)
#' @export
summary_var_importance <- function(performances, number.rows) {
  for (i in seq_along(performances$vips)) {
    performances$vips[[i]] <- performances$vips[[i]][1:number.rows, ]
  }
  performances$vips
}

# PCA PLOTS

#' Auto-exported function: pca_plot_3d
#'
#' Draw a 3D PCA scatter plot.
#'
#' @param dataset A dataset object containing metadata.
#' @param model A PCA result object containing component scores.
#' @param var.class The metadata variable used to define classes.
#' @param pcas A length-3 integer vector indicating which principal components to plot.
#' @param colors Optional vector of colors used for the classes.
#' @param legend.place Position of the legend.
#' @param ... Additional arguments passed to \code{legend()}.
#'
#' @return A 3D scatter plot of the selected principal components, drawn for its side effects.
#'
#' @examples
#' \dontrun{
#' datamat <- matrix(
#'   rnorm(24),
#'   nrow = 4,
#'   dimnames = list(paste0("v", 1:4), paste0("s", 1:6))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "A", "B", "B", "B")))
#' dataset <- list(data = datamat, metadata = metadata)
#' pca_model <- list(scores = prcomp(t(datamat))$x)
#' pca_plot_3d(dataset, pca_model, "class")
#' }
#' @keywords internal
#' @export
pca_plot_3d <- function(dataset, model, var.class, pcas = 1:3,
                        colors = NULL, legend.place = "topright", ...) {
  if (!requireNamespace("qdap", quietly = TRUE)) {
    if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
      stop("Packages qdap and scatterplot3d are needed for this function to work. Please install them: install.packages(c('qdap','scatterplot3d')).",
           call. = FALSE)
    } else {
      stop("Package qdap needed for this function to work. Please install it: install.packages('qdap').",
           call. = FALSE)
    }
  } else if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    stop("Package scatterplot3d needed for this function to work. Please install it: install.packages('scatterplot3d').",
         call. = FALSE)
  }
  
  if (length(pcas) != 3) stop("Wrong dimension in parameter pcas")
  if (ncol(model$scores) < 3) stop("Less than 3 components")
  
  classes <- dataset$metadata[, var.class]
  labs <- paste("Component", pcas)
  
  if (is.null(colors)) {
    colors <- seq_along(levels(classes))
  }
  
  colors_metadata <- qdap::mgsub(levels(classes), colors, classes)
  scatterplot3d::scatterplot3d(
    model$scores[, pcas],
    color = colors_metadata,
    pch = 17,
    xlab = labs[1],
    ylab = labs[2],
    zlab = labs[3]
  )
  legend(legend.place, levels(classes), col = colors, pch = 17, ...)
}