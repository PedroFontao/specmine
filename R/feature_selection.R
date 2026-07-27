###############################################################
#####################FEATURE SELECTION#########################
###############################################################

# feature selection
# dataset: data and metadata structure
# column.class: metadata column class
# method: "rfe" (recursive feature elimination) or "filter" (feature selection using univariate filters)
# functions:
#'
#' Auto-exported function: feature_selection
#'
#' @keywords internal
#' @return A feature selection object produced by \code{caret}. For \code{method = "rfe"}, this is an \code{rfe} object; for \code{method = "filter"}, this is an \code{sbf} object.
#' @examples
#' datamat <- matrix(
#'   rnorm(20),
#'   nrow = 4,
#'   dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#' )
#' metadata <- data.frame(class = factor(c("A", "A", "B", "B", "A")))
#' dataset <- list(data = datamat, metadata = metadata)
#' feature_selection(dataset, "class", method = "filter", functions = caret::rfSBF)
#' @export
feature_selection = function(dataset, column.class, method = "rfe", functions, validation = "cv",
                             repeats = 5, number = 10, subsets = 2^(2:4)){
  if (method == "rfe"){
    result = recursive_feature_elimination(dataset$data, dataset$metadata[,column.class], functions,
                                           validation, repeats, number, subsets)
  }
  else if (method == "filter"){
    result = filter_feature_selection(dataset$data, dataset$metadata[,column.class], functions,
                                      validation, repeats)
  }
  result
}


#Recursive Feature Elimination
#funcs list: lmFuncs, rfFuncs, treebagFuncs, ldaFuncs, nbFuncs, gamFuncs, lrFuncs
#' Recursive Feature Elimination
#'
#' Performs recursive feature elimination on a data matrix.
#'
#' @param datamat Data matrix with features in rows and samples in columns.
#' @param samples.class Sample class labels.
#' @param functions Caret RFE functions list.
#' @param method Resampling method passed to \code{caret::rfeControl()}.
#' @param repeats Number of repeats for resampling.
#' @param number Number of resampling folds or iterations.
#' @param subsets Subset sizes to evaluate.
#'
#' @return An \code{rfe} object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("randomForest", quietly = TRUE)) {
#'   datamat <- matrix(
#'     rnorm(20),
#'     nrow = 4,
#'     dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#'   )
#'   classes <- factor(c("A", "A", "B", "B", "A"))
#'   recursive_feature_elimination(datamat, classes, caret::rfFuncs)
#' }
#' }
#' @export
recursive_feature_elimination = function(datamat, samples.class, functions = caret::rfFuncs, method = "cv",
                                         repeats = 5, number = 10, subsets = 2^(2:4)){
  samples.df = data.frame(t(datamat))
  ctrl <- caret::rfeControl(functions = functions,
                            method = method,
                            repeats = repeats,
                            number = number,
                            verbose = FALSE)
  rfe.result <- caret::rfe(samples.df, samples.class,
                           sizes = subsets,
                           rfeControl = ctrl)
  
  rfe.result
}


#Feature Selection Using Univariate Filters
#functions list: lmSBF, rfSBF, treebagSBF, ldaSBF and nbSBF.
#' Feature Selection Using Univariate Filters
#'
#' Performs feature selection using univariate filters.
#'
#' @param datamat Data matrix with features in rows and samples in columns.
#' @param samples.class Sample class labels.
#' @param functions Caret SBF functions list.
#' @param method Resampling method passed to \code{caret::sbfControl()}.
#' @param repeats Number of repeats for resampling.
#'
#' @return An \code{sbf} object.
#'
#' @examples
#' datamat <- matrix(
#'   rnorm(20),
#'   nrow = 4,
#'   dimnames = list(paste0("x", 1:4), paste0("s", 1:5))
#' )
#' classes <- factor(c("A", "A", "B", "B", "A"))
#' filter_feature_selection(datamat, classes, caret::rfSBF)
#' @export
filter_feature_selection = function(datamat, samples.class, functions = caret::rfSBF, method = "cv",
                                    repeats = 5) {
  samples.df = data.frame(t(datamat))
  filterCtrl = caret::sbfControl(functions = functions, method = method, repeats = repeats)
  filter.result = caret::sbf(samples.df, samples.class, sbfControl = filterCtrl)
  filter.result
}