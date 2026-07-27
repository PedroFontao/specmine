#' Missing values imputation
#'
#' Impute missing values in a dataset using different methods.
#'
#' @param dataset A dataset object to process.
#' @param method Imputation method: `"value"`, `"mean"`, `"median"`, `"knn"`, or `"linapprox"`.
#' @param value If `method = "value"`, the value used to replace missing entries.
#' @param k If `method = "knn"`, the number of neighbors used for imputation.
#'
#' @return A dataset object with the same overall structure as the input, in which
#'   missing values in `dataset$data` have been imputed according to the selected
#'   method. The returned object preserves the dataset components and updates the
#'   description to record the imputation step.
#'
#' @examples
#' data <- matrix(
#'   c(1, NA, 3, 4, 5, NA),
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2", "s3"))
#' )
#' dataset <- list(data = data, description = "toy dataset")
#' missingvalues_imputation(dataset, method = "value", value = 0)
#'
#' @export
"missingvalues_imputation" = function(dataset, method = "value", value = 0.0005, k = 5){
  if (method == "value"){
    dataset = impute_nas_value(dataset, value)
  } 
  else if (method == "mean"){
    dataset = impute_nas_mean(dataset)
  } 
  else if (method == "median"){
    dataset = impute_nas_median(dataset)
  } 
  else if (method == "knn"){
    dataset = impute_nas_knn(dataset, k)
  } 
  else if (method == "linapprox"){
    dataset = impute_nas_linapprox(dataset)
  }
  add.desc = paste("Missing value imputation with method", method, sep=" ")
  dataset$description = paste(dataset$description, add.desc, sep="; ")
  dataset
}

impute_nas_linapprox <- function(dataset){
  dataset$data <- imputeTS::na_interpolation(dataset$data, option = "linear")
  dataset
}

#' Impute missing values with a constant
#'
#' Replace all `NA` values in `dataset$data` by a user-defined constant.
#'
#' @param dataset A dataset object to modify.
#' @param value A numeric or character value used to replace missing entries.
#'
#' @return A dataset object with the same structure as the input, where all
#'   missing values in `dataset$data` have been replaced by `value`. Other
#'   components of the dataset are preserved unchanged.
#'
#' @examples
#' data <- matrix(
#'   c(1, NA, 3, 4),
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2"))
#' )
#' dataset <- list(data = data)
#' impute_nas_value(dataset, 0)
#'
#' @export
"impute_nas_value" = function(dataset, value)
{
  dataset$data[is.na(dataset$data)] = value
  dataset
}

#' Impute missing values with mean
#'
#' Replace missing values in each variable by the mean of the observed values
#' for that variable.
#'
#' @param dataset A dataset object to modify.
#'
#' @return A dataset object with the same structure as the input, where missing
#'   values in `dataset$data` have been replaced by the mean of the corresponding
#'   variable calculated with `na.rm = TRUE`.
#'
#' @examples
#' data <- matrix(
#'   c(1, NA, 3, 4),
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2"))
#' )
#' dataset <- list(data = data)
#' impute_nas_mean(dataset)
#'
#' @export
"impute_nas_mean" = function(dataset){
  temp = apply(dataset$data, 1, function(x){
    if(sum(is.na(x))>0){
      x[is.na(x)] = mean(x, na.rm=TRUE)
    }
    x
  })
  dataset$data = t(temp)
  dataset
}

#' Impute missing values with median
#'
#' Replace missing values in each variable by the median of the observed values
#' for that variable.
#'
#' @param dataset A dataset object to modify.
#'
#' @return A dataset object with the same structure as the input, where missing
#'   values in `dataset$data` have been replaced by the median of the corresponding
#'   variable calculated with `na.rm = TRUE`.
#'
#' @examples
#' data <- matrix(
#'   c(1, NA, 5, 4),
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2"))
#' )
#' dataset <- list(data = data)
#' impute_nas_median(dataset)
#'
#' @export
"impute_nas_median" = function(dataset) {
  temp = apply(dataset$data, 1, function(x){
    if(sum(is.na(x))>0){
      x[is.na(x)] = median(x, na.rm=TRUE)
    }
    x
  })
  dataset$data = t(temp)
  dataset
}

#' Impute missing values with kNN
#'
#' Replace missing values using k-nearest neighbors imputation.
#'
#' @param dataset A dataset object to modify.
#' @param k Number of neighbors to use in the imputation procedure.
#' @param ... Additional arguments passed to `impute::impute.knn()`.
#'
#' @return A dataset object with the same structure as the input, where missing
#'   values in `dataset$data` have been imputed using the k-nearest neighbors
#'   method implemented in `impute::impute.knn()`. The returned object preserves
#'   the remaining dataset components unchanged.
#'
#' @examples
#' \donttest{
#' data <- matrix(
#'   c(1,  2, NA, 4,
#'     2,  3,  4, 5,
#'     3, NA,  5, 6,
#'     4,  5,  6, 7),
#'   nrow = 4,
#'   byrow = TRUE,
#'   dimnames = list(c("x1", "x2", "x3", "x4"), c("s1", "s2", "s3", "s4"))
#' )
#' dataset <- list(data = data)
#' impute_nas_knn(dataset, k = 2)
#' }
#'
#' @export
"impute_nas_knn" = function(dataset, k = 10, ...){
  dataset$data = impute::impute.knn(dataset$data, k = k, ...)$data
  dataset
}