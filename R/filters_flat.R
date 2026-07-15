# flat pattern filters

#' Flat pattern filter
#'
#' Performs a flat pattern filter over the dataset.
#'
#' @param dataset Dataset to filter.
#' @param filter.function Filtering function to use. One of
#'   `"iqr"`, `"rsd"`, `"rnsd"`, `"sd"`, `"mad"`, `"mean"`, or `"median"`.
#' @param by.percent Logical. If `TRUE`, the number of variables to filter
#'   will be defined as a percentage of the number of variables in the dataset;
#'   percentage is given by `red.value`.
#' @param by.threshold Logical. If `TRUE`, filtering selects variables where
#'   the filtering function is above or equal to a threshold.
#' @param red.value Reduction value. If `by.percent = TRUE`, this is the
#'   percentage of variables to remove, or `"auto"` for automatic calculation.
#'   If `by.threshold = TRUE`, this is the minimum value needed to keep
#'   the variable.
#'
#' @return Filtered dataset.
#'
#' @examples
#' if (requireNamespace("specmine.datasets", quietly = TRUE)) {
#'   data(propolis, package = "specmine.datasets")
#'   propolis_proc = missingvalues_imputation(propolis)
#'   propolis_proc = flat_pattern_filter(propolis_proc, "iqr",
#'     by.percent = TRUE, red.value = 75)
#' }
#'
#' @export
"flat_pattern_filter" = function(dataset, filter.function = "iqr", by.percent = TRUE, 
                                 by.threshold = FALSE, red.value = 0){
  
  if (by.percent & by.threshold) 
    warning("Both by.percent and by.threshold are TRUE; filtering by percentage")
  if (!by.percent & !by.threshold) 
    stop("Either by.percent or by.threshold need to be TRUE")
  
  filter.values = apply_filter_function(dataset$data, filter.function)
  
  if (by.percent)
    dataset$data = flat_pattern_filter_percentage(dataset$data, filter.values, red.value)
  else if (by.threshold)
    dataset$data = flat_pattern_filter_threshold(dataset$data, filter.values, red.value)
  
  add.desc = paste("Flat pattern filtering with function", filter.function, sep = " ")
  dataset$description = paste(dataset$description, add.desc, sep = "; ")
  dataset
}


# FLAT PATTERN FILTERING

# method: iqr, rsd, rnsd, sd, mad, mean, median
apply_filter_function = function(datamat, filter.fn = "iqr"){
  if (filter.fn == "iqr"){
    filter.values = apply(datamat, 1, IQR, na.rm = TRUE)
  }
  else if (filter.fn == "rsd"){
    sds = apply(datamat, 1, sd, na.rm = TRUE)
    mns = apply(datamat, 1, mean, na.rm = TRUE)
    filter.values = abs(sds / mns)
  }
  else if (filter.fn == "rnsd"){
    mads = apply(datamat, 1, mad, na.rm = TRUE)
    meds = apply(datamat, 1, median, na.rm = TRUE)
    filter.values = abs(mads / meds)
  }
  else if (filter.fn == "sd"){
    filter.values = apply(datamat, 1, sd, na.rm = TRUE)
  }
  else if (filter.fn == "mad"){
    filter.values = apply(datamat, 1, mad, na.rm = TRUE)
  }
  else if (filter.fn == "mean"){
    filter.values = apply(datamat, 1, mean, na.rm = TRUE)
  }
  else if (filter.fn == "median"){
    filter.values = apply(datamat, 1, median, na.rm = TRUE)
  }
  else stop("Invalid filter function")
  
  filter.values
}

flat_pattern_filter_percentage = function(datamat, filter.values, percentage = "auto") {
  
  rk = rank(-filter.values, ties.method = "first")
  var.num = nrow(datamat)
  
  if (identical(percentage, "auto")) {
    if (var.num < 250) datamat = datamat[rk < var.num * 0.95, , drop = FALSE]
    else if (var.num < 500) datamat = datamat[rk < var.num * 0.90, , drop = FALSE]
    else if (var.num < 1000) datamat = datamat[rk < var.num * 0.75, , drop = FALSE]
    else datamat = datamat[rk < var.num * 0.6, , drop = FALSE]
  } 
  else if (is.numeric(percentage)) {
    var.remain = rk < var.num * ((100 - percentage) / 100)
    datamat = datamat[var.remain, , drop = FALSE]
  }
  else stop("Invalid value for percentage parameter")
  
  datamat
}

flat_pattern_filter_threshold = function(datamat, filter.values, threshold = 0) {
  var.remain = filter.values >= threshold
  datamat[var.remain, , drop = FALSE]
}