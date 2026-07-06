"compare_regions_by_sample" = function(dataset1, dataset2, fn.to.apply, samples = NULL, ...){
  stats1 = apply_by_sample(dataset1, fn.to.apply, ...)
  stats2 = apply_by_sample(dataset2, fn.to.apply, ...)
  stats.total = data.frame(cbind(stats1,stats2))
  names(stats.total) = c(deparse(substitute(dataset1)), deparse(substitute(dataset2)))
  stats.total$ratio = stats1/stats2
  stats.total
}

# applies a function to the values of each variable
# fn.to.apply - function to apply (e.g. mean, max, min)
# variables - allows to define which variables to calculate the stats (if numbers, indexes are assumed)
# variable.bounds - allow to define an interval of variables (if numeric)
# samples - if defined restricts the application to a given set of samples 
"apply_by_variable" = function(dataset, fn.to.apply, variables = NULL, variable.bounds = NULL, 
                               samples = NULL, ...) {
  
  if (is.null(variables)) {
    if (is.null(variable.bounds)){
      variables = rownames(dataset$data)
    } 
    else {
      x.vars = get_x_values_as_num(dataset)
      variables = rownames(dataset$data)[x.vars > variable.bounds[1] & x.vars < variable.bounds[2]] 
    }  
  }  
  if (is.null(samples)) {
    samples = colnames(dataset$data)
  }
  
  apply(dataset$data[variables, samples, drop = FALSE], 1, fn.to.apply, ...)
}

"apply_by_sample" = function(dataset, fn.to.apply, samples = NULL, ...) {
  if (is.null(samples)) {
    samples = colnames(dataset$data)
  }
  apply(dataset$data[, samples, drop = FALSE], 2, fn.to.apply, ...)
}

"stats_by_variable" = function(dataset, variables = NULL, variable.bounds = NULL) {
  apply_by_variable(dataset, summary, variables, variable.bounds)
}

"stats_by_sample" = function(dataset, samples = NULL) {
  apply_by_sample(dataset, summary, samples)
}

#' Apply by group
#'
#' Applies a function to the variables of samples belonging to a given group.
#'
#' @param dataset Dataset to analyze.
#' @param fn.to.apply Function to apply.
#' @param metadata.var Metadata variable used to define the group.
#' @param var.value Value or values of the metadata variable to select.
#'
#' @return Result of applying the function by variable for the selected group.
#'
#' @examples
#' ## Example of applying a function to a group
#' library(specmine.datasets)
#' data(cachexia)
#' apply.group.result = apply_by_group(cachexia, mean, "Muscle.loss", "control")
#'
#' @keywords apply group
#' @export
"apply_by_group" = function(dataset, fn.to.apply, metadata.var, var.value) {
  indexes = which(dataset$metadata[, metadata.var] %in% var.value)
  apply_by_variable(dataset, fn.to.apply, samples = indexes)
}

#' Apply by groups
#'
#' Applies a function to groups defined by a metadata variable.
#'
#' @param dataset Dataset to analyze.
#' @param metadata.var Metadata variable used to define groups.
#' @param fn.to.apply Function to apply.
#' @param variables Variables to include.
#' @param variable.bounds Optional numeric bounds for variables.
#'
#' @return A table of grouped results.
#'
#' @examples
#' library(specmine.datasets)
#' data(cachexia)
#' apply.groups.result = apply_by_groups(cachexia, "Muscle.loss", mean)
#'
#' @keywords groups apply
#' @export
"apply_by_groups" = function(dataset, metadata.var, fn.to.apply = "mean",
                             variables = NULL, variable.bounds = NULL) {
  
  if (is.null(variables)) {
    if (is.null(variable.bounds)){
      variables = rownames(dataset$data)
    } 
    else {
      x.vars = get_x_values_as_num(dataset)
      variables = rownames(dataset$data)[x.vars > variable.bounds[1] & x.vars < variable.bounds[2]] 
    }  
  }
  df = NULL
  for (v in variables) {
    row = tapply(dataset$data[v,], dataset$metadata[, metadata.var], fn.to.apply)
    if (is.null(df)) df = row
    else df = rbind(df, row)
  }
  rownames(df) = variables
  df
}
