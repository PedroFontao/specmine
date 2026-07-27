# Functions that allow to filter the dataset by some criteria

# SUBSET functions - allow to define criteria for the information to keep

# returns dataset with selected set of samples
# samples - vector with indexes or names of the samples to select

#' Subset samples
#'
#' Returns a dataset with a selected set of samples.
#'
#' @param dataset Dataset to subset.
#' @param samples Vector with indexes or names of the samples to select.
#' @param rebuild.factors If TRUE, rebuild factors in metadata.
#'
#' @return A dataset object with the same overall structure as the input, containing
#'   only the selected samples in both `dataset$data` and `dataset$metadata`. If
#'   `rebuild.factors` is `TRUE`, unused factor levels in metadata are dropped.
#'
#' @examples
#' data <- matrix(
#'   c(1, 2, 3, 4, 5, 6),
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(
#'   class = factor(c("A", "B", "A")),
#'   row.names = c("s1", "s2", "s3")
#' )
#' dataset <- list(data = data, metadata = metadata)
#' subset_samples(dataset, c("s1", "s3"))
#'
#' @export
subset_samples = function(dataset, samples, rebuild.factors = TRUE) {
  
  dataset$metadata = dataset$metadata[samples,,drop= FALSE]
  if (rebuild.factors) dataset$metadata = rebuild_factors_df(dataset$metadata)
  
  dataset$data = dataset$data[,samples, drop=FALSE]
  dataset
}

# selects set of samples by the value of a metadata variable

#' Subset samples by metadata values
#'
#' Selects a set of samples by the value of a metadata variable.
#'
#' @param dataset Dataset to subset.
#' @param metadata.varname Metadata variable name.
#' @param values Values to keep.
#'
#' @return A dataset object with the same structure as the input, containing only
#'   samples whose `metadata.varname` value matches one of the requested `values`.
#'
#' @examples
#' data <- matrix(
#'   c(1, 2, 3, 4, 5, 6),
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(
#'   class = factor(c("A", "B", "A")),
#'   row.names = c("s1", "s2", "s3")
#' )
#' dataset <- list(data = data, metadata = metadata)
#' subset_samples_by_metadata_values(dataset, "class", "A")
#'
#' @export
subset_samples_by_metadata_values = function(dataset, metadata.varname, values)
{
  indexes = which(dataset$metadata[,metadata.varname] %in% values)
  subset_samples(dataset, indexes)
}

# selects a random subset of nsamples from the dataset
# returns a new dataset with the selected samples

#' Subset random samples
#'
#' Selects a random subset of samples from the dataset.
#'
#' @param dataset Dataset to subset.
#' @param nsamples Number of samples to select.
#'
#' @return A dataset object with the same structure as the input, containing a
#'   random subset of `nsamples` samples.
#'
#' @examples
#' set.seed(123)
#' data <- matrix(
#'   1:12,
#'   nrow = 3,
#'   dimnames = list(c("x1", "x2", "x3"), c("s1", "s2", "s3", "s4"))
#' )
#' metadata <- data.frame(
#'   class = factor(c("A", "A", "B", "B")),
#'   row.names = c("s1", "s2", "s3", "s4")
#' )
#' dataset <- list(data = data, metadata = metadata)
#' subset_random_samples(dataset, 2)
#'
#' @export
subset_random_samples = function(dataset, nsamples)
{
  indexes = sample(num_samples(dataset), nsamples)
  subset_samples(dataset, indexes)
}

#' Subset x values
#'
#' Selects rows of the dataset by x values or by row index.
#'
#' @param dataset Dataset to subset.
#' @param variables Variables to keep.
#' @param by.index Logical. If TRUE, variables are interpreted as row indexes.
#'
#' @return A dataset object with the same structure as the input, where
#'   `dataset$data` contains only the selected variables or row indexes.
#'
#' @examples
#' data <- matrix(
#'   1:9,
#'   nrow = 3,
#'   dimnames = list(c("10", "20", "30"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(group = c("A", "B", "A"), row.names = c("s1", "s2", "s3"))
#' dataset <- list(data = data, metadata = metadata)
#' subset_x_values(dataset, c("10", "30"))
#'
#' @export
subset_x_values = function(dataset, variables, by.index = FALSE) {
  if (!by.index) {
    variables = as.character(variables)
    indexes = which(rownames(dataset$data) %in% variables)
  }
  else indexes = variables
  dataset$data = dataset$data[indexes,,drop=FALSE]
  dataset
}

#' Subset x values by interval
#'
#' Selects x values within a numeric interval.
#'
#' @param dataset Dataset to subset.
#' @param min.value Minimum x value.
#' @param max.value Maximum x value.
#'
#' @return A dataset object with the same structure as the input, where
#'   `dataset$data` contains only variables whose x values fall between
#'   `min.value` and `max.value`, inclusive.
#'
#' @examples
#' data <- matrix(
#'   1:12,
#'   nrow = 4,
#'   dimnames = list(c("10", "20", "30", "40"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(group = c("A", "B", "A"), row.names = c("s1", "s2", "s3"))
#' dataset <- list(data = data, metadata = metadata)
#' subset_x_values_by_interval(dataset, 15, 35)
#'
#' @export
subset_x_values_by_interval = function(dataset, min.value, max.value)
{
  x.values = get_x_values_as_num(dataset)
  indexes = which(x.values >= min.value & x.values <= max.value)
  subset_x_values(dataset, indexes, by.index = TRUE)
}

#' Subset by samples and x values
#'
#' Selects samples and x values simultaneously.
#'
#' @param dataset Dataset to subset.
#' @param samples Sample indexes.
#' @param variables Variables to keep.
#' @param by.index Logical. If TRUE, variables are interpreted as indexes.
#' @param variable.bounds Optional numeric bounds for x values.
#' @param rebuild.factors If TRUE, rebuild factors in metadata.
#'
#' @return A dataset object with the same overall structure as the input,
#'   containing only the selected samples and selected x values. The returned
#'   object includes the corresponding subset of `dataset$data` and matching
#'   rows of `dataset$metadata`; factor levels are rebuilt when
#'   `rebuild.factors = TRUE`.
#'
#' @examples
#' data <- matrix(
#'   1:12,
#'   nrow = 4,
#'   dimnames = list(c("10", "20", "30", "40"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(
#'   class = factor(c("A", "B", "A")),
#'   row.names = c("s1", "s2", "s3")
#' )
#' dataset <- list(data = data, metadata = metadata)
#' subset_by_samples_and_xvalues(dataset, samples = c(1, 3), variables = c("10", "30"))
#'
#' @export
subset_by_samples_and_xvalues = function(dataset, samples, variables = NULL, by.index = FALSE, 
                                         variable.bounds = NULL, rebuild.factors = TRUE)
{
  if (!by.index) {
    if (is.null(variables)) {
      if (is.null(variable.bounds)) 
        stop("One of variables or variable.bounds parameters needs to be defined")
      else {
        x.values = get_x_values_as_num(dataset)
        x.indexes = which(x.values >= variable.bounds[1] & x.values <= variable.bounds[2])
      }
    }
    else {
      variables = as.character(variables)
      x.indexes = which(rownames(dataset$data) %in% variables)
    }
  }
  else x.indexes = variables
  
  dataset$metadata = dataset$metadata[samples,,drop= FALSE]
  if (rebuild.factors) dataset$metadata = rebuild_factors_df(dataset$metadata)
  
  dataset$data = dataset$data[x.indexes,samples, drop=FALSE]
  dataset
}

#' Subset metadata
#'
#' Selects metadata variables to keep.
#'
#' @param dataset Dataset to subset.
#' @param variables Metadata variables to keep.
#'
#' @return A dataset object with the same structure as the input, where
#'   `dataset$metadata` contains only the selected metadata variables.
#'
#' @examples
#' data <- matrix(
#'   1:6,
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(
#'   class = c("A", "B", "A"),
#'   batch = c(1, 1, 2),
#'   row.names = c("s1", "s2", "s3")
#' )
#' dataset <- list(data = data, metadata = metadata)
#' subset_metadata(dataset, "class")
#'
#' @export
subset_metadata = function(dataset, variables)
{
  dataset$metadata = dataset$metadata[,variables, drop = FALSE]
  dataset
}

# REMOVE functions
# Used to remove SAMPLES, DATA VARIABLES or METADATA

#' Remove data
#'
#' Removes selected samples, data variables, or metadata variables from a dataset.
#'
#' @param dataset Dataset to modify.
#' @param data.to.remove Data to remove.
#' @param type One of "sample", "data", or "metadata".
#' @param by.index Logical. If TRUE, data.to.remove has indexes.
#' @param rebuild.factors Logical. If TRUE, rebuild factors in metadata.
#'
#' @return A dataset object with the same overall structure as the input, where
#'   the requested samples, data variables, or metadata variables have been
#'   removed according to `type`.
#'
#' @examples
#' data <- matrix(
#'   1:9,
#'   nrow = 3,
#'   dimnames = list(c("x1", "x2", "x3"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(class = c("A", "B", "A"), row.names = c("s1", "s2", "s3"))
#' dataset <- list(data = data, metadata = metadata)
#' remove_data(dataset, "s2", type = "sample")
#'
#' @export
remove_data = function(dataset, data.to.remove, type = "sample", by.index = FALSE, rebuild.factors = TRUE) {
  if (type == "sample")
    dataset = remove_samples(dataset, data.to.remove, rebuild.factors)
  else if(type == "data")
    dataset = remove_data_variables(dataset, data.to.remove, by.index)
  else if(type == "metadata")
    dataset = remove_metadata_variables(dataset, data.to.remove)
  else stop("Type of data to remove is undefined")
  dataset
}   

#' Remove samples
#'
#' Removes samples from a dataset.
#'
#' @param dataset Dataset to modify.
#' @param samples.to.remove Samples to remove.
#' @param rebuild.factors Logical. If TRUE, rebuild factors in metadata.
#'
#' @return A dataset object with the same structure as the input, with the
#'   selected samples removed from both `dataset$data` and `dataset$metadata`.
#'   If `rebuild.factors` is `TRUE`, unused factor levels in metadata are dropped.
#'
#' @examples
#' data <- matrix(
#'   1:9,
#'   nrow = 3,
#'   dimnames = list(c("x1", "x2", "x3"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(
#'   class = factor(c("A", "B", "A")),
#'   row.names = c("s1", "s2", "s3")
#' )
#' dataset <- list(data = data, metadata = metadata)
#' remove_samples(dataset, "s2")
#'
#' @export
remove_samples = function(dataset, samples.to.remove, rebuild.factors = TRUE) {
  if (is.numeric(samples.to.remove))
    res = subset_samples(dataset, -samples.to.remove, rebuild.factors = rebuild.factors)
  else {
    indexes.to.remove = which(colnames(dataset$data) %in% samples.to.remove)
    res = subset_samples(dataset, -indexes.to.remove, rebuild.factors = rebuild.factors)
  }
  res
}

#' Remove data variables
#'
#' Removes x variables from a dataset.
#'
#' @param dataset Dataset to modify.
#' @param variables.to.remove Variables to remove.
#' @param by.index Logical. If TRUE, variables.to.remove are indexes.
#'
#' @return A dataset object with the same structure as the input, where the
#'   selected rows of `dataset$data` have been removed. If no matching variables
#'   are found, the input dataset is returned unchanged and a warning is issued.
#'
#' @examples
#' data <- matrix(
#'   1:9,
#'   nrow = 3,
#'   dimnames = list(c("10", "20", "30"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(class = c("A", "B", "A"), row.names = c("s1", "s2", "s3"))
#' dataset <- list(data = data, metadata = metadata)
#' remove_data_variables(dataset, "20")
#'
#' @export
remove_data_variables = function(dataset, variables.to.remove, by.index = FALSE) {
  if (length(variables.to.remove) == 0) {
    warning("No variables to remove")
    return (dataset)
  }
  
  if (!by.index) {
    variables.to.remove = as.character(variables.to.remove)
    indexes.to.remove = which(rownames(dataset$data) %in% variables.to.remove)
  }
  else indexes.to.remove = variables.to.remove
  subset_x_values(dataset, -indexes.to.remove, by.index = TRUE)
}

#' Remove x values by interval
#'
#' Removes x values inside an interval.
#'
#' @param dataset Dataset to modify.
#' @param min.value Minimum x value.
#' @param max.value Maximum x value.
#'
#' @return A dataset object with the same structure as the input, where all
#'   variables with x values between `min.value` and `max.value`, inclusive,
#'   have been removed from `dataset$data`.
#'
#' @examples
#' data <- matrix(
#'   1:12,
#'   nrow = 4,
#'   dimnames = list(c("10", "20", "30", "40"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(class = c("A", "B", "A"), row.names = c("s1", "s2", "s3"))
#' dataset <- list(data = data, metadata = metadata)
#' remove_x_values_by_interval(dataset, 15, 35)
#'
#' @export
remove_x_values_by_interval = function(dataset, min.value, max.value)
{
  x.values = get_x_values_as_num(dataset)
  indexes.to.remove = which(x.values >= min.value & x.values <= max.value)
  subset_x_values(dataset, -indexes.to.remove, by.index = TRUE)
}

#' Remove metadata variables
#'
#' Removes metadata variables from a dataset.
#'
#' @param dataset Dataset to modify.
#' @param variables.to.remove Metadata variables to remove.
#'
#' @return A dataset object with the same structure as the input, where the
#'   selected columns have been removed from `dataset$metadata`. If no matching
#'   metadata fields are found, the dataset is returned unchanged and a warning
#'   is issued.
#'
#' @examples
#' data <- matrix(
#'   1:6,
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(
#'   class = c("A", "B", "A"),
#'   batch = c(1, 1, 2),
#'   row.names = c("s1", "s2", "s3")
#' )
#' dataset <- list(data = data, metadata = metadata)
#' remove_metadata_variables(dataset, "batch")
#'
#' @export
remove_metadata_variables = function(dataset, variables.to.remove)
{
  if (!is.numeric(variables.to.remove))
    indexes.to.remove = which(colnames(dataset$metadata) %in% variables.to.remove)
  else indexes.to.remove = variables.to.remove
  
  if (!is.null(indexes.to.remove) & length(indexes.to.remove) > 0)
    dataset$metadata = dataset$metadata[,-indexes.to.remove, drop = FALSE]
  else warning("No metadata variables removed since no fields matched the criteria")
  dataset
}

# functions to remove samples / variables with NAs

#' Remove samples by NAs
#'
#' Removes samples with too many missing values.
#'
#' @param dataset Dataset to modify.
#' @param max.nas Maximum number of missing values allowed.
#' @param by.percent Logical. If TRUE, max.nas is treated as a percentage.
#'
#' @return A dataset object with the same structure as the input, where samples
#'   whose number of missing values exceeds `max.nas` have been removed from both
#'   `dataset$data` and `dataset$metadata`.
#'
#' @examples
#' data <- matrix(
#'   c(1, NA, 3,
#'     4,  5, NA,
#'     7,  8, 9),
#'   nrow = 3,
#'   dimnames = list(c("x1", "x2", "x3"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(class = c("A", "B", "A"), row.names = c("s1", "s2", "s3"))
#' dataset <- list(data = data, metadata = metadata)
#' remove_samples_by_nas(dataset, max.nas = 0)
#'
#' @export
remove_samples_by_nas = function(dataset, max.nas = 0, by.percent = FALSE)
{
  if (by.percent== TRUE) max.nas = max.nas * num_x_values(dataset) / 100
  res = apply(dataset$data, 2, function(x) sum(is.na(x)))
  to.remove = which(res > max.nas)
  remove_samples(dataset, to.remove)
}

#' Remove samples by NA metadata
#'
#' Removes samples with NA in a given metadata variable.
#'
#' @param dataset Dataset to modify.
#' @param metadata.var Metadata variable name.
#'
#' @return A dataset object with the same structure as the input, where samples
#'   with missing values in the selected metadata variable have been removed.
#'
#' @examples
#' data <- matrix(
#'   1:9,
#'   nrow = 3,
#'   dimnames = list(c("x1", "x2", "x3"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(
#'   class = c("A", NA, "B"),
#'   row.names = c("s1", "s2", "s3")
#' )
#' dataset <- list(data = data, metadata = metadata)
#' remove_samples_by_na_metadata(dataset, "class")
#'
#' @export
remove_samples_by_na_metadata = function(dataset, metadata.var)
{
  to.remove = which(is.na(dataset$metadata[,metadata.var]))
  remove_samples(dataset, to.remove)
}

#' Remove variables by NAs
#'
#' Removes variables with too many missing values.
#'
#' @param dataset Dataset to modify.
#' @param max.nas Maximum number of missing values allowed.
#' @param by.percent Logical. If TRUE, max.nas is treated as a percentage.
#'
#' @return A dataset object with the same structure as the input, where variables
#'   whose number of missing values exceeds `max.nas` have been removed from
#'   `dataset$data`.
#'
#' @examples
#' data <- matrix(
#'   c(1,  2, 3,
#'     NA, 5, NA,
#'     7,  8, 9),
#'   nrow = 3,
#'   dimnames = list(c("x1", "x2", "x3"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(class = c("A", "B", "A"), row.names = c("s1", "s2", "s3"))
#' dataset <- list(data = data, metadata = metadata)
#' remove_variables_by_nas(dataset, max.nas = 1)
#'
#' @export
remove_variables_by_nas = function(dataset, max.nas = 0, by.percent = FALSE)
{
  if (by.percent== TRUE) max.nas = max.nas * num_samples(dataset) / 100
  res = apply(dataset$data, 1, function(x) { sum(is.na(x)) } )
  to.remove = which(res > max.nas)
  remove_data_variables(dataset, to.remove, by.index = TRUE)
}

# aggregate samples 
#' Aggregate samples
#'
#' Aggregate samples according to an aggregate function like mean, median, etc.
#' This can be used to merge replicates.
#'
#' @param dataset List representing the dataset from a metabolomics experiment.
#' @param indexes Index vector with the samples that are going to be aggregated
#'   (e.g. c(1,1,2,2), this index vector will aggregate the first two samples and
#'   the last two samples).
#' @param aggreg.fn Aggregation function (e.g. "mean", "median", etc).
#' @param meta.to.remove Metadata variables to be removed.
#'
#' @return A dataset object in which samples have been aggregated according to
#'   `indexes` and `aggreg.fn`, with updated `data`, `metadata`, and preserved
#'   dataset-level fields such as labels, type, and description.
#'
#' @examples
#' data <- matrix(
#'   c(1, 2, 10, 12,
#'     3, 4, 14, 16),
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2", "s3", "s4"))
#' )
#' metadata <- data.frame(
#'   class = factor(c("A", "A", "B", "B")),
#'   batch = c(1, 1, 2, 2),
#'   row.names = c("s1", "s2", "s3", "s4")
#' )
#' dataset <- list(
#'   data = data,
#'   metadata = metadata,
#'   labels = list(),
#'   type = "example",
#'   description = "toy dataset"
#' )
#' aggregate_samples(dataset, c(1, 1, 2, 2), "mean")
#'
#' @keywords aggregation sample
#' @export
"aggregate_samples" = function(dataset, indexes, aggreg.fn = "mean", meta.to.remove = c()) {
  groups = unique(indexes)
  newdata = matrix(NA, nrow(dataset$data), length(groups))
  rownames(newdata) = rownames(dataset$data)
  colnames(newdata) = vector(mode = "character", length=length(groups) )
  
  newmeta = data.frame()
  for (i in 1:length(dataset$metadata)) {
    if (is.numeric(dataset$metadata[[i]]))
      newmeta[[i]] = vector(mode= "numeric", length=0 )
    else if (is.factor(dataset$metadata[[i]]))
      newmeta[[i]] = factor(levels=levels(dataset$metadata[[i]]) ) 
    else newmeta[[i]] = vector(mode="character", length=length(groups))
  }
  names(newmeta) = names(dataset$metadata)
  
  for (g in 1:length(groups)) {
    to.merge = which(indexes == groups[g])
    newdata[,g] = apply(dataset$data[,to.merge,drop=FALSE], 1, aggreg.fn)
    for (i in 1:length(dataset$metadata)) {
      if (is.numeric(dataset$metadata[[i]])) {
        func= match.fun(aggreg.fn)
        newmeta[g,i] = func(dataset$metadata[to.merge,i])
      }
      else if (is.factor(dataset$metadata[[i]]))
        newmeta[g,i] = levels(newmeta[[i]])[which.max(table(dataset$metadata[to.merge,i]))]
      else 
        newmeta[g,i] = dataset$metadata[to.merge[1],i]
    }
    first.index = to.merge[1]
    colnames(newdata)[g] = colnames(dataset$data)[first.index]
    rownames(newmeta)[g] = colnames(dataset$data)[first.index]
  }
  
  newdataset = list()
  newdataset$data = newdata
  newdataset$metadata = newmeta
  newdataset$labels = dataset$labels
  newdataset$type = dataset$type
  newdataset$description = dataset$description
  newdataset = remove_metadata_variables(newdataset, meta.to.remove)
  newdataset
}

#' Merge data and metadata
#'
#' Merges data and metadata into a data frame.
#'
#' @param dataset Dataset to merge.
#' @param samples Samples to keep.
#' @param metadata.vars Metadata variables to keep.
#' @param x.values X values to keep.
#' @param by.index Logical. If TRUE, x.values are indexes.
#'
#' @return A `data.frame` containing the selected data matrix, transposed so that
#'   samples are rows, combined with the selected metadata columns.
#'
#' @examples
#' data <- matrix(
#'   1:6,
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(
#'   class = c("A", "B", "A"),
#'   batch = c(1, 1, 2),
#'   row.names = c("s1", "s2", "s3")
#' )
#' dataset <- list(data = data, metadata = metadata)
#' merge_data_metadata(dataset, metadata.vars = "class")
#'
#' @export
merge_data_metadata = function(dataset, samples = NULL, metadata.vars = NULL, x.values = NULL, 
                               by.index = FALSE)
{
  if (!is.null(samples) )
    dataset = subset_samples(dataset, samples, rebuild.factors = TRUE)
  if (!is.null(x.values) )
    dataset = subset_x_values(dataset, x.values, by.index = by.index)
  if (!is.null(metadata.vars))
    dataset = subset_metadata(dataset, metadata.vars)
  
  df = as.data.frame(t(dataset$data))
  df = cbind(df, dataset$metadata)
  df
}

# diagnostics

#' Count missing values
#'
#' Returns the total number of missing values in the dataset.
#'
#' @param dataset Dataset to inspect.
#'
#' @return A single numeric value giving the total count of `NA` entries present
#'   in `dataset$data`.
#'
#' @examples
#' data <- matrix(
#'   c(1, NA, 3, 4),
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2"))
#' )
#' metadata <- data.frame(class = c("A", "B"), row.names = c("s1", "s2"))
#' dataset <- list(data = data, metadata = metadata)
#' count_missing_values(dataset)
#'
#' @export
count_missing_values = function(dataset)
{
  sum(is.na(dataset$data))
}

#' Count missing values per sample
#'
#' Returns the number of missing values per sample.
#'
#' @param dataset Dataset to inspect.
#' @param remove.zero If TRUE, removes zero counts.
#'
#' @return A named numeric vector giving the number of missing values for each
#'   sample in `dataset$data`. If `remove.zero` is `TRUE`, only samples with at
#'   least one missing value are returned.
#'
#' @examples
#' data <- matrix(
#'   c(1, NA, 3,
#'     4,  5, NA),
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(class = c("A", "B", "A"), row.names = c("s1", "s2", "s3"))
#' dataset <- list(data = data, metadata = metadata)
#' count_missing_values_per_sample(dataset)
#'
#' @export
count_missing_values_per_sample = function(dataset, remove.zero = TRUE) {
  res = apply(dataset$data, 2, function(x) sum(is.na(x)))
  if (remove.zero) res[res > 0]
  else res
}

#' Count missing values per variable
#'
#' Returns the number of missing values per variable.
#'
#' @param dataset Dataset to inspect.
#' @param remove.zero If TRUE, removes zero counts.
#'
#' @return A named numeric vector giving the number of missing values for each
#'   variable in `dataset$data`. If `remove.zero` is `TRUE`, only variables with
#'   at least one missing value are returned.
#'
#' @examples
#' data <- matrix(
#'   c(1, 2, NA,
#'     4, 5, NA),
#'   nrow = 2,
#'   dimnames = list(c("x1", "x2"), c("s1", "s2", "s3"))
#' )
#' metadata <- data.frame(class = c("A", "B", "A"), row.names = c("s1", "s2", "s3"))
#' dataset <- list(data = data, metadata = metadata)
#' count_missing_values_per_variable(dataset)
#'
#' @export
count_missing_values_per_variable = function(dataset, remove.zero = TRUE) {
  "count_na" = function(x) sum(is.na(x))
  res = apply(dataset$data, 1, count_na)
  if (remove.zero) res[res > 0]
  else res
}