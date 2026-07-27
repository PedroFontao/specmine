## FUNCTIONS TO DEFINE AND QUERY DATASET STRUCTURE

list.of.spectral.types = c("nmr-spectra","ir-spectra", "uvv-spectra", "raman-spectra", "fluor-spectra")

list.of.allowed.types = c(list.of.spectral.types, "lcms-spectra", "gcms-spectra", "nmr-peaks", "lcms-peaks", "gcms-peaks", "concentrations", "integrated-data", "undefined")

list.of.2d.spectral.types <- c("2d-nmr", "undefined")

# function to create a dataset from existing objects

#' Create dataset
#'
#' Creates a dataset from a numeric matrix.
#'
#' @param datamatrix Matrix with numerical data; rows are variables and columns are samples.
#' @param type Type of data, such as "nmr-spectra", "nmr-peaks", "ir-spectra", "uvv-spectra", "concentrations", or "undefined".
#' @param metadata Optional metadata as a data frame or matrix.
#' @param description Dataset description.
#' @param sample.names Optional sample names.
#' @param x.axis.values Optional x-axis values.
#' @param label.x Optional x-axis label.
#' @param label.values Optional value label.
#' @param xSet Optional xSet object.
#'
#' @return A list representing a specmine dataset. The returned object contains at least the elements
#'   \code{data}, a numeric matrix with variables in rows and samples in columns; \code{type}, a character
#'   string identifying the dataset type; \code{description}, a character string describing the dataset;
#'   \code{metadata}, a data frame with one row per sample when available; \code{labels}, a list with axis
#'   and value labels when provided; and \code{xSet}, an optional object associated with LC-MS processing.
#'   This object is the standard input structure used by downstream specmine analysis functions.
#'
#' @examples
#' datamatrix <- matrix(
#'   c(1.1, 2.2, 3.3, 4.4, 5.5, 6.6),
#'   nrow = 2,
#'   dimnames = list(NULL, NULL)
#' )
#' metadata <- data.frame(
#'   class = c("A", "B", "A"),
#'   row.names = c("s1", "s2", "s3")
#' )
#' create_dataset(
#'   datamatrix,
#'   type = "concentrations",
#'   metadata = metadata,
#'   sample.names = c("s1", "s2", "s3"),
#'   x.axis.values = c("v1", "v2"),
#'   label.x = "variable",
#'   label.values = "intensity"
#' )
#'
#' @importFrom utils capture.output
#' @export
"create_dataset" = function(datamatrix, type = "undefined", metadata = NULL, description = "",
                            sample.names = NULL, x.axis.values = NULL,
                            label.x = NULL, label.values = NULL, xSet = NULL) {
  
  if (is.null(datamatrix))
    stop("Invalid argument: datamatrix is null")
  
  if (!is.matrix(datamatrix)) {
    if (is.data.frame(datamatrix)) {
      warning("datamatrix is data.frame; converting to matrix")
      datamatrix = as.matrix(datamatrix)
    }
    else stop("Invalid argument: datamatrix is not matrix or data.frame")
  }
  if (!is.numeric(datamatrix))
    stop("datamatrix is not numeric")
  
  if (! type %in% list.of.allowed.types)
    stop("Type of data is not allowed")
  
  if (!is.null(metadata)) {
    if (nrow(metadata) != ncol(datamatrix))
      stop("Number of columns in data matrix not the same as number of rows in metadata")
    if (!is.data.frame(metadata)){
      if (is.matrix(metadata)) metadata = as.data.frame(metadata)
      else stop("metadata is not matrix or data.frame")
    }
  }
  else warning("Metadata is null; dataset will still be created with empty metadata")
  
  if (!is.null(label.x) | !is.null(label.values) )
    labels = list(x = label.x, val = label.values)
  else {
    labels = NULL
    warning("Labels are null")
  }
  
  if (!is.null(sample.names)) {
    if (length(sample.names) != ncol(datamatrix) )
      stop("Number of columns in data matrix not the same as length of sample names vector")
    colnames(datamatrix) = sample.names
    
    if (!is.null(metadata)) {
      if(length(sample.names) != nrow(metadata))
        stop("Number of rows in metadata not the same as length of sample names vector")
      rownames(metadata) = sample.names
    }
  }
  else {
    if (is.null(colnames(datamatrix))) {
      warning("Sample names not specified; will be assumed as sequential numbers")
      colnames(datamatrix) = as.character(1:ncol(datamatrix))
      if (!is.null(metadata)) rownames(metadata) = as.character(1:nrow(metadata))
    }
  }
  
  if (!is.null(x.axis.values)) {
    if (length(x.axis.values) != nrow(datamatrix))
      stop("Number of rows in data matrix not the same as length of x axis values vector")
    if (type %in% list.of.spectral.types & any(is.na(as.numeric(x.axis.values))) )
      stop("Invalid non numeric values for variable names in x.axis.values parameter (given spectral type)")
    rownames(datamatrix) = as.character(x.axis.values)
  }
  else {
    if (is.null(rownames(datamatrix))) {
      warning("data variable names not specified; will be assumed as sequential numbers")
      rownames(datamatrix) = as.character(1:nrow(datamatrix))
    }
    else if (type %in% list.of.spectral.types)
      if(any(is.na(as.numeric(rownames(datamatrix)))) )
        stop("Invalid non numeric values for variable names in rownames of matrix (given spectral type)")
  }
  
  if (!is.null(metadata)){
    if (!is.null(rownames(metadata))){
      metadata.ordered = data.frame(metadata[match(colnames(datamatrix),rownames(metadata)),])
      colnames(metadata.ordered) = colnames(metadata)
      rownames(metadata.ordered) = colnames(datamatrix)
      metadata = metadata.ordered
    } else {
      rownames(metadata) = colnames(datamatrix)
    }
  }
  
  dataset = list(data = datamatrix, type = type, description = description, metadata = metadata, labels = labels, xSet = xSet)
  
  dup.indexes = which(duplicated(rownames(dataset$data)))
  if (length(dup.indexes) != 0){
    dataset = remove_data_variables(dataset, dup.indexes, by.index = TRUE)
  }
  
  dataset
}

"create_2d_dataset" <- function(list_2d, type = "undefined", metadata = NULL, description = "",
                                sample.names = NULL, F1 = NULL, F2 = NULL, label.x = NULL,
                                label.y = NULL, label.values = NULL) {
  
  if (is.null(list_2d))
    stop("Invalid argument: list_2d is null")
  
  if (!is.list(list_2d)) {
    if (!any(unlist(lapply(list_2d, is.matrix)))) {
      indexes = which(unlist(lapply(list_2d, is.data.frame)), TRUE)
      warning("Some spectra are data frames, converting them to matrices")
      for (ind in indexes){
        list_2d[[ind]] <- as.matrix(list_2d[[ind]])
      }
    }
    else stop("Invalid argument: list_2d is not a list")
  }
  
  if (!any(unlist(lapply(list_2d, is.numeric))))
    stop("There is a non numeric spectra")
  
  if (! type %in% list.of.2d.spectral.types)
    stop("Type of data is not allowed")
  
  if (!is.null(metadata)) {
    if (nrow(metadata) != length(list_2d))
      stop("Number of 2D spectra samples is not the same as number of rows in metadata")
    if (!is.data.frame(metadata)){
      if (is.matrix(metadata)) metadata <- as.data.frame(metadata)
      else stop("metadata is not matrix or data.frame")
    }
  }
  else warning("Metadata is null; dataset will still be created with empty metadata")
  
  if (!is.null(label.x) | !is.null(label.values) | !is.null(label.y))
    labels = list(x = label.x, y = label.y, val = label.values)
  else {
    labels = NULL
    warning("Labels are null")
  }
  
  if (!is.null(sample.names)) {
    if (length(sample.names) != length(list_2d) )
      stop("Number of 2D spectra samples is not the same as length of sample names vector")
    names(list_2d) <- sample.names
    
    if (!is.null(metadata)) {
      if(length(sample.names) != nrow(metadata))
        stop("Number of rows in metadata not the same as length of sample names vector")
      rownames(metadata) <- sample.names
    }
  }
  else {
    if (is.null(names(list_2d))) {
      warning("Sample names not specified; will be assumed as sequential numbers")
      names(list_2d) <- as.character(1:length(list_2d))
      if (!is.null(metadata)) rownames(metadata) <- as.character(1:length(list_2d))
    }
  }
  
  if (!is.null(F1)) {
    if (length(F1) != nrow(list_2d[[1]]))
      stop("Number of rows in 2D spectra not the same as length of F1 dimension")
    if (type %in% list.of.spectral.types & any(is.na(as.numeric(F1))) )
      stop("Invalid non numeric values for variable names in F1 parameter (given spectral type)")
    for (i in 1:length(list_2d)){
      rownames(list_2d[[i]]) <- as.character(F1)
    }
  }
  else {
    if (is.null(rownames(list_2d[[1]]))) {
      warning("F1 dimension range not specified; will be assumed as sequential numbers")
      for (i in 1:length(list_2d)){
        rownames(list_2d[[i]]) <- as.character(1:nrow(list_2d[[i]]))
      }
    }
    else if (type %in% list.of.spectral.types)
      if(any(is.na(as.numeric(rownames(list_2d[[1]])))) )
        stop("Invalid non numeric values for variable names in rownames of 1st spectra (given spectral type)")
  }
  
  if (!is.null(F2)) {
    if (length(F2) != ncol(list_2d[[1]]))
      stop("Number of columns in 2D spectra not the same as length of F2 dimension")
    if (type %in% list.of.2d.spectral.types & any(is.na(as.numeric(F2))) )
      stop("Invalid non numeric values for variable names in F2 parameter (given spectral type)")
    for (i in 1:length(list_2d)){
      if (length(F2) == dim(list_2d[[i]])[2]){
        colnames(list_2d[[i]]) <- as.character(F2)
      }
    }
  }
  else {
    if (is.null(colnames(list_2d[[1]]))) {
      warning("F2 dimension range not specified; will be assumed as sequential numbers")
      for (i in 1:length(list_2d)){
        colnames(list_2d[[i]]) <- as.character(1:ncol(list_2d[[i]]))
      }
    }
    else if (type %in% list.of.spectral.types)
      if(any(is.na(as.numeric(colnames(list_2d[[1]])))) )
        stop("Invalid non numeric values for variable names in colnames of 1st spectra (given spectral type)")
  }
  
  if (!is.null(metadata)){
    if (!is.null(rownames(metadata))){
      metadata.ordered <- data.frame(metadata[match(names(list_2d),rownames(metadata)),])
      colnames(metadata.ordered) <- colnames(metadata)
      rownames(metadata.ordered) <- names(list_2d)
      metadata = metadata.ordered
    } else {
      rownames(metadata) <- names(list_2d)
    }
  }
  
  dataset <- list(data = list_2d, type = type, description = description, metadata = metadata, F1_ppm = F1, F2_ppm = F2, labels = labels)
  
  dataset
}

"check_dataset" = function(dataset, verbose = FALSE)
{
  if (is.null(dataset$data))
    stop("Invalid dataset: Data matrix is null")
  
  if (!is.null(dataset$metadata)) {
    if (nrow(dataset$metadata) != ncol(dataset$data) )
      stop("Invalid dataset: Number of columns in data matrix not the same as number of rows in metadata")
  }
  else warning("Metadata is null")
  
  if (!dataset$type %in% list.of.allowed.types) stop("Type of data is not allowed")
  
  if (dataset$type %in% list.of.spectral.types)
    if (any(is.na(as.numeric(rownames(dataset$data)))) )
      stop("Invalid non numeric values for variable names in rownames of matrix (given spectral type)")
  
  if (verbose) message("Valid dataset")
  invisible(TRUE)
}

"check_2d_dataset" <- function(dataset_2d, verbose = FALSE) {
  if (!is.null(dataset_2d$data)){
    if (any(unlist(lapply(dataset_2d$data, is.null)))){
      nulls <- which(unlist(lapply(dataset_2d$data, is.null)), TRUE)
      warning(paste("Spectra", nulls, "are null\n"))
    }
  }
  else stop("Invalid dataset: 2D Spectra List is null")
  
  if (!is.null(dataset_2d$metadata)) {
    if (nrow(dataset_2d$metadata) != length(dataset_2d$data) )
      stop("Invalid dataset: Number of 2D spectra samples in dataset not the same as number of rows in metadata")
  }
  else warning("Metadata is null")
  
  if (dataset_2d$type %in% list.of.2d.spectral.types) {
    if (any(is.na(as.numeric(rownames(dataset_2d$data[[1]]))))) {
      stop("Invalid non numeric values for variable names in rownames of 1st spectra (given spectral type)")
    }
  } else if (any(is.na(as.numeric(colnames(dataset_2d$data[[1]]))))) {
    stop("Invalid non numeric values for variable names in colnames of 1st spectra (given spectral type)")
  }
  
  if (verbose) message("Valid dataset")
  invisible(TRUE)
}

"sum_dataset" = function(dataset, stats = TRUE, verbose = TRUE)
{
  check_dataset(dataset, verbose = FALSE)
  
  out <- list(
    description = dataset$description,
    type = dataset$type,
    number.of.samples = ncol(dataset$data),
    number.of.data.points = nrow(dataset$data),
    number.of.metadata.variables = if (!is.null(dataset$metadata)) ncol(dataset$metadata) else NULL,
    label.x = if (!is.null(dataset$labels) && !is.null(dataset$labels$x)) as.character(dataset$labels$x) else NULL,
    label.values = if (!is.null(dataset$labels) && !is.null(dataset$labels$val)) as.character(dataset$labels$val) else NULL
  )
  
  if (stats) {
    out$statistics <- list(
      number.of.missing.values = sum(is.na(dataset$data)),
      mean = mean(dataset$data, na.rm = TRUE),
      median = median(dataset$data, na.rm = TRUE),
      standard.deviation = sd(dataset$data, na.rm = TRUE),
      range = range(dataset$data, na.rm = TRUE),
      quantiles = quantile(dataset$data, na.rm = TRUE)
    )
  }
  
  if (verbose) {
    out_lines <- c(
      "Dataset summary:",
      paste0("Description: ", out$description),
      paste0("Type of data: ", out$type),
      paste0("Number of samples: ", out$number.of.samples),
      paste0("Number of data points: ", out$number.of.data.points)
    )
    
    if (!is.null(out$number.of.metadata.variables)) {
      out_lines <- c(out_lines, paste0("Number of metadata variables: ", out$number.of.metadata.variables))
    }
    if (!is.null(out$label.x)) {
      out_lines <- c(out_lines, paste0("Label of x-axis values: ", out$label.x))
    }
    if (!is.null(out$label.values)) {
      out_lines <- c(out_lines, paste0("Label of data points: ", out$label.values))
    }
    
    if (stats) {
      out_lines <- c(
        out_lines,
        paste0("Number of missing values in data: ", out$statistics$number.of.missing.values),
        paste0("Mean of data values: ", out$statistics$mean),
        paste0("Median of data values: ", out$statistics$median),
        paste0("Standard deviation: ", out$statistics$standard.deviation),
        paste0("Range of values: ", paste(out$statistics$range, collapse = ", ")),
        "Quantiles:"
      )
      out_lines <- c(out_lines, capture.output(out$statistics$quantiles))
    }
    
    message(paste(out_lines, collapse = "\n"))
  }
  
  invisible(out)
}

"sum_2d_dataset" <- function(dataset_2d, stats = TRUE, verbose = TRUE)
{
  check_2d_dataset(dataset_2d, verbose = FALSE)
  
  out <- list(
    description = dataset_2d$description,
    type = dataset_2d$type,
    number.of.samples = length(dataset_2d$data),
    number.of.data.points = nrow(dataset_2d$data[[1]]) * ncol(dataset_2d$data[[1]]),
    number.of.metadata.variables = if (!is.null(dataset_2d$metadata)) ncol(dataset_2d$metadata) else NULL,
    label.x = if (!is.null(dataset_2d$labels) && !is.null(dataset_2d$labels$x)) as.character(dataset_2d$labels$x) else NULL,
    label.y = if (!is.null(dataset_2d$labels) && !is.null(dataset_2d$labels$y)) as.character(dataset_2d$labels$y) else NULL,
    label.values = if (!is.null(dataset_2d$labels) && !is.null(dataset_2d$labels$val)) as.character(dataset_2d$labels$val) else NULL
  )
  
  if (stats) {
    out$statistics <- list(
      number.of.missing.values = unlist(lapply(dataset_2d$data, function(x) sum(is.na(x)))),
      mean = unlist(lapply(dataset_2d$data, function(x) mean(x, na.rm = TRUE))),
      median = unlist(lapply(dataset_2d$data, function(x) median(x, na.rm = TRUE))),
      standard.deviation = unlist(lapply(dataset_2d$data, function(x) sd(x, na.rm = TRUE)))
    )
  }
  
  if (verbose) {
    out_lines <- c(
      "Dataset summary:",
      paste0("Description: ", out$description),
      paste0("Type of data: ", out$type),
      paste0("Number of samples: ", out$number.of.samples),
      paste0("Number of data points: ", out$number.of.data.points)
    )
    
    if (!is.null(out$number.of.metadata.variables)) {
      out_lines <- c(out_lines, paste0("Number of metadata variables: ", out$number.of.metadata.variables))
    }
    if (!is.null(out$label.x)) {
      out_lines <- c(out_lines, paste0("Label of x-axis values: ", out$label.x))
    }
    if (!is.null(out$label.y)) {
      out_lines <- c(out_lines, paste0("Label of y-axis values: ", out$label.y))
    }
    if (!is.null(out$label.values)) {
      out_lines <- c(out_lines, paste0("Label of pair'(x,y) values: ", out$label.values))
    }
    
    if (stats) {
      out_lines <- c(
        out_lines,
        "Number of missing values in data:",
        capture.output(out$statistics$number.of.missing.values),
        "Mean of data values:",
        capture.output(out$statistics$mean),
        "Median of data values:",
        capture.output(out$statistics$median),
        "Standard deviation:",
        capture.output(out$statistics$standard.deviation)
      )
    }
    
    message(paste(out_lines, collapse = "\n"))
  }
  
  invisible(out)
}

# QUERY functions

"get_data" = function(dataset)
{
  dataset$data
}

"get_data_as_df" = function(dataset)
{
  as.data.frame(dataset$data)
}

"get_sample_2d_data" <- function(dataset_2d, sample) {
  dataset_2d$data[[sample]]
}

"get_metadata" = function(dataset)
{
  dataset$metadata
}

"get_metadata_var" = function(dataset, var)
{
  dataset$metadata[,var]
}

"num_samples" = function(dataset)
{
  ncol(dataset$data)
}

"get_sample_names" = function(dataset)
{
  sample.names = colnames(dataset$data)
  sample.names
}

"num_x_values" = function(dataset)
{
  nrow(dataset$data)
}

#' Get x values as text
#'
#' Returns the x values of a dataset as text.
#'
#' @param dataset Dataset object.
#'
#' @return A character vector containing the variable identifiers stored in
#'   `rownames(dataset$data)`. Each element corresponds to one row of the data
#'   matrix and represents the x-axis value or variable label associated with
#'   that feature.
#'
#' @examples
#' datamatrix <- matrix(
#'   c(1, 2, 3, 4),
#'   nrow = 2,
#'   dimnames = list(c("10.5", "11.0"), c("s1", "s2"))
#' )
#' metadata <- data.frame(class = c("A", "B"), row.names = c("s1", "s2"))
#' dataset <- create_dataset(
#'   datamatrix,
#'   type = "concentrations",
#'   metadata = metadata
#' )
#' get_x_values_as_text(dataset)
#'
#' @export
"get_x_values_as_text" = function(dataset)
{
  x.values = rownames(dataset$data)
  as.character(x.values)
}

"get_x_values_as_num" = function(dataset)
{
  x.values = rownames(dataset$data)
  res = as.numeric(x.values)
  if (any(is.na(res))) stop("Variable labels are not all numeric")
  res
}

#' Get x label
#'
#' Returns the x-axis label associated with the dataset.
#'
#' @param dataset Dataset object.
#'
#' @return A character string giving the label of the x axis stored in
#'   `dataset$labels$x`. If no x-axis label has been defined, the function
#'   returns an empty string.
#'
#' @examples
#' datamatrix <- matrix(
#'   c(1, 2, 3, 4),
#'   nrow = 2,
#'   dimnames = list(c("10.5", "11.0"), c("s1", "s2"))
#' )
#' metadata <- data.frame(class = c("A", "B"), row.names = c("s1", "s2"))
#' dataset <- create_dataset(
#'   datamatrix,
#'   type = "concentrations",
#'   metadata = metadata,
#'   label.x = "ppm",
#'   label.values = "intensity"
#' )
#' get_x_label(dataset)
#'
#' @export
"get_x_label" = function(dataset) {
  if (is.null(dataset$labels) | is.null(dataset$labels$x)) return("")
  else return(dataset$labels$x)
}

"get_value_label" = function(dataset) {
  if (is.null(dataset$labels) | is.null(dataset$labels$val)) return("")
  else return(dataset$labels$val)
}

"get_type" = function(dataset) {
  dataset$type
}

"is_spectra" = function(dataset) {
  dataset$type %in% list.of.spectral.types
}

"get_data_value" = function(dataset, x.axis.val, sample, by.index = FALSE) {
  if (!by.index) {
    x.axis.val = as.character(x.axis.val)
    x.axis.index = which(rownames(dataset$data) == x.axis.val)
  }
  else {
    x.axis.index = x.axis.val
  }
  dataset$data[x.axis.index, sample]
}

"get_metadata_value" = function(dataset, variable, sample)
{
  dataset$metadata[sample, variable]
}

"get_data_values" = function(dataset, x.axis.val, by.index = FALSE)
{
  if (!by.index) {
    if (length(x.axis.val) >= 1) {
      x.axis.val = as.character(x.axis.val)
      x.axis.indexes = which(rownames(dataset$data) %in% x.axis.val)
    }
    else
      stop("Incorrect parameter x.axis.val: length not >= 1")
  }
  else {
    x.axis.indexes = x.axis.val
  }
  
  dataset$data[x.axis.indexes,]
}

"x_values_to_indexes" = function(dataset, x.values)
{
  x.values.ds = get_x_values_as_num(dataset)
  indexes = which(x.values.ds %in% x.values)
  indexes
}

"xvalue_interval_to_indexes" = function(dataset, min.value, max.value) {
  x.values = get_x_values_as_num(dataset)
  indexes = which(x.values >= min.value & x.values <= max.value)
  indexes
}

"indexes_to_xvalue_interval" = function(dataset, indexes) {
  x.values = get_x_values_as_num(dataset)
  x.val.inds = x.values[indexes]
  c(min(x.val.inds), max(x.val.inds))
}

variables_as_metadata = function(dataset, variables, by.index = FALSE){
  if (!by.index) {
    var.indexes = which(rownames(dataset$data) %in% variables)
  }
  else {
    var.indexes = variables
  }
  vars = t(dataset$data)[,var.indexes]
  
  if (!is.null(dataset$metadata)){
    metadata = dataset$metadata
    metadata.names = c(colnames(metadata), rownames(dataset$data)[var.indexes])
    metadata = cbind(metadata, vars)
  } else {
    metadata = vars
    metadata.names = rownames(dataset$data)[var.indexes]
  }
  metadata = as.data.frame(metadata)
  colnames(metadata) = metadata.names
  rownames(metadata) = colnames(dataset$data)
  dataset = set_metadata(dataset, metadata)
  dataset$data = dataset$data[-var.indexes,]
  dataset
}

metadata_as_variables = function(dataset, metadata.vars, by.index = FALSE){
  if (!by.index){
    metadata.indexes = which(colnames(dataset$metadata) %in% metadata.vars)
  } else {
    metadata.indexes = metadata.vars
  }
  metadata.variables = dataset$metadata[,metadata.indexes]
  var.names = colnames(dataset$metadata)[metadata.indexes]
  var.names2 = colnames(dataset$metadata)
  metadata.variables = t(as.matrix(metadata.variables))
  rownames(metadata.variables) = var.names
  dataset$data = rbind(dataset$data, metadata.variables)
  dataset$metadata = data.frame(dataset$metadata[,-metadata.indexes])
  colnames(dataset$metadata) = setdiff(var.names2, var.names)
  if (ncol(dataset$metadata) == 0) dataset$metadata = NULL
  dataset
}

"set_metadata" = function(dataset, new.metadata)
{
  if (nrow(new.metadata) != ncol(dataset$data))
    stop("Number of columns in data matrix not the same as number of rows in metadata")
  if (!is.data.frame(new.metadata)){
    if (is.matrix(new.metadata)) new.metadata = as.data.frame(new.metadata)
    else stop("metadata is not matrix or data.frame")
  }
  
  dataset$metadata = new.metadata
  dataset
}

"set_x_values" = function(dataset, new.x.values, new.x.label = NULL)
{
  if (length(new.x.values) != nrow(dataset$data) )
    stop("Length of new vector is not consistent with dataset")
  rownames(dataset$data) = as.character(new.x.values)
  if (!is.null(new.x.label))
    dataset = set_x_label(dataset, new.x.label)
  dataset
}

"set_x_label" = function(dataset, new.x.label)
{
  if (!is.null(dataset$label))
    dataset$labels$x = new.x.label
  else {
    dataset$labels = list()
    dataset$labels$x = new.x.label
  }
  dataset
}

"set_value_label" = function(dataset, new.val.label)
{
  if (!is.null(dataset$label))
    dataset$labels$val = new.val.label
  else {
    dataset$labels = list()
    dataset$labels$val = new.val.label
  }
  dataset
}

"set_sample_names" = function(dataset, new.sample.names)
{
  if (length(new.sample.names) != ncol(dataset$data))
    stop("Length of new sample names not consistent with dataset dimensions")
  
  colnames(dataset$data) = new.sample.names
  rownames(dataset$metadata) = new.sample.names
  
  dataset
}

"replace_data_value" = function(dataset, x.axis.val, sample, new.value, by.index = FALSE) {
  if (!by.index) {
    x.axis.val = as.character(x.axis.val)
    x.axis.index = which(rownames(dataset$data) == x.axis.val)
  }
  else {
    x.axis.index = x.axis.val
  }
  dataset$data[x.axis.index, sample] = new.value
  dataset
}

"replace_metadata_value" = function(dataset, variable, sample, new.value)
{
  dataset$metadata[sample, variable] = new.value
  dataset
}

"convert_to_factor" = function(dataset, metadata.var)
{
  dataset$metadata[,metadata.var] = factor(dataset$metadata[,metadata.var])
  dataset
}

"merge_datasets" = function(dataset1, dataset2)
{
  if (ncol(dataset1$metadata) != ncol(dataset2$metadata))
    stop("Different number of metadata variables")
  if (nrow(dataset1$data) != nrow(dataset2$data))
    stop("Different number of data variables")
  
  dataset1$data = cbind(dataset1$data, dataset2$data)
  dataset1$metadata = rbind(dataset1$metadata, dataset2$metadata)
  
  dataset1
}