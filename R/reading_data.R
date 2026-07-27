# reads a dataset from CSV files: one for data and (optionally) one for metadata

#' Reads a dataset from CSV files
#'
#' Reads a dataset from a CSV file containing the data matrix and, optionally,
#' a second CSV file containing metadata.
#'
#' @param filename.data Path to the CSV file containing the data matrix.
#' @param filename.meta Optional path to the CSV file containing metadata.
#' @param type Character string describing the dataset type.
#' @param description Character string with a dataset description.
#' @param label.x Optional x-axis label.
#' @param label.values Optional value labels.
#' @param sample.names Optional sample names.
#' @param format Data layout format, either \code{"row"} or \code{"col"}.
#' @param header.col Logical; whether the data file has a header row.
#' @param header.row Logical; whether the data file has a row names column.
#' @param sep Field separator used in the data file.
#' @param header.col.meta Logical; whether the metadata file has a header row.
#' @param header.row.meta Logical; whether the metadata file has a row names column.
#' @param sep.meta Field separator used in the metadata file.
#'
#' @return An object of class \code{dataset} created from the input files.
#'   This object contains the imported data matrix and, if provided, the sample
#'   metadata, together with dataset-level information such as type,
#'   description, labels, and sample names. The returned object is intended to
#'   be used as input to downstream processing and analysis functions in the
#'   package.
#'
#' @examples
#' data_file <- tempfile(fileext = ".csv")
#' meta_file <- tempfile(fileext = ".csv")
#'
#' data_in <- data.frame(
#'   x1 = c(1.1, 2.2, 3.3),
#'   x2 = c(4.4, 5.5, 6.6),
#'   row.names = c("s1", "s2", "s3")
#' )
#' metadata_in <- data.frame(
#'   class = c("A", "B", "A"),
#'   row.names = c("s1", "s2", "s3")
#' )
#'
#' utils::write.csv(data_in, data_file)
#' utils::write.csv(metadata_in, meta_file)
#'
#' dataset <- read_dataset_csv(
#'   filename.data = data_file,
#'   filename.meta = meta_file,
#'   format = "row",
#'   header.row = TRUE,
#'   header.row.meta = TRUE
#' )
#' class(dataset)
#'
#' @export
read_dataset_csv <- function(filename.data, filename.meta = NULL, type = "undefined",
                             description = "", label.x = NULL, label.values = NULL,
                             sample.names = NULL, format = "row", header.col = TRUE,
                             header.row = TRUE, sep = ",", header.col.meta = TRUE,
                             header.row.meta = TRUE, sep.meta = ",")
{
  if (!is.null(filename.meta)) {
    metadata <- read_metadata(
      filename.meta,
      header.col = header.col.meta,
      header.row = header.row.meta,
      sep = sep.meta
    )
  } else {
    metadata <- NULL
  }
  
  data <- read_data_csv(
    filename.data,
    format = format,
    header.col = header.col,
    header.row = header.row,
    sep = sep
  )
  
  if (!is.null(sample.names)) {
    dataset <- create_dataset(
      data,
      type = type,
      metadata = metadata,
      description = description,
      sample.names = sample.names,
      label.x = label.x,
      label.values = label.values
    )
  } else {
    dataset <- create_dataset(
      data,
      type = type,
      metadata = metadata,
      description = description,
      label.x = label.x,
      label.values = label.values
    )
  }
  
  dataset
}


# reads a CSV file creating a data matrix
# format: "row" -> samples are in rows (if header.row is TRUE, sample names are in the first column;
#               if header.col is TRUE, value labels are in the first row)
#       "col" -> samples are in columns (if header.col is TRUE, sample names are in the first row;
#               if header.row is TRUE, value labels are in the first column)

#' Reads a data matrix from a CSV file
#'
#' Imports a numeric data table from a CSV file and converts it to the internal
#' matrix structure used by the package.
#'
#' @param filename Path to the CSV file containing the data matrix.
#' @param format Data layout format, either \code{"row"} or \code{"col"}.
#'   If \code{"row"}, samples are assumed to be stored in rows and the matrix
#'   is transposed so that samples become columns in the returned matrix.
#'   If \code{"col"}, the matrix is kept as read.
#' @param header.col Logical; whether the file has a header row.
#' @param header.row Logical; whether the file has a row names column.
#' @param sep Field separator used in the file.
#'
#' @return A numeric matrix containing the imported data in the package's
#'   internal orientation, with variables in rows and samples in columns.
#'   If \code{format = "row"}, the input table is transposed before being
#'   returned. The function stops with an error if non-numeric values are found
#'   in the data table.
#'
#' @examples
#' data_file <- tempfile(fileext = ".csv")
#' data_in <- data.frame(
#'   x1 = c(1.1, 2.2, 3.3),
#'   x2 = c(4.4, 5.5, 6.6),
#'   row.names = c("s1", "s2", "s3")
#' )
#' utils::write.csv(data_in, data_file)
#' read_data_csv(data_file, format = "row", header.row = TRUE)
#'
#' @keywords internal
#' @export
read_data_csv <- function(filename, format = "row", header.col = TRUE,
                          header.row = TRUE, sep = ",")
{
  if (header.row) {
    rownames <- 1
  } else {
    rownames <- NULL
  }
  
  df <- read.table(
    filename,
    header = header.col,
    row.names = rownames,
    sep = sep,
    check.names = FALSE
  )
  
  if (sum(apply(df, c(1, 2), is.numeric)) != nrow(df) * ncol(df)) {
    stop("There are non-numerical values in data")
  }
  
  if (format == "row") {
    data <- as.matrix(t(df))
  } else if (format == "col") {
    data <- as.matrix(df)
  }
  
  data
}

#' Reads metadata from a CSV file
#'
#' Reads metadata from a CSV file and returns it as a data frame.
#'
#' @param filename Path to the metadata CSV file.
#' @param header.col Logical; whether the metadata file has a header row.
#' @param header.row Logical; whether the metadata file has a row names column.
#' @param sep Field separator used in the metadata file.
#'
#' @return A \code{data.frame} containing the imported metadata. Rows usually
#'   correspond to samples and columns correspond to metadata variables.
#'
#' @examples
#' meta_file <- tempfile(fileext = ".csv")
#' metadata_in <- data.frame(
#'   class = c("A", "B", "A"),
#'   batch = c("b1", "b1", "b2"),
#'   row.names = c("s1", "s2", "s3")
#' )
#' utils::write.csv(metadata_in, meta_file)
#' read_metadata(meta_file, header.row = TRUE)
#'
#' @export
read_metadata <- function(filename, header.col = TRUE, header.row = TRUE,
                          sep = ",")
{
  if (header.row) {
    rownames <- 1
  } else {
    rownames <- NULL
  }
  
  metadata <- read.table(
    filename,
    header = header.col,
    row.names = rownames,
    sep = sep,
    check.names = FALSE
  )
  
  metadata
}