## functions to convert ChemoSpec objects into our own
#' 
#' Auto-exported function: convert_from_chemospec
#'
#' @keywords internal
#' @return A \code{dataset} object converted from a ChemoSpec object, containing the transposed data matrix, metadata, sample names, and axis/unit labels.
#' @examples
#' csobj <- list(
#'   data = matrix(c(10, 20, 30, 40), nrow = 2, byrow = TRUE),
#'   freq = c(1, 2),
#'   groups = c("A", "B"),
#'   names = c("sample1", "sample2"),
#'   unit = c("ppm", "intensity")
#' )
#' dataset <- convert_from_chemospec(csobj)
#' class(dataset)
#' @export

convert_from_chemospec = function(csobj, type = "undefined", description = "") {
  datamatrix = t(csobj$data)
  x.values = csobj$freq
  metadata = data.frame(csobj$groups)
  samplenames = csobj$names
  label.x = csobj$unit[1]
  label.val = csobj$unit[2]
  dataset = create_dataset(datamatrix, type = type, metadata = metadata, description = description, 
                           x.axis.values = x.values, sample.names = samplenames, 
                           label.x = label.x, label.values = label.val)
  dataset
}