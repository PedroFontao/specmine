##########################
#####DETECT 2D NMR PEAKS#####
##########################

## Internal 2D peak picking function, retrieved from rNMR source code
## Finds points in a matrix that are larger than all surrounding points
## x  -  A numeric matrix containing the range of data to be peak picked
## thresh    - Numeric value specifying the minimum level to be included
## noiseFilt - Integer argument that can be set to 0, 1 or 2;
##              0 does not apply a noise filter, 1 applies a mild filter
##              (adjacent points in the direct dimension must be above the
##              noise threshold), 2 applies a strong filter (all adjacent points
##              must be above the noise threshold)
## Returns a vector of points defining the local maxima
localMax <- function(x, thresh, noiseFilt) {
  
  nC <- ncol(x)
  nR <- nrow(x)
  if (noiseFilt == 2) {
    x[x < thresh] <- NA
  }
  
  ## Find row/column local maxes
  if (noiseFilt == 1) {
    y <- x
    y[y < thresh] <- NA
    vMax <- intersect(which(c(NA, y) < c(y, NA)), which(c(NA, y) > c(y, NA)) - 1)
  } else {
    vMax <- intersect(which(c(NA, x) < c(x, NA)), which(c(NA, x) > c(x, NA)) - 1)
  }
  x <- t(x)
  hMax <- intersect(which(c(NA, x) < c(x, NA)), which(c(NA, x) > c(x, NA)) - 1) - 1
  hMax <- (hMax %% nC * nR) + hMax %/% nC + 1
  
  ## Find diagonal maxima
  x <- t(x)
  hvMax <- intersect(vMax, hMax)
  if (noiseFilt == 0) {
    hvMax <- hvMax[x[hvMax] > thresh]
  }
  dMax <- cbind(
    hvMax, hvMax - nR + 1, hvMax - nR - 1, hvMax + nR + 1,
    hvMax + nR - 1
  )
  dMax[dMax < 1 | dMax > nC * nR] <- NA
  dMax <- which(max.col(cbind(
    x[dMax[, 1]], x[dMax[, 2]], x[dMax[, 3]],
    x[dMax[, 4]], x[dMax[, 5]]
  )) == 1)
  
  hvMax[dMax]
}

# Function to change NA values in a data matrix according to a peak list
peaks_to_dataset <- function(empty_data, peaklst, reference) {
  for (i in seq_len(nrow(peaklst))) {
    if (!(is.na(peaklst$rows[i]) | is.na(peaklst$cols[i]))) {
      empty_data[peaklst$rows[i], peaklst$cols[i]] <- reference[peaklst$rows[i], peaklst$cols[i]]
    }
  }
  empty_data
}

## Function to give a set of peaklists found in 1 spectrum
## Finds pairs of ppms that have a peak
## spectrum  - Numeric matrix. A 2D NMR spectrum.
## threshold - Numeric value. Option to user establish a threshold defining a minimum value to be detected
## noise - Integer argument that can be set to 0, 1 or 2;
##              0 does not apply a noise filter, 1 applies a mild filter
##              (adjacent points in the direct dimension must be above the
##              noise threshold), 2 applies a strong filter (all adjacent points
##              must be above the noise threshold
## Returns a data frame containing the pairs (row/column) and the intensity of the peaks detected
peaklist <- function(spectrum, threshold = NULL, noise = 0) {
  
  if (is.null(spectrum)) {
    stop("Spectrum not provided")
  }
  
  if (!is.null(threshold)) {
    result <- localMax(spectrum, thresh = threshold, noiseFilt = noise)
    points <- spectrum[result]
    locations <- as.matrix(lapply(points, function(x) which(x == spectrum, arr.ind = TRUE)))
    peaks <- data.frame(
      rows = unlist(lapply(locations, function(x) rownames(spectrum)[x[1, ][1]])),
      cols = unlist(lapply(locations, function(x) colnames(spectrum)[x[1, ][2]]))
    )
  } else {
    threshold <- mean(base::Filter(isPositive, spectrum))
    result <- localMax(spectrum, thresh = threshold, noiseFilt = noise)
    points <- spectrum[result]
    locations <- as.matrix(sapply(points, function(x) which(x == spectrum, arr.ind = TRUE)))
    sample_rows <- sapply(locations[1, ], function(x) rownames(spectrum)[x])
    sample_cols <- sapply(locations[2, ], function(x) colnames(spectrum)[x])
    peaks <- data.frame(rows = sample_rows, cols = sample_cols)
  }
  
  peaks
}

## Function to detect peaks in a specmine dataset of 2D NMR data
## Finds peaks across sample, reducing dimensionality
## specmine_2d_dataset  - A list of variables containing at least a 2D matrix for all samples
## thresh - Numeric value. Option to user establish a threshold defining a minimum value to be detected
## noiseFilt - Integer argument that can be set to 0, 1 or 2;
##              0 does not apply a noise filter, 1 applies a mild filter
##              (adjacent points in the direct dimension must be above the
##              noise threshold), 2 applies a strong filter (all adjacent points
##              must be above the noise threshold
## negatives - Boolean value to decide if negative ppm values should be considered or not
## Returns a specmine dataset with only the variables that was found a peak for, a normal 1D
#' Detect peaks in 2D NMR spectra
#'
#' Detects local maxima in each spectrum of a 2D NMR dataset and converts the
#' detected peaks into a 1D specmine dataset suitable for downstream analysis.
#'
#' @param specmine_2d_dataset A specmine 2D dataset containing one numeric matrix per sample.
#' @param baseline_thresh Optional numeric threshold used as the minimum intensity for peak detection.
#'   If \code{NULL}, a threshold is estimated from the positive values of each spectrum.
#' @param noiseFilt Integer noise filter level: \code{0} for no filter, \code{1} for a mild
#'   filter, and \code{2} for a strong filter.
#' @param negatives Logical indicating whether variables with negative ppm values should be kept.
#' @param verbose Logical indicating whether progress messages should be shown.
#'
#' @return A specmine dataset object as a list. The \code{data} element is a numeric matrix where
#'   rows correspond to detected 2D peak coordinates encoded as combined F1/F2 positions and columns
#'   correspond to samples. Matrix entries contain peak intensities, with missing values indicating
#'   that a given peak was not detected in a sample. The returned object also includes the original
#'   sample metadata, the dataset description, the dataset type set to \code{"nmr-peaks"}, and
#'   labels describing the x-axis and data values. This output represents a reduced 1D peak table
#'   derived from the input 2D NMR spectra.
#'
#' @examples
#' m <- matrix(c(0, 5, 0, 0), nrow = 2)
#' rownames(m) <- c("1.0", "1.1")
#' colnames(m) <- c("2.0", "2.1")
#' x <- list(data = list(s1 = m), metadata = data.frame(), description = "Toy")
#' try(peak_detection2d(x, baseline_thresh = 1, negatives = TRUE, verbose = FALSE))
#'
#' @export
peak_detection2d <- function(specmine_2d_dataset,
                             baseline_thresh = NULL,
                             noiseFilt = 0,
                             negatives = FALSE,
                             verbose = TRUE) {
  
  data <- specmine_2d_dataset$data
  res_data <- list()
  
  for (i in seq_along(data)) {
    sample <- names(data)[[i]]
    
    rows <- nrow(data[[i]])
    cols <- ncol(data[[i]])
    
    res_data[[sample]] <- matrix(
      data = rep(NA, rows * cols),
      nrow = rows,
      ncol = cols,
      dimnames = list(rownames(data[[i]]), colnames(data[[i]]))
    )
    
    peaklist_res <- peaklist(
      spectrum = data[[i]],
      threshold = baseline_thresh,
      noise = noiseFilt
    )
    
    if (verbose) {
      message("Sample: ", sample, " has ", nrow(peaklist_res), " peaks")
    }
    
    res_data[[sample]] <- peaks_to_dataset(
      empty_data = res_data[[sample]],
      peaklst = peaklist_res,
      reference = data[[i]]
    )
  }
  
  dim_example <- dim(res_data[[1]])
  logical <- unlist(lapply(lapply(res_data, dim), function(x) identical(x, dim_example)))
  
  if (sum(logical) == length(res_data)) {
    res_data <- simplify2array(res_data)
  } else {
    res_data <- narray::stack(res_data)
  }
  
  if (!negatives) {
    col <- grep("-", colnames(res_data))
    row <- grep("-", rownames(res_data))
    
    if (length(col) > 0) {
      res_data <- res_data[, -col, ]
    }
    
    if (length(row) > 0) {
      res_data <- res_data[-row, , ]
    }
  }
  
  dimnames1 <- as.vector(t(outer(dimnames(res_data)[[1]], dimnames(res_data)[[2]], FUN = paste, sep = ".")))
  data_2d <- matrix(
    res_data,
    prod(dim(res_data)[1:2]),
    dim(res_data)[3],
    dimnames = list(dimnames1, dimnames(res_data)[[3]])
  )
  
  rownames(data_2d) <- base::make.names(rownames(data_2d), unique = TRUE)
  indexes <- which(rowSums(is.na(data_2d)) == ncol(data_2d))
  data_2d <- data_2d[-indexes, , drop = FALSE]
  
  dataset <- specmine::create_dataset(
    data_2d,
    type = "nmr-peaks",
    metadata = specmine_2d_dataset$metadata,
    description = specmine_2d_dataset$description,
    label.x = "F1 x F2 ppm",
    label.values = "intensity",
    sample.names = names(specmine_2d_dataset$data)
  )
  
  dataset
}