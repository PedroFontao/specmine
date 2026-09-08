#' Find peaks in Raman or SERS spectra
#'
#' Finds local maxima in one-dimensional spectra according to peak properties
#' such as height, threshold, distance, prominence, width and plateau size.
#'
#' The interface and peak-property terminology are inspired by
#' \code{scipy.signal.find_peaks()}. Prominence and width calculations use an
#' R implementation tailored to Raman/SERS spectra and are not intended to
#' reproduce SciPy results exactly.
#'
#' @param x Numeric vector containing the signal intensity values.
#' @param height Numeric value or length-two numeric vector specifying the
#'   required height of peaks. A single value is interpreted as a minimum height.
#'   A vector \code{c(min, max)} defines a closed interval.
#' @param threshold Numeric value or length-two numeric vector specifying the
#'   required vertical distance to neighbouring samples.
#' @param distance Positive integer specifying the minimum horizontal distance,
#'   in samples, between neighbouring peaks.
#' @param prominence Numeric value or length-two numeric vector specifying the
#'   required prominence of peaks.
#' @param width Numeric value or length-two numeric vector specifying the
#'   required width of peaks, in samples.
#' @param plateau_size Numeric value or length-two numeric vector specifying
#'   the required size of the flat top of peaks, in samples.
#' @param wlen Optional positive integer. It defines the maximum number of
#'   samples inspected on each side of a peak when calculating prominence.
#' @param rel_height Numeric value between zero and one used for peak-width
#'   calculations.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{peaks}}{Integer vector of one-based indices of peaks in
#'       \code{x} satisfying the requested conditions.}
#'     \item{\code{properties}}{Data frame containing properties of the
#'       returned peaks, including \code{peak_index} and \code{peak_height}.
#'       Additional columns are included when the respective properties are
#'       calculated.}
#'   }
#'
#' @examples
#' x = c(0, 1, 0, 2, 1, 0, 3, 2, 0, 4, 1)
#'
#' res = raman_find_peaks(
#'   x,
#'   height = 1,
#'   distance = 2
#' )
#'
#' res$peaks
#' res$properties
#'
#' @export
raman_find_peaks = function(
    x,
    height = NULL,
    threshold = NULL,
    distance = NULL,
    prominence = NULL,
    width = NULL,
    plateau_size = NULL,
    wlen = NULL,
    rel_height = 0.5
) {
  if (!is.numeric(x) ||
      !is.atomic(x) ||
      !is.null(dim(x)) ||
      anyNA(x) ||
      any(!is.finite(x))) {
    stop(
      "'x' must be a finite numeric vector without missing values.",
      call. = FALSE
    )
  }
  
  if (!is.null(distance) &&
      (!is.numeric(distance) ||
       length(distance) != 1L ||
       is.na(distance) ||
       !is.finite(distance) ||
       distance < 1 ||
       distance != as.integer(distance))) {
    stop(
      "'distance' must be a single positive integer.",
      call. = FALSE
    )
  }
  
  if (!is.null(wlen) &&
      (!is.numeric(wlen) ||
       length(wlen) != 1L ||
       is.na(wlen) ||
       !is.finite(wlen) ||
       wlen < 1 ||
       wlen != as.integer(wlen))) {
    stop(
      "'wlen' must be NULL or a single positive integer.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(rel_height) ||
      length(rel_height) != 1L ||
      is.na(rel_height) ||
      !is.finite(rel_height) ||
      rel_height < 0 ||
      rel_height > 1) {
    stop(
      "'rel_height' must be a single numeric value between 0 and 1.",
      call. = FALSE
    )
  }
  
  n = length(x)
  
  if (n < 3L) {
    return(
      list(
        peaks = integer(0),
        properties = data.frame()
      )
    )
  }
  
  all_peaks = .find_local_maxima(x)
  
  if (length(all_peaks) == 0L) {
    return(
      list(
        peaks = integer(0),
        properties = data.frame()
      )
    )
  }
  
  properties = data.frame(
    peak_index = all_peaks,
    peak_height = x[all_peaks],
    stringsAsFactors = FALSE
  )
  
  if (!is.null(plateau_size)) {
    plateau_size_interval = .as_interval(
      plateau_size,
      "plateau_size",
      lower = 1
    )
    
    properties = .filter_plateau_size(
      x = x,
      properties = properties,
      plateau_size_interval = plateau_size_interval
    )
  }
  
  if (!is.null(height)) {
    height_interval = .as_interval(
      height,
      "height"
    )
    
    keep = properties$peak_height >= height_interval[1L] &
      properties$peak_height <= height_interval[2L]
    
    properties = properties[keep, , drop = FALSE]
  }
  
  if (!is.null(threshold)) {
    threshold_interval = .as_interval(
      threshold,
      "threshold",
      lower = 0
    )
    
    left_diff = c(NA_real_, diff(x))
    right_diff = c(diff(x), NA_real_)
    
    left_threshold = ifelse(left_diff < 0, -left_diff, 0)
    right_threshold = ifelse(right_diff > 0, right_diff, 0)
    
    properties$left_threshold =
      left_threshold[properties$peak_index]
    
    properties$right_threshold =
      right_threshold[properties$peak_index]
    
    min_threshold = pmin(
      properties$left_threshold,
      properties$right_threshold
    )
    
    keep = min_threshold >= threshold_interval[1L] &
      min_threshold <= threshold_interval[2L]
    
    properties = properties[keep, , drop = FALSE]
  }
  
  if (!is.null(distance) && nrow(properties) > 1L) {
    properties = .filter_by_distance(
      properties = properties,
      distance = distance
    )
  }
  
  if (!is.null(prominence)) {
    prominence_interval = .as_interval(
      prominence,
      "prominence",
      lower = 0
    )
    
    prominence_result = .peak_prominences(
      x = x,
      peaks = properties$peak_index,
      wlen = wlen
    )
    
    properties$prominence = prominence_result$prominences
    properties$left_base = prominence_result$left_bases
    properties$right_base = prominence_result$right_bases
    
    keep = properties$prominence >= prominence_interval[1L] &
      properties$prominence <= prominence_interval[2L]
    
    properties = properties[keep, , drop = FALSE]
  }
  
  if (!is.null(width)) {
    width_interval = .as_interval(
      width,
      "width",
      lower = 0
    )
    
    if (!"prominence" %in% names(properties)) {
      prominence_result = .peak_prominences(
        x = x,
        peaks = properties$peak_index,
        wlen = wlen
      )
      
      properties$prominence = prominence_result$prominences
      properties$left_base = prominence_result$left_bases
      properties$right_base = prominence_result$right_bases
    }
    
    width_result = .peak_widths(
      x = x,
      peaks = properties$peak_index,
      rel_height = rel_height,
      prominence_data = properties
    )
    
    properties$width = width_result$widths
    properties$width_height = width_result$width_heights
    properties$left_ip = width_result$left_ips
    properties$right_ip = width_result$right_ips
    
    keep = properties$width >= width_interval[1L] &
      properties$width <= width_interval[2L]
    
    properties = properties[keep, , drop = FALSE]
  }
  
  properties = properties[
    order(properties$peak_index),
    ,
    drop = FALSE
  ]
  
  list(
    peaks = as.integer(properties$peak_index),
    properties = properties
  )
}


#' Find local maxima in a one-dimensional signal
#'
#' Detects strict maxima and flat-top maxima. Flat peaks are represented by the
#' centre index rounded down, matching the conventional handling of even-sized
#' plateaus.
#'
#' @param x Numeric vector.
#' @return Integer vector of local-maximum indices.
#' @noRd
.find_local_maxima = function(x) {
  n = length(x)
  
  if (n < 3L) {
    return(integer(0))
  }
  
  peaks = integer(0)
  i = 2L
  
  while (i < n) {
    if (x[i] > x[i - 1L]) {
      left_edge = i
      right_edge = i
      
      while (right_edge < n &&
             x[right_edge + 1L] == x[left_edge]) {
        right_edge = right_edge + 1L
      }
      
      if (right_edge < n &&
          x[right_edge] > x[right_edge + 1L]) {
        peaks = c(
          peaks,
          as.integer(floor((left_edge + right_edge) / 2))
        )
      }
      
      i = right_edge + 1L
    } else {
      i = i + 1L
    }
  }
  
  peaks
}


#' Convert user input to a lower/upper interval
#'
#' @param x Numeric value or length-two numeric vector.
#' @param name Property name used in error messages.
#' @param lower Optional lower bound permitted for the interval.
#' @return Numeric vector with lower and upper limits.
#' @noRd
.as_interval = function(
    x,
    name,
    lower = -Inf
) {
  if (is.null(x)) {
    return(c(-Inf, Inf))
  }
  
  if (!is.numeric(x) ||
      !is.atomic(x) ||
      length(x) < 1L ||
      length(x) > 2L ||
      anyNA(x) ||
      any(!is.finite(x))) {
    stop(
      sprintf(
        "'%s' must be a finite numeric value or a length-2 numeric vector.",
        name
      ),
      call. = FALSE
    )
  }
  
  interval = if (length(x) == 1L) c(x, Inf) else x
  
  if (interval[1L] > interval[2L]) {
    stop(
      sprintf(
        "The lower bound of '%s' must not exceed its upper bound.",
        name
      ),
      call. = FALSE
    )
  }
  
  if (interval[1L] < lower) {
    stop(
      sprintf(
        "'%s' must be greater than or equal to %s.",
        name,
        format(lower)
      ),
      call. = FALSE
    )
  }
  
  interval
}


#' Filter peaks by minimum distance
#'
#' @param properties Data frame containing \code{peak_index} and
#'   \code{peak_height}.
#' @param distance Minimum distance between retained peaks.
#' @return Data frame containing the retained peaks.
#' @noRd
.filter_by_distance = function(
    properties,
    distance
) {
  if (nrow(properties) <= 1L) {
    return(properties)
  }
  
  indices = properties$peak_index
  heights = properties$peak_height
  
  ord = order(
    -heights,
    indices
  )
  
  indices_sorted = indices[ord]
  keep = logical(length(indices_sorted))
  kept_indices = integer(0)
  
  for (i in seq_along(indices_sorted)) {
    idx = indices_sorted[i]
    
    if (length(kept_indices) == 0L ||
        min(abs(kept_indices - idx)) >= distance) {
      keep[i] = TRUE
      kept_indices = c(kept_indices, idx)
    }
  }
  
  properties[ord[keep], , drop = FALSE]
}


#' Filter peaks by plateau size
#'
#' @param x Numeric vector.
#' @param properties Data frame containing peak properties.
#' @param plateau_size_interval Length-two numeric vector.
#' @return Data frame containing the retained peaks.
#' @noRd
.filter_plateau_size = function(
    x,
    properties,
    plateau_size_interval
) {
  n = length(x)
  plateau_sizes = numeric(nrow(properties))
  
  for (i in seq_len(nrow(properties))) {
    idx = properties$peak_index[i]
    value = x[idx]
    
    left = idx
    while (left > 1L && x[left - 1L] == value) {
      left = left - 1L
    }
    
    right = idx
    while (right < n && x[right + 1L] == value) {
      right = right + 1L
    }
    
    plateau_sizes[i] = right - left + 1L
  }
  
  properties$plateau_size = plateau_sizes
  
  keep = plateau_sizes >= plateau_size_interval[1L] &
    plateau_sizes <= plateau_size_interval[2L]
  
  properties[keep, , drop = FALSE]
}


#' Calculate approximate peak prominences
#'
#' @param x Numeric vector.
#' @param peaks Integer vector of peak indices.
#' @param wlen Optional positive integer defining the maximum search distance
#'   on each side of a peak.
#' @return List with prominence values and indices of left and right bases.
#' @noRd
.peak_prominences = function(
    x,
    peaks,
    wlen = NULL
) {
  n = length(x)
  m = length(peaks)
  
  prominences = numeric(m)
  left_bases = integer(m)
  right_bases = integer(m)
  
  for (i in seq_len(m)) {
    p = peaks[i]
    h = x[p]
    
    left_offset = if (is.null(wlen)) p - 1L else wlen
    right_offset = if (is.null(wlen)) n - p else wlen
    
    left_start = max(1L, p - left_offset)
    right_end = min(n, p + right_offset)
    
    left_segment = x[left_start:p]
    right_segment = x[p:right_end]
    
    left_base_idx = which.min(left_segment)
    right_base_idx = which.min(right_segment)
    
    left_base_val = left_segment[left_base_idx]
    right_base_val = right_segment[right_base_idx]
    
    left_bases[i] = left_start + left_base_idx - 1L
    right_bases[i] = p + right_base_idx - 1L
    
    prominences[i] = h - max(
      left_base_val,
      right_base_val
    )
  }
  
  list(
    prominences = prominences,
    left_bases = left_bases,
    right_bases = right_bases
  )
}


#' Calculate approximate peak widths
#'
#' Peak widths are calculated from discrete crossings of the width-height level
#' between the left and right prominence bases. Values are expressed in samples.
#'
#' @param x Numeric vector.
#' @param peaks Integer vector of peak indices.
#' @param rel_height Relative height used to calculate widths.
#' @param prominence_data Data frame containing \code{prominence},
#'   \code{left_base} and \code{right_base}.
#' @return List with widths, width heights and crossing indices.
#' @noRd
.peak_widths = function(
    x,
    peaks,
    rel_height,
    prominence_data
) {
  m = length(peaks)
  
  widths = numeric(m)
  width_heights = numeric(m)
  left_ips = numeric(m)
  right_ips = numeric(m)
  
  for (i in seq_len(m)) {
    p = peaks[i]
    h = x[p]
    prom = prominence_data$prominence[i]
    left_base = prominence_data$left_base[i]
    right_base = prominence_data$right_base[i]
    
    width_height = h - prom * rel_height
    width_heights[i] = width_height
    
    left_idx = p
    while (left_idx > left_base &&
           x[left_idx] > width_height) {
      left_idx = left_idx - 1L
    }
    
    right_idx = p
    while (right_idx < right_base &&
           x[right_idx] > width_height) {
      right_idx = right_idx + 1L
    }
    
    left_ips[i] = left_idx
    right_ips[i] = right_idx
    widths[i] = right_idx - left_idx
  }
  
  list(
    widths = widths,
    width_heights = width_heights,
    left_ips = left_ips,
    right_ips = right_ips
  )
}


#' Align detected peaks across spectra into wavenumber bins
#'
#' Groups peaks detected in multiple spectra into wavenumber bins defined by a
#' tolerance. The resulting feature matrix can be used in downstream analyses
#' such as PCA, clustering or classification.
#'
#' Missing peaks are represented by zero in the returned feature matrix.
#'
#' @param peak_results A non-empty list of peak-detection results, one per
#'   sample. Each element must contain \code{sample_id}, \code{wavenumbers},
#'   \code{peaks} and \code{properties}.
#' @param tolerance Positive numeric value defining the half-width, in
#'   wavenumbers, around each reference position.
#' @param reference_wavenumbers Optional non-empty numeric vector of reference
#'   wavenumbers. If \code{NULL}, reference positions are derived from detected
#'   peaks and clustered according to \code{tolerance}.
#' @param feature_method Character string defining the feature summarised within
#'   each bin. One of \code{"max_height"}, \code{"max_prominence"},
#'   \code{"presence"} or \code{"mean_height"}.
#'
#' @return Data frame with one row per reference wavenumber and one column per
#'   sample. The first column, \code{wavenumber}, contains bin reference
#'   positions.
#'
#' @examples
#' wn = seq(400, 1800, by = 1)
#'
#' x1 = dnorm(wn, mean = 1000, sd = 20) +
#'   0.01 * rnorm(length(wn))
#'
#' x2 = dnorm(wn, mean = 1002, sd = 20) +
#'   0.01 * rnorm(length(wn))
#'
#' p1 = raman_find_peaks(x1, height = 0.01, distance = 5)
#' p2 = raman_find_peaks(x2, height = 0.01, distance = 5)
#'
#' peak_results = list(
#'   list(
#'     sample_id = "sample1",
#'     wavenumbers = wn,
#'     peaks = p1$peaks,
#'     properties = p1$properties
#'   ),
#'   list(
#'     sample_id = "sample2",
#'     wavenumbers = wn,
#'     peaks = p2$peaks,
#'     properties = p2$properties
#'   )
#' )
#'
#' features = raman_align_peaks(
#'   peak_results,
#'   tolerance = 4,
#'   feature_method = "max_height"
#' )
#'
#' head(features)
#'
#' @export
raman_align_peaks = function(
    peak_results,
    tolerance = 4,
    reference_wavenumbers = NULL,
    feature_method = "max_height"
) {
  if (!is.list(peak_results) ||
      length(peak_results) == 0L) {
    stop(
      "'peak_results' must be a non-empty list.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(tolerance) ||
      length(tolerance) != 1L ||
      is.na(tolerance) ||
      !is.finite(tolerance) ||
      tolerance <= 0) {
    stop(
      "'tolerance' must be a single positive finite numeric value.",
      call. = FALSE
    )
  }
  
  valid_methods = c(
    "max_height",
    "max_prominence",
    "presence",
    "mean_height"
  )
  
  if (!is.character(feature_method) ||
      length(feature_method) != 1L ||
      is.na(feature_method) ||
      !feature_method %in% valid_methods) {
    stop(
      sprintf(
        "'feature_method' must be one of: %s.",
        paste(valid_methods, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  sample_ids = character(length(peak_results))
  
  for (i in seq_along(peak_results)) {
    res = peak_results[[i]]
    
    if (!is.list(res)) {
      stop(
        sprintf(
          "Element %d of 'peak_results' must be a list.",
          i
        ),
        call. = FALSE
      )
    }
    
    required = c(
      "sample_id",
      "wavenumbers",
      "peaks",
      "properties"
    )
    
    missing = setdiff(required, names(res))
    
    if (length(missing) > 0L) {
      stop(
        sprintf(
          "Element %d of 'peak_results' is missing: %s.",
          i,
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    
    if (!is.character(res$sample_id) ||
        length(res$sample_id) != 1L ||
        is.na(res$sample_id) ||
        !nzchar(res$sample_id)) {
      stop(
        sprintf(
          "'sample_id' in element %d must be one non-empty character string.",
          i
        ),
        call. = FALSE
      )
    }
    
    if (!is.numeric(res$wavenumbers) ||
        anyNA(res$wavenumbers) ||
        any(!is.finite(res$wavenumbers))) {
      stop(
        sprintf(
          "'wavenumbers' in element %d must be a finite numeric vector.",
          i
        ),
        call. = FALSE
      )
    }
    
    if (!is.numeric(res$peaks) ||
        anyNA(res$peaks) ||
        any(res$peaks < 1L) ||
        any(res$peaks > length(res$wavenumbers)) ||
        any(res$peaks != as.integer(res$peaks))) {
      stop(
        sprintf(
          "'peaks' in element %d must contain valid integer indices.",
          i
        ),
        call. = FALSE
      )
    }
    
    if (!is.data.frame(res$properties) ||
        !"peak_height" %in% names(res$properties) ||
        !is.numeric(res$properties$peak_height) ||
        nrow(res$properties) != length(res$peaks)) {
      stop(
        sprintf(
          paste0(
            "'properties' in element %d must be a data frame containing ",
            "one numeric 'peak_height' value per peak."
          ),
          i
        ),
        call. = FALSE
      )
    }
    
    if (feature_method == "max_prominence" &&
        (!"prominence" %in% names(res$properties) ||
         !is.numeric(res$properties$prominence))) {
      stop(
        paste0(
          "'properties' must contain a numeric 'prominence' column when ",
          "feature_method = 'max_prominence'."
        ),
        call. = FALSE
      )
    }
    
    sample_ids[i] = res$sample_id
  }
  
  if (anyDuplicated(sample_ids)) {
    stop(
      "'sample_id' values in 'peak_results' must be unique.",
      call. = FALSE
    )
  }
  
  if (!is.null(reference_wavenumbers)) {
    if (!is.numeric(reference_wavenumbers) ||
        length(reference_wavenumbers) == 0L ||
        anyNA(reference_wavenumbers) ||
        any(!is.finite(reference_wavenumbers))) {
      stop(
        paste0(
          "'reference_wavenumbers' must be NULL or a non-empty finite ",
          "numeric vector."
        ),
        call. = FALSE
      )
    }
    
    reference_wavenumbers = sort(unique(reference_wavenumbers))
  }
  
  all_peak_wn = numeric(0)
  
  for (res in peak_results) {
    all_peak_wn = c(
      all_peak_wn,
      res$wavenumbers[as.integer(res$peaks)]
    )
  }
  
  if (length(all_peak_wn) == 0L &&
      is.null(reference_wavenumbers)) {
    result = data.frame(
      wavenumber = numeric(0),
      check.names = FALSE
    )
    
    for (sample_id in sample_ids) {
      result[[sample_id]] = numeric(0)
    }
    
    return(result)
  }
  
  if (is.null(reference_wavenumbers)) {
    reference_wavenumbers = .cluster_wavenumbers(
      wavenumbers = all_peak_wn,
      tolerance = tolerance
    )
  }
  
  n_ref = length(reference_wavenumbers)
  n_samples = length(peak_results)
  
  feature_matrix = matrix(
    0,
    nrow = n_ref,
    ncol = n_samples,
    dimnames = list(
      as.character(reference_wavenumbers),
      sample_ids
    )
  )
  
  for (j in seq_along(peak_results)) {
    res = peak_results[[j]]
    
    wn_peaks = res$wavenumbers[as.integer(res$peaks)]
    props = res$properties
    
    for (i in seq_len(n_ref)) {
      in_window = abs(
        wn_peaks - reference_wavenumbers[i]
      ) <= tolerance
      
      if (!any(in_window)) {
        next
      }
      
      if (feature_method == "presence") {
        feature_matrix[i, j] = 1
      } else if (feature_method == "max_height") {
        feature_matrix[i, j] = max(
          props$peak_height[in_window]
        )
      } else if (feature_method == "mean_height") {
        feature_matrix[i, j] = mean(
          props$peak_height[in_window]
        )
      } else {
        feature_matrix[i, j] = max(
          props$prominence[in_window]
        )
      }
    }
  }
  
  result = as.data.frame(
    feature_matrix,
    check.names = FALSE
  )
  
  result$wavenumber = reference_wavenumbers
  
  result[
    ,
    c(
      "wavenumber",
      setdiff(names(result), "wavenumber")
    ),
    drop = FALSE
  ]
}


#' Cluster wavenumbers into reference positions
#'
#' @param wavenumbers Numeric vector of detected peak wavenumbers.
#' @param tolerance Positive numeric value specifying maximum separation between
#'   neighbouring positions in the same cluster.
#' @return Numeric vector with the median wavenumber of each cluster.
#' @noRd
.cluster_wavenumbers = function(
    wavenumbers,
    tolerance
) {
  if (!is.numeric(wavenumbers) ||
      length(wavenumbers) == 0L) {
    return(numeric(0))
  }
  
  if (anyNA(wavenumbers) ||
      any(!is.finite(wavenumbers))) {
    stop(
      "'wavenumbers' must contain only finite numeric values.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(tolerance) ||
      length(tolerance) != 1L ||
      is.na(tolerance) ||
      !is.finite(tolerance) ||
      tolerance <= 0) {
    stop(
      "'tolerance' must be a single positive finite numeric value.",
      call. = FALSE
    )
  }
  
  wavenumbers = sort(wavenumbers)
  
  cluster_ids = cumsum(
    c(
      TRUE,
      diff(wavenumbers) > tolerance
    )
  )
  
  cluster_centers = vapply(
    split(wavenumbers, cluster_ids),
    stats::median,
    numeric(1L)
  )
  
  as.numeric(cluster_centers)
}


#' Normalize peak feature matrix
#'
#' Normalizes a peak-feature matrix returned by \code{raman_align_peaks()}.
#'
#' With \code{method = "min-max"}, normalization is carried out row-wise across
#' samples. Rows containing identical values for every sample are set to zero.
#' With \code{method = "L2"}, normalization is performed column-wise for each
#' sample. Columns with zero L2 norm are retained as zero.
#'
#' @param features Data frame with a numeric \code{wavenumber} column and one
#'   or more numeric sample-feature columns.
#' @param method Normalization method: \code{"min-max"} or \code{"L2"}.
#'
#' @return Data frame with the same structure as \code{features}, containing
#'   normalized sample-feature values.
#'
#' @examples
#' features = data.frame(
#'   wavenumber = c(1000, 1200, 1600),
#'   sample1 = c(1, 2, 3),
#'   sample2 = c(2, 4, 6)
#' )
#'
#' raman_normalize_peak_features(features, method = "L2")
#'
#' @export
raman_normalize_peak_features = function(
    features,
    method = "min-max"
) {
  if (!is.data.frame(features) ||
      !"wavenumber" %in% names(features)) {
    stop(
      "'features' must be a data frame containing a 'wavenumber' column.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(features$wavenumber) ||
      anyNA(features$wavenumber) ||
      any(!is.finite(features$wavenumber))) {
    stop(
      "'wavenumber' must be a finite numeric column.",
      call. = FALSE
    )
  }
  
  if (!is.character(method) ||
      length(method) != 1L ||
      is.na(method)) {
    stop(
      "'method' must be a single character string.",
      call. = FALSE
    )
  }
  
  valid_methods = c(
    "min-max",
    "L2"
  )
  
  if (!method %in% valid_methods) {
    stop(
      sprintf(
        "'method' must be one of: %s.",
        paste(valid_methods, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  feat_cols = setdiff(
    names(features),
    "wavenumber"
  )
  
  if (length(feat_cols) == 0L) {
    return(features)
  }
  
  numeric_columns = vapply(
    features[feat_cols],
    is.numeric,
    logical(1L)
  )
  
  if (any(!numeric_columns)) {
    stop(
      sprintf(
        "All feature columns must be numeric; invalid columns: %s.",
        paste(
          feat_cols[!numeric_columns],
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
  
  feature_matrix = as.matrix(
    features[, feat_cols, drop = FALSE]
  )
  
  if (anyNA(feature_matrix) ||
      any(!is.finite(feature_matrix))) {
    stop(
      "Feature columns must not contain missing or non-finite values.",
      call. = FALSE
    )
  }
  
  if (method == "min-max") {
    row_min = apply(
      feature_matrix,
      1L,
      min
    )
    
    row_max = apply(
      feature_matrix,
      1L,
      max
    )
    
    row_range = row_max - row_min
    
    for (i in seq_len(nrow(feature_matrix))) {
      if (row_range[i] > 0) {
        feature_matrix[i, ] =
          (feature_matrix[i, ] - row_min[i]) / row_range[i]
      } else {
        feature_matrix[i, ] = 0
      }
    }
  } else {
    for (j in seq_len(ncol(feature_matrix))) {
      col_norm = sqrt(
        sum(feature_matrix[, j]^2)
      )
      
      if (col_norm > 0) {
        feature_matrix[, j] =
          feature_matrix[, j] / col_norm
      }
    }
  }
  
  features[, feat_cols] = feature_matrix
  
  features
}