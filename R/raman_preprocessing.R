#' Crop Raman or SERS spectra to a wavenumber interval
#'
#' Restricts a Raman or surface-enhanced Raman spectroscopy (SERS) dataset to
#' an inclusive wavenumber interval.
#'
#' The dataset must contain wavenumbers in rows and samples in columns.
#' Wavenumbers are obtained from \code{dataset$wavenumbers} when available.
#' Otherwise, numeric row names in \code{dataset$data} are used.
#'
#' @param dataset A dataset list containing a numeric matrix in
#'   \code{dataset$data}, with wavenumbers in rows and samples in columns.
#'   Optionally, it may contain a numeric vector in
#'   \code{dataset$wavenumbers}.
#' @param min_wn Single finite numeric value specifying the lower inclusive
#'   wavenumber limit.
#' @param max_wn Single finite numeric value specifying the upper inclusive
#'   wavenumber limit.
#'
#' @return A dataset equivalent to \code{dataset}, with \code{dataset$data}
#'   restricted to wavenumbers between \code{min_wn} and \code{max_wn}.
#'   The \code{dataset$wavenumbers} element is updated or created to match the
#'   cropped data.
#'
#' @examples
#' dataset = list(
#'   data = matrix(
#'     seq_len(30),
#'     nrow = 6,
#'     dimnames = list(
#'       as.character(seq(400, 900, by = 100)),
#'       paste0("sample", 1:5)
#'     )
#'   ),
#'   wavenumbers = seq(400, 900, by = 100)
#' )
#'
#' cropped = raman_crop_spectra(
#'   dataset,
#'   min_wn = 500,
#'   max_wn = 800
#' )
#'
#' cropped$wavenumbers
#'
#' @export
raman_crop_spectra = function(
    dataset,
    min_wn,
    max_wn
) {
  if (!is.list(dataset) || is.null(dataset$data)) {
    stop(
      "'dataset' must be a list containing a 'data' element.",
      call. = FALSE
    )
  }
  
  if (!is.matrix(dataset$data) && !is.data.frame(dataset$data)) {
    stop(
      "'dataset$data' must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }
  
  data = as.matrix(dataset$data)
  
  if (!is.numeric(data)) {
    stop(
      "'dataset$data' must contain numeric values only.",
      call. = FALSE
    )
  }
  
  if (nrow(data) < 1L || ncol(data) < 1L) {
    stop(
      "'dataset$data' must contain at least one wavenumber and one sample.",
      call. = FALSE
    )
  }
  
  if (anyNA(data) || any(!is.finite(data))) {
    stop(
      "'dataset$data' must not contain missing or non-finite values.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(min_wn) ||
      length(min_wn) != 1L ||
      is.na(min_wn) ||
      !is.finite(min_wn)) {
    stop(
      "'min_wn' must be a single finite numeric value.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(max_wn) ||
      length(max_wn) != 1L ||
      is.na(max_wn) ||
      !is.finite(max_wn)) {
    stop(
      "'max_wn' must be a single finite numeric value.",
      call. = FALSE
    )
  }
  
  if (min_wn > max_wn) {
    stop(
      "'min_wn' must be less than or equal to 'max_wn'.",
      call. = FALSE
    )
  }
  
  if (!is.null(dataset$wavenumbers)) {
    wavenumbers = dataset$wavenumbers
    
    if (!is.numeric(wavenumbers) ||
        length(wavenumbers) != nrow(data) ||
        anyNA(wavenumbers) ||
        any(!is.finite(wavenumbers))) {
      stop(
        "'dataset$wavenumbers' must be a finite numeric vector with ",
        "one value per row in 'dataset$data'.",
        call. = FALSE
      )
    }
  } else {
    if (is.null(rownames(data))) {
      stop(
        "Provide 'dataset$wavenumbers' or numeric row names in ",
        "'dataset$data'.",
        call. = FALSE
      )
    }
    
    wavenumbers = suppressWarnings(
      as.numeric(rownames(data))
    )
    
    if (anyNA(wavenumbers) || any(!is.finite(wavenumbers))) {
      stop(
        "Row names in 'dataset$data' must be numeric wavenumbers.",
        call. = FALSE
      )
    }
  }
  
  indices_to_keep = which(
    wavenumbers >= min_wn &
      wavenumbers <= max_wn
  )
  
  if (length(indices_to_keep) == 0L) {
    limits = range(wavenumbers)
    
    stop(
      sprintf(
        paste0(
          "The selected interval is outside the available ",
          "wavenumber range [%s, %s]."
        ),
        format(limits[1L], trim = TRUE),
        format(limits[2L], trim = TRUE)
      ),
      call. = FALSE
    )
  }
  
  dataset$data = data[indices_to_keep, , drop = FALSE]
  dataset$wavenumbers = wavenumbers[indices_to_keep]
  
  dataset
}


#' Remove spikes from Raman or SERS spectra
#'
#' Detects and replaces intensity spikes in Raman or surface-enhanced Raman
#' spectroscopy (SERS) spectra using modified Z-scores based on the median
#' absolute deviation (MAD).
#'
#' The dataset must contain wavenumbers in rows and samples in columns. Each
#' sample spectrum is processed independently. Detected spikes are replaced by
#' the mean intensity of neighbouring non-spike points.
#'
#' @param dataset A dataset list containing a numeric matrix in
#'   \code{dataset$data}, with wavenumbers in rows and samples in columns.
#' @param ma Non-negative integer specifying the number of neighbouring points
#'   considered on each side of a detected spike.
#' @param threshold Positive finite numeric value specifying the modified
#'   Z-score threshold used to identify spikes.
#'
#' @return A dataset equivalent to \code{dataset}, with despiked spectra stored
#'   in \code{dataset$data}. Information about the operation is stored in
#'   \code{dataset$despike}.
#'
#' @examples
#' dataset = list(
#'   data = matrix(
#'     c(
#'       1, 1.1, 1.2, 20, 1.1, 1.2,
#'       2, 2.1, 2.2, 30, 2.1, 2.2
#'     ),
#'     nrow = 6,
#'     ncol = 2,
#'     dimnames = list(
#'       paste0("wavenumber", 1:6),
#'       paste0("sample", 1:2)
#'     )
#'   )
#' )
#'
#' despiked = raman_despike(
#'   dataset,
#'   ma = 1,
#'   threshold = 3
#' )
#'
#' despiked$despike$n_spikes
#'
#' @export
raman_despike = function(
    dataset,
    ma = 10,
    threshold = 7
) {
  if (!is.list(dataset) || is.null(dataset$data)) {
    stop(
      "'dataset' must be a list containing a 'data' element.",
      call. = FALSE
    )
  }
  
  if (!is.matrix(dataset$data) && !is.data.frame(dataset$data)) {
    stop(
      "'dataset$data' must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }
  
  data = as.matrix(dataset$data)
  
  if (!is.numeric(data)) {
    stop(
      "'dataset$data' must contain numeric values only.",
      call. = FALSE
    )
  }
  
  if (nrow(data) < 1L || ncol(data) < 1L) {
    stop(
      "'dataset$data' must contain at least one wavenumber and one sample.",
      call. = FALSE
    )
  }
  
  if (anyNA(data) || any(!is.finite(data))) {
    stop(
      "'dataset$data' must not contain missing or non-finite values.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(ma) ||
      length(ma) != 1L ||
      is.na(ma) ||
      !is.finite(ma) ||
      ma < 0 ||
      ma != as.integer(ma)) {
    stop(
      "'ma' must be a single non-negative integer.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(threshold) ||
      length(threshold) != 1L ||
      is.na(threshold) ||
      !is.finite(threshold) ||
      threshold <= 0) {
    stop(
      "'threshold' must be a single positive finite numeric value.",
      call. = FALSE
    )
  }
  
  ma = as.integer(ma)
  
  clean_single_spectrum = function(intensity) {
    n_wavenumbers = length(intensity)
    
    median_intensity = stats::median(intensity)
    
    median_absolute_deviation = stats::median(
      abs(intensity - median_intensity)
    )
    
    if (median_absolute_deviation == 0) {
      return(
        list(
          spectrum = intensity,
          n_spikes = 0L
        )
      )
    }
    
    modified_z_score = 0.6745 *
      (intensity - median_intensity) /
      median_absolute_deviation
    
    is_spike = abs(modified_z_score) > threshold
    spike_indices = which(is_spike)
    
    corrected_intensity = intensity
    
    for (spike_index in spike_indices) {
      neighbourhood = seq.int(
        from = max(1L, spike_index - ma),
        to = min(n_wavenumbers, spike_index + ma)
      )
      
      valid_neighbours = neighbourhood[
        !is_spike[neighbourhood]
      ]
      
      if (length(valid_neighbours) > 0L) {
        corrected_intensity[spike_index] = mean(
          intensity[valid_neighbours]
        )
      }
    }
    
    list(
      spectrum = corrected_intensity,
      n_spikes = length(spike_indices)
    )
  }
  
  corrected_data = matrix(
    0,
    nrow = nrow(data),
    ncol = ncol(data),
    dimnames = dimnames(data)
  )
  
  spikes_per_sample = integer(ncol(data))
  
  for (sample_index in seq_len(ncol(data))) {
    result = clean_single_spectrum(
      data[, sample_index]
    )
    
    corrected_data[, sample_index] = result$spectrum
    spikes_per_sample[sample_index] = result$n_spikes
  }
  
  sample_names = colnames(data)
  
  if (is.null(sample_names)) {
    sample_names = paste0(
      "sample",
      seq_len(ncol(data))
    )
  }
  
  dataset$data = corrected_data
  
  dataset$despike = list(
    ma = ma,
    threshold = threshold,
    n_spikes = sum(spikes_per_sample),
    spikes_per_sample = stats::setNames(
      spikes_per_sample,
      sample_names
    )
  )
  
  dataset
}


#' Baseline correction of Raman and SERS spectra using airPLS
#'
#' Applies adaptive iteratively reweighted penalized least squares (airPLS)
#' baseline correction to Raman or surface-enhanced Raman spectroscopy (SERS)
#' spectra.
#'
#' The dataset is expected to contain spectral variables or wavenumbers in rows
#' and samples in columns. Each sample is baseline-corrected independently.
#' Corrected spectra are shifted so that their minimum intensity is zero.
#'
#' @param dataset A dataset list containing a numeric matrix in
#'   \code{dataset$data}, with wavenumbers in rows and samples in columns.
#' @param lambda Positive numeric smoothing parameter. Larger values produce a
#'   smoother estimated baseline.
#' @param max_iter Positive integer giving the maximum number of airPLS
#'   iterations.
#'
#' @return A dataset equivalent to \code{dataset}, with baseline-corrected
#'   values stored in \code{dataset$data}. Details of the correction are stored
#'   in \code{dataset$baseline}.
#'
#' @references
#' Zhang, Z.-M., Chen, S., and Liang, Y.-Z. (2010). Baseline correction using
#' adaptive iteratively reweighted penalized least squares. \emph{Analyst},
#' 135(5), 1138-1146. \doi{10.1039/B922045C}.
#'
#' @examples
#' set.seed(123)
#'
#' dataset = list(
#'   data = matrix(
#'     rnorm(100, mean = 10, sd = 2),
#'     nrow = 10,
#'     dimnames = list(
#'       paste0("wavenumber", 1:10),
#'       paste0("sample", 1:10)
#'     )
#'   )
#' )
#'
#' corrected = airPLS_fast_dataset(
#'   dataset,
#'   lambda = 1e4,
#'   max_iter = 5
#' )
#'
#' range(corrected$data)
#'
#' @export
airPLS_fast_dataset = function(
    dataset,
    lambda = 1e5,
    max_iter = 50
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "Package 'Matrix' is required for this function. ",
      "Install it with install.packages('Matrix').",
      call. = FALSE
    )
  }
  
  if (!is.list(dataset) || is.null(dataset$data)) {
    stop(
      "'dataset' must be a list containing a 'data' element.",
      call. = FALSE
    )
  }
  
  if (!is.matrix(dataset$data) && !is.data.frame(dataset$data)) {
    stop(
      "'dataset$data' must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }
  
  data = as.matrix(dataset$data)
  
  if (!is.numeric(data)) {
    stop(
      "'dataset$data' must contain numeric values only.",
      call. = FALSE
    )
  }
  
  if (nrow(data) < 3L || ncol(data) < 1L) {
    stop(
      "'dataset$data' must contain at least three wavenumbers ",
      "and one sample.",
      call. = FALSE
    )
  }
  
  if (anyNA(data) || any(!is.finite(data))) {
    stop(
      "'dataset$data' must not contain missing or non-finite values.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(lambda) ||
      length(lambda) != 1L ||
      is.na(lambda) ||
      !is.finite(lambda) ||
      lambda <= 0) {
    stop(
      "'lambda' must be a single positive finite numeric value.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(max_iter) ||
      length(max_iter) != 1L ||
      is.na(max_iter) ||
      !is.finite(max_iter) ||
      max_iter < 1 ||
      max_iter != as.integer(max_iter)) {
    stop(
      "'max_iter' must be a single positive integer.",
      call. = FALSE
    )
  }
  
  max_iter = as.integer(max_iter)
  
  n_wavenumbers = nrow(data)
  n_samples = ncol(data)
  
  difference_matrix = Matrix::bandSparse(
    n = n_wavenumbers - 2L,
    m = n_wavenumbers,
    k = 0:2,
    diagonals = list(
      rep(1, n_wavenumbers - 2L),
      rep(-2, n_wavenumbers - 2L),
      rep(1, n_wavenumbers - 2L)
    )
  )
  
  penalty_matrix = lambda *
    Matrix::crossprod(difference_matrix)
  
  corrected_data = matrix(
    0,
    nrow = n_wavenumbers,
    ncol = n_samples,
    dimnames = dimnames(data)
  )
  
  tolerance = 0.001
  
  for (sample_index in seq_len(n_samples)) {
    intensity = data[, sample_index]
    weights = rep(1, n_wavenumbers)
    baseline = rep(0, n_wavenumbers)
    
    for (iteration in seq_len(max_iter)) {
      system_matrix = Matrix::Diagonal(x = weights) + penalty_matrix
      
      baseline = as.numeric(
        Matrix::solve(
          system_matrix,
          weights * intensity
        )
      )
      
      residuals = intensity - baseline
      negative_residuals = residuals[residuals < 0]
      
      if (length(negative_residuals) == 0L) {
        break
      }
      
      negative_sum = sum(abs(negative_residuals))
      intensity_sum = sum(abs(intensity))
      
      if (intensity_sum > 0 &&
          iteration > 1L &&
          negative_sum / intensity_sum < tolerance) {
        break
      }
      
      weights = numeric(n_wavenumbers)
      
      weights[residuals < 0] = exp(
        iteration * residuals[residuals < 0] / negative_sum
      )
    }
    
    corrected_spectrum = intensity - baseline
    
    corrected_data[, sample_index] = corrected_spectrum -
      min(corrected_spectrum)
  }
  
  dataset$data = corrected_data
  
  dataset$baseline = list(
    method = "airPLS",
    lambda = lambda,
    max_iter = max_iter,
    tolerance = tolerance
  )
  
  dataset
}


#' Normalize Raman or SERS spectra
#'
#' Normalizes Raman or surface-enhanced Raman spectroscopy (SERS) spectra
#' using min-max normalization, L2 vector normalization, or normalization by
#' the maximum intensity near a target silicon peak.
#'
#' The dataset must contain wavenumbers in rows and samples in columns. Each
#' sample spectrum is normalized independently.
#'
#' @param dataset A dataset list containing a numeric matrix in
#'   \code{dataset$data}, with wavenumbers in rows and samples in columns.
#'   For \code{method = "silicon"}, the dataset must also contain a numeric
#'   vector in \code{dataset$wavenumbers}, with one value per row of
#'   \code{dataset$data}.
#' @param method Character string specifying the normalization method. One of
#'   \code{"min-max"}, \code{"L2"}, or \code{"silicon"}.
#' @param target_peak Single finite numeric value specifying the target
#'   wavenumber of the silicon peak when \code{method = "silicon"}.
#' @param search_window Positive finite numeric value specifying the half-width
#'   of the wavenumber search interval around \code{target_peak} when
#'   \code{method = "silicon"}.
#'
#' @return A dataset equivalent to \code{dataset}, with normalized spectra
#'   stored in \code{dataset$data}. Information about the selected
#'   normalization method and its parameters is stored in
#'   \code{dataset$normalization}.
#'
#' @examples
#' dataset = list(
#'   data = matrix(
#'     c(
#'       1, 2, 3, 4, 5,
#'       2, 4, 6, 8, 10
#'     ),
#'     nrow = 5,
#'     ncol = 2,
#'     dimnames = list(
#'       as.character(seq(500, 900, by = 100)),
#'       paste0("sample", 1:2)
#'     )
#'   ),
#'   wavenumbers = seq(500, 900, by = 100)
#' )
#'
#' normalized = raman_normalize(
#'   dataset,
#'   method = "min-max"
#' )
#'
#' range(normalized$data[, 1])
#'
#' @export
raman_normalize = function(
    dataset,
    method = "min-max",
    target_peak = 520,
    search_window = 20
) {
  if (!is.list(dataset) || is.null(dataset$data)) {
    stop(
      "'dataset' must be a list containing a 'data' element.",
      call. = FALSE
    )
  }
  
  if (!is.matrix(dataset$data) && !is.data.frame(dataset$data)) {
    stop(
      "'dataset$data' must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }
  
  data = as.matrix(dataset$data)
  
  if (!is.numeric(data)) {
    stop(
      "'dataset$data' must contain numeric values only.",
      call. = FALSE
    )
  }
  
  if (nrow(data) < 1L || ncol(data) < 1L) {
    stop(
      "'dataset$data' must contain at least one wavenumber and one sample.",
      call. = FALSE
    )
  }
  
  if (anyNA(data) || any(!is.finite(data))) {
    stop(
      "'dataset$data' must not contain missing or non-finite values.",
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
    "L2",
    "silicon"
  )
  
  if (!method %in% valid_methods) {
    stop(
      "'method' must be one of: 'min-max', 'L2', or 'silicon'.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(target_peak) ||
      length(target_peak) != 1L ||
      is.na(target_peak) ||
      !is.finite(target_peak)) {
    stop(
      "'target_peak' must be a single finite numeric value.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(search_window) ||
      length(search_window) != 1L ||
      is.na(search_window) ||
      !is.finite(search_window) ||
      search_window <= 0) {
    stop(
      "'search_window' must be a single positive finite numeric value.",
      call. = FALSE
    )
  }
  
  normalized_data = matrix(
    0,
    nrow = nrow(data),
    ncol = ncol(data),
    dimnames = dimnames(data)
  )
  
  if (method == "min-max") {
    for (sample_index in seq_len(ncol(data))) {
      spectrum = data[, sample_index]
      spectrum_minimum = min(spectrum)
      spectrum_range = max(spectrum) - spectrum_minimum
      
      if (spectrum_range == 0) {
        normalized_data[, sample_index] = spectrum
      } else {
        normalized_data[, sample_index] = (
          spectrum - spectrum_minimum
        ) / spectrum_range
      }
    }
  }
  
  if (method == "L2") {
    for (sample_index in seq_len(ncol(data))) {
      spectrum = data[, sample_index]
      l2_norm = sqrt(sum(spectrum^2))
      
      if (l2_norm == 0) {
        normalized_data[, sample_index] = spectrum
      } else {
        normalized_data[, sample_index] = spectrum / l2_norm
      }
    }
  }
  
  if (method == "silicon") {
    if (is.null(dataset$wavenumbers) ||
        !is.numeric(dataset$wavenumbers) ||
        length(dataset$wavenumbers) != nrow(data) ||
        anyNA(dataset$wavenumbers) ||
        any(!is.finite(dataset$wavenumbers))) {
      stop(
        "'dataset$wavenumbers' must be a finite numeric vector with ",
        "one value per row in 'dataset$data' when method is 'silicon'.",
        call. = FALSE
      )
    }
    
    lower_bound = target_peak - search_window
    upper_bound = target_peak + search_window
    
    peak_indices = which(
      dataset$wavenumbers >= lower_bound &
        dataset$wavenumbers <= upper_bound
    )
    
    if (length(peak_indices) == 0L) {
      stop(
        "No wavenumbers were found within the selected silicon peak window.",
        call. = FALSE
      )
    }
    
    for (sample_index in seq_len(ncol(data))) {
      spectrum = data[, sample_index]
      
      silicon_peak_intensity = max(
        spectrum[peak_indices]
      )
      
      if (silicon_peak_intensity == 0) {
        normalized_data[, sample_index] = spectrum
      } else {
        normalized_data[, sample_index] = spectrum /
          silicon_peak_intensity
      }
    }
  }
  
  dataset$data = normalized_data
  
  dataset$normalization = list(
    method = method,
    target_peak = target_peak,
    search_window = search_window
  )
  
  dataset
}



#' Compute Savitzky-Golay second derivative of Raman or SERS spectra
#'
#' Applies a Savitzky-Golay smoothing filter and computes the second derivative
#' of Raman or surface-enhanced Raman spectroscopy (SERS) spectra.
#'
#' This function is a wrapper around \code{savitzky_golay()} specialized for
#' Raman/SERS datasets. It fixes the differentiation order to 2 and provides
#' argument names and validation tailored to spectroscopic workflows.
#'
#' The dataset must contain wavenumbers in rows and samples in columns.
#' Each sample spectrum is processed independently. The second derivative is
#' computed using polynomial fitting within a moving window.
#'
#' @param dataset A dataset list containing a numeric matrix in
#'   \code{dataset$data}, with wavenumbers in rows and samples in columns.
#'   The dataset must also contain a numeric vector in
#'   \code{dataset$wavenumbers}, with one value per row of
#'   \code{dataset$data}.
#' @param window_size Odd positive integer specifying the size of the moving
#'   window (number of points). Must be less than or equal to the number of
#'   wavenumbers.
#' @param polynomial_order Positive integer specifying the order of the
#'   polynomial used for local fitting. Must be less than \code{window_size}.
#'
#' @return A dataset equivalent to \code{dataset}, with second derivative
#'   spectra stored in \code{dataset$data}. The original wavenumbers are stored
#'   in \code{dataset$original_wavenumbers}, when available. Information about
#'   the transform is stored in \code{dataset$derivative}.
#'
#' @examples
#' wavenumbers = seq(400, 1000, length.out = 21)
#'
#' dataset = list(
#'   data = cbind(
#'     sample1 = sin(seq(0, 2 * pi, length.out = 21)),
#'     sample2 = cos(seq(0, 2 * pi, length.out = 21))
#'   ),
#'   wavenumbers = wavenumbers
#' )
#'
#' deriv2 = raman_sgolay_derivative(
#'   dataset,
#'   window_size = 7,
#'   polynomial_order = 2
#' )
#'
#' dim(deriv2$data)
#'
#' @export
raman_sgolay_derivative = function(
    dataset,
    window_size = 7,
    polynomial_order = 2
) {
  if (!is.list(dataset) || is.null(dataset$data)) {
    stop(
      "'dataset' must be a list containing a 'data' element.",
      call. = FALSE
    )
  }
  
  if (!is.matrix(dataset$data) && !is.data.frame(dataset$data)) {
    stop(
      "'dataset$data' must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }
  
  data = as.matrix(dataset$data)
  
  if (!is.numeric(data)) {
    stop(
      "'dataset$data' must contain numeric values only.",
      call. = FALSE
    )
  }
  
  if (nrow(data) < 3L || ncol(data) < 1L) {
    stop(
      "'dataset$data' must contain at least three wavenumbers and one sample.",
      call. = FALSE
    )
  }
  
  if (anyNA(data) || any(!is.finite(data))) {
    stop(
      "'dataset$data' must not contain missing or non-finite values.",
      call. = FALSE
    )
  }
  
  if (is.null(dataset$wavenumbers) ||
      !is.numeric(dataset$wavenumbers) ||
      length(dataset$wavenumbers) != nrow(data) ||
      anyNA(dataset$wavenumbers) ||
      any(!is.finite(dataset$wavenumbers))) {
    stop(
      "'dataset$wavenumbers' must be a finite numeric vector with ",
      "one value per row in 'dataset$data'.",
      call. = FALSE
    )
  }
  
  wavenumbers = dataset$wavenumbers
  
  if (!is.numeric(window_size) ||
      length(window_size) != 1L ||
      is.na(window_size) ||
      !is.finite(window_size) ||
      window_size < 3 ||
      window_size != as.integer(window_size) ||
      window_size %% 2 == 0) {
    stop(
      "'window_size' must be a single odd integer >= 3.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(polynomial_order) ||
      length(polynomial_order) != 1L ||
      is.na(polynomial_order) ||
      !is.finite(polynomial_order) ||
      polynomial_order < 1 ||
      polynomial_order != as.integer(polynomial_order)) {
    stop(
      "'polynomial_order' must be a single positive integer.",
      call. = FALSE
    )
  }
  
  window_size = as.integer(window_size)
  polynomial_order = as.integer(polynomial_order)
  
  if (window_size > nrow(data)) {
    stop(
      "'window_size' must not exceed the number of wavenumbers.",
      call. = FALSE
    )
  }
  
  if (polynomial_order >= window_size) {
    stop(
      "'polynomial_order' must be less than 'window_size'.",
      call. = FALSE
    )
  }
  
  n_wavenumbers = nrow(data)
  
  result = savitzky_golay(
    dataset = dataset,
    p.order = polynomial_order,
    window = window_size,
    deriv = 2
  )
  
  
  result$derivative = list(
    method = "savitzky_golay_second_derivative",
    window_size = window_size,
    polynomial_order = polynomial_order,
    n_wavenumbers = n_wavenumbers
  )
  
  result
}

