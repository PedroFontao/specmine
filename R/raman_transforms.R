#' Extract Fourier power features from Raman or SERS spectra
#'
#' Computes one-sided Fourier power spectra for each Raman or
#' surface-enhanced Raman spectroscopy (SERS) sample.
#'
#' This function follows the \code{specmine} Raman workflow by preserving the
#' original spectrum in \code{dataset$data}. Fourier features are stored in
#' \code{dataset$fourier$data}; their normalized frequencies are stored in
#' \code{dataset$fourier$frequencies}.
#'
#' The input data matrix must have wavenumbers in rows and samples in columns.
#' Each sample is transformed independently.
#'
#' @param dataset A dataset list containing a numeric matrix in
#'   \code{dataset$data}, with wavenumbers in rows and samples in columns.
#' @param remove_dc Logical value indicating whether the zero-frequency
#'   direct-current component should be excluded from the stored output.
#'
#' @return The input \code{dataset}, with a \code{fourier} element containing:
#'   \describe{
#'     \item{\code{data}}{Matrix of one-sided Fourier power features, with
#'       frequencies in rows and samples in columns.}
#'     \item{\code{frequencies}}{Numeric vector of normalized frequencies.}
#'     \item{\code{remove_dc}}{Logical value indicating whether the DC
#'       component was removed.}
#'     \item{\code{n_wavenumbers}}{Number of original spectral variables.}
#'   }
#'
#' @examples
#' wavenumbers = seq(400, 1000, length.out = 16)
#'
#' dataset = list(
#'   data = cbind(
#'     sample1 = sin(seq(0, 2 * pi, length.out = 16)),
#'     sample2 = cos(seq(0, 2 * pi, length.out = 16))
#'   ),
#'   wavenumbers = wavenumbers
#' )
#'
#' transformed = raman_transform_fourier(dataset)
#'
#' dim(transformed$fourier$data)
#'
#' @export
raman_transform_fourier = function(
    dataset,
    remove_dc = FALSE
) {
  .validate_raman_dataset(dataset)
  
  if (!is.logical(remove_dc) ||
      length(remove_dc) != 1L ||
      is.na(remove_dc)) {
    stop(
      "'remove_dc' must be a single non-missing logical value.",
      call. = FALSE
    )
  }
  
  data = as.matrix(dataset$data)
  n_wavenumbers = nrow(data)
  n_samples = ncol(data)
  
  n_frequencies = floor(n_wavenumbers / 2) + 1L
  frequency_indices = seq.int(0L, n_frequencies - 1L)
  frequencies = frequency_indices / n_wavenumbers
  
  power_spectrum = matrix(
    0,
    nrow = n_frequencies,
    ncol = n_samples,
    dimnames = list(
      as.character(frequencies),
      colnames(data)
    )
  )
  
  for (sample_index in seq_len(n_samples)) {
    fourier_coefficients = stats::fft(data[, sample_index])
    
    power_values = Mod(
      fourier_coefficients[seq_len(n_frequencies)]
    )^2 / n_wavenumbers
    
    if (n_wavenumbers %% 2L == 0L) {
      if (n_frequencies > 2L) {
        power_values[2L:(n_frequencies - 1L)] =
          2 * power_values[2L:(n_frequencies - 1L)]
      }
    } else if (n_frequencies > 1L) {
      power_values[2L:n_frequencies] =
        2 * power_values[2L:n_frequencies]
    }
    
    power_spectrum[, sample_index] = power_values
  }
  
  if (remove_dc) {
    power_spectrum = power_spectrum[-1L, , drop = FALSE]
    frequencies = frequencies[-1L]
  }
  
  rownames(power_spectrum) = as.character(frequencies)
  
  dataset$fourier = list(
    data = power_spectrum,
    frequencies = frequencies,
    remove_dc = remove_dc,
    n_wavenumbers = n_wavenumbers
  )
  
  dataset
}


#' Extract Haar wavelet features from Raman or SERS spectra
#'
#' Computes a multilevel discrete Haar wavelet transform for each Raman or
#' surface-enhanced Raman spectroscopy (SERS) sample.
#'
#' This function preserves the original data in \code{dataset$data}, following
#' the \code{specmine} pipeline. Wavelet coefficients are stored in
#' \code{dataset$wavelet$data}.
#'
#' The input data matrix must have wavenumbers in rows and samples in columns.
#' If needed, each spectrum is extended by repeating its final intensity value
#' until its length is a power of two.
#'
#' @param dataset A dataset list containing a numeric matrix in
#'   \code{dataset$data}, with wavenumbers in rows and samples in columns.
#' @param level Positive integer specifying the number of Haar decomposition
#'   levels.
#'
#' @return The input \code{dataset}, with a \code{wavelet} element containing:
#'   \describe{
#'     \item{\code{data}}{Matrix of Haar wavelet coefficients, with
#'       coefficients in rows and samples in columns.}
#'     \item{\code{filter}}{The character string \code{"haar"}.}
#'     \item{\code{level}}{Number of decomposition levels.}
#'     \item{\code{original_length}}{Number of original wavenumbers.}
#'     \item{\code{padded_length}}{Length used after optional padding.}
#'   }
#'
#' @examples
#' wavenumbers = seq(400, 1100, length.out = 8)
#'
#' dataset = list(
#'   data = cbind(
#'     sample1 = sin(seq(0, 2 * pi, length.out = 8)),
#'     sample2 = cos(seq(0, 2 * pi, length.out = 8))
#'   ),
#'   wavenumbers = wavenumbers
#' )
#'
#' transformed = raman_transform_wavelet(
#'   dataset,
#'   level = 2
#' )
#'
#' dim(transformed$wavelet$data)
#'
#' @export
raman_transform_wavelet = function(
    dataset,
    level = 1
) {
  .validate_raman_dataset(dataset)
  
  if (!is.numeric(level) ||
      length(level) != 1L ||
      is.na(level) ||
      !is.finite(level) ||
      level < 1L ||
      level != as.integer(level)) {
    stop(
      "'level' must be a single positive integer.",
      call. = FALSE
    )
  }
  
  level = as.integer(level)
  
  data = as.matrix(dataset$data)
  n_wavenumbers = nrow(data)
  n_samples = ncol(data)
  
  padded_length = 2^ceiling(log2(n_wavenumbers))
  max_level = as.integer(log2(padded_length))
  
  if (level > max_level) {
    stop(
      sprintf(
        "'level' must not exceed %d for spectra with %d wavenumbers.",
        max_level,
        n_wavenumbers
      ),
      call. = FALSE
    )
  }
  
  padded_data = data
  
  if (padded_length > n_wavenumbers) {
    padding = matrix(
      rep(
        data[n_wavenumbers, ],
        each = padded_length - n_wavenumbers
      ),
      nrow = padded_length - n_wavenumbers,
      ncol = n_samples,
      dimnames = list(
        NULL,
        colnames(data)
      )
    )
    
    padded_data = rbind(
      data,
      padding
    )
  }
  
  wavelet_coefficients = matrix(
    0,
    nrow = padded_length,
    ncol = n_samples,
    dimnames = list(
      paste0(
        "wavelet_coefficient_",
        seq_len(padded_length)
      ),
      colnames(data)
    )
  )
  
  for (sample_index in seq_len(n_samples)) {
    wavelet_coefficients[, sample_index] = .haar_transform(
      spectrum = padded_data[, sample_index],
      level = level
    )
  }
  
  dataset$wavelet = list(
    data = wavelet_coefficients,
    filter = "haar",
    level = level,
    original_length = n_wavenumbers,
    padded_length = padded_length,
    padding = "repeat_last_value"
  )
  
  dataset
}


#' Validate an input Raman dataset
#'
#' @param dataset Dataset list containing a \code{data} element.
#' @return Invisibly returns \code{TRUE} if validation succeeds.
#' @noRd
.validate_raman_dataset = function(dataset) {
  if (!is.list(dataset) ||
      is.null(dataset$data)) {
    stop(
      "'dataset' must be a list containing a 'data' element.",
      call. = FALSE
    )
  }
  
  if (!is.matrix(dataset$data) &&
      !is.data.frame(dataset$data)) {
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
  
  if (nrow(data) < 2L ||
      ncol(data) < 1L) {
    stop(
      paste0(
        "'dataset$data' must contain at least two wavenumbers ",
        "and one sample."
      ),
      call. = FALSE
    )
  }
  
  if (anyNA(data) ||
      any(!is.finite(data))) {
    stop(
      "'dataset$data' must not contain missing or non-finite values.",
      call. = FALSE
    )
  }
  
  if (!is.null(dataset$wavenumbers)) {
    if (!is.numeric(dataset$wavenumbers) ||
        length(dataset$wavenumbers) != nrow(data) ||
        anyNA(dataset$wavenumbers) ||
        any(!is.finite(dataset$wavenumbers))) {
      stop(
        paste0(
          "'dataset$wavenumbers' must be a finite numeric vector ",
          "with one value per data row."
        ),
        call. = FALSE
      )
    }
  }
  
  invisible(TRUE)
}


#' Perform a multilevel discrete Haar wavelet transform
#'
#' @param spectrum Numeric vector whose length is a power of two.
#' @param level Number of Haar decomposition levels.
#' @return Numeric vector of Haar wavelet coefficients.
#' @noRd
.haar_transform = function(
    spectrum,
    level
) {
  coefficients = spectrum
  active_length = length(coefficients)
  
  for (current_level in seq_len(level)) {
    active_values = coefficients[seq_len(active_length)]
    
    odd_indices = seq.int(
      from = 1L,
      to = active_length,
      by = 2L
    )
    
    even_indices = seq.int(
      from = 2L,
      to = active_length,
      by = 2L
    )
    
    approximation = (
      active_values[odd_indices] +
        active_values[even_indices]
    ) / sqrt(2)
    
    detail = (
      active_values[odd_indices] -
        active_values[even_indices]
    ) / sqrt(2)
    
    half_length = active_length %/% 2L
    
    coefficients[seq_len(half_length)] = approximation
    
    coefficients[
      seq.int(
        from = half_length + 1L,
        to = active_length
      )
    ] = detail
    
    active_length = half_length
  }
  
  coefficients
}