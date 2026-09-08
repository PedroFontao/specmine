# specmine 4.0.0

## New features

- Added Raman/SERS spectral preprocessing functions:
  `raman_crop_spectra()`, `raman_despike()`,
  `airPLS_fast_dataset()`, `raman_normalize()` and
  `raman_sgolay_derivative()`.

- Added Raman/SERS peak-analysis functions:
  `raman_find_peaks()`, `raman_align_peaks()` and
  `raman_normalize_peak_features()`.

- Added optional Fourier and Haar wavelet feature extraction through
  `raman_transform_fourier()` and `raman_transform_wavelet()`.

- Added interactive PCA 3D score plots with `pca_scoresplot3D()` and
  improved support for non-consecutive principal-component selections
  in 3D PCA visualizations.

## Improvements

- Improved input validation, documentation and reproducible examples for
  Raman/SERS and PCA workflows.

- Improved compatibility of PCA 3D group colouring for metadata factors with
  unused levels.
# specmine 3.1.8

- Fixed malformed roxygen example blocks.
- Added and improved `\value{}` documentation in `.Rd` files.
- Removed commented example code and improved example robustness for checks.
- Improved documentation consistency.
- Updated package metadata.

# specmine 3.1.7

- Fixed CRAN resubmission issues related to optional dependencies.
- Updated repository links and package metadata.
- Improved package documentation.