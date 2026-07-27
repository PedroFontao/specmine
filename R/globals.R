#' @noRd
"_PACKAGE"

# R/globals.R
utils::globalVariables(c("profmethod"))

if (getRversion() >= '2.15.1') {
  utils::globalVariables(c(
    '%>%', '.data', 'defaultSummary', 'pca',
    '.spc.error', '.spc.filehdr', '.spc.log',
    '.spc.read.x', '.spc.read.y', '.spc.subfiledir', '.spc.subhdr'
  ))
}
