#' List public MetaboLights studies
#'
#' Returns the identifiers of public studies available in MetaboLights.
#'
#' @return A \code{character} vector with public MetaboLights study identifiers.
#'
#' @examples
#' \dontrun{
#' head(metabolights_studies_list())
#' }
#'
#' @export
metabolights_studies_list = function() {
  if (!requireNamespace("RCurl", quietly = TRUE)) {
    stop(
      "Package 'RCurl' needed for this function to work. Please install it with install.packages('RCurl').",
      call. = FALSE
    )
  }
  
  ftp_base = "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/"
  listStudies_character = RCurl::getURL(ftp_base, dirlistonly = TRUE)
  listStudies_vector = strsplit(listStudies_character, "\n")[[1]]
  listStudies_vector
}

#'
#' Auto-exported function: get_files_list_per_assay
#'
#' @param studyID Character string with the MetaboLights study identifier.
#'
#' @return A named \code{list} of \code{data.frame} objects, one per assay, each containing
#'   the sample names and corresponding file names for the specified study. Each data frame has
#'   two columns: \code{Samples}, with sample identifiers, and \code{Files}, with the associated
#'   raw or spectral data file names available in MetaboLights.
#'
#' @examples
#' \dontrun{
#' assays <- get_files_list_per_assay("MTBLS1")
#' names(assays)
#' }
#'
#' @keywords internal
#' @export
get_files_list_per_assay = function(studyID) {
  if (!requireNamespace("curl", quietly = TRUE)) {
    stop(
      "Package 'curl' needed for this function to work. Please install it with install.packages('curl').",
      call. = FALSE
    )
  }
  
  if (!studyID %in% metabolights_studies_list()) {
    stop("Invalid studyID.", call. = FALSE)
  }
  
  if (!curl::has_internet()) {
    stop("No internet connection.", call. = FALSE)
  }
  
  ftp_base = paste0("ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/", studyID, "/")
  assays_files = list()
  
  i_file = readLines(paste0(ftp_base, "i_Investigation.txt"))
  
  assays_names_line = grep("^Study Assay File Name", i_file, value = TRUE)
  assays_names = grep("^a_", strsplit(assays_names_line, "\"")[[1]], value = TRUE)
  assay_files = paste0(ftp_base, assays_names)
  
  for (i in seq_along(assay_files)) {
    assay = read.table(assay_files[i], header = TRUE)
    
    if ("Free.Induction.Decay.Data.File" %in% colnames(assay)) {
      files_column = "Free.Induction.Decay.Data.File"
    } else {
      files_column = "Raw.Spectral.Data.File"
    }
    
    assay_info = assay[, c("Sample.Name", files_column)]
    colnames(assay_info) = c("Samples", "Files")
    
    assays_files[[paste("Assay", i)]] = assay_info
  }
  
  assays_files
}

#' Get sample-file mapping for one MetaboLights assay
#'
#' Returns the sample-to-file mapping for one assay in a MetaboLights study and can optionally
#' save it as a CSV file.
#'
#' @param studyID Character string with the MetaboLights study identifier.
#' @param assay Numeric or character index selecting the assay in the study.
#' @param directory Optional output directory where \code{samples_files.csv} will be written.
#'
#' @return Invisibly returns a two-column object with samples and corresponding files.
#'
#' @examples
#' \dontrun{
#' get_metabolights_study_samples_files("MTBLS1", assay = 1)
#' }
#'
#' @export
get_metabolights_study_samples_files = function(studyID, assay, directory = NULL) {
  assays_in_study = get_files_list_per_assay(studyID)
  
  samples_files = assays_in_study[[assay]]
  colnames(samples_files) = NULL
  
  if (!is.null(directory)) {
    if (!nzchar(directory)) {
      stop("Please provide a valid 'directory'.", call. = FALSE)
    }
    if (!dir.exists(directory)) {
      dir.create(directory, recursive = TRUE)
    }
    utils::write.csv(
      samples_files,
      file = file.path(directory, "samples_files.csv"),
      row.names = FALSE
    )
  }
  
  invisible(samples_files)
}

#' Download files for one MetaboLights assay
#'
#' Downloads the data files associated with one assay of a MetaboLights study.
#'
#' @param studyID Character string with the MetaboLights study identifier.
#' @param assay Numeric or character index selecting the assay in the study.
#' @param directory Output directory where assay files will be downloaded.
#' @param verbose Logical indicating whether progress messages should be shown.
#'
#' @return Invisibly returns a \code{character} vector with downloaded file paths.
#'
#' @examples
#' \dontrun{
#' tmpdir <- tempdir()
#' get_metabolights_study_files_assay("MTBLS1", assay = 1, directory = tmpdir, verbose = FALSE)
#' }
#'
#' @export
get_metabolights_study_files_assay = function(studyID, assay, directory, verbose = TRUE) {
  if (!requireNamespace("curl", quietly = TRUE)) {
    stop(
      "Package 'curl' needed for this function to work. Please install it with install.packages('curl').",
      call. = FALSE
    )
  }
  
  if (missing(directory) || is.null(directory) || !nzchar(directory)) {
    stop("Please provide a valid 'directory'.", call. = FALSE)
  }
  
  files_per_assay = get_files_list_per_assay(studyID)
  
  files_to_download = as.character(files_per_assay[[assay]][, "Files"])
  ftp_base = paste0("ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/", studyID, "/")
  files_to_download_paths = paste0(ftp_base, files_to_download)
  
  assay_dir = file.path(directory, paste0("assay_", assay))
  data_dir = file.path(assay_dir, "data")
  files_dest = file.path(data_dir, files_to_download)
  
  if (!dir.exists(assay_dir)) {
    dir.create(assay_dir, recursive = TRUE)
  }
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  
  for (i in seq_along(files_to_download)) {
    if (!curl::has_internet()) {
      stop("No internet connection.", call. = FALSE)
    } else {
      if (verbose) message("Downloading file ", i, " of ", length(files_to_download), ": ", files_to_download[i])
      curl::curl_download(files_to_download_paths[i], files_dest[i], quiet = !verbose)
    }
  }
  
  invisible(files_dest)
}

#' Get metadata for one MetaboLights assay
#'
#' Returns the metadata associated with the samples in one assay of a MetaboLights study and can
#' optionally save it as a CSV file.
#'
#' @param studyID Character string with the MetaboLights study identifier.
#' @param assay Numeric or character index selecting the assay in the study.
#' @param directory Optional output directory where the metadata CSV will be written.
#'
#' @return Invisibly returns a \code{data.frame} with sample metadata for the selected assay.
#'
#' @examples
#' \dontrun{
#' md <- get_metabolights_study_metadata_assay("MTBLS1", assay = 1)
#' head(md)
#' }
#'
#' @export
get_metabolights_study_metadata_assay = function(studyID, assay, directory = NULL) {
  if (!requireNamespace("curl", quietly = TRUE)) {
    stop(
      "Package 'curl' needed for this function to work. Please install it with install.packages('curl').",
      call. = FALSE
    )
  }
  
  if (!studyID %in% metabolights_studies_list()) {
    stop("Invalid studyID.", call. = FALSE)
  }
  
  ftp_base = paste0("ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/", studyID, "/")
  
  if (!curl::has_internet()) {
    stop("No internet connection.", call. = FALSE)
  }
  
  i_file = readLines(paste0(ftp_base, "i_Investigation.txt"))
  
  factors_names = strsplit(grep("^Study Factor Name", i_file, value = TRUE), "\"")[[1]]
  factors_names_vec = grep("\t", factors_names, value = TRUE, invert = TRUE)
  
  factors = c()
  for (factor in factors_names_vec) {
    factors = c(
      factors,
      paste0("Factor.Value.", paste(strsplit(factor, " ")[[1]], collapse = "."), ".")
    )
  }
  
  files_per_assay = get_files_list_per_assay(studyID)
  samples_in_assay = as.character(files_per_assay[[assay]][, 1])
  
  sample_file_o = strsplit(grep("^Study File Name", i_file, value = TRUE), "\"")[[1]][2]
  sample_file = paste0(ftp_base, sample_file_o)
  
  sample_metadata = read.table(sample_file, header = TRUE)
  
  metadata = as.data.frame(sample_metadata[, factors])
  rownames(metadata) = as.character(sample_metadata$Sample.Name)
  metadata = as.data.frame(metadata[samples_in_assay, , drop = FALSE])
  colnames(metadata) = factors_names_vec
  rownames(metadata) = as.character(samples_in_assay)
  
  if (!is.null(directory)) {
    if (!nzchar(directory)) {
      stop("Please provide a valid 'directory'.", call. = FALSE)
    }
    if (!dir.exists(directory)) {
      dir.create(directory, recursive = TRUE)
    }
    metadata_filepath = file.path(directory, paste0("metadata", assay, ".csv"))
    utils::write.csv(metadata, metadata_filepath)
  }
  
  invisible(metadata)
}

#' Download a complete MetaboLights study
#'
#' Downloads files and metadata for every assay in a public MetaboLights study.
#'
#' @param studyID Character string with the MetaboLights study identifier.
#' @param directory Output directory where study files and metadata will be written.
#' @param verbose Logical indicating whether progress messages should be shown.
#'
#' @return Invisibly returns a \code{list} with one entry per assay containing downloaded files,
#'   metadata, and sample-file mappings.
#'
#' @examples
#' \dontrun{
#' tmpdir <- tempdir()
#' get_metabolights_study("MTBLS1", directory = tmpdir, verbose = FALSE)
#' }
#'
#' @export
get_metabolights_study = function(studyID, directory, verbose = TRUE) {
  if (!requireNamespace("curl", quietly = TRUE)) {
    stop(
      "Package 'curl' needed for this function to work. Please install it with install.packages('curl').",
      call. = FALSE
    )
  }
  
  if (missing(directory) || is.null(directory) || !nzchar(directory)) {
    stop("Please provide a valid 'directory'.", call. = FALSE)
  }
  
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  if (verbose) message("Getting assays...")
  assays_in_study = get_files_list_per_assay(studyID)
  if (verbose) message("Done.")
  
  results = vector("list", length(assays_in_study))
  
  for (assay in seq_along(assays_in_study)) {
    assay_dir = file.path(directory, paste0("assay_", assay))
    
    if (!dir.exists(assay_dir)) {
      dir.create(assay_dir, recursive = TRUE)
    }
    
    if (verbose) message("Assay ", assay)
    
    if (verbose) message("Getting files from assay ", assay, "...")
    files_downloaded = get_metabolights_study_files_assay(
      studyID = studyID,
      assay = assay,
      directory = directory,
      verbose = verbose
    )
    if (verbose) message("Done.")
    
    if (verbose) message("Getting metadata from assay ", assay, "...")
    metadata = get_metabolights_study_metadata_assay(
      studyID = studyID,
      assay = assay,
      directory = assay_dir
    )
    if (verbose) message("Done.")
    
    if (verbose) message("Getting samples_files file from assay ", assay, "...")
    samples_files = assays_in_study[[assay]]
    colnames(samples_files) = NULL
    utils::write.csv(
      samples_files,
      file = file.path(assay_dir, "samples_files.csv"),
      row.names = FALSE
    )
    if (verbose) message("Done.")
    
    results[[assay]] = list(
      files = files_downloaded,
      metadata = metadata,
      samples_files = samples_files
    )
  }
  
  invisible(results)
}