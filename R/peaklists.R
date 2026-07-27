## Functions to read, create and handle peak lists

# Structure: list with one field per sample
# each sample is a data frame representing a peak list (with two columns)

# returns a dataset in standard format from all peaks in all samples
#' Create a dataset from peak lists
#'
#' Merges equivalent peaks across samples, extracts intensities, and creates a
#' dataset object in standard package format.
#'
#' @param sample.list A named list of peak tables, one per sample.
#' @param metadata Optional sample metadata.
#' @param description Character string with a dataset description.
#' @param type Character string describing the dataset type.
#'
#' @return A dataset object created from the input peak lists.
#'
#' @examples
#' sample.list <- list(
#'   s1 = data.frame(ppm = c(1.0, 2.0), int = c(10, 20)),
#'   s2 = data.frame(ppm = c(1.0, 3.0), int = c(15, 30))
#' )
#' metadata <- data.frame(class = c("A", "B"), row.names = c("s1", "s2"))
#' dataset <- dataset_from_peaks(sample.list, metadata = metadata)
#' class(dataset)
#'
#' @export
"dataset_from_peaks" = function(sample.list, metadata = NULL, description = "", type = "nmr-peaks") {
  
  merged.peaks = merge_eq_peaks_samplelist(sample.list)
  samples.df = get_all_intensities(merged.peaks)
  create_dataset(as.matrix(samples.df), metadata = metadata, type = type,
                 description = description)
  
}

# READING

# reads csv files, each with a sample;
# filenames - list of file names of the files to read
# returns list of data frames
#' Read multiple CSV peak files
#'
#' Reads multiple CSV files, one per sample, and returns them as a named list.
#'
#' @param filenames Character vector with file names or file paths.
#' @param ext File extension appended to each file name.
#' @param verbose Logical; whether progress messages should be shown.
#' @param ... Additional arguments passed to \code{read.csv()}.
#'
#' @return A named list of data frames, one per input file.
#'
#' @examples
#' f1 <- tempfile(fileext = ".csv")
#' f2 <- tempfile(fileext = ".csv")
#' utils::write.csv(data.frame(ppm = c(1, 2), int = c(10, 20)), f1, row.names = FALSE)
#' utils::write.csv(data.frame(ppm = c(1, 3), int = c(15, 30)), f2, row.names = FALSE)
#' read_multiple_csvs(c(f1, f2), ext = "", verbose = FALSE)
#'
#' @export
"read_multiple_csvs" = function(filenames, ext = ".csv", verbose = TRUE, ...)
{
  sampleList = list()
  sampleNames = c()
  snames <- gsub("\\.[^.]*$", "", basename(filenames))
  for (i in 1:length(filenames)) {
    if (verbose) message("Reading sample ", filenames[i])
    sampleList[[i]] = read.csv(paste(filenames[i], ext, sep=""), ...)
  }
  sampleNames = snames
  names(sampleList) = sampleNames
  sampleList
}

# reads list of CSV files from a given folder
# returns list of data frames
#' Read all CSV peak files in a folder
#'
#' Lists CSV files in a folder and reads them into a named list of data frames.
#'
#' @param foldername Path to the folder containing CSV files.
#' @param verbose Logical; whether progress messages should be shown.
#' @param ... Additional arguments passed to \code{read.csv()}.
#'
#' @return A named list of data frames, one per CSV file found in the folder.
#'
#' @examples
#' folder <- file.path(tempdir(), "peak_csvs")
#' dir.create(folder, showWarnings = FALSE)
#' utils::write.csv(data.frame(ppm = c(1, 2), int = c(10, 20)),
#'                  file.path(folder, "sample1.csv"), row.names = FALSE)
#' utils::write.csv(data.frame(ppm = c(1, 3), int = c(15, 30)),
#'                  file.path(folder, "sample2.csv"), row.names = FALSE)
#' read_csvs_folder(folder, verbose = FALSE)
#'
#' @export
"read_csvs_folder" = function(foldername, verbose = TRUE, ...)
{
  files <- dir(foldername, pattern=".[Cc][Ss][Vv]$", recursive=TRUE, full.names=TRUE)
  sampleList = read_multiple_csvs(files, ext = "", verbose = verbose, ...)
  sampleList
}

# GETTING INFO

# counts number of peaks in a sample (given its index)
"peaks_per_sample" = function(sample.list, sample.index)
{
  nrow(sample.list[[sample.index]])
}

# counts the number of peaks in each sample in the peaklist
"peaks_per_samples" = function(sample.list)
{
  res = c()
  for(i in 1:length(sample.list)) res[i] = peaks_per_sample(sample.list, i)
  res
}

# finds samples that have the same peak values- x and y (equal data frames)
"find_equal_samples" = function(sample.list)
{
  eq1 = c()
  eq2 = c()
  for (i in 1:(length(sample.list)-1))
  {
    for(j in (i+1):length(sample.list))
    {
      res = compare::compare(sample.list[[i]], sample.list[[j]])
      if(res$result == TRUE)
      {
        eq1 = c(eq1, names(sample.list)[i])
        eq2 = c(eq2, names(sample.list)[j])
      }
    }
  }
  data.frame(eq1, eq2)
}

# get the full list of frequencies from all samples in a list
"get_overall_freq_list" = function(sample.list)
{
  res = sample.list[[1]][[1]]
  for(i in 2:length(sample.list))
    res = union(res, sample.list[[i]][[1]])
  sort(res)
}

# get the value of an internsity given the sample (data frame) and the frequency
"get_intensity" = function(sample.df, freq, tolerance = 0.001)
{
  cond = (sample.df[[1]] > freq - tolerance) & (sample.df[[1]] < freq + tolerance)
  if (any(cond)) res = sample.df[cond,][[2]]
  else res = NA
  res
}

# gets all intensities for a sample list; result is a data frame
"get_all_intensities" = function(sample.list, tol = 0.001)
{
  all.freqs = get_overall_freq_list(sample.list)
  
  for(i in 1:length(sample.list))
  {
    intens.vals = c()
    for(k in 1:length(all.freqs))
    {
      intens.vals[k] = get_intensity(sample.list[[i]], all.freqs[k], tol)
    }
    if (i==1) res.df = data.frame(intens.vals)
    else res.df = cbind(res.df, intens.vals)
  }
  names(res.df) = names(sample.list)
  rownames(res.df) = all.freqs
  res.df
}

"remove_peaks_interval" = function(sample.df, peak.val.min, peak.val.max) {
  sample.df[sample.df$ppm < peak.val.min | sample.df$ppm > peak.val.max,]
}

"remove_peaks_interval_sample_list" = function(sample.list, peak.val.min, peak.val.max) {
  sample.list.res = list()
  for(i in 1:length(sample.list)) {
    sample.list.res[[i]] = remove_peaks_interval(sample.list[[i]], peak.val.min, peak.val.max)
  }
  names(sample.list.res) = names(sample.list)
  sample.list.res
}

# functions working over the list of all intensities
#'
#' Auto-exported function: get_peak_values
#'
#' @param samples.df A data frame of peak intensities with peaks in rows and samples in columns.
#' @param peak.val Peak value identifying the row to extract.
#'
#' @return A numeric vector with the intensity values for the peak specified by \code{peak.val} across all samples in \code{samples.df}.
#'
#' @examples
#' samples.df <- data.frame(
#'   s1 = c(10, 20),
#'   s2 = c(15, NA),
#'   row.names = c("1", "2")
#' )
#' get_peak_values(samples.df, "1")
#'
#' @keywords internal
#' @export
"get_peak_values" = function(samples.df, peak.val)
{
  index = which(rownames(samples.df) == peak.val)
  values = c()
  for (s in (samples.df[index,])) {
    values = c(values, s)
  }
  values
}

"values_per_peak" = function(samples.df)
{
  res = c()
  for(i in 1:nrow(samples.df)) res[i] = sum(!is.na(samples.df[i,]))
  res
}

"values_per_sample" = function(samples.df)
{
  res = c()
  for(i in 1:ncol(samples.df)) res[i] = sum(!is.na(samples.df[,i]))
  res
}

# PROCESSING

# merge peaks with equal frequencies (ppm) in a sample given by a data frame
# intensity values are summed
"merge_equal_peaks" = function(sample.df, tolerance = 0.0)
{
  d = diff(sample.df[[1]])
  indexes = which(d <= tolerance)
  if (length(indexes) != 0){
    new.sample.df = sample.df[-(indexes+1),]
    new.sample.df[[2]] = sum_vec(sample.df[[2]], indexes)
  }
  else {
    new.sample.df = sample.df
  }
  
  new.sample.df
}

"sum_vec" = function(orig.vec, indexes)
{
  newvec = c()
  for (i in 1:length(orig.vec)) newvec[i] = orig.vec[i]
  
  for(k in length(indexes):1)
  {
    newvec[indexes[k]] = newvec[indexes[k]] + newvec[indexes[k]+1]
    newvec = newvec[-(indexes[k]+1)]
  }
  newvec
}

# merge peaks with equal freqs (ppm) in all samples of a list
"merge_eq_peaks_samplelist" = function(sample.list, tolerance = 0.0)
{
  newlist = list()
  for(i in 1:length(sample.list))
  {
    newlist[[i]] = merge_equal_peaks(sample.list[[i]], tolerance)
  }
  names(newlist) = names(sample.list)
  newlist
}

# create a matrix with all values for all samples (as used by MetaboAnalyst)
"create_full_matrix" = function(sample.list)
{
  resmat = NULL
  for(i in 1:length(sample.list))
  {
    for(j in 1:nrow(sample.list[[i]]))
      resmat = rbind(resmat, cbind(sample.list[[i]][j,], i))
  }
  colnames(resmat) = c("ppm","int","sample")
  resmat
}