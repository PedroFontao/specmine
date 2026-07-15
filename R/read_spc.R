read_dataset_spc <- function(folder.data, filename.meta = NULL,
                             type = "undefined", description = "",
                             nosubhdr = FALSE, label.x = NULL,
                             label.values = NULL,
                             header.col.meta = TRUE,
                             header.row.meta = TRUE,
                             sep.meta = ",") {
  
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
  
  data.spc <- read_data_spc(folder.data, nosubhdr = nosubhdr)
  
  freqs <- data.spc[[1]]$wavelength
  datamat <- matrix(NA, nrow = length(freqs), ncol = length(data.spc))
  
  for (i in seq_along(data.spc)) {
    datamat[, i] <- data.spc[[i]]$spc
  }
  
  rownames(datamat) <- as.character(freqs)
  colnames(datamat) <- names(data.spc)
  
  if (is.null(label.x)) {
    label.x <- data.spc[[1]]$labels$.wavelength
  }
  
  if (is.null(label.values)) {
    label.values <- data.spc[[1]]$labels$spc
  }
  
  dataset <- create_dataset(
    datamat,
    type = type,
    metadata = metadata,
    description = description,
    label.x = label.x,
    label.values = label.values
  )
  
  dataset
}


read_data_spc <- function(foldername, nosubhdr = FALSE) {
  filenames <- dir(foldername, pattern = "\\.[Ss][Pp][Cc]$", full.names = TRUE)
  sampleList <- vector("list", length(filenames))
  sampleNames <- gsub("\\.[^.]*$", "", basename(filenames))
  
  for (i in seq_along(filenames)) {
    print(paste("Reading sample", filenames[i]))
    sampleList[[i]] <- read_spc_nosubhdr(
      filenames[i],
      no.object = TRUE,
      nosubhdr = nosubhdr
    )
  }
  
  names(sampleList) <- sampleNames
  sampleList
}


get_samples_names_spc <- function(foldername) {
  files <- list.files(
    foldername,
    pattern = "\\.[Ss][Pp][Cc]$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  samples.names <- gsub("\\.[^.]*$", "", basename(files))
  samples.names
}
