##' Import for Thermo Galactic's spc file format
##' These functions allow to import .spc files.
##' A detailed description of the .spc file format is available at
##'
##' @param filename The complete file name of the .spc file.
##' @param keys.hdr2data,keys.hdr2log,keys.log2data,keys.log2log character
##'   vectors with the names of parameters in the .spc file's log block
##'   (log2xxx) or header (hdr2xxx) that should go into the extra data
##'   (yyy2data) or into the \code{long.description} field of the returned
##'   object's log (yyy2log).
##' @param log.txt Should the text part of the .spc file's log block be read?
##' @param log.bin,log.disk Should the normal and on-disk binary parts of the
##'   .spc file's log block be read?
##' @param hdr A list with fileheader fields that overwrite the settings of
##'   actual file's header.
##' @param nosubhdr Boolean value to decide if the header should be read or not.
##' @param no.object If \code{TRUE}, a list with wavelengths, spectra, labels,
##'   log and data are returned instead of a hyperSpec object.
##' @return A list with imported spectra information.
##' @author C. Beleites
##' @rdname read-spc
##' @references Reference information for the SPC file format.
##' @export
read_spc_nosubhdr <- function(filename,
                              keys.hdr2data = c("fexper", "fres", "fsource"),
                              keys.hdr2log = c("fdate", "fpeakpt"),
                              keys.log2data = FALSE, keys.log2log = TRUE,
                              log.txt = TRUE, log.bin = FALSE, log.disk = FALSE,
                              hdr = list(), nosubhdr = TRUE,
                              no.object = FALSE) {
  f <- readBin(filename, "raw", file.info(filename)$size, 1)
  
  hdr <- modifyList(.spc.filehdr(f), hdr)
  fpos <- hdr$.last.read
  
  if (!hdr$ftflgs["TXYXYS"]) {
    if (!hdr$ftflgs["TXVALS"]) {
      wavelength <- seq(hdr$ffirst, hdr$flast, length.out = hdr$fnpts)
    } else {
      if (!hdr$ftflgs["TMULTI"]) {
        tmp <- .spc.read.x(f, fpos, hdr$fnpts)
        wavelength <- tmp$x
        fpos <- tmp$.last.read
      }
    }
  }
  
  label <- list(.wavelength = hdr$fxtype, spc = hdr$fytype,
                z = hdr$fztype, z.end = hdr$fztype)
  
  data <- list(z = NA, z.end = NA)
  if (hdr$fwplanes > 0) {
    data <- c(data, w = NA)
  }
  
  tmp <- .spc.log(f, hdr$flogoff, log.bin, log.disk, log.txt, keys.log2data, keys.log2log)
  
  log <- list(
    short = "read.spc",
    long = list(
      call = match.call(),
      log = tmp$log.long,
      header = getbynames(hdr, keys.hdr2log)
    )
  )
  
  data <- c(data, tmp$extra.data, getbynames(hdr, keys.hdr2data))
  
  if (hdr$ftflgs["TXYXYS"] && hdr$ftflgs["TMULTI"]) {
    spc <- list()
    data <- as.data.frame(data)
  } else {
    spc <- matrix(NA, nrow = hdr$fnsub, ncol = hdr$fnpts)
    data <- as.data.frame(lapply(data, rep, hdr$fnsub))
  }
  
  if (hdr$subfiledir) {
    hdr$subfiledir <- .spc.subfiledir(f, hdr$subfiledir, hdr$fnsub)
    
    for (s in seq_len(hdr$fnsub)) {
      if (!nosubhdr) {
        hdr <- .spc.subhdr(f, hdr$subfiledir$ssfposn[s], hdr)
      }
      
      fpos <- hdr$.last.read
      wavelength <- .spc.read.x(f, fpos, hdr$fnpts)
      fpos <- wavelength$.last.read
      
      y <- .spc.read.y(f, fpos, npts = hdr$fnpts, exponent = hdr$fexp,
                       word = hdr$ftflgs["TSPREC"])
      fpos <- y$.last.read
      
      if (!nosubhdr) {
        data$z <- hdr$subhdr$subtime
        data$z.end <- hdr$subhdr$subnext
        if (hdr$fwplanes > 0) {
          data$w <- hdr$subhdr$w
        }
      }
      
      if (!exists("wavelength")) {
        .spc.error(
          "read.spc", list(hdr = hdr),
          "wavelength not read. This may be caused by wrong header information."
        )
      }
      
      spc[[s]] <- list(
        spc = y$y,
        wavelength = wavelength$x,
        data = data,
        log = log,
        labels = label
      )
    }
  } else {
    for (s in seq_len(hdr$fnsub)) {
      if (!nosubhdr) {
        hdr <- .spc.subhdr(f, fpos, hdr)
      }
      
      fpos <- hdr$.last.read
      tmp <- .spc.read.y(f, fpos, npts = hdr$fnpts, exponent = hdr$fexp,
                         word = hdr$ftflgs["TSPREC"])
      fpos <- tmp$.last.read
      
      spc[s, ] <- tmp$y
      
      if (!nosubhdr) {
        data[s, c("z", "z.end")] <- unlist(hdr$subhdr[c("subtime", "subnext")])
        if (hdr$fwplanes > 0) {
          data[s, "w"] <- hdr$subhdr$w
        }
      }
    }
  }
  
  if (hdr$ftflgs["TXYXYS"] && hdr$ftflgs["TMULTI"]) {
    spc
  } else if (no.object) {
    list(spc = spc, wavelength = wavelength, data = data, log = log, labels = label)
  } else {
    list(
      spc = spc,
      wavelength = wavelength,
      data = data[rep(1, hdr$fnsub), ],
      log = log,
      labels = label
    )
  }
}

getbynames <- function(x, e) {
  x <- x[e]
  if (length(x) > 0) {
    if (is.character(e)) {
      names(x) <- e
    }
    x[sapply(x, is.null)] <- NA
    x
  } else {
    list()
  }
}

split_line <- function(x, separator, trim.blank = TRUE) {
  tmp <- regexpr(separator, x)
  
  key <- substr(x, 1, tmp - 1)
  value <- substr(x, tmp + 1, nchar(x))
  
  if (trim.blank) {
    blank.pattern <- "^[[:blank:]]*([^[:blank:]]+.*[^[:blank:]]+)[[:blank:]]*$"
    key <- sub(blank.pattern, "\\1", key)
    value <- sub(blank.pattern, "\\1", value)
  }
  
  value <- as.list(value)
  names(value) <- key
  value
}

split_string <- function(x, separator, trim.blank = TRUE, remove.empty = TRUE) {
  pos <- gregexpr(separator, x)
  if (length(pos) == 1 && pos[[1]] == -1) {
    return(x)
  }
  
  pos <- pos[[1]]
  
  pos <- matrix(
    c(1, pos + attr(pos, "match.length"),
      pos - 1, nchar(x)),
    ncol = 2
  )
  
  if (pos[nrow(pos), 1] > nchar(x)) {
    pos <- pos[-nrow(pos), ]
  }
  
  x <- apply(pos, 1, function(p, x) substr(x, p[1], p[2]), x)
  
  if (trim.blank) {
    blank.pattern <- "^[[:blank:]]*([^[:blank:]]+.*[^[:blank:]]+)[[:blank:]]*$"
    x <- sub(blank.pattern, "\\1", x)
  }
  
  if (remove.empty) {
    x <- x[sapply(x, nchar) > 0]
  }
  
  x
}
