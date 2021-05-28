# function to read in cal report files from SciexOS software

getDat <- function(folderN) {
  #browser()
  # cnames <- c("Target_mz", "MZ_found", "Intensity_cps", "Width_Da", "Mass_shift_ppm")
  cnames <- c("used", "mz_expected", "min_search_mz", "max_search_mz", 
              "mz_before_cal", "delta_mz_before_cal", "ppm_error_before_cal", 
              "mz_after_cal", "delta_mz_after_cal", "ppm_error_after_cal", 
              "intensity_cps", "resolution", "charge_state", "monoisotopic", "X")
  readFile <- function(file) {
    if (!grepl("pos|neg", readLines(file, n = 2)[2]))
      return(NULL)
    #browser()
    posneg <- stringr::str_match(readLines(file, n = 2)[2], "([posneg]{3})")[,2]
    #browser()
    if (posneg == "pos") {
      dataMs1 <- read.csv(file, skip = 17, nrows = 11, col.names = cnames)
      dataMs1$msLevel <- 1
      dataMs2 <- read.csv(file, skip = 44, nrows = 7, col.names = cnames)
      dataMs2$msLevel <- 2
      dats <- rbind(dataMs1, dataMs2)
      dats$pol <- posneg
    } else if (posneg == "neg") {
      dataMs1 <- read.csv(file, skip = 17, nrows = 21, col.names = cnames)
      dataMs1$msLevel <- 1
      dataMs2 <- read.csv(file, skip = 54, nrows = 9, col.names = cnames)
      dataMs2$msLevel <- 2
      dats <- rbind(dataMs1, dataMs2)
      dats$pol <- posneg
    } else {
      stop("polarity not found")
    }
    dats <- dats[dats$used == "Yes", c("mz_expected", "ppm_error_after_cal", "intensity_cps", "resolution", "pol", "msLevel")]
    dats[, c("mz_expected", "ppm_error_after_cal", "intensity_cps", "resolution", "msLevel")] <-
      lapply(dats[, c("mz_expected", "ppm_error_after_cal", "intensity_cps", "resolution", "msLevel")], as.numeric)
    dats  
  }
  allcsvs <- list.files(folderN, full.names = T, pattern = "\\.csv$")
  dat <- lapply(allcsvs, readFile)
  info <- file.info(allcsvs)
  keep <- vapply(dat, Negate(is.null), logical(1))
  dat <- dat[keep]
  info <- info[keep,]
  
  # add time to dataframes and rbind
  dat <- Map(function(x, y) transform(
    x, 
    target = stringr::str_pad(as.character(round(x[, 1])), 3, pad = "0"), 
    time = y
  ), dat, info$mtime)
  
  # apply weighting for intensities
  applyWeight <- function(x) {
    x$intensity_corr_cps <- NA
    if (exists("weighting")) {
      for (r in seq_len(nrow(weighting))) {
        tmass <- weighting[r, "target"]
        w <- weighting[r, "weight"]
        if (tmass %in% x$target) {
          x[x$target == tmass & x$msLevel == 1, "intensity_corr_cps"] <- x[x$target == tmass & x$msLevel == 1, "intensity_cps"] * w
        }
      }
    }
    x$intensity_corr_cps <- ifelse(is.na(x$intensity_corr_cps),
                                   x$intensity_cps, x$intensity_corr_cps)
    x
  }
  dat <- lapply(dat, applyWeight)
  
  # Compute average
  calcAv <- function(x) {
    newRow <- as.data.frame(t(apply(x[, c("mz_expected", "intensity_cps", "resolution", "intensity_corr_cps")], 2, mean)))
    newRow$ppm_error_after_cal <- mean(abs(x$ppm_error_after_cal))
    newRow$target <- "mean"
    newRow$time <- x$time[1]
    newRow$pol <- x$pol[1]
    newRow$msLevel <- x$msLevel[1]
    rbind(x, newRow)
  }
  
  datms1 <- lapply(dat, function(y) y[y$msLevel == 1,])
  datms2 <- lapply(dat, function(y) y[y$msLevel == 2,])
  datms1 <- lapply(datms1, calcAv)
  datms2 <- lapply(datms2, calcAv)
  datms1 <- do.call("rbind", datms1)
  datms2 <- do.call("rbind", datms2)
  dat <- rbind(datms1,datms2)
  
  # calculate mDa shift in mass
  dat <- transform(dat, mDa_error_after_cal = abs(ppm_error_after_cal / 1e6 * mz_expected * 1e3))
  dat <- transform(dat, ppm_error_after_cal = abs(ppm_error_after_cal))

  dat <- dat[order(dat$time), ]
  dat
}