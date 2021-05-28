# function to read in cal report files from Analyst software


# get all files and info
getDat <- function(folderN) {
  
  cnames <- c("selected", "used", "mz_expected", "min_search_mz", "max_search_mz", 
              "mz_before_cal", "delta_mz_before_cal", "ppm_error_before_cal", 
              "mz_after_cal", "delta_mz_after_cal", "ppm_error_after_cal", 
              "intensity_cps", "resolution", "charge_state", "monoisotopic", "X")
  readFile <- function(filePath) {
    #browser()
    suppressWarnings({
      fileLines <- readLines(filePath, skipNul = T)
    })
    
    if (any(grepl("Calibration Failed!", fileLines))) 
      return(NULL) 
    
    if (!any(grepl("(Pos)|(Neg)", fileLines))) 
      return(NULL)
    
    posnegPosition <- which(grepl("(Pos)|(Neg)", fileLines))[1]
    posneg <- stringr::str_match(fileLines[posnegPosition], "([PosNeg]{3})")[,2]
    posneg <- stringr::str_to_lower(posneg)
    if (!(posneg %in% c("pos", "neg")))
      stop("error in file:", filePath)
    if (posneg == "pos") {
      startMS1 <- which(grepl("Selected", fileLines))[1] - 1 
      endMS1 <- which(sapply(fileLines, nchar) == 0)
      endMS1 <- endMS1[endMS1 > startMS1][1]
      numMS1 <- endMS1 - startMS1 - 2
      suppressWarnings({
        dataMs1 <- read.delim(filePath, skip = startMS1, nrows = numMS1, col.names = cnames, skipNul = T)
      })
      
      dataMs1$msLevel <- 1
      startMS2 <- which(grepl("Selected", fileLines))[2] - 1
      endMS2 <- which(sapply(fileLines, nchar) == 0)
      endMS2 <- endMS2[endMS2 > startMS2][1]
      numMS2 <- endMS2 - startMS2 - 2
      suppressWarnings({
        dataMs2 <- read.delim(filePath, skip = startMS2, nrows = numMS2, col.names = cnames, skipNul = T)
      })
      
      dataMs2$msLevel <- 2
      dats <- rbind(dataMs1, dataMs2)
      dats$pol <- posneg
    } else if (posneg == "neg") {
      startMS1 <- which(grepl("Selected", fileLines))[1] - 1 
      endMS1 <- which(sapply(fileLines, nchar) == 0)
      endMS1 <- endMS1[endMS1 > startMS1][1]
      numMS1 <- endMS1 - startMS1 - 2
      suppressWarnings({
        dataMs1 <- read.delim(filePath, skip = startMS1, nrows = numMS1, col.names = cnames, skipNul = T)
      })
      
      dataMs1$msLevel <- 1
      startMS2 <- which(grepl("Selected", fileLines))[2] - 1
      endMS2 <- which(sapply(fileLines, nchar) == 0)
      endMS2 <- endMS2[endMS2 > startMS2][1]
      numMS2 <- endMS2 - startMS2 - 2
      suppressWarnings({
        dataMs2 <- read.delim(filePath, skip = startMS2, nrows = numMS2, col.names = cnames, skipNul = T)
      })
      
      dataMs2$msLevel <- 2
      dats <- rbind(dataMs1, dataMs2)
      dats$pol <- posneg
    } else {
      stop("polarity not found")
    }
    
    dats <- dats[dats$used == "Yes", c("mz_expected", "ppm_error_after_cal", "intensity_cps", "resolution", "pol", "msLevel")]
    numericCols <- c("mz_expected", "ppm_error_after_cal", "intensity_cps", "resolution", "msLevel")
    dats[, numericCols] <- lapply(dats[, numericCols], as.numeric)
    dats  
  }
  allcsvs <- list.files(folderN, full.names = T, pattern = "\\.txt$")
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
  
  # calc means
  calcAv <- function(x) {
    newRow <- as.data.frame(t(apply(x[, c("mz_expected", "intensity_cps", "resolution", "intensity_corr_cps")], 2, mean)))
    newRow$ppm_error_after_cal <- mean(abs(x$ppm_error_after_cal))
    newRow$target <- "mean"
    newRow$time <- x$time[1]
    newRow$pol <- x$pol[1]
    newRow$msLevel <- x$msLevel[1]
    rbind(x, newRow)
  }
  #browser()
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
  #browser()
  # sort by time
  dat <- dat[order(dat$time), ]
  dat
}