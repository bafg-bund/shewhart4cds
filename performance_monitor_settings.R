# Settings file for performance monitor app

instrumentName <- "6600 Sediment"

# Cal report type used can be either "analyst" or "sciexos" 
calReportType <- "analyst"

# Location of calibration report files
# CAL_REPORTS <- c("D:\\Analyst Data\\Projects\\Plastik\\Data\\Cal Data",
#                  "D:\\Analyst Data\\Projects\\Ute\\Data\\Cal Data")
CAL_REPORTS <- c("../example_files_analyst")

# Warning levels
POS_INT_RANGE <- c(2e5, 4e5)
POS_INT_RANGE_MS2 <- c(4.5e4, 6.0e4)
NEG_INT_RANGE <- c(9.05e5, 1.65e6)
NEG_INT_RANGE_MS2 <- c(1.5e2, 3.0e2)

MAX_MASS_SHIFT_PPM <- 1.3
MAX_MASS_SHIFT_MS2_PPM <- 5

MIN_RESOLUTION_POS <- 30000
MIN_RESOLUTION_MS2_POS <- 30000
MIN_RESOLUTION_NEG <- 29000
MIN_RESOLUTION_MS2_NEG <- 27000

# Weighting MS1 (optional, remove or comment-out if not needed)
weighting <- rbind(
  data.frame(target = "609", weight = 0.5),
  data.frame(target = "618", weight = 10),
  data.frame(target = "922", weight = 20)
)

# Intensity warning levels after weighting
POS_INT_RANGE_CORR <- c(2e5, 4e5)
POS_INT_RANGE_MS2_CORR <- c(4.5e4, 6.0e4)
NEG_INT_RANGE_CORR <- c(1e6, 1.75e6)
NEG_INT_RANGE_MS2_CORR <- c(1.5e2, 3.0e2)

