
library(shiny)
library(ggplot2)

# location of cal data
PATH_TO_CAL_FILES <- "..\\calibration_data\\ms1\\pos"

# Warning levels
POS_INT_RANGE <- c(1.4e6, 2.1e6)
POS_INT_RANGE_MS2 <- c(4.5e4, 6.0e4)
NEG_INT_RANGE <- c(1.15e6, 1.85e6)
NEG_INT_RANGE_MS2 <- c(1.5e2, 3.0e2)

MAX_MASS_SHIFT_PPM <- 1.3
MAX_MASS_SHIFT_MS2_PPM <- 5
MAX_PEAKWIDTH_MDA <- 17
MAX_PEAKWIDTH_MS2_MDA <- 15

# Weighting
weighting <- rbind(
  data.frame(target = "609", weight = 0.5),
  data.frame(target = "618", weight = 10),
  data.frame(target = "922", weight = 20)
)

# Intensity warning levels after weighting
POS_INT_RANGE_CORR <- c(1.0e6, 1.7e6)
POS_INT_RANGE_MS2_CORR <- c(4.5e4, 6.0e4)
NEG_INT_RANGE_CORR <- c(1e6, 1.75e6)
NEG_INT_RANGE_MS2_CORR <- c(1.5e2, 3.0e2)

# get all files and info
getDat <- function(folderN) {
  #browser()
  cnames <- c("Target_mz", "MZ_found", "Intensity_cps", "Width_Da", "Mass_shift_ppm")
  dat <- lapply(list.files(folderN, full.names = T), read.delim, col.names = cnames)
  info <- file.info(list.files(folderN, full.names = T))
  
  # add time to dataframes and rbind
  dat <- Map(function(x, y) transform(
    x, 
    target = stringr::str_pad(as.character(round(x[, 1])), 3, pad = "0"), 
    time = y
    ), dat, info$mtime)
  
  # apply weighting for intensities
  applyWeight <- function(x) {
    x$Intensity_corr_cps <- NA
    for (r in seq_len(nrow(weighting))) {
      tmass <- weighting[r, "target"]
      w <- weighting[r, "weight"]
      if (tmass %in% x$target) {
        x[x$target == tmass, "Intensity_corr_cps"] <- x[x$target == tmass, "Intensity_cps"] * w
      } 
    }
    x$Intensity_corr_cps <- ifelse(is.na(x$Intensity_corr_cps),
                                     x$Intensity_cps, x$Intensity_corr_cps)
    x
  }
  dat <- lapply(dat, applyWeight)
  
  # calc means
  calcAv <- function(x) {
    newRow <- as.data.frame(t(apply(x[, c("Target_mz", "MZ_found", "Intensity_cps", "Width_Da", "Intensity_corr_cps")], 2, mean)))
    newRow$Mass_shift_ppm <- mean(abs(x$Mass_shift_ppm))
    newRow$target <- "mean"
    newRow$time <- x$time[1]
    rbind(x, newRow)
  }
  dat <- lapply(dat, calcAv)
  dat <- Reduce(rbind, dat)
  
  
  
  
  # calculate mDa shift in mass
  dat <- transform(dat, Mass_shift_mDa = abs(Mass_shift_ppm / 1e6 * MZ_found * 1e3))
  dat <- transform(dat, Mass_shift_ppm = abs(Mass_shift_ppm))
  dat <- transform(dat, Width_mDa = Width_Da * 1e3)
  
  dat
}

dat1 <- getDat(PATH_TO_CAL_FILES) 

fixDates <- function(dates) {
  attr(dates, "tzone") <- "Europe/Paris"
  c(dates[1]-10, dates[2]+10)
}

addRelIntensity <- function(df, dates) {
  calcRel <- function(partDf) {
    startRow <- which.min(abs(partDf$time - dates[1]))
    endRow <- which.min(abs(partDf$time - dates[2]))
    while (partDf$time[startRow] < dates[1])  # start date must come after dates[1]
      startRow <- startRow + 1
    intenStart <- partDf[startRow, "Intensity_cps"]
    stopifnot(length(intenStart) == 1)
    partDf$Intensity_rel <- partDf$Intensity_cps / intenStart
    partDf <- partDf[startRow:endRow, ]
    partDf
  }
  l <- by(df, df$target, calcRel, simplify = F)
  df <- Reduce(rbind, l)
  df
}

makeTrendPlot <- function(df, dates, type) {
  #browser()
  
  # Calculate relative intensity according to chosen start date
  if (type == "Intensity_rel") {
    df <- addRelIntensity(df, dates)
    # recalculate mean based on rel intensities
    df[df$target == "mean", "Intensity_rel"] <- as.vector(by(df, df$time, function(x) {
      mean(x[x$target != "mean", "Intensity_rel"], na.rm = T)
    }))
  }
  ploti <- ggplot(df, aes_(quote(time), as.name(type))) + geom_point() + 
    geom_line() + facet_wrap(~ target, scales = "free") +
    xlab("Time") +
    scale_x_datetime(limits = dates) + theme_grey(18)
  if (type == "Intensity_cps") 
    ploti <- ploti + scale_y_continuous(labels = function(x) format(x, digits = 1, scientific = TRUE))
  ploti
}
makeBellCurve <- function(df, dates, type, pol_i, ms_i) {
  # Calculate relative intensity according to chosen start date
  if (type == "Intensity_rel") {
    df <- addRelIntensity(df, dates)
    # recalculate mean based on relative intensities
    df[df$target == "mean", "Intensity_rel"] <- as.vector(by(df, df$time, function(x) {
      mean(x[x$target != "mean", "Intensity_rel"], na.rm = T)
    }))
  } else {
    df$Intensity_rel <- NA  # otherwise mean calculation will not work
  }
  
  
  newDat <- df[df$target == "mean", ]
  newDat <- newDat[newDat$time >= dates[1] &
                     newDat$time <= dates[2], ]
  types <- c("Intensity_cps", "Mass_shift_ppm", "Mass_shift_mDa", "Width_mDa", "Intensity_rel", "Intensity_corr_cps")
  means <- vapply(types, function(x) mean(newDat[, x], na.rm = T, trim = .15), numeric(1))
  stdevs <- vapply(types, function(x) sd(newDat[, x], na.rm = T), numeric(1))
  ploti <- ggplot(newDat, aes_(as.name(type))) + geom_density() + 
    geom_vline(xintercept = newDat[which.max(newDat$time), type]) + 
    geom_vline(xintercept = means[type], color = "blue", alpha = .2) +
    geom_vline(xintercept = means[type] + stdevs[type], color = "red", alpha = .2) +
    geom_vline(xintercept = means[type] - stdevs[type], color = "red", alpha = .2) +
    ylab("Density") +
    annotate("text", means[type], y = Inf, label = "mean", color = "blue", alpha = .2, hjust = -.1,
             vjust = 1) +
    annotate("text", means[type] + stdevs[type], y = Inf, label = "1 sd", color = "red", alpha = .2, 
             hjust = -.1, vjust = 1) +
    annotate("text", means[type] - stdevs[type], y = Inf, label = "1 sd", color = "red", alpha = .2,
             hjust = -.1, vjust = 1) +
    annotate("text", newDat[which.max(newDat$time), type], y = Inf, label = "latest", 
             color = "black", hjust = -.1, vjust = 1) +
    theme_grey(18)
  # Warning levels for intensity
  if (type == "Intensity_cps") {
    switch(paste0(pol_i, ms_i),
           posMS1 = {
             ploti <- ploti + 
               geom_vline(xintercept = POS_INT_RANGE[1], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE[1], Inf, label = "low warning", 
                      color = "darkred", alpha = .3, hjust = -.1, vjust = 1) + 
               geom_vline(xintercept = POS_INT_RANGE[2], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE[2], Inf, label = "high warning", 
                        color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
           },
           posMS2 = {
             ploti <- ploti + 
               geom_vline(xintercept = POS_INT_RANGE_MS2[1], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE_MS2[1], Inf, label = "low warning", 
                        color = "darkred", alpha = .3, hjust = -.1, vjust = 1) + 
               geom_vline(xintercept = POS_INT_RANGE_MS2[2], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE_MS2[2], Inf, label = "high warning", 
                        color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
           },
           negMS1 = {
             ploti <- ploti + geom_vline(xintercept = NEG_INT_RANGE[1], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE[1], Inf, label = "low warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1) + geom_vline(xintercept = NEG_INT_RANGE[2], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE[2], Inf, label = "high warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1)
           },
           negMS1 = {
             ploti <- ploti + geom_vline(xintercept = NEG_INT_RANGE_MS2[1], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_MS2[1], Inf, label = "low warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1) + geom_vline(xintercept = NEG_INT_RANGE_MS2[2], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_MS2[2], Inf, label = "high warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1)
           })
    
  }
  
  if (type == "Intensity_corr_cps") {
    switch(paste0(pol_i, ms_i),
           posMS1 = {
             ploti <- ploti + 
               geom_vline(xintercept = POS_INT_RANGE_CORR[1], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE_CORR[1], Inf, label = "low warning", 
                        color = "darkred", alpha = .3, hjust = -.1, vjust = 1) + 
               geom_vline(xintercept = POS_INT_RANGE_CORR[2], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE_CORR[2], Inf, label = "high warning", 
                        color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
           },
           posMS2 = {
             ploti <- ploti + 
               geom_vline(xintercept = POS_INT_RANGE_MS2_CORR[1], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE_MS2_CORR[1], Inf, label = "low warning", 
                        color = "darkred", alpha = .3, hjust = -.1, vjust = 1) + 
               geom_vline(xintercept = POS_INT_RANGE_MS2_CORR[2], color = "darkred", 
                          alpha = .3) +
               annotate("text", POS_INT_RANGE_MS2_CORR[2], Inf, label = "high warning", 
                        color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
           },
           negMS1 = {
             ploti <- ploti + geom_vline(xintercept = NEG_INT_RANGE_CORR[1], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_CORR[1], Inf, label = "low warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1) + geom_vline(xintercept = NEG_INT_RANGE_CORR[2], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_CORR[2], Inf, label = "high warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1)
           },
           negMS1 = {
             ploti <- ploti + geom_vline(xintercept = NEG_INT_RANGE_MS2_CORR[1], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_MS2_CORR[1], Inf, label = "low warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1) + geom_vline(xintercept = NEG_INT_RANGE_MS2_CORR[2], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_MS2_CORR[2], Inf, label = "high warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1)
           })
    
  }
  
  # warning level for mass shift in ppm (mDa is ignored)
  if (type == "Mass_shift_ppm" && ms_i == "MS1") {
    ploti <- ploti + 
      geom_vline(xintercept = MAX_MASS_SHIFT_PPM, color = "darkred", 
                 alpha = .3) +
      annotate("text", MAX_MASS_SHIFT_PPM, Inf, label = "high warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  }
  if (type == "Mass_shift_ppm" && ms_i == "MS2") {
    ploti <- ploti + 
      geom_vline(xintercept = MAX_MASS_SHIFT_MS2_PPM, color = "darkred", 
                 alpha = .3) +
      annotate("text", MAX_MASS_SHIFT_MS2_PPM, Inf, label = "high warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  }
  # warning level for peak width
  if (type == "Width_mDa" && ms_i == "MS1") {
    ploti <- ploti + 
      geom_vline(xintercept = MAX_PEAKWIDTH_MDA, color = "darkred", 
                 alpha = .3) +
      annotate("text", MAX_PEAKWIDTH_MDA, Inf, label = "high warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  }
  if (type == "Width_mDa" && ms_i == "MS2") {
    ploti <- ploti + 
      geom_vline(xintercept = MAX_PEAKWIDTH_MS2_MDA, color = "darkred", 
                 alpha = .3) +
      annotate("text", MAX_PEAKWIDTH_MS2_MDA, Inf, label = "high warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  }
  
  
  ploti
}



ui <- fluidPage(
  fluidRow(
    column(3, h2("Performance 6600 'Wasser'")),
    column(7, sliderInput("dates", "", timeFormat = "%F %R", value = c(min(dat1$time), max(dat1$time)),
                     min = min(dat1$time), max = max(dat1$time), width = "90%")),
    column(1, radioButtons("pol", "", choices = c("pos", "neg"), selected = "pos", inline = F)),
    column(1, radioButtons("ms", "", choices = c("MS1", "MS2"), selected = "MS1", inline = F))
  ),
  fluidRow(column(12,
      tabsetPanel(
        tabPanel("Intensity_cps", 
                 fluidRow(plotOutput("inten")), fluidRow(plotOutput("intenHist"))),
        tabPanel("Intensity_rel",
                 fluidRow(plotOutput("intenRel")), fluidRow(plotOutput("intenRelHist"))),
        tabPanel("Intensity_corr",
                 fluidRow(plotOutput("intenCorr")), fluidRow(plotOutput("intenCorrHist"))),
        tabPanel("m/z shift mDa", 
                 fluidRow(plotOutput("mDa")), fluidRow(plotOutput("mDaHist"))
                 ),
        tabPanel("m/z shift ppm", 
                 fluidRow(plotOutput("ppm")), fluidRow(plotOutput("ppmHist"))
                 ),
        tabPanel("Spec peak width mDa", 
                 fluidRow(plotOutput("width")), fluidRow(plotOutput("widthHist"))
                 )
      )
    )
  )
)

server <- function(input, output, session) {
   dat <- reactive({
     fol <- switch(paste0(input$pol, input$ms), 
                   posMS1 = "..\\calibration_data\\ms1\\pos",
                   posMS2 = "..\\calibration_data\\ms2\\pos",
                   negMS1 = "..\\calibration_data\\ms1\\neg",
                   negMS2 = "..\\calibration_data\\ms2\\neg")
     getDat(fol)
   })
   
   observe({
     updateSliderInput(session, "dates", timeFormat = "%F %R", min = min(dat()$time), max = max(dat()$time),
                       value = c(min(dat()$time), max(dat()$time)))  
   })
   
   output$inten <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "Intensity_cps"))
   output$intenHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "Intensity_cps", 
                                                input$pol, input$ms))
   
   output$intenRel <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "Intensity_rel"))
   output$intenRelHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "Intensity_rel", 
                                                   input$pol, input$ms))
   
   output$intenCorr <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "Intensity_corr_cps"))
   output$intenCorrHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "Intensity_corr_cps", 
                                                   input$pol, input$ms))
   
   output$mDa <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "Mass_shift_mDa"))
   output$mDaHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "Mass_shift_mDa", 
                                              input$pol, input$ms))
   
   output$ppm <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "Mass_shift_ppm"))
   output$ppmHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "Mass_shift_ppm", 
                                              input$pol, input$ms))
   
   output$width <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "Width_mDa"))
   output$widthHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "Width_mDa", 
                                                input$pol, input$ms))
}

shinyApp(ui = ui, server = server)

