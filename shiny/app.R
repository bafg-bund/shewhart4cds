
# performance monitor viewer for 6600 wasser using automatically produced result files


library(shiny)
library(ggplot2)

source("../performance_monitor_settings.R")

stopifnot(calReportType %in% c("analyst" ,"sciexos"))

source(paste0("../get_cal_data_", calReportType, ".R"))

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
    intenStart <- partDf[startRow, "intensity_cps"]
    stopifnot(length(intenStart) == 1)
    partDf$intensity_rel <- partDf$intensity_cps / intenStart
    partDf <- partDf[startRow:endRow, ]
    partDf
  }
  l <- by(df, df$target, calcRel, simplify = F)
  df <- Reduce(rbind, l)
  df
}

makeTrendPlot <- function(df, dates, type) {
  #browser(expr = type == "intensity_rel")
  
  # Calculate relative intensity according to chosen start date
  if (type == "intensity_rel") {
    df <- addRelIntensity(df, dates)
    # recalculate mean based on rel intensities
    df[df$target == "mean", "intensity_rel"] <- as.vector(by(df, df$time, function(x) {
      mean(x[x$target != "mean", "intensity_rel"], na.rm = T)
    }))
  }
  ploti <- ggplot(df, aes_(quote(time), as.name(type))) + geom_point() + 
    geom_line() + facet_wrap(~ target, scales = "free") +
    xlab("Time") +
    scale_x_datetime(limits = dates) + theme_grey(18)
  if (type == "intensity_cps") 
    ploti <- ploti + scale_y_continuous(labels = function(x) format(x, digits = 1, scientific = TRUE))
  ploti
}

makeBellCurve <- function(df, dates, type, pol_i, ms_i) {
  # Calculate relative intensity according to chosen start date
  if (type == "intensity_rel") {
    df <- addRelIntensity(df, dates)
    # recalculate mean based on relative intensities
    df[df$target == "mean", "intensity_rel"] <- as.vector(by(df, df$time, function(x) {
      mean(x[x$target != "mean", "intensity_rel"], na.rm = T)
    }))
  } else {
    df$intensity_rel <- NA  # otherwise mean calculation will not work
  }
  
  newDat <- df[df$target == "mean", ]
  newDat <- newDat[newDat$time >= dates[1] &
                     newDat$time <= dates[2], ]
  types <- c("intensity_cps", "intensity_corr_cps", "ppm_error_after_cal", "mDa_error_after_cal", "resolution", "intensity_rel")
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
  
  if (type == "intensity_cps") {
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
           negMS2 = {
             ploti <- ploti + geom_vline(xintercept = NEG_INT_RANGE_MS2[1], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_MS2[1], Inf, label = "low warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1) + geom_vline(xintercept = NEG_INT_RANGE_MS2[2], color = "darkred", alpha = .3) +
               annotate("text", NEG_INT_RANGE_MS2[2], Inf, label = "high warning", color = "darkred", alpha = .3, 
                        hjust = -.1, vjust = 1)
           })
    
  }
  
   if (type == "intensity_corr_cps") {
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
  if (type == "ppm_error_after_cal" && ms_i == "MS1") {
    ploti <- ploti + 
      geom_vline(xintercept = MAX_MASS_SHIFT_PPM, color = "darkred", 
                 alpha = .3) +
      annotate("text", MAX_MASS_SHIFT_PPM, Inf, label = "high warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  }
  if (type == "ppm_error_after_cal" && ms_i == "MS2") {
    ploti <- ploti + 
      geom_vline(xintercept = MAX_MASS_SHIFT_MS2_PPM, color = "darkred", 
                 alpha = .3) +
      annotate("text", MAX_MASS_SHIFT_MS2_PPM, Inf, label = "high warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  }
  # warning level for resolution
  if (type == "resolution" && ms_i == "MS1" && pol_i == "pos") {
    ploti <- ploti + 
      geom_vline(xintercept = MIN_RESOLUTION_POS, color = "darkred", 
                 alpha = .3) +
      annotate("text", MIN_RESOLUTION_POS, Inf, label = "low warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  }
  if (type == "resolution" && ms_i == "MS2" && pol_i == "pos") {
    ploti <- ploti + 
      geom_vline(xintercept = MIN_RESOLUTION_MS2_POS, color = "darkred", 
                 alpha = .3) +
      annotate("text", MIN_RESOLUTION_MS2_POS, Inf, label = "low warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  }
  if (type == "resolution" && ms_i == "MS1" && pol_i == "neg") {
    ploti <- ploti + 
      geom_vline(xintercept = MIN_RESOLUTION_NEG, color = "darkred", 
                 alpha = .3) +
      annotate("text", MIN_RESOLUTION_NEG, Inf, label = "low warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  }
  if (type == "resolution" && ms_i == "MS2" && pol_i == "neg") {
    ploti <- ploti + 
      geom_vline(xintercept = MIN_RESOLUTION_MS2_NEG, color = "darkred", 
                 alpha = .3) +
      annotate("text", MIN_RESOLUTION_MS2_NEG, Inf, label = "low warning", 
               color = "darkred", alpha = .3, hjust = -.1, vjust = 1)
  }
  
  ploti
}

# Get measurement times for setting up the ui
getTimes <- function(folderN) {
  allcsvs <- list.files(folderN, full.names = T, pattern = "\\.txt$")
  info <- file.info(allcsvs)
  info$mtime
}
datTimes <- getTimes(CAL_REPORTS) 

ui <- fluidPage(
  fluidRow(
    column(3, h2(paste("Performance", instrumentName))),
    column(7, sliderInput("dates", "", timeFormat = "%F %R", value = c(min(datTimes, na.rm =T), max(datTimes, na.rm =T)),
                     min = min(datTimes, na.rm =T), max = max(datTimes, na.rm =T), width = "90%")),
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
        tabPanel("Resolution", 
                 fluidRow(plotOutput("width")), fluidRow(plotOutput("widthHist"))
                 )
      )
    )
  )
)

server <- function(input, output, session) {
  withProgress({
    dfAll <- getDat(CAL_REPORTS)
  }, message = "Loading cal files...", session = session, min = 1, value = 1)
  
   
  dat <- reactive({
     newdat <- switch(paste0(input$pol, input$ms), 
                   posMS1 = subset(dfAll, pol == "pos" & msLevel == 1),
                   posMS2 = subset(dfAll, pol == "pos" & msLevel == 2),
                   negMS1 = subset(dfAll, pol == "neg" & msLevel == 1),
                   negMS2 = subset(dfAll, pol == "neg" & msLevel == 2))
     newdat
   })

   output$inten <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "intensity_cps"))
   output$intenHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "intensity_cps", 
                                                input$pol, input$ms))
   
   output$intenRel <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "intensity_rel"))
   output$intenRelHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "intensity_rel",
                                                   input$pol, input$ms))
   
   output$intenCorr <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "intensity_corr_cps"))
   output$intenCorrHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "intensity_corr_cps",
                                                   input$pol, input$ms))
   
   output$mDa <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "mDa_error_after_cal"))
   output$mDaHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "mDa_error_after_cal", 
                                              input$pol, input$ms))
   
   output$ppm <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "ppm_error_after_cal"))
   output$ppmHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "ppm_error_after_cal", 
                                              input$pol, input$ms))
   
   output$width <- renderPlot(makeTrendPlot(dat(), fixDates(input$dates), "resolution"))
   output$widthHist <- renderPlot(makeBellCurve(dat(), fixDates(input$dates), "resolution", 
                                                input$pol, input$ms))
}

shinyApp(ui = ui, server = server)

