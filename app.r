# Sinclair Z Shiny Tester
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
packages = c("shiny",       # interactive components
             "shinyBS",     # pop up help boxes
             "ggplot2",     # nice graphics
             "dplyr",       # data handling
             "reshape2",    # data handling
             "tidyr",       # data handling
             "Hmisc")       # error bar plotting

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# dimensions 
nyears <- 40
nages <- 30

#-------------------------------------------------------------------------------------
# function calc_Sinclair_Z
# estimates catch curve Z using moving window with cohort-specific coefficients
# matrix has year in rows and age in columns
# window.years determines the number of years in each moving window (defaults to 4)
calc_Sinclair_Z <- function(mat, window.years=4){
  res <- list()
  mat <- as.matrix(mat)
  if (all(is.na(mat))){
    res$error <- TRUE
    return(res)
  }
  year <- as.numeric(rownames(mat))
  age <- as.numeric(colnames(mat))
  ny <- length(year)
  nw <- ny - window.years + 1 # number of windows to be evaluated
  est.Sinclair.Z <- matrix(NA, nrow=nw, ncol=3)
  plot.year <- rep(NA, nw)
  for (i in 1:nw){
    submat <- mat[year %in% year[i:(i+window.years-1)],]
    data <- melt(submat)
    colnames(data) <- c("Year","Age","Value")
    data$cohort <- data$Year - data$Age
    data$Value[data$Value == 0] <- NA
    data <- data[!is.na(data$Value) ,]
    data$lnVal <- log(data$Value)
    can.calc <- FALSE
    if (length(data[,1]) >= 2){
      if (max(table(data$cohort)) >= 3){
        can.calc <- TRUE
      }
    }
    if (can.calc == TRUE){
      my.lm <- lm(data$lnVal ~ as.factor(data$cohort) + data$Age)
      data$pred <- predict(my.lm)
      data$resid <- residuals(my.lm)
      res[[i]] <- data
      est.Sinclair.Z[i,1] <- -1 * my.lm$coefficients[names(my.lm$coefficients) == "data$Age"]
      est.Sinclair.Z[i,2:3] <- -1 * rev(confint(my.lm, "data$Age", level=0.90))
    }
    else{  
      res[[i]] <- data
      est.Sinclair.Z[i,] <- rep(NA, 3)
    }
    plot.year[i] <- year[i] + (window.years - 1) / 2 
  }
  res$est.Sinclair.Z <- est.Sinclair.Z
  res$plot.year <- plot.year
  res$error <- FALSE
  res$window.years <- window.years
  return(res)
}

#-------------------------------------------------------------------------------------
give.n <- function(x){
  return(c(y = median(x)+0.10, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Define UI for application that draws a histogram
ui <- navbarPage("Sinclair Z Shiny Tester",  
                 
   tabPanel("Introduction",
            sidebarLayout(
              sidebarPanel(
                fluidRow(
                  column(12,
                         h2("Welcome"),
                         br(),
                         p("This Shiny app is designed to let you explore the performance of the Sinclair total mortality (Z) estimator using simulated data that you create. The Sinclair Z approach is similar to a catch curve, except that it includes multiple cohorts from adjacent years and estimates a single total mortality value for the group of cohorts (a common slope with different intercepts for each cohort). A moving window of 3 to 6 years is used and you decide which ages to include within each window. Diagnostics are provided to allow you to see the consequences of breaking the assumptions of constant Z within a window and full selectivity for all selected ages, as well as the influence of sampling variability. The interactive nature of the app allows quick feedback as you make a change in the data or options. Pop up windows provide information about the inputs and plots as you roll the mouse over them."),
                         br(),
                         p("Have fun!")
                         )
                )
              ),
              mainPanel(
                h3("Data Creator"),
                p("This is where you create the data by moving sliders and selecting options that define the population, fishery, and survey characteristics. You have the ability to make the fishery and survey selectivities logistic, domed, or full at all ages. You also determine how the fishing mortality changes over time (M is assumed equal to 0.2 for all years and ages, of course). You can add noise about the trend in F over time, including autocorrelation in the noise if you wish. The randomg number generator (RNG) seed value choice allows you to return to a single setting or explore different age ranges for estimating Z on a fixed data set."),
                br(),
                h3("Estimate Z"),
                p("This is where you make the choice of which ages to use in the estimation of Z. Use the Diagnostic tab to see if the residuals are centered on zero. A partially selected age can produce a distribution of residuals that does not overlap the zero line. The number of years in the moving windows can be changed from the default value of 4 to examine the impact on estimates and diagnostics. Unlike the real world, the Z estimate plot shows the actual Z trend in red, so you can see how well the estimator performed. This should help you understand how much to pay attention to individual years using the Sinclair Z estimator. The individual fits can be paged through using the four buttons (display first, next, previous, and last fit), so that you can look at unusual or typical year windows in the time series of Z estimates.")
              )
            )
   ),
   
   tabPanel("Data Creator",
            sidebarLayout(
               sidebarPanel(
                 fluidRow(
                   column(6,
                          sliderInput("rngSeed",
                                      label = "RNG seed",
                                      min = 1,
                                      max = 1000,
                                      step = 1,
                                      value = 100),
                          
                          bsTooltip("rngSeed",
                                    "Changes recruitment and survey observations (multinomial)",
                                    "right"),
                          
                          sliderInput("FisheryFullSelectivity",
                                      label = "Fishery Full Selectivity",
                                      min = 1,
                                      max = nages,
                                      step = 1,
                                      value = c(9, 25)),
                          
                          bsTooltip("FisheryFullSelectivity",
                                    "Use sliders to determine ages that are and are not fully selected",
                                    "right"),
                          
                          sliderInput("FisherySelectivityWidth",
                                      label = "Fishery Sel Width",
                                      min = 0.1,
                                      max = 5,
                                      step = 0.1,
                                      value = 2.0),
                          
                          bsTooltip("FisherySelectivityWidth",
                                    "Determines how steep logistic curve is",
                                    "right"),
                          
                          sliderInput("StartingF",
                                      label = "F in first year",
                                      min = 0,
                                      max = 1,
                                      step = 0.1,
                                      value = 0.2),
                          
                          bsTooltip("StartingF",
                                    "Change F over time starting at this value",
                                    "right"),
                          
                          sliderInput("EndingF",
                                      label = "F in last year",
                                      min = 0,
                                      max = 1,
                                      step = 0.1,
                                      value = 0.6),
                          
                          bsTooltip("EndingF",
                                    "Change F over time ending at this value",
                                    "right"),
                          
                          radioButtons("changeF",
                                       label = "How F changes",
                                       choices = c("step","linear"),
                                       selected = "step",
                                       inline = FALSE),
                          
                          bsTooltip("changeF",
                                    "Does F change all at once halfway through time period or slowly over time",
                                    "right")
                          
                          ),
                   
                   column(6,
                          selectInput("surveyESS",
                                      label = "survey ESS",
                                      choices = c(50, 100, 1000, 10000),
                                      selected = 1000),
                          
                          bsTooltip("surveyESS",
                                    "ESS = Effective Sample Size. Creates observation error in surveys at age. Cannot see any response to changes here, need to look in the Estimate Z tab Data plot or Fits plots to see the effect of changing this value",
                                    "right"),
                          
                          sliderInput("SurveyFullSelectivity",
                                      label = "Survey Full Selectivity",
                                      min = 1,
                                      max = nages,
                                      step = 1,
                                      value = c(1, 15)),
                          
                          bsTooltip("SurveyFullSelectivity",
                                    "Use sliders to determine ages that are and are not fully selected",
                                    "right"),
                          
                          sliderInput("SurveySelectivityWidth",
                                      label = "Survey Sel Width",
                                      min = 0.1,
                                      max = 5,
                                      step = 0.1,
                                      value = 4.0),
                          
                          bsTooltip("SurveySelectivityWidth",
                                    "Determines how steep logistic curve is",
                                    "right"),
                          
                          sliderInput("sigmaF",
                                      "sigmaF",
                                      min = 0.0,
                                      max = 0.9,
                                      step = 0.1,
                                      value = 0.0),
                          
                          bsTooltip("sigmaF",
                                    "Standard deviation of lognormal distribution to add noise to annual F values",
                                    "right"),
                          
                          sliderInput("phiF",
                                      "Autocorrelation",
                                      min = -1.0,
                                      max = 1.0,
                                      step = 0.1,
                                      value = 0.0),
                          
                          bsTooltip("phiF",
                                    "Corellation in random deviates for F from one year to the next",
                                    "right")
                          )
                   )
               ),
                   
               mainPanel(
                  plotOutput("selectivityPlot"),
                  plotOutput("ZvaluesPlot"),
                  bsTooltip("selectivityPlot",
                            "Shows selectivity at age for the fishery and survey, note the oldest age is not a plus group",
                            "left"),
                  bsTooltip("ZvaluesPlot",
                            "Shows the true Z values at age over time as blue-black dots for the age range used in the Z estimation only, the red line is the average value over this age range (also shown in the Z estimates plot in the Estimate Z tab). Note the y-axis adjusts to the range of values and that M is always 0.2",
                            "left")
               )
            )
   ),
                 
   tabPanel("Estimate Z",
            sidebarLayout(
               sidebarPanel(
                  sliderInput("ageInput", 
                              label = "Age Range", 
                              min = 1, 
                              max = 30, 
                              value = c(4,10)),
                  
                  bsTooltip("ageInput",
                            "Move these sliders to determine which ages are used in the Sinclair Z estimation",
                            "right"),
                  
                  sliderInput("nwSelection",
                              label = "Number of Years in Moving Window",
                              min = 3,
                              max = 6,
                              step = 1,
                              value = 4),
                  
                  bsTooltip("nwSelection",
                            "Width of moving window used in Sinclair Z estimation. Default value used in publication is 4",
                            "right"),
                  
                  uiOutput("fitslider"),
                  
                  bsTooltip("fitslider",
                            "Determines which regression is shown in the Fits tab. Can use the play button to loop through all the fits. Can also select specific fit by moving slider as ususal.",
                            "right"),
                  
                  br(),
                  br(),
                  h4("Reference"),
                  p("Sinclair, A.F. 2001. Natural mortality of cod (Gadus morhua) in the Southern Gulf of St Lawrence. ICES Journal of Marine Science. 58: 1-10. https://doi.org/10.1006/jmsc.1999.0490")
                  
               ),
               
               # Main panel (plots)
               mainPanel(
                  tabsetPanel(
                     tabPanel("Data", plotOutput("dataPlot")),
                     tabPanel("Z estimates", plotOutput("ZPlot")),
                     tabPanel("Diagnostic", plotOutput("diagnosticPlot")),
                     tabPanel("Fits", plotOutput("fitsPlot")),
                     bsTooltip("dataPlot",
                               "Survey catch at age and year for the ages selected",
                               "left"),
                     bsTooltip("ZPlot",
                               "Estimates of total mortality (points) with 90% confidence intervals (vertical lines) and true value from Data Creator (red line)",
                               "left"),
                     bsTooltip("diagnosticPlot",
                               "Box and whisker plot of residuals from all regressions (shown in Fits plots) for each age. Numbers denote the number of residuals at age (can differ due to zero catch at age for some cohorts)",
                               "left"),
                     bsTooltip("fitsPlot",
                               "One moving window regression with common total mortality estimate (negative of the slope) and different intercepts for each cohort. Colors and symbols denote different cohorts). The y-axis for these plots is the natural logarithm of the survey catch at age for a year in the moving window. Use the buttons to cycle through different year windows",
                               "left")
                     
                     
                  )
               )
        )
   )
)


#-------------------------------------------------------------------------------------
# Define server logic required to draw a histogram
server <- function(input, output) {

  # selectivity is based on logistic function Sel = 1 / (1 + exp(-slope*(age-A50)))
  # slope is determined by the  inter-quartile width (distance between 25% and 75% selectivity)
  # boundary of full selectivity used to solve for Sel = 0.999 and then set equal 1.0 for the slider age
  # logistic down at older ages uses Sel = 1 - (1 / (1 + exp(-slope*(age-A50))))
  # when sliders overlap, that age has full selectivity 
  # (cannot move increasing logistic to right of decreasing logistic)
  Selectivitydf <- reactive({
    fleetSel <- rep(1.0, nages)
    FS1 <- input$FisheryFullSelectivity[1]
    FS2 <- input$FisheryFullSelectivity[2]
    if (FS1 > 1) {
      selSlope <- 2.0 * log(3) / input$FisherySelectivityWidth
      selA50 <- (log((1.0 / 0.999) - 1.0) / selSlope) + FS1
      for (i in 1:(FS1-1)){
        fleetSel[i] <- 1.0 / (1.0 + exp(-selSlope * (i - selA50)))
      }
    } 
    if (FS2 < nages){
      selSlope <- 2.0 * log(3.0) / input$FisherySelectivityWidth
      selA50 <- (log((1.0 / 0.001) - 1.0) / selSlope) + FS2
      for (i in (FS2 + 1):nages){
        fleetSel[i] <- 1.0 - (1.0 / (1.0 + exp(-selSlope * (i - selA50))))
      }
    }
    surveySel <- rep(1.0, nages)
    SS1 <- input$SurveyFullSelectivity[1]
    SS2 <- input$SurveyFullSelectivity[2]
    if (SS1 > 1){
      survSlope <- 2.0 * log(3) / input$SurveySelectivityWidth
      survA50 <- (log((1.0 / 0.999) - 1.0) / survSlope) + SS1
      for (i in 1:(SS1-1)){
        surveySel[i] <- 1.0 / (1.0 + exp(-survSlope * (i - survA50)))
      }
    }
    if (SS2 < nages){
      survSlope <- 2.0 * log(3.0) / input$SurveySelectivityWidth
      survA50 <- (log((1.0 / 0.001) - 1.0) / survSlope) + SS2
      for (i in (SS2+1):nages){
        surveySel[i] <- 1.0 - (1.0 / (1.0 + exp(-survSlope * (i - survA50))))
      }
    }
    data.frame(Age = c(1:nages, 1:nages),
                        Source = c(rep("Fishery", nages), rep("Survey", nages)),
                        Selectivity = c(fleetSel, surveySel))
  })
  
  ZAA <- reactive({
    set.seed(input$rngSeed)
    Fdev <- rep(NA, nyears)
    FdevAR1 <- rep(NA, nyears)
    Fdev <- rnorm(nyears, mean = 0.0, sd = input$sigmaF)
    FdevAR1[1] <- Fdev[1]
    for (i in 2:nyears){
      FdevAR1[i] <- input$phiF * FdevAR1[i-1] + Fdev[i] 
    }
    Fmult <- rep(input$StartingF, nyears)
    if (input$EndingF != Fmult[1]){
      if (input$changeF == "step"){
        Fmult[ceiling(nyears / 2):nyears] <- input$EndingF
      }
      if (input$changeF == "linear"){
        Fmult <- seq(input$StartingF, input$EndingF, length.out = nyears)
      }
    }
    FmultNoisy <- Fmult * exp(FdevAR1)
    seldf <- Selectivitydf()
    fleetSel <- filter(seldf, Source == "Fishery") 
    FAA <- matrix(NA, nrow = nyears, ncol = nages)
    for (i in 1:nyears){
      for (j in 1:nages){
        FAA[i,j] <- FmultNoisy[i] * fleetSel$Selectivity[j]
      }
    }
    MAA <- matrix(0.2, nrow = nyears, ncol = nages)
    MAA + FAA
  })
  
  dat <- reactive({
    ZAA_use <- ZAA()
    sigmaR <- 0.6
    NAA <- matrix(NA,  nrow = nyears, ncol = nages)
    NAA[,1] <- exp(rnorm(nyears, mean = log(10000), sd = sigmaR))
    for (j in 2:nages){
      NAA[1,j] <- NAA[1,(j-1)] * exp(-ZAA_use[1,(j-1)])
    }
    for (i in 2:nyears){
      for(j in 2:nages){
        NAA[i,j] <- NAA[(i-1),(j-1)] * exp(-ZAA_use[(i-1),(j-1)])
      }
    }
    # calculate survey catch at age
    seldf <- Selectivitydf()
    surveySel <- filter(seldf, Source == "Survey") 
    surveyq <- 0.1
    surveyCAA <- matrix(NA, nrow = nyears, ncol = nages)
    for (i in 1:nyears){
      for (j in 1:nages){
        surveyCAA[i,j] <- NAA[i,j] * surveySel$Selectivity[j] * surveyq
      }
    }
    
    # apply observation error to survey catch at age using multinomial distribution and input ESS
    surveyObs <- matrix(NA, nrow = nyears, ncol = nages)
    surveyProp <- rep(NA, nages)
    for (i in 1:nyears){
      surveyProp <- surveyCAA[i,] / sum(surveyCAA[i,]) 
      surveyObs[i,] <- rmultinom(1, input$surveyESS, surveyProp)
    }
    
    # create data frame for survey data
    surveyDat <- data.frame(YEAR      = integer(),
                            AGE       = integer(),
                            NO_AT_AGE = double() )
    for (i in 1:nyears){
      for (j in 1:nages){
        thisDat <- data.frame(YEAR = i,
                              AGE  = j,
                              NO_AT_AGE = surveyObs[i,j])
        surveyDat <- rbind(surveyDat,thisDat)
      }
    }
    
    filter(surveyDat, 
           AGE %in% seq(input$ageInput[1], input$ageInput[2])) 
  })
  
  res <- reactive({
    mat <- select(dat(), c("YEAR", "AGE", "NO_AT_AGE")) %>%
      spread(key=AGE, value=NO_AT_AGE, fill=0)
    years <- mat$YEAR
    mat <- select(mat, -c(YEAR))
    rownames(mat) <- years
    calc_Sinclair_Z(mat, input$nwSelection) 
  })
  
  output$ZvaluesPlot <- renderPlot({
    ZAA_use <- ZAA()
    ZAAdf <- data.frame(YEAR = integer(),
                        AGE  = integer(),
                        Z    = double())
    for (i in 1:nyears){
      for (j in 1:nages){
        thisZAA <- data.frame(YEAR = i,
                              AGE  = j,
                              Z    = ZAA_use[i,j])
        ZAAdf <- rbind(ZAAdf, thisZAA)
      }
    }
    ZAAdf <- filter(ZAAdf, AGE %in% seq(input$ageInput[1], input$ageInput[2]))
    Ztrue <- apply(ZAA_use[,input$ageInput[1]: input$ageInput[2]], 1, mean)
    Ztruedf <- data.frame(Year = 1:nyears,
                          Z = Ztrue)
    ggplot(ZAAdf, aes(x=YEAR, y=Z, color=AGE)) +
      geom_point() +
      geom_line(data=Ztruedf, aes(x=Year, y=Z), color="tomato") + 
      theme_bw()
  })
  
  output$selectivityPlot <- renderPlot({
    seldf <- Selectivitydf()
    ggplot(seldf, aes(x=Age, y=Selectivity, color=Source)) +
      geom_point() +
      geom_line() +
      expand_limits(y = 0) +
      theme_bw()
  })
  
  output$dataPlot <- renderPlot({
    ggplot(filter(dat(), NO_AT_AGE > 0), aes(x=AGE, y=YEAR, size=NO_AT_AGE)) +
      geom_point() +
      scale_y_reverse() +
      theme_bw()
  })
  
  output$ZPlot <- renderPlot({
    # add Ztrue
    ZAA_use <- ZAA()
    Ztrue <- apply(ZAA_use[,input$ageInput[1]: input$ageInput[2]], 1, mean)
    Ztruedf <- data.frame(Year = 1:nyears,
                          Z = Ztrue)
    # end add Ztrue section
    est.Z <- data.frame(Year = res()$plot.year,
                        Z = res()$est.Sinclair.Z[,1],
                        low = res()$est.Sinclair.Z[,2],
                        high = res()$est.Sinclair.Z[,3])
    ggplot(est.Z, aes(x=Year, y=Z)) +
      geom_errorbar(aes(ymin=low, ymax=high), na.rm = TRUE) +
      geom_point(na.rm = TRUE) +
      geom_line(data=Ztruedf, aes(x=Year, y=Z), color="tomato") + # add Ztrue line
      geom_hline(yintercept = 0) +
      theme_bw()
  })
  
  output$diagnosticPlot <- renderPlot({
    ny <- length(res()$plot.year)
    if (!all(is.na(res()$est.Sinclair.Z))){
      resids <- data.frame(Age = numeric(),
                           Resid = numeric())
      for (i in 1:ny){
        if (!is.na(res()$est.Sinclair.Z[i,1])){
          data <- res()[[i]]
          resids <- rbind(resids,data.frame(Age = data$Age,
                                            Resid = data$resid))
        }  
      }
      ggplot(resids, aes(x=as.factor(Age), y=Resid)) +
        geom_boxplot() +
        geom_hline(yintercept = 0, color = "tomato") +
        stat_summary_bin(fun.data=give.n, geom="text", fun.y=median) +
        xlab("Age") +
        theme_bw()
    }
  })
  
  output$fitslider <- renderUI({
    sliderInput("whichfit",
                label = "Which fit to view",
                min = 1,
                max = length(res()$plot.year),
                step = 1,
                value = 1,
                animate = TRUE)
  })
  
  output$fitsPlot <- renderPlot({
    ny <- length(res()$plot.year)
    my.col <- 1:100
    my.col[7] <- "blue"
    data <- res()[[input$whichfit]]
    n.coh <- length(unique(data$cohort))
    i.coh <- (data$cohort-min(data$cohort, na.rm=T)+1)
    plot(data$Age,data$lnVal,pch=i.coh,col=my.col[i.coh],xlab="Age",ylab="ln Val")
    title(main=paste0("Years ",min(data$Year)," to ",max(data$Year),
                      "\nZ = ",round(res()$est.Sinclair.Z[input$whichfit,1],3), 
                      "  (",round(res()$est.Sinclair.Z[input$whichfit,2],3),", ",
                      round(res()$est.Sinclair.Z[input$whichfit,3],3),")"), outer=F)
    for (j in 1:n.coh){
      subcoh <- data[i.coh == j,]
      if (length(subcoh$lnVal) >= 2){
        lines(subcoh$Age,subcoh$pred,col=my.col[j],lty=j) 
      }
    }
  })
  
}

#-------------------------------------------------------------------------------------
# Run the application 
shinyApp(ui = ui, server = server)

