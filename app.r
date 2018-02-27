# Sinclair Z Shiny Tester
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape)
library(Hmisc)

# dimensions 
nyears <- 40
nages <- 30
selectivityOption <- c("All Ages Fully Selected", "Asymptotic", "Domed")

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
ui <- navbarPage("Sinclair Z Shiny Tester",  #fluidPage(
                 
                 # Application title
                 #  titlePanel("Sinclair Z Shiny Tester"),
                 
                 # Sidebar  
                 tabPanel("Data Creator",
                          sidebarLayout(
                            sidebarPanel(
                              sliderInput("rngSeed",
                                          label = "RNG seed value",
                                          min = 1,
                                          max = 1000,
                                          step = 1,
                                          value = 14),
                              
                              selectInput("surveyESS",
                                          label = "survey ESS",
                                          choices = c(50, 100, 1000, 10000),
                                          selected = 1000),
                              
                              sliderInput("FisheryFullSelectivity",
                                          label = "Fishery Full Selectivity",
                                          min = 1,
                                          max = nages,
                                          step = 1,
                                          value = c(5, nages)),

                              sliderInput("FisherySelectivityWidth",
                                          label = "Fish Sel Inner-Quartile Width",
                                          min = 0.1,
                                          max = 5,
                                          step = 0.1,
                                          value = 2),

                              sliderInput("SurveyFullSelectivity",
                                          label = "Survey Full Selectivity",
                                          min = 1,
                                          max = nages,
                                          step = 1,
                                          value = c(2, nages)),
                              
                              sliderInput("SurveySelectivityWidth",
                                          label = "Survey Sel Inner-Quartile Width",
                                          min = 0.1,
                                          max = 5,
                                          step = 0.1,
                                          value = 1)
                            ),
                            
                            mainPanel(
                              plotOutput("creatorPlot"),
                              plotOutput("selectivityPlot")
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
                              
                              sliderInput("nwSelection",
                                          label = "Number of Years in Moving Window",
                                          min = 3,
                                          max = 6,
                                          step = 1,
                                          value = 4),
                              
                              actionButton("firstfit",
                                           label = "Display First Fit"),
                              
                              actionButton("previousfit",
                                           label = "Display Previous Fit"),
                              
                              actionButton("nextfit",
                                           label = "Display Next Fit"),
                              
                              actionButton("lastfit",
                                           label = "Display Last Fit"),
                              
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
                                tabPanel("Fits", plotOutput("fitsPlot"))
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
    Fmult <- rep(0.4, nyears)
    seldf <- Selectivitydf()
    fleetSel <- filter(seldf, Source == "Fishery") 
    FAA <- matrix(NA, nrow = nyears, ncol = nages)
    for (i in 1:nyears){
      for (j in 1:nages){
        FAA[i,j] <- Fmult[i] * fleetSel$Selectivity[j]
      }
    }
    MAA <- matrix(0.2, nrow = nyears, ncol = nages)
    MAA + FAA
  })
  
  dat <- reactive({
    set.seed(input$rngSeed)
    #dropvals <- rnorm(nyears * nages)  # randM
    #dropvals <- rnorm(nyears * nages)  # randF
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
  
  output$creatorPlot <- renderPlot({
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
    ggplot(ZAAdf, aes(x=YEAR, y=Z, color=AGE)) +
      geom_point() +
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
    ggplot(dat(), aes(x=AGE, y=YEAR, size=NO_AT_AGE)) +
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
    mat <- select(dat(), c("YEAR", "AGE", "NO_AT_AGE")) %>%
      spread(key=AGE, value=NO_AT_AGE, fill=0)
    years <- mat$YEAR
    mat <- select(mat, -c(YEAR))
    rownames(mat) <- years
    res <- calc_Sinclair_Z(mat, input$nwSelection) 
    est.Z <- data.frame(Year = res$plot.year,
                        Z = res$est.Sinclair.Z[,1],
                        low = res$est.Sinclair.Z[,2],
                        high = res$est.Sinclair.Z[,3])
    ggplot(est.Z, aes(x=Year, y=Z)) +
      geom_errorbar(aes(ymin=low, ymax=high), na.rm = TRUE) +
      geom_point(na.rm = TRUE) +
      geom_line(data=Ztruedf, aes(x=Year, y=Z), color="tomato") + # add Ztrue line
      geom_hline(yintercept = 0) +
      theme_bw()
  })
  
  output$diagnosticPlot <- renderPlot({
    mat <- select(dat(), c("YEAR", "AGE", "NO_AT_AGE")) %>%
      spread(key=AGE, value=NO_AT_AGE, fill=0)
    years <- mat$YEAR
    mat <- select(mat, -c(YEAR))
    rownames(mat) <- years
    res <- calc_Sinclair_Z(mat, input$nwSelection)
    ny <- length(res$plot.year)
    if (!all(is.na(res$est.Sinclair.Z))){
      resids <- data.frame(Age = numeric(),
                           Resid = numeric())
      for (i in 1:ny){
        if (!is.na(res$est.Sinclair.Z[i,1])){
          data <- res[[i]]
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
  
  icount <- reactiveValues(fit = 1)
  observeEvent(input$firstfit, {icount$fit <- 1})
  observeEvent(input$previousfit, {
    icount$fit <- icount$fit - 1
    icount$fit <- max(icount$fit, 1)
    icount$fit <- min(icount$fit, ny)
  })
  observeEvent(input$nextfit, {
    icount$fit <- icount$fit + 1
    icount$fit <- max(icount$fit, 1)
    icount$fit <- min(icount$fit, ny)
  })
  observeEvent(input$lastfit, {icount$fit <- ny})
  
  output$fitsPlot <- renderPlot({
    mat <- select(dat(), c("YEAR", "AGE", "NO_AT_AGE")) %>%
      spread(key=AGE, value=NO_AT_AGE, fill=0)
    years <- mat$YEAR
    mat <- select(mat, -c(YEAR))
    rownames(mat) <- years
    res <- calc_Sinclair_Z(mat, input$nwSelection)
    ny <<- length(res$plot.year)
    my.col <- 1:100
    my.col[7] <- "blue"
    data <- res[[icount$fit]]
    n.coh <- length(unique(data$cohort))
    i.coh <- (data$cohort-min(data$cohort, na.rm=T)+1)
    plot(data$Age,data$lnVal,pch=i.coh,col=my.col[i.coh],xlab="Age",ylab="ln Val")
    title(main=paste0("Years ",min(data$Year)," to ",max(data$Year),
                      "\nZ = ",round(res$est.Sinclair.Z[icount$fit,1],3), 
                      "  (",round(res$est.Sinclair.Z[icount$fit,2],3),", ",
                      round(res$est.Sinclair.Z[icount$fit,3],3),")"), outer=F)
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

