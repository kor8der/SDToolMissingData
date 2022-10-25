#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library (ggplot2)
library(plotly)


# Define UI for application that draws a histogram
ui <- fluidPage(
      
    # Application title
    titlePanel("Incomplete data signal detection"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput("TruePositivesNum", "Define True Positive observations:", 0, min = 0, max = 1000000, step = 1),
            numericInput("FalsePositivesNum", "Define False Positive observations:", 0, min = 0, max = 1000000, step = 1),
            numericInput("FalseNegativesNum", "Define False Negative observations:", 0, min = 0, max = 1000000, step = 1),
            numericInput("TrueNegativesNum", "Define True Negative observations:", 0, min = 0, max = 1000000, step = 1),
            numericInput("PositiveRate", "Define positive percentage:", 0, min = 0, max = 100, step = 0.001),
            numericInput("DesiredPrecision", "Define Desired Precision (note that changing this may alter performance by orders of magnitude):", 0.001, min = 0.00000000001, max = 1, step = 0.00001),
            numericInput("dMinimum", "Define minimum d'", 0.5, min = 0.001, max = 15, step = 0.001),
            numericInput("dStepSize", "Define d' stepsize", 0.5, min = 0.001, max = 15, step = 0.001),
            numericInput("dNumber", "Define number of d'", 3, min = 1, max = 14, step = 1),
            sliderInput("brRange", "Define br range:", c(0.2,0.8), min = 0.01, max = 0.999, step = 0.001),
            selectInput("SelectedOutput", "Pick output:", c("payoff ratio", "bias", "result", "distance")),
            numericInput("yMin", "(optional) define minimum displayed y:", NA, min = 0, max = 1000000, step = 1),
            numericInput("yMax", "(optional) define maximum displayed y:", NA, min = 0, max = 1000000, step = 1),
            numericInput("yAnchor", "(optional) define y reference:", NA, min = 0, max = 1000000, step = 1),
            textInput("yAnchorName", "(optional) define name for y reference:", ""),
            textInput("yTitle", "(optional) define y-axis title:", ""),
            numericInput("xMin", "(optional) define minimum displayed x:", NA, min = 0, max = 1000000, step = 1),
            numericInput("xMax", "(optional) define maximum displayed x:", NA, min = 0, max = 1000000, step = 1),
            numericInput("xAnchor", "(optional) define x reference:", NA, min = 0, max = 1000000, step = 1),
            textInput("xAnchorName", "(optional) define name x reference:", ""),
            textInput("xTitle", "(optional) define x-axis title:", ""),
            numericInput("dAnchor", "(optional) define d' reference:", NA, min = 0, max = 1000000, step = 1),
            numericInput("dMin", "(optional) define d' error margin minimum.", NA, min = 0, max = 1000000, step = 1),
            numericInput("dMax", "(optional) define d' error margin maximum:", NA, min = 0, max = 1000000, step = 1),
            textInput("graphTitle", "(optional) define graph title:", ""),
            width = 3
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          h3(textOutput("chosenAnalysis")),
            plotlyOutput("distPlot"),
            width = 9
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    ChosenAnalysis <- reactiveVal(0)
    TP <- reactive({input$TruePositivesNum})
    FP <- reactive({input$FalsePositivesNum})
    FN <- reactive({input$FalseNegativesNum})
    TN <- reactive({input$TrueNegativesNum})
    PRate <- reactive({input$PositiveRate})
    dPrimes <- reactive({seq(from = input$dMinimum, length.out = input$dNumber, by = input$dStepSize)})
    dPrimeErrors <- reactive(c(input$dMin, input$dMax))
    dPrimeReference <- reactive({input$dAnchor})
    baseRates <- reactive({seq(input$brRange[1],input$brRange[2], 0.001)})
    desiredPrecision <- reactive({input$DesiredPrecision})
    selectedOutputValue <-reactive({input$SelectedOutput})
    yAxisMin <- reactive({input$yMin})
    yAxisMax <- reactive({input$yMax})
    xAxisMin <- reactive({input$xMin})
    xAxisMax <- reactive({input$xMax})
    yAnchorVal <- reactive({input$yAnchor})
    xAnchorVal <- reactive({input$xAnchor})
    yAnchorName <- reactive({input$yAnchorName})
    xAnchorName <- reactive({input$xAnchorName})
    customGraphTitle <- reactive({input$graphTitle})
    xAxisTitle <- reactive({input$xTitle})
    yAxisTitle <- reactive({input$yTitle})
    
    font <- list(family = "Times New Roman",size = 24,color = "#7f7f7f")
    initialStepsize <- 0.1
    
    output$distPlot <- renderPlotly({
      resultMatrixes <- list()
      resultMatrixCount <- 1
      newbr = baseRates()
      realbrLength <- length(newbr)
      resultNamesList <- list()
      realdPrimeVector <- dPrimes()
      
      if(!is.na(dPrimeReference()) && dPrimeReference() > 0){
        realdPrimeVector <- append(realdPrimeVector, dPrimeReference())  
      }
      
      if(length(dPrimeErrors()) == 2&&!is.na(dPrimeErrors()[1]) && !is.na(dPrimeErrors()[2]) && dPrimeErrors()[1] > 0 && dPrimeErrors()[2] > 0){
        realdPrimeVector <- append(realdPrimeVector, dPrimeErrors())  
      }
      
        if(1%in%ChosenAnalysis()){
            AnalysisName <- "FP/TP"
            Goal <- FP()/TP()
            for (i in realdPrimeVector){
              imat <- matrix(nrow = 0, ncol = 10)
                for (j in newbr){
                    if((1-j)/j < Goal){
                        SKIPPED <- TRUE
                    }
                  else{
                    stepSize <- initialStepsize
                    optimalBias <- 0
                    nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                    nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                    FPdivTP <- nFalsePositives / nTruePositives
                    
                    distance <- Goal - FPdivTP
                    
                    if(distance > 0){
                        MovingUp <- TRUE
                    }
                    else{
                        MovingUp <- FALSE
                    }
                    Halved <- FALSE
                    
                    while (abs(distance) > desiredPrecision()){
                        
                        SKIPPED <- FALSE
                        if(distance > 0){
                            if (MovingUp == FALSE || Halved == TRUE){
                                stepSize <- stepSize/2
                                Halved = TRUE
                            }
                            optimalBias <- optimalBias - stepSize 
                            MovingUp <- TRUE
                        }
                        else{
                            if (MovingUp == TRUE || Halved == TRUE){
                                stepSize <- stepSize/2
                                Halved = TRUE
                            }
                            optimalBias <- optimalBias + stepSize  
                            MovingUp <- FALSE
                        }
                        nFalsePositives <- (1-pnorm((optimalBias+i/2))) * ((1-j))
                        nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                        FPdivTP <- nFalsePositives / nTruePositives
                        
                        distance <- Goal - FPdivTP
                        
                        if(nTruePositives == 0  || stepSize < 0.000000000001){
                          SKIPPED <- TRUE
                          break 
                        }
                      }  
                    }
                    if (SKIPPED == FALSE){
                        FNPayOff <- log10(1/(exp(optimalBias*i)*j/(1-j)))
                        nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                        nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                        
                        imat <- rbind(imat, c(i, j, FNPayOff, optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FPdivTP, distance))
                    }
                    else{
                      
                      imat <- rbind(imat, c(i, j, NA, NA, NA, NA, NA, NA, NA, NA)) ##10 colums
                    }
                }
                
              resultMatrixes[[resultMatrixCount]] <- imat
              
              resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
              resultMatrixCount = resultMatrixCount + 1
            }
        }   
        if(2%in%ChosenAnalysis()){
            AnalysisName <- "TP/FN"
            for (i in realdPrimeVector)
            {
                imat <- matrix(nrow = 0, ncol = 10)
                Goal <- FN()/(TP()+FN())
                optimalBias = qnorm(Goal) + i/2
                
                for (j in newbr){
                    nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                    nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                    nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                    nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                    FNPayOff <- log10(1/(exp(optimalBias*i)*j/(1-j)))
                    FNdivTP <- nFalseNegatives / nTruePositives
                    distance <- FNdivTP - FP()/TN()
                    imat <- rbind(imat, c(i, j, FNPayOff, optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FNdivTP, distance))
                }
                
                resultMatrixes[[resultMatrixCount]] <- imat
                
                resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
                
                resultMatrixCount = resultMatrixCount + 1
            }
        }
        if(3%in%ChosenAnalysis()){
            AnalysisName <- "TN/TP"
            Goal <- TN()/TP()
            for (i in realdPrimeVector)
            {
              imat <- matrix(nrow = 0, ncol = 10)
              for (j in baseRates()){
                    stepSize <- initialStepsize
                    optimalBias <- 0
                    nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                    nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                    TNdivTP <- nTrueNegatives / nTruePositives
                    distance <- Goal - TNdivTP
                    
                    if(distance < 0){
                        MovingUp <- TRUE
                    }
                    else{
                        MovingUp <- FALSE
                    }
                    
                    Halved <- FALSE
                    
                    while (abs(distance) > desiredPrecision()){
                        if(distance < 0)
                        {
                            if (MovingUp == FALSE || Halved == TRUE)
                            {
                                stepSize <- stepSize/2
                                Halved = TRUE
                            }
                            
                            optimalBias <- optimalBias - stepSize 
                            MovingUp <- TRUE
                        }
                        else
                        {
                            if (MovingUp == TRUE || Halved == TRUE)
                            {
                                stepSize <- stepSize/2
                                Halved = TRUE
                            }
                            
                            optimalBias <- optimalBias + stepSize  
                            MovingUp <- FALSE
                        }
                        
                        nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                        nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                        TNdivTP <- nTrueNegatives / nTruePositives
                        distance <- Goal - TNdivTP
                        
                    }
                    
                    nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                    nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                    FNPayOff <- log10(1/(exp(optimalBias*i)*j/(1-j)))
                    
                    imat <- rbind(imat, c(i, j, FNPayOff, optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, TNdivTP, distance))
              }
              resultMatrixes[[resultMatrixCount]] <- imat
              
              resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
              
              resultMatrixCount = resultMatrixCount + 1
            }
            
       }
        if(4%in%ChosenAnalysis()){
            AnalysisName <- "FP/FN"
            Goal <- FP()/FN()
            for (i in realdPrimeVector)
            {
              imat <- matrix(nrow = 0, ncol = 10)
                for (j in baseRates()){
                    stepSize <- initialStepsize
                    optimalBias <- 0
                    nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                    nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                    FPdivFN <- nFalsePositives / nFalseNegatives
                    
                    distance <- Goal - FPdivFN
                    
                    if(distance > 0){
                        MovingUp <- TRUE
                    }
                    else{
                        MovingUp <- FALSE
                    }
                    
                    Halved <- FALSE
                    
                    while (abs(distance) > desiredPrecision()){
                        if(distance > 0)
                        {
                            if (MovingUp == FALSE || Halved == TRUE)
                            {
                                stepSize <- stepSize/2
                                Halved = TRUE
                            }
                            optimalBias <- optimalBias - stepSize 
                            MovingUp <- TRUE
                        }
                        else
                        {
                            if (MovingUp == TRUE || Halved == TRUE)
                            {
                                stepSize <- stepSize/2
                                Halved = TRUE
                            }
                            optimalBias <- optimalBias + stepSize  
                            MovingUp <- FALSE
                        }
                        
                        nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                        nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                        FPdivFN <- nFalsePositives / nFalseNegatives
                        distance <- Goal - FPdivFN
                    }
                    
                    nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                    nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                    FNPayOff <- log10(1/(exp(optimalBias*i)*j/(1-j)))
                    
                    imat <- rbind(imat, c(i, j, FNPayOff, optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FPdivFN, distance))
                }
              resultMatrixes[[resultMatrixCount]] <- imat
              
              resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
              
              resultMatrixCount = resultMatrixCount + 1
            }
        }
        if(5%in%ChosenAnalysis()){
          AnalysisName <- "TN/FP"
          Goal <- FP()/(TN()+FP())
            for (i in realdPrimeVector)
            {
                optimalBias = qnorm(Goal) - i/2
                imat <- matrix(nrow = 0, ncol = 10)
                
                for (j in newbr){
                    nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                    nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                    nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                    nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                    FNPayOff <- log10(1/(exp(optimalBias*i)*j/(1-j)))
                    FPdivTN <- nFalsePositives / nTrueNegatives
                    distance <- FPdivTN - FP()/TN()
                    imat <- rbind(imat, c(i, j, FNPayOff, optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FPdivTN, distance))
                }
                
                resultMatrixes[[resultMatrixCount]] <- imat
                
                resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
                
                resultMatrixCount = resultMatrixCount + 1
            }
        }
        if(6%in%ChosenAnalysis()){
          AnalysisName <- "FN/TN"
            Goal <- FN()/TN()
            for (i in realdPrimeVector)
            {
              imat <- matrix(nrow = 0, ncol = 10)
              for (j in newbr){
                  SKIPPED <- FALSE
                  if((1-j)/j > Goal)
                  {
                    SKIPPED <- TRUE
                  }else{
                  
                    stepSize <- initialStepsize
                    optimalBias <- 0
                    nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                    nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                    FNdivTN <- nFalseNegatives / nTrueNegatives
                    distance <- Goal - FNdivTN
                    if(distance < 0){
                        MovingUp <- TRUE
                    }
                    else{
                        MovingUp <- FALSE
                    }
                    Halved <- FALSE
                    
                    while (abs(distance) > desiredPrecision()){
                        if(distance < 0)
                        {
                            if (MovingUp == FALSE || Halved == TRUE)
                            {
                                stepSize <- stepSize/2
                                Halved = TRUE
                            }
                            optimalBias <- optimalBias - stepSize 
                            MovingUp <- TRUE
                        } 
                        else
                        {
                            if (MovingUp == TRUE || Halved == TRUE)
                            {
                                stepSize <- stepSize/2
                                Halved = TRUE
                            }
                            optimalBias <- optimalBias + stepSize  
                            MovingUp <- FALSE
                        }
                        
                        nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                        nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                        FNdivTN <- nFalseNegatives / nTrueNegatives
                        distance <- Goal - FNdivTN
                        if(nTruePositives == 0  || stepSize < 0.000000000001)
                        {
                            SKIPPED <- TRUE
                            break 
                        }
                    }
                    
                  }
                if (SKIPPED == FALSE){
                  nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                  nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                  FNPayOff <- log10(1/(exp(optimalBias*i)*j/(1-j)))
                  imat <- rbind(imat, c(i, j, FNPayOff, optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FNdivTN, distance))
                }
                else{
                  imat <- rbind(imat, c(i, j, NA, NA, NA, NA, NA, NA, NA, NA))
                }
              }
              
              resultMatrixes[[resultMatrixCount]] <- imat
              
              resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
              
              resultMatrixCount = resultMatrixCount + 1
          }
        }
        if(7%in%ChosenAnalysis()){
          AnalysisName <- "PositiveRate"
          Goal <- PRate()
          for (i in realdPrimeVector){
            imat <- matrix(nrow = 0, ncol = 10)
            for (j in newbr){
              print(i)
              print(j)
              SKIPPED <- FALSE
              if(0 > Goal){
                SKIPPED <- TRUE
              }
              else{
                
                stepSize <- initialStepsize
                optimalBias <- 0
                nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                PosRate <- (nFalsePositives+nTruePositives) / (nFalsePositives+nTruePositives+nFalseNegatives+nTruePositives)
                distance <- Goal - PosRate
                if(distance > 0){
                  MovingUp <- TRUE
                }
                else{
                  MovingUp <- FALSE
                }
                Halved <- FALSE
                
                while (abs(distance) > desiredPrecision()){
                  if(distance > 0){
                    if (MovingUp == FALSE || Halved == TRUE){
                      stepSize <- stepSize/2
                      Halved = TRUE
                    }
                    optimalBias <- optimalBias - stepSize 
                    MovingUp <- TRUE
                  } 
                  else{
                    if (MovingUp == TRUE || Halved == TRUE){
                      stepSize <- stepSize/2
                      Halved = TRUE
                    }
                    optimalBias <- optimalBias + stepSize  
                    MovingUp <- FALSE
                  }
                  
                  nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                  nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                  nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                  nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                  PosRate <- (nFalsePositives+nTruePositives) / (nFalsePositives+nTruePositives+nFalseNegatives+nTruePositives)
                  distance <- Goal - PosRate
                  if(stepSize < 0.000000000001){
                    SKIPPED <- TRUE
                    break 
                  }
                }
                
              }
              if (SKIPPED == FALSE){
                FNPayOff <- log10(1/(exp(optimalBias*i)*j/(1-j)))
                imat <- rbind(imat, c(i, j, FNPayOff, optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, PosRate, distance))
              }
              else{
                imat <- rbind(imat, c(i, j, NA, NA, NA, NA, NA, NA, NA, NA))
              }
            }
            
            resultMatrixes[[resultMatrixCount]] <- imat
            
            resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
            
            resultMatrixCount = resultMatrixCount + 1
          }
        }
      ## Restriction area 
      if (length(resultNamesList) > 0){
        endV <- vector()     
        startV <- vector()      
        for (n in 1:length(resultMatrixes)){
          test = matrix(resultMatrixes[[n]],nrow = length(newbr), ncol = 10)
  
          startN <- 1
          endN <- length(newbr)
  
          count <- 0
          reachedStart <- FALSE
          reachedEnd <- FALSE
          
          for (m in 1:nrow(test)){
            if (reachedStart == FALSE){
              if (is.na(test[m,3])){
                startN <- startN+1
              }
              else{
                reachedStart <- TRUE
              }
            }
            else if (reachedEnd == FALSE){
              if (is.na(test[m,3])){
                reachedEnd  <- TRUE
                endN <- count
              }
            }
            count <- count+1
          }
          startV <- append (startV, startN)
          endV <- append (endV, endN)
        }
        startV
        endV
        
        realStart <- max(startV)
        realEnd <- min(endV)
        
        if (selectedOutputValue() == "bias"){
          selectedColumn <- 4    
        }
        else if (selectedOutputValue() == "payoff ratio") {
          selectedColumn <- 3    
        }
        else if (selectedOutputValue() == "result") {
          selectedColumn <- 9
        }
        else if (selectedOutputValue() == "distance") {
          selectedColumn <- 10    
        }
        
        displayMatrix <- matrix(resultMatrixes[[1]],nrow = length(newbr), ncol = 10)
        
        
        data <- data.frame(displayMatrix[realStart:realEnd,2], displayMatrix[realStart:realEnd,selectedColumn])
        k <- 2
        
        while (k <= length(resultMatrixes)){
          displayMatrix <- matrix(resultMatrixes[[k]],nrow = length(newbr), ncol = 10)
          newdata <- data.frame(displayMatrix[realStart:realEnd,selectedColumn])
          data <- cbind(data, newdata)
          k <- k+1
        }
        
        xTitle <- "Base Rate"
        yTitle <- selectedOutputValue()
        graphTitle <- "Base Rate"
        
        if (nchar(xAxisTitle()) > 0)
        {
          xTitle<- xAxisTitle()
        }
        if (nchar(yAxisTitle()) > 0)
        {
          yTitle<- yAxisTitle()
        }
        if (nchar(customGraphTitle()) > 0)
        {
          graphTitle<- customGraphTitle()
        }
        
        x <- list(title = xTitle,titlefont = font)
        
        y <- list(title = yTitle,titlefont = font)
        
        theColors <- getColors(length(dPrimes()))
        
        fig <- plot_ly(data = data, type ='scatter', mode = 'lines') %>% layout(title = graphTitle, xaxis = x, yaxis = y)
        
        if(!is.na(yAnchorVal())){
          fig <- fig %>% add_lines(y = yAnchorVal(), x = baseRates(), name = yAnchorName(),  line = list(color = "grey")) 
        }
        
        if(!is.na(xAnchorVal())){
          fig <- fig %>% add_lines(y = (c(max(data), min(data))), x = xAnchorVal(), name = xAnchorName(),  line = list(color = "grey")) 
        }
        
        if(length(dPrimeErrors()) == 2&&!is.na(dPrimeErrors()[1]) && !is.na(dPrimeErrors()[2]) && dPrimeErrors()[1] > 0 && dPrimeErrors()[2] > 0){
          lastNum <- length(realdPrimeVector)
          fig <- fig %>% add_polygons(x=c(data[,1], rev(data[,1])), y= c(data[,(ncol(data)-1)],rev(data[,(ncol(data))])), mode='lines', name="Error", color = I("dark grey"), opacity =0.8, line = list (width = 0))
        }
        
        
        
        dPrime <- 1
        roundsofDprime <- 1
        
        for (l in 2:length(dPrimes())){
          if (roundsofDprime == 1){
            style <- "line"
          } 
          else if(roundsofDprime == 2){
            style <- "dash"
          }
          else{
            style <- "dot"
          }
          fig <- fig %>% add_trace(x=data[,1], y= data[,l], mode='lines', name=resultNamesList[[l-1]], line = list(color= theColors[dPrime], dash = style))
          l <- l+1
          dPrime <- dPrime+1
          if(length(realdPrimeVector)<dPrime){
            dPrime <- 1
            roundsofDprime <- roundsofDprime + 1
          }
        }
        
        if (!is.na(dPrimeReference()) && dPrimeReference() > 0){
          fig <- fig %>% add_trace(x=data[,1], y= data[,(length(dPrimes())+2)], mode='lines', name=resultNamesList[[(length(dPrimes())+1)]], line = list(color= "000000", dash = "dot"))
        }
        
        
        fig
      }
    })
    
    output$chosenAnalysis <- renderText({
      def <- "No Chosen Analysis"
      ret <- ""
      vec <- vector()
        if (TP() > 0){
            if (FP() > 0){
                vec <- append(vec, 1)
                if (nchar(ret) > 1){
                  ret <- paste(ret, "TP/FP", sep = " & ")  
                }
                else{
                  ret <- "TP/FP"
                }
                
            }
            else if (FN()> 0){
              vec <- append(vec, 2)
                if (nchar(ret) > 1){
                  ret <- paste(ret, "TP/FN", sep = " & ")  
                }
                else{
                  ret <- "TP/FN"
                }
            }
            else if (TN() > 0){
              vec <- append(vec, 3)
                if (nchar(ret) > 1){
                  ret <- paste(ret, "TP/TN", sep = " & ")  
                }
                else{
                  ret <- "TP/TN"
                }
            }
        }
        if (FP()> 0){
            if (FN() > 0){
              vec <- append(vec, 4)
                if (nchar(ret) > 1){
                  ret <- paste(ret, "FP/FN", sep = " & ")
                }
                else{
                  ret <- "FP/FN"
                }
            }
            else if (TN() > 0){
              vec <- append(vec, 5)
                if (nchar(ret) > 1){
                  ret <- paste(ret, "FP/TN", sep = " & ")
                }
                else{
                  ret <- "FP/TN"
                }
            }
        }
        if (FN() > 0){
            if (TN() > 0){
              vec <- append(vec, 6)
                if (nchar(ret) > 1){
                  ret <- paste(ret, "FN/TN", sep = " & ") 
                }
                else{
                  ret <- "FN/TN"
                }
            }
        }
        if (PRate() >0){
          vec <- append(vec, 7)
          if (nchar(ret) > 1){
            ret <- paste(ret, "PositiveRate", sep = " & ") 
          }
          else{
            ret <- "PositiveRate"
          }
        }

      ChosenAnalysis(vec)
      if (nchar(ret) > 1){
        return(ret)
      }
      else {
        return(def)
      }
    })
}

getColors <- function(num){
  theColors <- c("CC0000","00CC00","0000CC", "606060", "FFFF00", "00FFFF", "FF00FF", "994C00", "00994C", "4C0099", "99FF33", "3399FF", "FF3399", "C0C0C0", "000000")
  if(num < 16){
    return(theColors[1:num])
  }
  else{
    return(c(theColors,theColors[1:num-15]))
  }
}

# Run the application 
shinyApp(ui = ui, server = server)
