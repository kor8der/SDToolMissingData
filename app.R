##Before running this code for the first time, run the commented command below:
#install.packages("shiny", "ggplot2", "plotly", "rsconnect")

library(shiny)
library (ggplot2)
library(plotly)
library(rsconnect)


# Define UI
ui <- fluidPage(
      
    # Application
    titlePanel("Incomplete data signal detection"),
    
    # Sidebar
    sidebarLayout(
        sidebarPanel(
            numericInput("TruePositivesNum", "Define True Positive observations:", 0, min = 0, max = 1000000, step = 1),
            numericInput("FalsePositivesNum", "Define False Positive observations:", 0, min = 0, max = 1000000, step = 1),
            numericInput("FalseNegativesNum", "Define False Negative observations:", 0, min = 0, max = 1000000, step = 1),
            numericInput("TrueNegativesNum", "Define True Negative observations:", 0, min = 0, max = 1000000, step = 1),
            numericInput("PositiveRate", "Define positive percentage:", 0, min = 0, max = 0.999, step = 0.001),
            numericInput("DesiredPrecision", "Define Desired Precision (note that changing this may alter performance by orders of magnitude):", 0.0001, min = 0.00000000001, max = 1, step = 0.00001),
            numericInput("dMinimum", "Define minimum d'", 0.5, min = 0.001, max = 15, step = 0.001),
            numericInput("dStepSize", "Define d' stepsize", 0.5, min = 0.001, max = 15, step = 0.001),
            numericInput("dNumber", "Define number of d' to chart", 3, min = 1, max = 14, step = 1),
            sliderInput("brRange", "Define br range:", c(0.2,0.8), min = 0.01, max = 0.999, step = 0.001),
            selectInput("SelectedOutput", "Pick output:", c("payoff ratio","log 10 of payoff ratio", "bias", "result", "distance"), selected = "log 10 of payoff ratio"),
            textInput("graphTitle", "(optional) define graph title:", ""),
            textInput("yTitle", "(optional) define y-axis title:", ""),
            numericInput("yMin", "(optional) define minimum displayed y:", NA, min = 0, max = 1000000, step = 1),
            numericInput("yMax", "(optional) define maximum displayed y:", NA, min = 0, max = 1000000, step = 1),
            textInput("xTitle", "(optional) define x-axis title:", ""),
            numericInput("xMin", "(optional) define minimum displayed x:", NA, min = 0, max = 1000000, step = 1),
            numericInput("xMax", "(optional) define maximum displayed x:", NA, min = 0, max = 1000000, step = 1),
            numericInput("yAnchor", "(optional) define y reference:", NA, min = 0, max = 1000000, step = 1),
            textInput("yAnchorName", "(optional) define name for y reference:", ""),
            numericInput("xAnchor", "(optional) define x reference:", NA, min = 0, max = 1000000, step = 1),
            textInput("xAnchorName", "(optional) define name x reference:", ""),
            numericInput("dAnchor", "(optional) define d' reference:", NA, min = 0, max = 1000000, step = 1),
            numericInput("dMin", "(optional) define d' error margin minimum.", NA, min = 0, max = 1000000, step = 1),
            numericInput("dMax", "(optional) define d' error margin maximum:", NA, min = 0, max = 1000000, step = 1),
            width = 3
        ),
        
        # Show the header for the chosen analysis, and the plot
        mainPanel(
          h3(textOutput("chosenAnalysis")),
            plotlyOutput("distPlot"),
            width = 9
        )
    )
)

# Define server logic
server <- function(input, output) {
  #defining reactive values that communicate with the UI
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
  
  #Defining global variables. 
  font <- list(family = "Times New Roman",size = 24,color = "#7f7f7f")
  initialStepsize <- 0.1
  analysisNames <- c("TP/FP", "TP/FN", "TP/TN", "FP/TN", "FP/FN", "FN/TN", "Positive Rate")
  
  output$distPlot <- renderPlotly({
    resultMatrixes <- list()  #List of results from 
    resultMatrixCount <- 0
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
            imat <- matrix(nrow = 0, ncol = 11)
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
                      FNPayOff <- (1/(exp(optimalBias*i)*j/(1-j)))
                      nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                      nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                      
                      imat <- rbind(imat, c(i, j, FNPayOff, log10(FNPayOff), optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FPdivTP, distance))
                  }
                  else{
                    
                    imat <- rbind(imat, c(i, j, NA, NA, NA, NA, NA, NA, NA, NA, NA)) ##11 colums
                  }
              }
            resultMatrixCount = resultMatrixCount + 1  
            resultMatrixes[[resultMatrixCount]] <- imat
            
            resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
            
          }
      }   
    if(2%in%ChosenAnalysis()){
          AnalysisName <- "TP/FN"
          for (i in realdPrimeVector)
          {
              imat <- matrix(nrow = 0, ncol = 11)
              Goal <- FN()/(TP()+FN())
              optimalBias = qnorm(Goal) + i/2
              
              for (j in newbr){
                  nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                  nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                  nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                  nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                  FNPayOff <- (1/(exp(optimalBias*i)*j/(1-j)))
                  FNdivTP <- nFalseNegatives / nTruePositives
                  distance <- FNdivTP - FP()/TN()
                  imat <- rbind(imat, c(i, j, FNPayOff, log10(FNPayOff), optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FNdivTP, distance))
              }
              resultMatrixCount = resultMatrixCount + 1
              resultMatrixes[[resultMatrixCount]] <- imat
              
              resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
              
              
          }
      }
    if(3%in%ChosenAnalysis()){
          AnalysisName <- "TN/TP"
          Goal <- TN()/TP()
          for (i in realdPrimeVector)
          {
            imat <- matrix(nrow = 0, ncol = 11)
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
                  FNPayOff <- (1/(exp(optimalBias*i)*j/(1-j)))
                  
                  imat <- rbind(imat, c(i, j, FNPayOff, log10(FNPayOff), optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, TNdivTP, distance))
            }
            resultMatrixCount = resultMatrixCount + 1
            resultMatrixes[[resultMatrixCount]] <- imat
            
            resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
          }
          
     }
    if(4%in%ChosenAnalysis()){
          AnalysisName <- "FP/FN"
          Goal <- FP()/FN()
          for (i in realdPrimeVector)
          {
            imat <- matrix(nrow = 0, ncol = 11)
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
                  FNPayOff <- (1/(exp(optimalBias*i)*j/(1-j)))
                  
                  imat <- rbind(imat, c(i, j, FNPayOff, log10(FNPayOff), optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FPdivFN, distance))
              }
            resultMatrixCount = resultMatrixCount + 1
            resultMatrixes[[resultMatrixCount]] <- imat
            
            resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
          }
      }
    if(5%in%ChosenAnalysis()){
        AnalysisName <- "TN/FP"
        Goal <- FP()/(TN()+FP())
          for (i in realdPrimeVector)
          {
              optimalBias = qnorm(Goal) - i/2
              imat <- matrix(nrow = 0, ncol = 11)
              
              for (j in newbr){
                  nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                  nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                  nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                  nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                  FNPayOff <- (1/(exp(optimalBias*i)*j/(1-j)))
                  FPdivTN <- nFalsePositives / nTrueNegatives
                  distance <- FPdivTN - FP()/TN()
                  imat <- rbind(imat, c(i, j, FNPayOff, log10(FNPayOff), optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FPdivTN, distance))
              }
              resultMatrixCount = resultMatrixCount + 1
              resultMatrixes[[resultMatrixCount]] <- imat
              resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
          }
      }
    if(6%in%ChosenAnalysis()){
        AnalysisName <- "FN/TN"
          Goal <- FN()/TN()
          for (i in realdPrimeVector)
          {
            imat <- matrix(nrow = 0, ncol = 11)
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
                FNPayOff <- (1/(exp(optimalBias*i)*j/(1-j)))
                imat <- rbind(imat, c(i, j, FNPayOff, log10(FNPayOff), optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FNdivTN, distance))
              }
              else{
                imat <- rbind(imat, c(i, j, NA, NA, NA, NA, NA, NA, NA, NA, NA)) ##11 colums
              }
            }
            resultMatrixCount = resultMatrixCount + 1
            resultMatrixes[[resultMatrixCount]] <- imat
            resultNamesList <- append(resultNamesList, paste(AnalysisName, i, sep = " "))
        }
      }
    if(7%in%ChosenAnalysis()){
        AnalysisName <- "PositiveRate"
        Goal <- PRate()
        for (i in realdPrimeVector){
          imat <- matrix(nrow = 0, ncol = 11)
          for (j in newbr){
            SKIPPED <- FALSE
            if(0 >= Goal){
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
                PosRate <- (nFalsePositives+nTruePositives) / (nFalsePositives+nTruePositives+nFalseNegatives+nTrueNegatives)
                distance <- Goal - PosRate
                
                
                if(stepSize < 0.000000000001){
                  SKIPPED <- TRUE
                  break 
                }
              }
              
            }
            if (SKIPPED == FALSE){
              FNPayOff <- (1/(exp(optimalBias*i)*j/(1-j)))
              imat <- rbind(imat, c(i, j, FNPayOff, log10(FNPayOff), optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, PosRate, distance))
            }
            else{
              imat <- rbind(imat, c(i, j, NA, NA, NA, NA, NA, NA, NA, NA, NA)) ##11 colums
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
        test = matrix(resultMatrixes[[n]],nrow = length(newbr), ncol = 11)
  
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
      else if (selectedOutputValue() == "payoff ratio"){
        selectedColumn <- 3    
      }
      else if (selectedOutputValue() == "log 10 of payoff ratio"){
        selectedColumn <- 4    
      }
      else if (selectedOutputValue() == "result"){
        selectedColumn <- 10
      }
      else if (selectedOutputValue() == "distance"){
        selectedColumn <- 11    
      }
      
      displayMatrix <- matrix(resultMatrixes[[1]],nrow = length(newbr), ncol = 11)
      
      
      data <- data.frame(displayMatrix[realStart:realEnd,2], displayMatrix[realStart:realEnd,selectedColumn])
      k <- 2
      
      while (k <= length(resultMatrixes)){
        displayMatrix <- matrix(resultMatrixes[[k]],nrow = length(newbr), ncol = 11)
        newdata <- data.frame(displayMatrix[realStart:realEnd,selectedColumn])
        data <- cbind(data, newdata)
        k <- k+1
      }
      
      print(resultMatrixCount)
      
      ###location columns overview.
      mainColumns <- 0
      anchorColumns <- 0
      errorColumns <- 0
      
      if(!is.na(dPrimeReference()) && dPrimeReference() > 0){
        if(length(dPrimeErrors()) == 2&&!is.na(dPrimeErrors()[1]) && !is.na(dPrimeErrors()[2]) && dPrimeErrors()[1] > 0 && dPrimeErrors()[2] > 0){
          if(length(ChosenAnalysis()) > 1){
            mainColumns <- c((2:(length(dPrimes())+1)),((length(dPrimes())+5):(length(dPrimes())*2+4)))
            anchorColumns <- c((length(dPrimes())+2), (length(dPrimes())*2+5))
            errorColumns <- c((length(dPrimes())+3),(length(dPrimes())+4), (length(dPrimes())*2+6), (length(dPrimes())*2+7))
          }else{
            mainColumns <- (2:(length(dPrimes())+1))
            anchorColumns <- (length(dPrimes())+2)
            errorColumns <- c((length(dPrimes())+3),(length(dPrimes())+4))
          }
        }else{
          if(length(ChosenAnalysis()) > 1){
            mainColumns <- c((2:(length(dPrimes())+1)),((length(dPrimes())+3):(length(dPrimes())*2+2)))
            anchorColumns <- c((length(dPrimes())+2), (length(dPrimes())*2+3))
          }else{
            mainColumns <- (2:(length(dPrimes())+1))
            anchorColumns <- (length(dPrimes())+2)
          }
        }
      }else{
        if(length(dPrimeErrors()) == 2&&!is.na(dPrimeErrors()[1]) && !is.na(dPrimeErrors()[2]) && dPrimeErrors()[1] > 0 && dPrimeErrors()[2] > 0){
          if(length(ChosenAnalysis()) > 1){
            mainColumns <- c((2:(length(dPrimes())+1)),((length(dPrimes())+4):(length(dPrimes())*2+3)))
            errorColumns <- c((length(dPrimes())+2),(length(dPrimes())+3), (length(dPrimes())*2+4), (length(dPrimes())*2+5))
          }else{
            mainColumns <- (2:(length(dPrimes())+1))
            errorColumns <- c((length(dPrimes())+2),(length(dPrimes())+3))
          }
        }else{
          if(length(ChosenAnalysis()) > 1){
            mainColumns <- c((2:(length(dPrimes())+1)),((length(dPrimes())+2):(length(dPrimes())*2+1)))
          }else{
            mainColumns <- (2:(length(dPrimes())+1))
          }
        }
      }
      
      xTitle <- "Base Rate"
      yTitle <- selectedOutputValue()
      graphTitle <- "Base Rate"
      
      if (nchar(xAxisTitle()) > 0){
        xTitle<- xAxisTitle()
      }
      if (nchar(yAxisTitle()) > 0){
        yTitle<- yAxisTitle()
      }
      if (nchar(customGraphTitle()) > 0){
        graphTitle<- customGraphTitle()
      }
      
      if (!is.na(xAxisMin()) && !is.na(xAxisMax())){
        
        x <- list(title = xTitle,titlefont = font, range = c(xAxisMin(), xAxisMax()))
      }else{
        x <- list(title = xTitle,titlefont = font)  
      }
      
      if (!is.na(yAxisMin()) && !is.na(yAxisMax())){
        y <- list(title = yTitle,titlefont = font, range = c(yAxisMin(), yAxisMax()))
      }else{
        y <- list(title = yTitle,titlefont = font)  
      }
      
      theColors <- getColors(length(dPrimes()))
      
      fig <- plot_ly(data = data, type ='scatter', mode = 'lines') %>% layout(title = graphTitle, xaxis = x, yaxis = y)
      
      if(!is.na(yAnchorVal())){
        fig <- fig %>% add_lines(y = yAnchorVal(), x = baseRates(), name = yAnchorName(),  line = list(color = "grey")) 
      }
      
      if(!is.na(xAnchorVal())){
        fig <- fig %>% add_lines(y = (c(max(data), min(data))), x = xAnchorVal(), name = xAnchorName(),  line = list(color = "grey")) 
      }
      
      if(length(errorColumns) > 1){
        fig <- fig %>% add_polygons(x=c(data[,1], rev(data[,1])), y= c(data[,(errorColumns[1])],rev(data[,errorColumns[2]])), mode='lines', name=paste(analysisNames[ChosenAnalysis()[1]],"Error"), color = I("dark grey"), opacity =0.8, line = list (width = 0))
        if(length(errorColumns)>2){
          fig <- fig %>% add_polygons(x=c(data[,1], rev(data[,1])), y= c(data[,errorColumns[3]],rev(data[,errorColumns[4]])), mode='lines', name=paste(analysisNames[ChosenAnalysis()[2]],"Error"), color = I("dark grey"), opacity =0.8, line = list (width = 0))
        }
      }
      dPrime <- 1
      roundsofDprime <- 1
      
      for (l in 1:length(mainColumns)){
          if (roundsofDprime == 1){
            style <- "line"
          } 
          else if(roundsofDprime == 2){
            style <- "dash"
          }
          else{
            style <- "dot"
          }
          fig <- fig %>% add_trace(x=data[,1], y= data[,mainColumns[l]], mode='lines', name=resultNamesList[(mainColumns[l]-1)], line = list(color= theColors[dPrime], dash = style))
          l <- l+1
          dPrime <- dPrime+1
          if(length(dPrimes())<dPrime){
            dPrime <- 1
            roundsofDprime <- roundsofDprime + 1
          }
      }
      
      
      if (length(anchorColumns) > 1 || (length(anchorColumns) ==1 && anchorColumns != 0)){
        fig <- fig %>% add_trace(x=data[,1], y= data[,anchorColumns[1]], mode='lines', name=resultNamesList[(anchorColumns[1]-1)], line = list(color= "000000", dash = "dot"))
        if(length(anchorColumns) > 1){
            fig <- fig %>% add_trace(x=data[,1], y= data[,anchorColumns[2]], mode='lines', name=resultNamesList[(anchorColumns[2]-1)], line = list(color= "000000", dash = "dot"))
        }
      }
      
      fig
    }
  })
  
  output$chosenAnalysis <- renderText({
    def <- "No Chosen Analysis"
    ret <- ""
    vec <- vector()
      if (!is.na(TP()) && TP() > 0){
          if (!is.na(FP()) && FP() > 0){
              vec <- append(vec, 1)
              if (nchar(ret) > 1){
                ret <- paste(ret, "TP/FP", sep = " & ")  
              }
              else{
                ret <- "TP/FP"
              }
            }
          else if (!is.na(FN()) && FN()> 0){
            vec <- append(vec, 2)
              if (nchar(ret) > 1){
                ret <- paste(ret, "TP/FN", sep = " & ")  
              }
              else{
                ret <- "TP/FN"
              }
          }
          else if (!is.na(TN()) && TN() > 0){
            vec <- append(vec, 3)
              if (nchar(ret) > 1){
                ret <- paste(ret, "TP/TN", sep = " & ")  
              }
              else{
                ret <- "TP/TN"
              }
          }
      }
      if (!is.na(FP()) && FP()> 0){
          if (!is.na(FN()) && FN() > 0){
            vec <- append(vec, 4)
              if (nchar(ret) > 1){
                ret <- paste(ret, "FP/FN", sep = " & ")
              }
              else{
                ret <- "FP/FN"
              }
          }
          else if (!is.na(TN()) && TN() > 0){
            vec <- append(vec, 5)
              if (nchar(ret) > 1){
                ret <- paste(ret, "FP/TN", sep = " & ")
              }
              else{
                ret <- "FP/TN"
              }
          }
      }
      if (!is.na(FN()) && FN() > 0){
          if(!is.na(TN()) && TN() > 0){
            vec <- append(vec, 6)
              if (nchar(ret) > 1){
                ret <- paste(ret, "FN/TN", sep = " & ") 
              }
              else{
                ret <- "FN/TN"
              }
          }
      }
      if (!is.na(PRate()) &&PRate() >0 && PRate() < 1){
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
