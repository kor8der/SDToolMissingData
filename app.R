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
            numericInput("DesiredPrecision", "Define Desired Precision (note that changing this may alter performance by orders of magnitude):", 0.001, min = 0.00000000001, max = 1, step = 0.00001),
            sliderInput("dRange", "Define D' range:", c(0.5,3.0), min = 0.5, max = 8, step = 0.5),
            sliderInput("brRange", "Define br range:", c(0.049,0.801), min = 0.01, max = 0.999, step = 0.001),
            selectInput("SelectedOutput", "Pick output:", c("payoff ratio", "bias", "result", "distance")),
            ##numericInput("TotalPositives", "Define Total Positive observations:", 0, min = 0, max = 1000000, step = 1),
            ##numericInput("TotalNegatives", "Define Total Negative rate:", 0, min = 0, max = 1000000, step = 1),
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
    dPrimes <- reactive({seq(input$dRange[1],input$dRange[2], 0.5)})
    baseRates <- reactive({seq(input$brRange[1],input$brRange[2], 0.001)})
    desiredPrecision <- reactive({input$DesiredPrecision})
    selectedOutputValue <-reactive({input$SelectedOutput})
    font <- list(
        family = "Courier New, monospace",
        size = 24,
        color = "#7f7f7f"
    )
    initialStepsize <- 0.1
    
    output$distPlot <- renderPlotly({
      print (ChosenAnalysis())
      resultMatrixes <- list()
      resultMatrixCount <- 1
      newbr = baseRates()
      realbrLength <- length(newbr)
      
        if(1%in%ChosenAnalysis()){
            AnalysisName <- "FP/TP"
            Goal <- FP()/TP()
            print(dPrimes())
            for (i in dPrimes())
            {
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
              
              resultMatrixCount = resultMatrixCount + 1
            }
        }   
        if(2%in%ChosenAnalysis()){
            AnalysisName <- "TP/FN"
            for (i in dPrimes())
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
                
                resultMatrixCount = resultMatrixCount + 1
            }
        }
        if(3%in%ChosenAnalysis()){
            AnalysisName <- "TN/TP"
            Goal <- TN()/TP()
            for (i in dPrimes())
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
              
              resultMatrixCount = resultMatrixCount + 1
            }
            
       }
        if(4%in%ChosenAnalysis()){
            AnalysisName <- "FP/FN"
            Goal <- FP()/FN()
            for (i in dPrimes())
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
              
              resultMatrixCount = resultMatrixCount + 1
            }
        }
        if(5%in%ChosenAnalysis()){
          AnalysisName <- "TN/FP"
          Goal <- FP()/(TN()+FP())
            for (i in dPrimes())
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
                
                resultMatrixCount = resultMatrixCount + 1
            }
        }
        if(6%in%ChosenAnalysis()){
          AnalysisName <- "FN/TN"
            Goal <- FN()/TN()
            for (i in dPrimes())
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
              
              resultMatrixCount = resultMatrixCount + 1
          }
        }
      ## Restriction area 

      endV <- vector()     
      startV <- vector()      
      for (n in 1:length(resultMatrixes))
      {
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
      
      
      ##
      
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
      
      
      data <- data.frame(displayMatrix[realStart:realEnd,2], imat[realStart:realEnd,selectedColumn])
      ##colnames <- c("br", dPrimes()[1])
      ##colnames(data) <- colnames
      k <- 2
      
      while (k <= length(resultMatrixes))
      {
        displayMatrix <- matrix(resultMatrixes[[k]],nrow = length(newbr), ncol = 10)
        newdata <- data.frame(displayMatrix[realStart:realEnd,selectedColumn])
        data <- cbind(data, newdata)
        ##colnames <- c(colnames, dPrimes()[k+1])             
        k <- k+1
      }
      
      #colnames(data) <- colnames
      
      x <- list(
        title = "Base Rate",
        titlefont = font
      )
      
      y <- list(
        title = "log10 nFP/nFN",
        titlefont = font
      )
      
      #print(data)
      
      fig <- plot_ly(x= data[,1],y= data[,2], type ='scatter', mode='lines')
      #fig <- plot_ly(x= data[,1],y= data[,2], name = paste("d' ", colnames[2]), type ='scatter', mode='lines')
      
      l <- 3 
      
      ##print(resultMatrixes)
      
      while (l <= ncol(data))
      {
        fig <- fig %>% add_trace(y= data[,l], mode='lines')
        l <- l+1
        
      }
      
      fig
    })
    
    output$chosenAnalysis <- renderText({
      def <- "No Chosen Analysis"
      ret <- ""
      vec <- vector()
        if (TP() > 0)
        {
            if (FP() > 0){
                vec <- append(vec, 1)
                ret <- paste(ret, "TP/FP", sep = " ")
            }
            else if (FN()> 0){
              vec <- append(vec, 2)
                ret <- "TP/FN"
                ret <- paste(ret, "TP/FN", sep = " ")
            }
            else if (TN() > 0){
              vec <- append(vec, 3)
                ret <- "TP/TN"
                ret <- paste(ret, "TP/TN", sep = " ")
            }
        }
        if (FP()> 0)
        {
            if (FN() > 0){
              vec <- append(vec, 4)
                ret <- "FP/FN"
                ret <- paste(ret, "FP/FN", sep = " ")
            }
            else if (TN() > 0){
              vec <- append(vec, 5)
                ret <- "FP/TN"
                ret <- paste(ret, "FP/TN", sep = " ")
            }
        }
        if (FN() > 0)
        {
            if (TN() > 0){
              vec <- append(vec, 6)
                ret <- "FN/TN"
                ret <- paste(ret, "FN/TN", sep = " ")
            }
        }
      
      
      ChosenAnalysis(vec)
      if (length(ret) > 1){
        return(ret)
      }
      else {
        return(def)
      }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
