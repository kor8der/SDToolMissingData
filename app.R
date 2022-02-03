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
            sliderInput("brRange", "Define br range:", c(0.01,0.8), min = 0.01, max = 0.999, step = 0.001),
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

analysisLogic <- function(TP, FP, FN, TN, BR, CR) {
    
    if (TP > 0){
        if (FP > 0){
            return(1)
        }
        else if (FN > 0){
            return(2)
        }
        else if (TN > 0){
            return(3)
        }
        else{
            return(0)
        }
        
    }
    else if (FP > 0){
        if (FN > 0){
            return(4)
        }
        else if (TN > 0){
            return(5)
        }
        else{
            return(0)
        }
    }
    if (FN > 0)
    {
        if (TN > 0){
            return(6)
        }
        else{
            return(0)
        }
    }
    else{
        return(0)
    }
    
}

showaplot <- function(teepee, efpee){
    plot(c(teepee, efpee))
    
}

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
    imatCount = 0
    font <- list(
        family = "Courier New, monospace",
        size = 24,
        color = "#7f7f7f"
    )
    initialStepsize <- 0.1
    
    output$distPlot <- renderPlotly({
        if(ChosenAnalysis() == 1){
            imat <- matrix(nrow = 0, ncol = 10)
            Goal <- FP()/TP()
            realbrLength <- 0
            newbr = baseRates()
            for (i in dPrimes())
            {
                for (j in newbr){
                    if((1-j)/j < Goal){
                        break
                    }
                    
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
                            print("SKIPPED")
                            break 
                        }
                        
                    }
                    if (SKIPPED == FALSE){
                        FNPayOff <- log10(1/(exp(optimalBias*i)*j/(1-j)))
                        nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
                        nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
                        
                        imat[imatCount] <- rbind(imat, c(i, j, FNPayOff, optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FPdivTP, distance))
                    }
                }
                
                imatCount = imatCount + 1
                
                if (realbrLength == 0){
                    realbrLength <- nrow(imat)
                    newbr = imat[1:realbrLength,2]
                }
            }
            
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
            
            data <- data.frame(imat[1:realbrLength,2], imat[1:realbrLength,selectedColumn])
            colnames <- c("br", dPrimes()[1])
            colnames(data) <- colnames
            k <- 1
            
            while (k < length(dPrimes()))
            {
                ##print(c("pre data: ", k, length(baseRates()), nrow(imat)))
                ##print(length(baseRates())*k+1)
                ##print(length(baseRates())*(k+1))
                ##test <- imat[length(baseRates())*k+1:length(baseRates()*(k+1)),2]
                ##print(c("test1 success: ", k))
                ##test2 <- imat[length(baseRates())*k+1:length(baseRates()*(k+1)),3]
                ##print(c("test2 success: ", k))
                ##print(c("start",realbrLength*k+1))
                ##print(c("end",realbrLength*(k+1)))
                ##print(c("length",realbrLength*(k+1)-realbrLength*k))
                ##print(c("total", nrow(imat)))
                ##START HERE
                ##o <- 0
                ##while (o < 1000)
                ##{
                ##    print(c(o, imat[o,selectedColumn]))    
                ##    o <- o+1
                ##}    
                
                ##print(imat[realbrLength*k+1, selectedColumn])
                ##print(imat[realbrLength*(k+1), selectedColumn])
                ##print(imat[402,1:5])
                ##print(imat[402:802,selectedColumn])
                ##print(realbrLength*k+1)
                ##alpha <- 401
                ##newdata <- data.frame(imat[(alpha*k+1):802,selectedColumn])
                newdata <- data.frame(imat[(realbrLength*k+1):(realbrLength*(k+1)),selectedColumn])
                ##data <- data.frame(imat[1:length(baseRates()),2], imat[1:length(baseRates()),selectedColumn])
                ##print("pre newdata created")
                ##print(newdata[1:2,])
                
                data <- cbind(data, newdata)
                
                ##print(data[1:2,])
                
                ##print("newdata created")
                
                ##print ("row bound")
                colnames <- c(colnames, dPrimes()[k+1])             
                ##print ("col named")
                ##print(c("post data: ", k))
                k <- k+1
                ##print(k)
            }
            
            colnames(data) <- colnames
            
            x <- list(
                title = "Base Rate",
                titlefont = font
            )
            
            y <- list(
                title = "log10 nFP/nFN",
                titlefont = font
            )
            
            fig <- plot_ly(x= data[,1],y= data[,2], name = paste("d' ", colnames[2]), type ='scatter', mode='lines')
            l <- 3 
            while (l <= length(colnames))
            {
                fig <- fig %>% add_trace(y= data[,l], name = paste("d' ", colnames[l]), mode='lines')
                l <- l+1
            }
            fig
        }   
        else if(ChosenAnalysis() == 2){
            imat <- matrix(nrow = 0, ncol = 10)
            newbr = baseRates()
            realbrLength <- length(newbr)
            for (i in dPrimes())
            {
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
                
            }
            
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
            
            data <- data.frame(imat[1:realbrLength,2], imat[1:realbrLength,selectedColumn])
            
            colnames <- c("br", dPrimes()[1])
            
            colnames(data) <- colnames
            
            k <- 1
            
            while (k < length(dPrimes()))
            {
                newdata <- data.frame(imat[(realbrLength*k+1):(realbrLength*(k+1)),selectedColumn])
                data <- cbind(data, newdata)
                
                colnames <- c(colnames, dPrimes()[k+1])             
                k <- k+1
            }
            
            colnames(data) <- colnames
            
            x <- list(
                title = "Base Rate",
                titlefont = font
            )
            
            y <- list(
                title = "log10 nFP/nFN",
                titlefont = font
            )
            fig <- plot_ly(x= data[,1],y= data[,2], name = paste("d' ", colnames[2]), type ='scatter', mode='lines')
            l <- 3 
            
            while (l <= length(colnames))
            {
                fig <- fig %>% add_trace(y= data[,l], name = paste("d' ", colnames[l]), mode='lines')
                l <- l+1
            }
            
            fig
        }
        else if(ChosenAnalysis() == 3){
            imat <- matrix(nrow = 0, ncol = 10)
            Goal <- TN()/TP()
            
            for (i in dPrimes())
            {
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
            }
            
            print(c("Start Plot!", "5"))
            
            
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
            
            data <- data.frame(imat[1:length(baseRates()),2], imat[1:length(baseRates()),selectedColumn])
            colnames <- c("br", dPrimes()[1])
            colnames(data) <- colnames
            k <- 1
            
            while (k < length(dPrimes()))
            {
                newdata <- data.frame(imat[length(baseRates())*k+1:length(baseRates()*(k+1)),selectedColumn])
                
                data <- cbind(data, newdata)
                
                colnames <- c(colnames, dPrimes()[k+1])             
                
                k <- k+1
            }
            
            colnames(data) <- colnames
            
            x <- list(
                title = "Base Rate",
                titlefont = font
            )
            
            y <- list(
                title = "log10 nFP/nFN",
                titlefont = font
            )
            
            fig <- plot_ly(x= data[,1],y= data[,2], name = paste("d' ", colnames[2]), type ='scatter', mode='lines')
            
            l <- 3 
            
            while (l <= length(colnames))
            {
                fig <- fig %>% add_trace(y= data[,l], name = paste("d' ", colnames[l]), mode='lines')
                l <- l+1
            }
            
            fig
        }
        else if(ChosenAnalysis() == 4){
            imat <- matrix(nrow = 0, ncol = 10)
            Goal <- FP()/FN()
            for (i in dPrimes())
            {
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
            }
            
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
            
            data <- data.frame(imat[1:length(baseRates()),2], imat[1:length(baseRates()),selectedColumn])
            
            colnames <- c("br", dPrimes()[1])
            
            colnames(data) <- colnames
            
            k <- 1
            
            while (k < length(dPrimes()))
            {
                newdata <- data.frame(imat[length(baseRates())*k+1:length(baseRates()*(k+1)),selectedColumn])
                
                data <- cbind(data, newdata)
                
                colnames <- c(colnames, dPrimes()[k+1])             
                k <- k+1
            }
            
            colnames(data) <- colnames
            
            x <- list(
                title = "Base Rate",
                titlefont = font
            )
            
            y <- list(
                title = "log10 nFP/nFN",
                titlefont = font
            )
            
            fig <- plot_ly(x= data[,1],y= data[,2], name = paste("d' ", colnames[2]), type ='scatter', mode='lines')
            l <- 3 
            
            while (l <= length(colnames))
            {
                fig <- fig %>% add_trace(y= data[,l], name = paste("d' ", colnames[l]), mode='lines')
                l <- l+1
            }
            
            fig
        }
        else if(ChosenAnalysis() == 5){
            imat <- matrix(nrow = 0, ncol = 10)
            newbr = baseRates()
            realbrLength <- length(newbr)
            for (i in dPrimes())
            {
                Goal <- FP()/(TN()+FP())
                optimalBias = qnorm(Goal) - i/2
                
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
                
            }
            
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
            
            data <- data.frame(imat[1:realbrLength,2], imat[1:realbrLength,selectedColumn])
            colnames <- c("br", dPrimes()[1])
            colnames(data) <- colnames
            k <- 1
            while (k < length(dPrimes()))
            {
                newdata <- data.frame(imat[(realbrLength*k+1):(realbrLength*(k+1)),selectedColumn])
                data <- cbind(data, newdata)
                colnames <- c(colnames, dPrimes()[k+1])             
                k <- k+1
            }
            
            colnames(data) <- colnames
            
            x <- list(
                title = "Base Rate",
                titlefont = font
            )
            
            y <- list(
                title = "log10 nFP/nFN",
                titlefont = font
            )
            
            fig <- plot_ly(x= data[,1],y= data[,2], name = paste("d' ", colnames[2]), type ='scatter', mode='lines')
            l <- 3 
            
            while (l <= length(colnames))
            {
                fig <- fig %>% add_trace(y= data[,l], name = paste("d' ", colnames[l]), mode='lines')
                l <- l+1
            }
            
            fig
        }
        else if(ChosenAnalysis() == 6){
            imat <- matrix(nrow = 0, ncol = 10)
            Goal <- FN()/TN()
            realbrLength <- 0
            newbr = baseRates()
            for (i in dPrimes())
            {
                for (j in newbr){
                    if((1-j)/j > Goal)
                    {
                        next
                    }
                    
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
                            print("SKIPPED")
                            break 
                        }
                        
                    }
                    nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
                    nTruePositives <- (1- pnorm((optimalBias-i/2))) * (j)
                    FNPayOff <- log10(1/(exp(optimalBias*i)*j/(1-j)))
                    imat <- rbind(imat, c(i, j, FNPayOff, optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FNdivTN, distance))
                }
                
                if (realbrLength == 0)
                {
                    realbrLength <- nrow(imat)
                    newbr = imat[1:realbrLength,2]
                }
            }
            
            
            
            print(c("Start Plot!", "5"))
            
            
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
            
            data <- data.frame(imat[1:realbrLength,2], imat[1:realbrLength,selectedColumn])
            
            colnames <- c("br", dPrimes()[1])
            
            colnames(data) <- colnames
            
            k <- 1
            
            while (k < length(dPrimes()))
            {
                ##print(c("pre data: ", k, length(baseRates()), nrow(imat)))
                ##print(length(baseRates())*k+1)
                ##print(length(baseRates())*(k+1))
                ##test <- imat[length(baseRates())*k+1:length(baseRates()*(k+1)),2]
                ##print(c("test1 success: ", k))
                ##test2 <- imat[length(baseRates())*k+1:length(baseRates()*(k+1)),3]
                ##print(c("test2 success: ", k))
                ##print(c("start",realbrLength*k+1))
                ##print(c("end",realbrLength*(k+1)))
                ##print(c("length",realbrLength*(k+1)-realbrLength*k))
                ##print(c("total", nrow(imat)))
                ##START HERE
                ##o <- 0
                ##while (o < 1000)
                ##{
                ##    print(c(o, imat[o,selectedColumn]))    
                ##    o <- o+1
                ##}    
                
                ##print(imat[realbrLength*k+1, selectedColumn])
                ##print(imat[realbrLength*(k+1), selectedColumn])
                ##print(imat[402,1:5])
                ##print(imat[402:802,selectedColumn])
                ##print(realbrLength*k+1)
                ##alpha <- 401
                ##newdata <- data.frame(imat[(alpha*k+1):802,selectedColumn])
                newdata <- data.frame(imat[(realbrLength*k+1):(realbrLength*(k+1)),selectedColumn])
                ##data <- data.frame(imat[1:length(baseRates()),2], imat[1:length(baseRates()),selectedColumn])
                ##print("pre newdata created")
                ##print(newdata[1:2,])
                
                data <- cbind(data, newdata)
                
                ##print(data[1:2,])
                
                ##print("newdata created")
                
                ##print ("row bound")
                colnames <- c(colnames, dPrimes()[k+1])             
                ##print ("col named")
                ##print(c("post data: ", k))
                k <- k+1
                ##print(k)
            }
            
            ##print(colnames)
            
            ##print ("pre name binding")
            
            colnames(data) <- colnames
            
            ##print ("post name binding")
            
            ##print(colnames)
            
            ##print ("printing first two rows")
            ##print(data[1:2,])
            
            
            ##colnames(data) <- c("d", "br", "payOff", "optimalBias", "nFalsePositives", "nFalseNegatives", "nTruePositives", "nTrueNegatives", "FPdivFN", "distance")
            
            x <- list(
                title = "Base Rate",
                titlefont = font
            )
            
            ##print("x created")
            y <- list(
                title = "log10 nFP/nFN",
                titlefont = font
            )
            ##print("y created")
            
            ##print(dPrimes()[1])
            
            ##print ("printing first two rows")
            ##print(data[800:802,])
            
            fig <- plot_ly(x= data[,1],y= data[,2], name = paste("d' ", colnames[2]), type ='scatter', mode='lines')
            ##print(nrow(data[,1]))
            ##print ("printing first two columns and rows")
            ##print(data[1:2,1:2])
            
            ##print("fig created")
            
            l <- 3 
            ##print(c("l defined", l))
            
            ##print("entering while")
            while (l <= length(colnames))
            {
                fig <- fig %>% add_trace(y= data[,l], name = paste("d' ", colnames[l]), mode='lines')
                l <- l+1
            }
            
            
            
            
            ##print("returning fig")
            fig
        }
    })
    
    output$chosenAnalysis <- renderText({
        ret <- "No Chosen Analysis"
        if (TP() > 0)
        {
            if (FP() > 0){
                ChosenAnalysis(1)
                ret <- "TP/FP"
            }
            else if (FN()> 0){
                ChosenAnalysis(2)
                ret <- "TP/FN"
            }
            else if (TN() > 0){
                ChosenAnalysis(3)
                ret <- "TP/TN"
            }
        }
        else if (FP()> 0)
        {
            if (FN() > 0){
                ChosenAnalysis(4)
                ret <- "FP/FN"
            }
            else if (TN() > 0){
                ChosenAnalysis(5)
                ret <- "FP/TN"
            }
        }
        else if (FN() > 0)
        {
            if (TN() > 0){
                ChosenAnalysis(6)
                ret <- "FN/TN"
            }
        }
        ret
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
