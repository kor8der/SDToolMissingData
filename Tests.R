

getThing <- function(Goal){
  print("before first for")
  imat <- matrix(nrow = 0, ncol = 10)
  dPrimes <- c(0.5)
  baseRates <- c(0.49, 0.51)
  Goal <- 1
}
  

  
  findOptimalRange <- function (matrixes){
    
    for (matrix in matrixes){
      length <- matrix.length
      
    }
    
  }
  

  
  megamats = list()
  for (i in 1:2){
    imats = array(dim= c(0,2,3))  
    for (j in 1:3){
      imat <- matrix(nrow =0, ncol=2)          
      for (k in 1:5){
        imats[,,] <- rbind(imats, c(j, k))    
      } 
      imats[[,,j]] <-  imat
    }
    megamats[[i]] <-  imats
  }
  
  lil <- matrix(megamats[1], nrow=3, ncol=3)
  lil
  
  
  
  
  # Create two vectors 
  data1 <- c(1,2,3,4,5,6)
  data2 <- c(60, 18, 12, 13, 14, 19)

  
  # pass these vectors as input to the array.
  #  3 rows,3 columns and 3 arrays
  result <- array(c(data1, data2), dim = c(3,3,3))
  print(result)
  
  
  imats = array(dim= c(3,3,3))
  
  imats[1,1,1] <- 22
  
  
  imats
  
  
  
  
  megamats <- list()
  for (methods in 1:5){
    imats <- matrix(nrow=0, ncol=4)
    for (j in 1:3){
      imats <- rbind(imats, c(33,55,11,00))
    }
    megamats[[methods]] <- imats
  }
    
  megamats
  
  
  test = matrix(megamats[[1]],nrow = 3, ncol = 4)
  is.na(test)
  
  startN <- 1
  endN <- 12
  count <- 0
  reachedStart <- FALSE
  reachedEnd <- FALSE
  
  for (x in is.na(test)){
    
    print(x)
    if (reachedStart == FALSE){
      if (x == FALSE){
        reachedStart <- TRUE
      }
      else{
        startN <- startN+1
      }
    }
    else{
      if (reachedEnd == FALSE){
        if (x == TRUE){
          
          reachedEnd  <- TRUE
          endN <- count
          
        }
      }
    }
    count <- count+1
  }
  
startN  
endN



resultMatrixes <- list()
resultMatrixCount <- 1
newbr = c(0.112,0.113,0.114,0.115,0.116,0.117)
realbrLength <- length(newbr)
resultNamesList <- list()
AnalysisName <- "PositiveRate"
Goal <- 0.93
initialStepsize <- 0.1
for (i in c(0.5,1.0)){
  imat <- matrix(nrow = 0, ncol = 10)
  for (j in c(0.112,0.113,0.114,0.115,0.116,0.117)){
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
      
      while (abs(distance) > 0.001){
        print(distance)
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
imat



