

getThing <- function(Goal){
  print("before first for")
  imat <- matrix(nrow = 0, ncol = 10)
  dPrimes <- c(0.5)
  baseRates <- c(0.49, 0.51)
  Goal <- 1
  
for (i in dPrimes)
  {
  print("before second for")
    for (j in baseRates){
      FNPayOff <- 0
      stepSize <- 1
      optimalBeta <-(1-j)/j*(1)/(exp(FNPayOff))
      optimalBias <- log(optimalBeta)/i
      nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
      nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
      nTruePositives <- 1 - pnorm((optimalBias-i/2)) * (j)
      nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
      FPdivFN <- nFalsePositives / nFalseNegatives
      
      distance <- Goal - FPdivFN
      print(c(FNPayOff, stepSize, Goal, FPdivFN, distance))
      
      
      if(distance < 0){
        MovingUp <- TRUE
      }
      else{
        MovingUp <- FALSE
      }
      
      while (abs(distance) > 0.01){
        if(distance < 0)
        {
          if (MovingUp == FALSE)
          {
            stepSize <- stepSize/2
          }
          FNPayOff <- FNPayOff - stepSize 
          MovingUP <- TRUE
        }
        else
        {
          if (MovingUp == TRUE)
          {
            stepSize <- stepSize/2
          }
          FNPayOff <- FNPayOff + stepSize  
          MovingUP <- FALSE
        }
        
        optimalBeta <-(1-j)/j*(1)/(exp(FNPayOff))
        optimalBias <- log(optimalBeta)/i
        nFalsePositives <- (1-pnorm((optimalBias+i/2))) * (1-j)
        nFalseNegatives <- pnorm((optimalBias-i/2)) * (j)
        nTruePositives <- 1 - pnorm((optimalBias-i/2)) * (j)
        nTrueNegatives <- pnorm((optimalBias+i/2)) * (1-j)
        FPdivFN <- nFalsePositives / nFalseNegatives
        
        distance <- Goal - FPdivFN
        print(c(FNPayOff, stepSize, Goal, FPdivFN, distance))
        
      }
      
      imat <- rbind(imat, c(i, j, FNPayOff, optimalBias, nFalsePositives, nFalseNegatives, nTruePositives, nTrueNegatives, FPdivFN, distance))
    }
  }
  
 plot(imat[1:length(baseRates),2], imat[1:length(baseRates),3])
}

getThing(1)
