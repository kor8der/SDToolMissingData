

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