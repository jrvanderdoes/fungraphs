###########
#
## Libraries
#
###########
library(fda)
library(fda.usc)
library(gSeg)
library(ade4)
library(lattice)
library(plotly)
library(tidyr)
library(nbpMatching)

###########
#
## Functions
#
###########

# constrFunc <- function(x, minVal, maxVal){
#   max(minVal,min(maxVal,x))
# }

createMDT <- function(distance, K){
  # Setup Data
  dropLast <- F
  distance_tmp <- distance
  if(nrow(distance_tmp)%%2){
    # If odd rows, add blank
    distance_tmp <- rbind(distance_tmp,0)
    distance_tmp <- cbind(distance_tmp,0)
    dropLast <- T
  }

  MDT <- matrix(0, 0, 2)

  for (i in 1:K) {
    # Format distance data
    data_distmat <- distancematrix(distance_tmp)

    # Get MDT with non-directional edges
    res <- nonbimatch(data_distmat)
    edge_new <- unique(res$matches[,c('Group1.Row','Group2.Row')])

    MDT <- rbind(MDT, edge_new)

    # Reduce data (If have another tree, this will make it perpendicular)
    for(j in 1:nrow(edge_new)){
      distance_tmp[edge_new$Group1.Row[j], edge_new$Group2.Row[j]] <-
        distance_tmp[edge_new$Group2.Row[j], edge_new$Group1.Row[j]] <- Inf
    }
  }

  # Drop fake pair
  if(dropLast){
    MDT <- MDT[MDT[,1]!=nrow(distance_tmp),]
    MDT <- MDT[MDT[,2]!=nrow(distance_tmp),]
  }

  MDT
}

calculateDistanceMatrix <- function(data, silent=F,
                                    errType='L2', dataUse='Orig', ln_n=3){

  n <- ncol(data)
  # Generate distance matrix with distances
  lenSol <- n*(n-1)/2
  distData <- matrix(0, nrow = n, ncol = n)

  iter <- 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      iter <- iter+1

      distData[i,j] <- distData[j,i] <-
        calcDistance(data[,i], data[,j], errType, dataUse, ln_n)
    }
    if(!silent)
      cat(paste0('Completed: ',round(iter/lenSol,4)*100,'%\n'))
  }

  distData
}

calcDistance <- function(data1, data2, type, dataUse, n=3, M=100){
  ## Function
  .getDist <- function(data1, data2, type, n=3){
    if(type == 'Max'){
      result <- max(data1-data2)
    }else if(type == 'L1'){
      result <- sum(abs(data1-data2))
    }else if(type == 'L2'){
      result <- sum(sqrt(abs((data1-data2)^2)))
    }else if(type == 'LN'){
      result <- sum((data1-data2)^n)^(1/n)
    } else{
      stop('Error: Use L1, L2, LN (with n), or Max for type')
    }

    result
  }

  ## Code
  result <- 0

  ## Use CF or Orig Data as needed
  if(dataUse=='CE'){
    # W <- as.data.frame(sapply(rep(0,M),sde::BM, N=length(data1)-1))
    # ce1s <- exp(complex(imaginary = 1)*(t(data1) %*% as.matrix(W)))
    # ce2s <- exp(complex(imaginary = 1)*(t(data2) %*% as.matrix(W)))
    #
    # result_n <- sum(sapply(rbind(ce1s,ce2s), .getDist1,type=type, n=n)) / M
    # # dist1 just changes dist to take a vector of two values for the first value

    # This loop is faster than array above
    for(i in 1:M){
      v <- sde::BM(N=length(data1)-1)
      ce1 <- exp(complex(imaginary = 1)*(data1 %*% v))
      ce2 <- exp(complex(imaginary = 1)*(data2 %*% v))

      result <- result + .getDist(ce1, ce2, type, n)
    }
    result <- result / M
  }else if(dataUse=='Orig'){
    result <- .getDist(data1, data2, type, n)
  }else{
    stop('Error: Use CE or Orig for dataUse')
  }

  result
}

