msft <- readRDS('C:/Users/jerem/OneDrive/Documents/msft.rds')
msft_use <- linear_imputatation(as.data.frame(msft[-1]))

mm <- mean_change(msft_use)
cv <- cov_change(msft_use)
gc <- graph_change(msft_use)

mm_bs <- bs_func(msft_use,mean_change)
cv_bs <- bs_func(msft_use,cov_change)
gc_bs <- bs_func(msft_use,graph_change)

plot_fd(msft_use,CPs=mm_bs)
plot_fd(msft_use,CPs=cv_bs)
plot_fd(msft_use,CPs=gc_bs)

graph_change <- function(data, treeType='MST', kTrees=5, 
                         error='L2', dataUse='Orig'){
  
  data_dist <- calculateDistanceMatrix(
    data=data, silent = T, 
    errType=error, dataUse = dataUse)
  
  if(treeType=='MST'){
    tree <- mstree(as.dist(data_dist), ngmax=kTrees)
  }else{
    stop('Sorry no other trees')
  }
  
  gc <- tryCatch({
    tmp <- gseg1(ncol(data), tree, pval.perm = T, B=100, pval.appr = F)
    ifelse(tmp$pval.perm$ori$pval<0.05,tmp$scanZ$ori$tauhat,NA)
  },error= function(x){NA})
   
  gc
}

bs_func <- function(data,method, addAmt=0, silent=F){
  
  potential_cp <- method(data)
  
  # No Change Point Detected
  if(is.na(potential_cp)) return()
  
  # Display progress
  if(!silent)
    cat(paste0('ChangePoint Detected (',1+addAmt,'-' ,addAmt+ncol(data),' at ',
               addAmt+potential_cp,'): Segment Data and Re-Search\n'))
  
  return(c(
    bs_func(data=data[,1:potential_cp], 
            method=method,
            addAmt = addAmt,
            silent=silent),
    potential_cp + addAmt,
    bs_func(data=data[,(potential_cp+1):ncol(data)], 
            method=method,
            addAmt = addAmt+potential_cp,
            silent=silent)
  ))
}

msft_cidr <- msft_use
for(i in 1:nrow(msft_use)){
  msft_cidr[i,] <- log(msft_use[i,]) - log(msft_use[1, ])
}

mm_bs1 <- bs_func(msft_cidr,mean_change)
cv_bs1 <- bs_func(msft_cidr,cov_change)
gc_bs1 <- bs_func(msft_cidr,graph_change)

plot_fd(msft_cidr,CPs=mm_bs1)
plot_fd(msft_cidr,CPs=cv_bs1)
plot_fd(msft_cidr,CPs=gc_bs1)
