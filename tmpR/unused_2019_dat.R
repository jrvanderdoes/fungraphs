library(tidyverse)

# 2019 - val2

# filenames <- list.files("C:\\Users\\jerem\\Downloads\\marginalpdbcpt_2019",
#                         pattern="*.*", full.names=TRUE)
# filenames <- c(filenames,
#                list.files("C:\\Users\\jerem\\Downloads\\marginalpdbcpt_2020",
#                pattern="*.*", full.names=TRUE) )
filenames <- list.files("C:\\Users\\jerem\\Downloads\\marginalpdbcpt_2020",
                        pattern="*.*", full.names=TRUE)
ldf <- lapply(filenames, read.csv,sep = ";",skip = 1,header = F)
df <- do.call("rbind", ldf)
df <- df[,-7]
colnames(df) <- c('year','month','day','hour','val1','val2')
df$date <- as.Date(paste0(df$year,'-',df$month,'-',df$day),
                   format = '%Y-%m-%d')

df_plot <- df %>%
  na.omit %>%
  pivot_wider(id_cols = hour,names_from = date,values_from = val2)
df_plot <- df_plot[-25,-1]
df_plot <- as.data.frame(df_plot)

#plot_fd(df_plot)

df_plot <- linear_imputatation(df_plot,use.prev.curve = T)
data_dist <- calculateDistanceMatrix(
  df_plot, silent = T, errType='L2', dataUse = 'Orig')

# MDT5 <- createMDT(data_dist,5)
# gSeg_MDT5 <- gseg1_new(ncol(df_plot), MDT5, pval.perm = T, pval.appr = F)

###########################
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

mm_bs1 <- bs_func(data_dist,graph_change)
mm_bs1

plot_fd(df_plot,mm_bs1)
