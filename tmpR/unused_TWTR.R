library(tidyverse)

TWTR <- read.csv(paste0("C:/Users/jerem/OneDrive/",
                       "Documents/School/Waterloo/",
                       "Research/RPackages/fungraphs/",
                       "tmpR/TWTR20192022.csv"))
TWTR$dateTime <- as.POSIXct(TWTR$Local.time,format = "%d.%m.%Y %H:%M:%S")
TWTR$hour <- format(TWTR$dateTime, "%H")
TWTR$dayYear <- format(TWTR$dateTime, "%d.%m.%Y")
TWTR_format <- TWTR %>%
  select(Close,hour,dayYear, Volume) %>%
  #filter(as.numeric(hour)>=9 & as.numeric(hour)<17) %>%
  filter(Volume>0) %>%
  pivot_wider(id_cols = dayYear,
              names_from = hour,
              values_from = Close,
              values_fn = list) %>%
  unnest(cols = everything() ) %>%
  t() %>%
  as.matrix()
colnames(TWTR_format) <- TWTR_format[1,]
TWTR_format <- TWTR_format[-1,]
TWTR_format <- apply(TWTR_format, 2 ,as.numeric)
TWTR_final <- TWTR_format

for(i in 1:nrow(TWTR_format)){
  TWTR_final[i,] <- log(TWTR_format[i,]) - log(TWTR_format[1, ])
}

TWTR_final <- linear_imputatation(TWTR_final,use.prev.curve = T)

plot_fd(TWTR_final)

########################################
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
    n = ncol(data)
    tmp <- gseg1(n, tree,
                 n0=max(0.05*n,20), n1=min(0.95*n,n-20),
                 pval.perm = T, B=100, pval.appr = F)
    ifelse(tmp$pval.perm$ori$pval<0.05,tmp$scanZ$ori$tauhat,NA)
  },error= function(x){NA})

  gc
}

bs_func <- function(data, method, addAmt=0, silent=F){
  if(ncol(data)<=30) return()

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

ver_func <- function(data, cps){
  cps_new <- c()
  cps_ver <- c(0,cps,ncol(data))

  for(i in 2:length(cps_ver)){
    dat <- data[,(cps_ver[i-1]+1):cps_ver[i]]


    data_dist <- calculateDistanceMatrix(
      data=dat, silent = T,
      errType='L2', dataUse = 'Orig')
    tree <- mstree(as.dist(data_dist), ngmax=5)

    gc <- tryCatch({
      n = ncol(dat)
      tmp <- gseg1(n, tree,
                   n0=max(0.05*n,20), n1=min(0.95*n,n-20),
                   pval.perm = T, B=100, pval.appr = F)
      ifelse(tmp$pval.perm$ori$pval<0.05,tmp$scanZ$ori$tauhat,NA)
    })
    cps_new <- c(cps,gc)

  }

  cps_new[!is.na(cps_new)]
}

set.seed(1234)
gc <- bs_func(data=TWTR_final, method=graph_change)
gc_ver <- ver_func(TWTR_final, gc)

plot_fd(TWTR_final,gc_ver)

