base_path <- 'C:/Users/jerem/Downloads/Spanish-DayAhead/marginalpdbc_'

data <- data.frame()
year_data <- list()

years <- 2019:2020

for(i in 1:length(years)){
  files <- list.files(paste0(base_path,years[i]))

  day_data <- list()
  for(j in 1:length(files)){
    file <- files[j]
    day_data[[j]] <- read.csv2(paste0(base_path,years[i],'/',file),skip=1)[,1:6]
    colnames(day_data[[j]]) <- c('Year','Month','Day','Hour','P1','P2')
    day_data[[j]] <- na.omit(day_data[[j]])
  }
  year_data[[i]] <- do.call('rbind',day_data)
}

data <- do.call('rbind', year_data)
data$date <- as.Date(paste0(data$Year,'-',data$Month,'-',data$Day),'%Y-%m-%d')
data <- data[order(data$date,data$Hour),]
dates <- unique(data$date)

data_daily <- matrix(nrow=23,ncol = length(dates))
data_daily <- cbind(data.frame('Hours'=1:23),data_daily)

# data_daily <- data %>%
#   na.omit() %>%
#   pivot_wider(id_cols = Hour, names_from = date,values_from = P1)

for(i in 1:length(dates)){
  tmp <- data[data$date==dates[i],c('Hour',"P1")]
  if(nrow(tmp)==23){
    data_daily[,i+1] <- as.numeric(tmp[,2])
  }else{
    tmp <- merge(data.frame('Hour'=1:23),
                 unique(tmp),
                 all.x=T)
    data_daily[,i+1] <- as.numeric(tmp[,2])
  }
}


#plot_fd(as.matrix(data_daily[,-1]))

##################################################

graph_change <- function(data, treeType='MST', kTrees=5,
                         error='L2', dataUse='Orig'){

  n = ncol(data)
  n0 = round(max(0.05*n,30))
  n1 = round(min(0.95*n,n-30))
  if(n0>=n1) return(NA)

  data_dist <- calculateDistanceMatrix(
    data=data, silent = T,
    errType=error, dataUse = dataUse)

  if(treeType=='MST'){
    tree <- mstree(as.dist(data_dist), ngmax=kTrees)
  }else{
    stop('Sorry no other trees')
  }

  gc <- tryCatch({
    tmp <- gseg1(n, tree,
                 n0=n0, n1=n1,
                 pval.perm = T, B=300, pval.appr = F)
    ifelse(tmp$pval.perm$generalized$pval<0.05,
           tmp$scanZ$generalized$tauhat,
           NA)
  },
  error = function(x){
    NA
  })

  gc
}

bs_func <- function(data, addAmt=0, silent=F){
  potential_cp <- graph_change(data)

  # No Change Point Detected
  if(is.na(potential_cp)) return()

  # Display progress
  if(!silent)
    cat(paste0('ChangePoint Detected (',1+addAmt,'-' ,addAmt+ncol(data),' at ',
               addAmt+potential_cp,'): Segment Data and Re-Search\n'))

  return(c(
    bs_func(data=data[,1:potential_cp],
            addAmt = addAmt,
            silent=silent),
    potential_cp + addAmt,
    bs_func(data=data[,(potential_cp+1):ncol(data)],
            addAmt = addAmt+potential_cp,
            silent=silent)
  ))
}

ver_func <- function(data, cps){
  cps_new <- c()
  cps_ver <- na.omit(c(0,cps,ncol(data)))

  for(i in 2:length(cps_ver)){
    gc <- graph_change( data[, (cps_ver[i-1]+1):cps_ver[i]] )


    cps_new <- c(cps_new,cps_ver[i-1]+gc)

  }

  cps_new[!is.na(cps_new)]
}


use_data <- linear_imputatation(as.matrix(data_daily[,-1]),use.prev.curve = T)
set.seed(12345)
gc <- bs_func(data=use_data)
gc_ver <- ver_func(data=use_data, cps=gc)



plot_fd(use_data, gc)
