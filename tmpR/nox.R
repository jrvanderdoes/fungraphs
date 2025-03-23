nox <- read.csv("C:/Users/jerem/Downloads/nox.csv")[,-1]
nox$date <- as.Date(nox$date,'%Y-%m-%d')

nox1 <- t(nox)
colnames(nox1) <- nox1[1,]
nox1 <- as.matrix(nox1[-1,])

plot_fd(nox1)

# 2020-2022
nox0 <- nox[year(nox$date) %in% c(2018:2023),]
nox0 <- t(nox0)
colnames(nox0) <- nox0[1,]
nox0 <- as.matrix(nox0[-1,])

plot_fd(nox0)

### Month as FD
nox2 <- data.frame('date'=nox$date,
                   'mean'=rowMeans(nox[,-1]))
nox2$year = year(nox2$date)
nox2$month = month(nox2$date)
nox2$day = day(nox2$date)
nox2<- nox2[!(nox2$month==2 & nox2$day==29),]
nox2$percent = (day(nox2$date)-1) /
  ifelse(month(nox2$date)==2,27,
         ifelse(month(nox2$date)%in%c(1,3,5,7,8,10,12),30,29))

nox2 <- nox2 %>% pivot_wider(id_cols = percent,
                     names_from = c(year,month),
                     values_from = mean)
nox2 <- nox2[order(nox2$percent),]

plot_fd(nox2)

######################################
use_dat <- apply(nox0,MARGIN = 2, as.numeric)
nox_int <- linear_imputatation(use_dat,use.prev.curve = T)

nox_fd <- fda::Data2fd(1:nrow(nox_int), nox_int)
# Play with more components
nPCs <- 10
nox_fpca <- fda::pca.fd(nox_fd, nharm = nPCs)
nox_fpca_comp <- nox_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for (i in 1:nPCs) {
  ts_dat[[i]] <- ts(nox_fpca_comp[, i], freq = 7)
  comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps_resids[[i]] <- resid(comps[[i]])
}
nox_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
#   Want: 24 x 365
#     forecast: 365 x 3
#    coefs: 26  x 3
#      coefs %*% comps: 26 x 365
#    Eval: 24 x 26
#       eval %*% orig: 24 x 365
orig_coefs <- nox_fpca$harmonics$coefs %*% t(nox_fpca_forecast)
eval_fd_vals <- eval.basis(1:nrow(nox_int), nox_fd$basis) %*% orig_coefs

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

  gc_ver[!is.na(gc_ver)]
}

gc <- bs_func(data=eval_fd_vals, method=graph_change)
gc_ver <- ver_func(eval_fd_vals, gc)


