library(devtools)
load_all()
library(fda)
library(ade4)
library(gSeg)
library(ggplot2)

path <- 'C://Users//jerem//OneDrive//Documents//School//Waterloo//Research//RPackages//fungraphs//paper_trials_known'

trials <- 1000
evalPts <- seq(0,1,0.05)
ns <- c(15, 25, 50, 100, 200)
break_location <- 0.5

results <- data.frame('MDT1ori'=1:trials, 'MDT1wei'=NA,
                      'MDT1max'=NA, 'MDT1gen'=NA,
                      'MDT3ori'=NA, 'MDT3wei'=NA,
                      'MDT3max'=NA, 'MDT3gen'=NA,
                      'MDT5ori'=NA, 'MDT5wei'=NA,
                      'MDT5max'=NA, 'MDT5gen'=NA,
                      'MDT7ori'=NA, 'MDT7wei'=NA,
                      'MDT7max'=NA, 'MDT7gen'=NA,
                      'MDT15ori'=NA, 'MDT15wei'=NA,
                      'MDT15max'=NA, 'MDT15gen'=NA,

                      'NNL1ori'=NA, 'NNL1wei'=NA,
                      'NNL1max'=NA, 'NNL1gen'=NA,
                      'NNL3ori'=NA, 'NNL3wei'=NA,
                      'NNL3max'=NA, 'NNL3gen'=NA,
                      'NNL5ori'=NA, 'NNL5wei'=NA,
                      'NNL5max'=NA, 'NNL5gen'=NA,
                      'NNL7ori'=NA, 'NNL7wei'=NA,
                      'NNL7max'=NA, 'NNL7gen'=NA,
                      'NNL15ori'=NA, 'NNL15wei'=NA,
                      'NNL15max'=NA, 'NNL15gen'=NA,

                      'MST1ori'=NA, 'MST1wei'=NA,
                      'MST1max'=NA, 'MST1gen'=NA,
                      'MST3ori'=NA, 'MST3wei'=NA,
                      'MST3max'=NA, 'MST3gen'=NA,
                      'MST5ori'=NA, 'MST5wei'=NA,
                      'MST5max'=NA, 'MST5gen'=NA,
                      'MST7ori'=NA, 'MST7wei'=NA,
                      'MST7max'=NA, 'MST7gen'=NA,
                      'MST15ori'=NA, 'MST15wei'=NA,
                      'MST15max'=NA, 'MST15gen'=NA
)

results_save <- list()
results_simplify <- data.frame('Length'=rep(NA,length(ns)),
                               'MDT1ori'=NA, 'MDT1wei'=NA,
                               'MDT1max'=NA, 'MDT1gen'=NA,
                               'MDT3ori'=NA, 'MDT3wei'=NA,
                               'MDT3max'=NA, 'MDT3gen'=NA,
                               'MDT5ori'=NA, 'MDT5wei'=NA,
                               'MDT5max'=NA, 'MDT5gen'=NA,
                               'MDT7ori'=NA, 'MDT7wei'=NA,
                               'MDT7max'=NA, 'MDT7gen'=NA,
                               'MDT15ori'=NA, 'MDT15wei'=NA,
                               'MDT15max'=NA, 'MDT15gen'=NA,

                               'NNL1ori'=NA, 'NNL1wei'=NA,
                               'NNL1max'=NA, 'NNL1gen'=NA,
                               'NNL3ori'=NA, 'NNL3wei'=NA,
                               'NNL3max'=NA, 'NNL3gen'=NA,
                               'NNL5ori'=NA, 'NNL5wei'=NA,
                               'NNL5max'=NA, 'NNL5gen'=NA,
                               'NNL7ori'=NA, 'NNL7wei'=NA,
                               'NNL7max'=NA, 'NNL7gen'=NA,
                               'NNL15ori'=NA, 'NNL15wei'=NA,
                               'NNL15max'=NA, 'NNL15gen'=NA,

                               'MST1ori'=NA, 'MST1wei'=NA,
                               'MST1max'=NA, 'MST1gen'=NA,
                               'MST3ori'=NA, 'MST3wei'=NA,
                               'MST3max'=NA, 'MST3gen'=NA,
                               'MST5ori'=NA, 'MST5wei'=NA,
                               'MST5max'=NA, 'MST5gen'=NA,
                               'MST7ori'=NA, 'MST7wei'=NA,
                               'MST7max'=NA, 'MST7gen'=NA,
                               'MST15ori'=NA, 'MST15wei'=NA,
                               'MST15max'=NA, 'MST15gen'=NA
)

for(l in 2:length(ns)){

  cat(paste0('Length, ',ns[l],' (',trials,'): '))

  for(i in 1:trials){
    set.seed(123456 * (l+1) + i)
    cat(paste0(i,', '))
    sink(nullfile())    # now suppresses

    ## Generate Data
    data_cp <- gen_FD_KL_Expansion(
      ns = c(ns[l],1),
      eigsList = list(c(3,2,1,0.5),
                      c(3,2,1,0.5)),
      basesList = list(create.bspline.basis(nbasis=4, norder=4),
                       create.bspline.basis(nbasis=4, norder=4)),
      meansList = c(0,0),
      distsArray = c('Normal','Normal'),
      evals = evalPts,
      kappasArray = 0,
      silent = T)
    data_cp <- data_cp[,-ncol(data_cp)]
    data_dist <- calculateDistanceMatrix(data_cp, silent = T,
                                         errType='L2', dataUse = 'Orig')

    MDT1 <- createMDT(data_dist,1)
    MDT3 <- createMDT(data_dist,3)
    MDT5 <- createMDT(data_dist,5)
    MDT7 <- createMDT(data_dist,7)
    MDT15 <- createMDT(data_dist,15)
    #gSeg_MDT_tmp <- gseg1(ns[l], MDT3, pval.appr = T)
    gSeg_MDT1 <- gseg1_new(ns[l], as.matrix(MDT1), pval.perm = T, pval.appr = F)
    gSeg_MDT3 <- gseg1_new(ns[l], MDT3, pval.perm = T, pval.appr = F)
    gSeg_MDT5 <- gseg1_new(ns[l], MDT5, pval.perm = T, pval.appr = F)
    gSeg_MDT7 <- gseg1_new(ns[l], MDT7, pval.perm = T, pval.appr = F)
    if(ns[l]>30) {
      gSeg_MDT15 <- gseg1_new(ns[l], MDT15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_MDT15 <- list()
      gSeg_MDT15$pval.perm$ori$pval <- NA
      gSeg_MDT15$pval.perm$weighted$pval <- NA
      gSeg_MDT15$pval.perm$max.type$pval <- NA
      gSeg_MDT15$pval.perm$generalized$pval <- NA
    }


    NNL1 <- nnl1(data_dist,1)
    NNL3 <- nnl1(data_dist,3)
    NNL5 <- nnl1(data_dist,5)
    NNL7 <- nnl1(data_dist,7)
    NNL15 <- nnl1(data_dist,15)
    gSeg_NNL1 <- gseg1_new(ns[l], NNL1, pval.perm = T, pval.appr = F)
    gSeg_NNL3 <- gseg1_new(ns[l], NNL3, pval.perm = T, pval.appr = F)
    gSeg_NNL5 <- gseg1_new(ns[l], NNL5, pval.perm = T, pval.appr = F)
    gSeg_NNL7 <- gseg1_new(ns[l], NNL7, pval.perm = T, pval.appr = F)
    if(ns[l]>30) {
      gSeg_NNL15 <- gseg1_new(ns[l], NNL15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_NNL15 <- list()
      gSeg_NNL15$pval.perm$ori$pval <- NA
      gSeg_NNL15$pval.perm$weighted$pval <- NA
      gSeg_NNL15$pval.perm$max.type$pval <- NA
      gSeg_NNL15$pval.perm$generalized$pval <- NA
    }

    MST1 <- mstree(as.dist(data_dist), ngmax=1)
    MST3 <- mstree(as.dist(data_dist), ngmax=3)
    MST5 <- mstree(as.dist(data_dist), ngmax=5)
    MST7 <- mstree(as.dist(data_dist), ngmax=7)
    MST15 <- mstree(as.dist(data_dist), ngmax=15)
    gSeg_MST1 <- gseg1_new(ns[l], MST1, pval.perm = T, pval.appr = F)
    gSeg_MST3 <- gseg1_new(ns[l], MST3, pval.perm = T, pval.appr = F)
    gSeg_MST5 <- gseg1_new(ns[l], MST5, pval.perm = T, pval.appr = F)
    gSeg_MST7 <- gseg1_new(ns[l], MST7, pval.perm = T, pval.appr = F)
    if(ns[l]>30) {
      gSeg_MST15 <- gseg1_new(ns[l], MST15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_MST15 <- list()
      gSeg_MST15$pval.perm$ori$pval <- NA
      gSeg_MST15$pval.perm$weighted$pval <- NA
      gSeg_MST15$pval.perm$max.type$pval <- NA
      gSeg_MST15$pval.perm$generalized$pval <- NA
    }

    results[i, ] <- c(
      # MDT 1
      gSeg_MDT1$pval.perm$ori$pval,
      gSeg_MDT1$pval.perm$weighted$pval,
      gSeg_MDT1$pval.perm$max.type$pval,
      gSeg_MDT1$pval.perm$generalized$pval,
      # MDT 3
      gSeg_MDT3$pval.perm$ori$pval,
      gSeg_MDT3$pval.perm$weighted$pval,
      gSeg_MDT3$pval.perm$max.type$pval,
      gSeg_MDT3$pval.perm$generalized$pval,
      # MDT 5
      gSeg_MDT5$pval.perm$ori$pval,
      gSeg_MDT5$pval.perm$weighted$pval,
      gSeg_MDT5$pval.perm$max.type$pval,
      gSeg_MDT5$pval.perm$generalized$pval,
      # MDT 7
      gSeg_MDT7$pval.perm$ori$pval,
      gSeg_MDT7$pval.perm$weighted$pval,
      gSeg_MDT7$pval.perm$max.type$pval,
      gSeg_MDT7$pval.perm$generalized$pval,
      # MDT 15
      gSeg_MDT15$pval.perm$ori$pval,
      gSeg_MDT15$pval.perm$weighted$pval,
      gSeg_MDT15$pval.perm$max.type$pval,
      gSeg_MDT15$pval.perm$generalized$pval,

      # NNL 1
      gSeg_NNL1$pval.perm$ori$pval,
      gSeg_NNL1$pval.perm$weighted$pval,
      gSeg_NNL1$pval.perm$max.type$pval,
      gSeg_NNL1$pval.perm$generalized$pval,
      # NNL 3
      gSeg_NNL3$pval.perm$ori$pval,
      gSeg_NNL3$pval.perm$weighted$pval,
      gSeg_NNL3$pval.perm$max.type$pval,
      gSeg_NNL3$pval.perm$generalized$pval,
      # NNL 5
      gSeg_NNL5$pval.perm$ori$pval,
      gSeg_NNL5$pval.perm$weighted$pval,
      gSeg_NNL5$pval.perm$max.type$pval,
      gSeg_NNL5$pval.perm$generalized$pval,
      # NNL 7
      gSeg_NNL7$pval.perm$ori$pval,
      gSeg_NNL7$pval.perm$weighted$pval,
      gSeg_NNL7$pval.perm$max.type$pval,
      gSeg_NNL7$pval.perm$generalized$pval,
      # NNL 15
      gSeg_NNL15$pval.perm$ori$pval,
      gSeg_NNL15$pval.perm$weighted$pval,
      gSeg_NNL15$pval.perm$max.type$pval,
      gSeg_NNL15$pval.perm$generalized$pval,

      # MST 1
      gSeg_MST1$pval.perm$ori$pval,
      gSeg_MST1$pval.perm$weighted$pval,
      gSeg_MST1$pval.perm$max.type$pval,
      gSeg_MST1$pval.perm$generalized$pval,
      # MST 3
      gSeg_MST3$pval.perm$ori$pval,
      gSeg_MST3$pval.perm$weighted$pval,
      gSeg_MST3$pval.perm$max.type$pval,
      gSeg_MST3$pval.perm$generalized$pval,
      # MST 5
      gSeg_MST5$pval.perm$ori$pval,
      gSeg_MST5$pval.perm$weighted$pval,
      gSeg_MST5$pval.perm$max.type$pval,
      gSeg_MST5$pval.perm$generalized$pval,
      # MST 7
      gSeg_MST7$pval.perm$ori$pval,
      gSeg_MST7$pval.perm$weighted$pval,
      gSeg_MST7$pval.perm$max.type$pval,
      gSeg_MST7$pval.perm$generalized$pval,
      # MST 15
      gSeg_MST15$pval.perm$ori$pval,
      gSeg_MST15$pval.perm$weighted$pval,
      gSeg_MST15$pval.perm$max.type$pval,
      gSeg_MST15$pval.perm$generalized$pval
    )
    sink()
    # saveRDS(results, paste0('C:/Users/jerem/OneDrive/Documents/',
    #                         'School/Waterloo/Research/',
    #                         'GraphCP/Data/Simulations/',
    #                         'tmp_t',trials,'_l',l,'.rds'))
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(ns[l], as.numeric(colSums(results>0.05,na.rm = T)/colSums(!is.na(results))))

  # saveRDS(results_save, paste0('C:/Users/jerem/OneDrive/Documents/',
  #                              'School/Waterloo/Research/',
  #                              'GraphCP/Data/Simulations/',
  #                              'long_results.rds'))
  # saveRDS(results_simplify, paste0('C:/Users/jerem/OneDrive/Documents/',
  #                                  'School/Waterloo/Research/',
  #                                  'GraphCP/Data/Simulations/',
  #                                  'long_simplify.rds'))
}

saveRDS(results_save,
        file=paste0(path,"\\results\\null_full.rds") )
saveRDS(results_simplify,
        file=paste0(path,"\\results\\null_simple.rds") )

lens <- unique(results_simplify$Length)
for(i in 1:length(lens)){
  tmp <- t(results_simplify[results_simplify$Length==lens[i],])[-1,,drop=FALSE]
  tmp <-
    cbind(tmp,
          data.frame('tree'=substr(rownames(tmp),1,3),
                     'stat'=stringr::str_replace(substr(rownames(tmp),5,8),'[0-9]',''),
                     'treenum'=paste0(substr(rownames(tmp),1,3),stringr::str_extract(rownames(tmp),'[0-9]?[0-9]')),
                     'all'=rownames(tmp)
          )
    )
  colnames(tmp) <- c('value','tree','stat','treenum','all')
  tmp$value <- 1-as.numeric(tmp$value)
  res0 <- tmp %>% pivot_wider(names_from = stat,values_from = value,id_cols = treenum)
  if(i==1){
    res <- res0
  }else if(i==4){
    res1 <- res0
  }else if(i<4){
    res <- cbind(res, res0[,-1])
  }else{
    res1 <- cbind(res1, res0[,-1])
  }
}


for(i in 1:nrow(res)){
  tmp_prt <- res[i,-1]
  for(j in 1:length(tmp_prt)){
    if(is.nan(tmp_prt[[j]]))
      tmp_prt[[j]] <- 999
  }

  cat(res[i,1],' & ',paste(round(tmp_prt,3), collapse=" & "),'\\\\ \n')
  cat('% \\hline \n')

}

for(i in 1:nrow(res1)){
  tmp_prt <- res1[i,-1]
  for(j in 1:length(tmp_prt)){
    if(is.nan(tmp_prt[[j]]))
      tmp_prt[[j]] <- 999
  }

  cat(res1[i,1],' & ',paste(round(tmp_prt,3), collapse=" & "),'\\\\ \n')
  cat('% \\hline \n')

}
