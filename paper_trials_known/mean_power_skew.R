library(devtools)
load_all()
library(fda)
library(ade4)
library(gSeg)
library(ggplot2)

path <- 'C://Users//jerem//OneDrive//Documents//School//Waterloo//Research//RPackages//fungraphs//paper_trials_known'

trials <- 1000
evalPts <- seq(0,1,0.05)
n <- 50
break_location <- 0.5

results <- data.frame(# 'MDT1ori'=1:trials, 'MDT1wei'=NA,
                      # 'MDT1max'=NA, 'MDT1gen'=NA,
                      # 'MDT3ori'=NA, 'MDT3wei'=NA,
                      # 'MDT3max'=NA, 'MDT3gen'=NA,
                      # 'MDT5ori'=NA, 'MDT5wei'=NA,
                      # 'MDT5max'=NA, 'MDT5gen'=NA,
                      # 'MDT7ori'=NA, 'MDT7wei'=NA,
                      # 'MDT7max'=NA, 'MDT7gen'=NA,
                      'MDT15ori'=1:trials, 'MDT15wei'=NA,
                      'MDT15max'=NA, 'MDT15gen'=NA,

                      # 'NNL1ori'=NA, 'NNL1wei'=NA,
                      # 'NNL1max'=NA, 'NNL1gen'=NA,
                      # 'NNL3ori'=NA, 'NNL3wei'=NA,
                      # 'NNL3max'=NA, 'NNL3gen'=NA,
                      # 'NNL5ori'=NA, 'NNL5wei'=NA,
                      # 'NNL5max'=NA, 'NNL5gen'=NA,
                      # 'NNL7ori'=NA, 'NNL7wei'=NA,
                      # 'NNL7max'=NA, 'NNL7gen'=NA,
                      'NNL15ori'=NA, 'NNL15wei'=NA,
                      'NNL15max'=NA, 'NNL15gen'=NA,

                      # 'MST1ori'=NA, 'MST1wei'=NA,
                      # 'MST1max'=NA, 'MST1gen'=NA,
                      # 'MST3ori'=NA, 'MST3wei'=NA,
                      # 'MST3max'=NA, 'MST3gen'=NA,
                      # 'MST5ori'=NA, 'MST5wei'=NA,
                      # 'MST5max'=NA, 'MST5gen'=NA,
                      # 'MST7ori'=NA, 'MST7wei'=NA,
                      # 'MST7max'=NA, 'MST7gen'=NA,
                      'MST15ori'=NA, 'MST15wei'=NA,
                      'MST15max'=NA, 'MST15gen'=NA,

                      'mean'=NA
)

results_save <- list()
means = c(0, 0.05,0.1,0.15,0.25,0.5,0.75,1)
results_simplify <- data.frame('mean'=means,
                               # 'MDT1ori'=NA, 'MDT1wei'=NA,
                               # 'MDT1max'=NA, 'MDT1gen'=NA,
                               # 'MDT3ori'=NA, 'MDT3wei'=NA,
                               # 'MDT3max'=NA, 'MDT3gen'=NA,
                               # 'MDT5ori'=NA, 'MDT5wei'=NA,
                               # 'MDT5max'=NA, 'MDT5gen'=NA,
                               # 'MDT7ori'=NA, 'MDT7wei'=NA,
                               # 'MDT7max'=NA, 'MDT7gen'=NA,
                               'MDT15ori'=NA, 'MDT15wei'=NA,
                               'MDT15max'=NA, 'MDT15gen'=NA,

                               # 'NNL1ori'=NA, 'NNL1wei'=NA,
                               # 'NNL1max'=NA, 'NNL1gen'=NA,
                               # 'NNL3ori'=NA, 'NNL3wei'=NA,
                               # 'NNL3max'=NA, 'NNL3gen'=NA,
                               # 'NNL5ori'=NA, 'NNL5wei'=NA,
                               # 'NNL5max'=NA, 'NNL5gen'=NA,
                               # 'NNL7ori'=NA, 'NNL7wei'=NA,
                               # 'NNL7max'=NA, 'NNL7gen'=NA,
                               'NNL15ori'=NA, 'NNL15wei'=NA,
                               'NNL15max'=NA, 'NNL15gen'=NA,

                               # 'MST1ori'=NA, 'MST1wei'=NA,
                               # 'MST1max'=NA, 'MST1gen'=NA,
                               # 'MST3ori'=NA, 'MST3wei'=NA,
                               # 'MST3max'=NA, 'MST3gen'=NA,
                               # 'MST5ori'=NA, 'MST5wei'=NA,
                               # 'MST5max'=NA, 'MST5gen'=NA,
                               # 'MST7ori'=NA, 'MST7wei'=NA,
                               # 'MST7max'=NA, 'MST7gen'=NA,
                               'MST15ori'=NA, 'MST15wei'=NA,
                               'MST15max'=NA, 'MST15gen'=NA,

                               'mean'=NA
)

data <- data_dist <- list()


for(l in 1:length(means)){

  cat(paste0('Mean, ',means[l],' (',trials,'): '))

  for(i in 1:trials){
    set.seed(123456 * (l+1) + i + 1)
    cat(paste0(i,', '))
    sink(nullfile())    # now suppresses

    ## Generate Data
    data[[i]] <- gen_FD_KL_Expansion(
      ns = c(n*break_location,n*(1-break_location)),
      eigsList = list(c(3,2,1,0.5),
                      c(3,2,1,0.5)),
      basesList = list(create.bspline.basis(nbasis=4, norder=4),
                       create.bspline.basis(nbasis=4, norder=4)),
      meansList = c(0,means[l]),
      distsArray = c('sn','sn'),
      evals = evalPts,
      kappasArray = 0,
      silent = T)
    data_dist[[i]] <- calculateDistanceMatrix(data[[i]], silent = T,
                                              errType='L2', dataUse = 'Orig')

    # MDT1 <- createMDT(data_dist[[i]],1)
    # MDT3 <- createMDT(data_dist[[i]],3)
    # MDT5 <- createMDT(data_dist[[i]],5)
    # MDT7 <- createMDT(data_dist[[i]],7)
    MDT15 <- createMDT(data_dist[[i]],15)
    # gSeg_MDT1 <- gseg1_new(n, MDT1, pval.perm = T, pval.appr = F)
    # gSeg_MDT3 <- gseg1_new(n, MDT3, pval.perm = T, pval.appr = F)
    # gSeg_MDT5 <- gseg1_new(n, MDT5, pval.perm = T, pval.appr = F)
    # gSeg_MDT7 <- gseg1_new(n, MDT7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_MDT15 <- list()
      gSeg_MDT15$pval.perm$ori$pval <- NA
      gSeg_MDT15$pval.perm$weighted$pval <- NA
      gSeg_MDT15$pval.perm$max.type$pval <- NA
      gSeg_MDT15$pval.perm$generalized$pval <- NA
    }


    # NNL1 <- nnl1(data_dist[[i]],1)
    # NNL3 <- nnl1(data_dist[[i]],3)
    # NNL5 <- nnl1(data_dist[[i]],5)
    # NNL7 <- nnl1(data_dist[[i]],7)
    NNL15 <- nnl1(data_dist[[i]],15)
    # gSeg_NNL1 <- gseg1_new(n, NNL1, pval.perm = T, pval.appr = F)
    # gSeg_NNL3 <- gseg1_new(n, NNL3, pval.perm = T, pval.appr = F)
    # gSeg_NNL5 <- gseg1_new(n, NNL5, pval.perm = T, pval.appr = F)
    # gSeg_NNL7 <- gseg1_new(n, NNL7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_NNL15 <- gseg1_new(n, NNL15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_NNL15 <- list()
      gSeg_NNL15$pval.perm$ori$pval <- NA
      gSeg_NNL15$pval.perm$weighted$pval <- NA
      gSeg_NNL15$pval.perm$max.type$pval <- NA
      gSeg_NNL15$pval.perm$generalized$pval <- NA
    }

    # MST1 <- mstree(as.dist(data_dist[[i]]), ngmax=1)
    # MST3 <- mstree(as.dist(data_dist[[i]]), ngmax=3)
    # MST5 <- mstree(as.dist(data_dist[[i]]), ngmax=5)
    # MST7 <- mstree(as.dist(data_dist[[i]]), ngmax=7)
    MST15 <- mstree(as.dist(data_dist[[i]]), ngmax=15)
    # gSeg_MST1 <- gseg1_new(n, MST1, pval.perm = T, pval.appr = F)
    # gSeg_MST3 <- gseg1_new(n, MST3, pval.perm = T, pval.appr = F)
    # gSeg_MST5 <- gseg1_new(n, MST5, pval.perm = T, pval.appr = F)
    # gSeg_MST7 <- gseg1_new(n, MST7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_MST15 <- gseg1_new(n, MST15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_MST15 <- list()
      gSeg_MST15$pval.perm$ori$pval <- NA
      gSeg_MST15$pval.perm$weighted$pval <- NA
      gSeg_MST15$pval.perm$max.type$pval <- NA
      gSeg_MST15$pval.perm$generalized$pval <- NA
    }

    results[i, ] <- c(
      # # MDT 1
      # gSeg_MDT1$pval.perm$ori$pval,
      # gSeg_MDT1$pval.perm$weighted$pval,
      # gSeg_MDT1$pval.perm$max.type$pval,
      # gSeg_MDT1$pval.perm$generalized$pval,
      # # MDT 3
      # gSeg_MDT3$pval.perm$ori$pval,
      # gSeg_MDT3$pval.perm$weighted$pval,
      # gSeg_MDT3$pval.perm$max.type$pval,
      # gSeg_MDT3$pval.perm$generalized$pval,
      # # MDT 5
      # gSeg_MDT5$pval.perm$ori$pval,
      # gSeg_MDT5$pval.perm$weighted$pval,
      # gSeg_MDT5$pval.perm$max.type$pval,
      # gSeg_MDT5$pval.perm$generalized$pval,
      # # MDT 7
      # gSeg_MDT7$pval.perm$ori$pval,
      # gSeg_MDT7$pval.perm$weighted$pval,
      # gSeg_MDT7$pval.perm$max.type$pval,
      # gSeg_MDT7$pval.perm$generalized$pval,
      # MDT 15
      gSeg_MDT15$pval.perm$ori$pval,
      gSeg_MDT15$pval.perm$weighted$pval,
      gSeg_MDT15$pval.perm$max.type$pval,
      gSeg_MDT15$pval.perm$generalized$pval,

      # # NNL 1
      # gSeg_NNL1$pval.perm$ori$pval,
      # gSeg_NNL1$pval.perm$weighted$pval,
      # gSeg_NNL1$pval.perm$max.type$pval,
      # gSeg_NNL1$pval.perm$generalized$pval,
      # # NNL 3
      # gSeg_NNL3$pval.perm$ori$pval,
      # gSeg_NNL3$pval.perm$weighted$pval,
      # gSeg_NNL3$pval.perm$max.type$pval,
      # gSeg_NNL3$pval.perm$generalized$pval,
      # # NNL 5
      # gSeg_NNL5$pval.perm$ori$pval,
      # gSeg_NNL5$pval.perm$weighted$pval,
      # gSeg_NNL5$pval.perm$max.type$pval,
      # gSeg_NNL5$pval.perm$generalized$pval,
      # # NNL 7
      # gSeg_NNL7$pval.perm$ori$pval,
      # gSeg_NNL7$pval.perm$weighted$pval,
      # gSeg_NNL7$pval.perm$max.type$pval,
      # gSeg_NNL7$pval.perm$generalized$pval,
      # NNL 15
      gSeg_NNL15$pval.perm$ori$pval,
      gSeg_NNL15$pval.perm$weighted$pval,
      gSeg_NNL15$pval.perm$max.type$pval,
      gSeg_NNL15$pval.perm$generalized$pval,

      # # MST 1
      # gSeg_MST1$pval.perm$ori$pval,
      # gSeg_MST1$pval.perm$weighted$pval,
      # gSeg_MST1$pval.perm$max.type$pval,
      # gSeg_MST1$pval.perm$generalized$pval,
      # # MST 3
      # gSeg_MST3$pval.perm$ori$pval,
      # gSeg_MST3$pval.perm$weighted$pval,
      # gSeg_MST3$pval.perm$max.type$pval,
      # gSeg_MST3$pval.perm$generalized$pval,
      # # MST 5
      # gSeg_MST5$pval.perm$ori$pval,
      # gSeg_MST5$pval.perm$weighted$pval,
      # gSeg_MST5$pval.perm$max.type$pval,
      # gSeg_MST5$pval.perm$generalized$pval,
      # # MST 7
      # gSeg_MST7$pval.perm$ori$pval,
      # gSeg_MST7$pval.perm$weighted$pval,
      # gSeg_MST7$pval.perm$max.type$pval,
      # gSeg_MST7$pval.perm$generalized$pval,
      # MST 15
      gSeg_MST15$pval.perm$ori$pval,
      gSeg_MST15$pval.perm$weighted$pval,
      gSeg_MST15$pval.perm$max.type$pval,
      gSeg_MST15$pval.perm$generalized$pval,

      ifelse(is.na(mean_change(data[[i]])),0,1)
    )
    sink()
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(means[l], as.numeric(colSums(results<=0.05,na.rm = T)/colSums(!is.na(results))))

  saveRDS(results_save,
          file=paste0(path,"\\results\\skew_change_full.rds") )
  saveRDS(results_simplify,
          file=paste0(path,"\\results\\skew_change_simple.rds") )
}

######################

results <- data.frame(# 'MDT1ori'=1:trials, 'MDT1wei'=NA,
  # 'MDT1max'=NA, 'MDT1gen'=NA,
  # 'MDT3ori'=NA, 'MDT3wei'=NA,
  # 'MDT3max'=NA, 'MDT3gen'=NA,
  # 'MDT5ori'=NA, 'MDT5wei'=NA,
  # 'MDT5max'=NA, 'MDT5gen'=NA,
  # 'MDT7ori'=NA, 'MDT7wei'=NA,
  # 'MDT7max'=NA, 'MDT7gen'=NA,
  'MDT15ori'=1:trials, 'MDT15wei'=NA,
  'MDT15max'=NA, 'MDT15gen'=NA,

  # 'NNL1ori'=NA, 'NNL1wei'=NA,
  # 'NNL1max'=NA, 'NNL1gen'=NA,
  # 'NNL3ori'=NA, 'NNL3wei'=NA,
  # 'NNL3max'=NA, 'NNL3gen'=NA,
  # 'NNL5ori'=NA, 'NNL5wei'=NA,
  # 'NNL5max'=NA, 'NNL5gen'=NA,
  # 'NNL7ori'=NA, 'NNL7wei'=NA,
  # 'NNL7max'=NA, 'NNL7gen'=NA,
  'NNL15ori'=NA, 'NNL15wei'=NA,
  'NNL15max'=NA, 'NNL15gen'=NA,

  # 'MST1ori'=NA, 'MST1wei'=NA,
  # 'MST1max'=NA, 'MST1gen'=NA,
  # 'MST3ori'=NA, 'MST3wei'=NA,
  # 'MST3max'=NA, 'MST3gen'=NA,
  # 'MST5ori'=NA, 'MST5wei'=NA,
  # 'MST5max'=NA, 'MST5gen'=NA,
  # 'MST7ori'=NA, 'MST7wei'=NA,
  # 'MST7max'=NA, 'MST7gen'=NA,
  'MST15ori'=NA, 'MST15wei'=NA,
  'MST15max'=NA, 'MST15gen'=NA,

  'mean'=NA
)

results_save <- list()
means = c(0, 0.05,0.1,0.15,0.25,0.5,0.75,1)
results_simplify <- data.frame('mean'=means,
                               # 'MDT1ori'=NA, 'MDT1wei'=NA,
                               # 'MDT1max'=NA, 'MDT1gen'=NA,
                               # 'MDT3ori'=NA, 'MDT3wei'=NA,
                               # 'MDT3max'=NA, 'MDT3gen'=NA,
                               # 'MDT5ori'=NA, 'MDT5wei'=NA,
                               # 'MDT5max'=NA, 'MDT5gen'=NA,
                               # 'MDT7ori'=NA, 'MDT7wei'=NA,
                               # 'MDT7max'=NA, 'MDT7gen'=NA,
                               'MDT15ori'=NA, 'MDT15wei'=NA,
                               'MDT15max'=NA, 'MDT15gen'=NA,

                               # 'NNL1ori'=NA, 'NNL1wei'=NA,
                               # 'NNL1max'=NA, 'NNL1gen'=NA,
                               # 'NNL3ori'=NA, 'NNL3wei'=NA,
                               # 'NNL3max'=NA, 'NNL3gen'=NA,
                               # 'NNL5ori'=NA, 'NNL5wei'=NA,
                               # 'NNL5max'=NA, 'NNL5gen'=NA,
                               # 'NNL7ori'=NA, 'NNL7wei'=NA,
                               # 'NNL7max'=NA, 'NNL7gen'=NA,
                               'NNL15ori'=NA, 'NNL15wei'=NA,
                               'NNL15max'=NA, 'NNL15gen'=NA,

                               # 'MST1ori'=NA, 'MST1wei'=NA,
                               # 'MST1max'=NA, 'MST1gen'=NA,
                               # 'MST3ori'=NA, 'MST3wei'=NA,
                               # 'MST3max'=NA, 'MST3gen'=NA,
                               # 'MST5ori'=NA, 'MST5wei'=NA,
                               # 'MST5max'=NA, 'MST5gen'=NA,
                               # 'MST7ori'=NA, 'MST7wei'=NA,
                               # 'MST7max'=NA, 'MST7gen'=NA,
                               'MST15ori'=NA, 'MST15wei'=NA,
                               'MST15max'=NA, 'MST15gen'=NA,

                               'mean'=NA
)

data <- data_dist <- list()


for(l in 1:length(means)){

  cat(paste0('Mean, ',means[l],' (',trials,'): '))

  for(i in 1:trials){
    set.seed(123456 * (l+1) + i + 1)
    cat(paste0(i,', '))
    sink(nullfile())    # now suppresses

    ## Generate Data
    data[[i]] <- gen_FD_KL_Expansion(
      ns = c(n*break_location,n*(1-break_location)),
      eigsList = list(c(3,2,1,0.5),
                      c(3,2,1,0.5)),
      basesList = list(create.bspline.basis(nbasis=4, norder=4),
                       create.bspline.basis(nbasis=4, norder=4)),
      meansList = c(0,means[l]),
      distsArray = c('Normal','Normal'),
      evals = evalPts,
      silent = T)
    data_dist[[i]] <- calculateDistanceMatrix(data[[i]], silent = T,
                                              errType='L2', dataUse = 'Orig')

    # MDT1 <- createMDT(data_dist[[i]],1)
    # MDT3 <- createMDT(data_dist[[i]],3)
    # MDT5 <- createMDT(data_dist[[i]],5)
    # MDT7 <- createMDT(data_dist[[i]],7)
    MDT15 <- createMDT(data_dist[[i]],15)
    # gSeg_MDT1 <- gseg1_new(n, MDT1, pval.perm = T, pval.appr = F)
    # gSeg_MDT3 <- gseg1_new(n, MDT3, pval.perm = T, pval.appr = F)
    # gSeg_MDT5 <- gseg1_new(n, MDT5, pval.perm = T, pval.appr = F)
    # gSeg_MDT7 <- gseg1_new(n, MDT7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_MDT15 <- list()
      gSeg_MDT15$pval.perm$ori$pval <- NA
      gSeg_MDT15$pval.perm$weighted$pval <- NA
      gSeg_MDT15$pval.perm$max.type$pval <- NA
      gSeg_MDT15$pval.perm$generalized$pval <- NA
    }


    # NNL1 <- nnl1(data_dist[[i]],1)
    # NNL3 <- nnl1(data_dist[[i]],3)
    # NNL5 <- nnl1(data_dist[[i]],5)
    # NNL7 <- nnl1(data_dist[[i]],7)
    NNL15 <- nnl1(data_dist[[i]],15)
    # gSeg_NNL1 <- gseg1_new(n, NNL1, pval.perm = T, pval.appr = F)
    # gSeg_NNL3 <- gseg1_new(n, NNL3, pval.perm = T, pval.appr = F)
    # gSeg_NNL5 <- gseg1_new(n, NNL5, pval.perm = T, pval.appr = F)
    # gSeg_NNL7 <- gseg1_new(n, NNL7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_NNL15 <- gseg1_new(n, NNL15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_NNL15 <- list()
      gSeg_NNL15$pval.perm$ori$pval <- NA
      gSeg_NNL15$pval.perm$weighted$pval <- NA
      gSeg_NNL15$pval.perm$max.type$pval <- NA
      gSeg_NNL15$pval.perm$generalized$pval <- NA
    }

    # MST1 <- mstree(as.dist(data_dist[[i]]), ngmax=1)
    # MST3 <- mstree(as.dist(data_dist[[i]]), ngmax=3)
    # MST5 <- mstree(as.dist(data_dist[[i]]), ngmax=5)
    # MST7 <- mstree(as.dist(data_dist[[i]]), ngmax=7)
    MST15 <- mstree(as.dist(data_dist[[i]]), ngmax=15)
    # gSeg_MST1 <- gseg1_new(n, MST1, pval.perm = T, pval.appr = F)
    # gSeg_MST3 <- gseg1_new(n, MST3, pval.perm = T, pval.appr = F)
    # gSeg_MST5 <- gseg1_new(n, MST5, pval.perm = T, pval.appr = F)
    # gSeg_MST7 <- gseg1_new(n, MST7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_MST15 <- gseg1_new(n, MST15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_MST15 <- list()
      gSeg_MST15$pval.perm$ori$pval <- NA
      gSeg_MST15$pval.perm$weighted$pval <- NA
      gSeg_MST15$pval.perm$max.type$pval <- NA
      gSeg_MST15$pval.perm$generalized$pval <- NA
    }

    results[i, ] <- c(
      # # MDT 1
      # gSeg_MDT1$pval.perm$ori$pval,
      # gSeg_MDT1$pval.perm$weighted$pval,
      # gSeg_MDT1$pval.perm$max.type$pval,
      # gSeg_MDT1$pval.perm$generalized$pval,
      # # MDT 3
      # gSeg_MDT3$pval.perm$ori$pval,
      # gSeg_MDT3$pval.perm$weighted$pval,
      # gSeg_MDT3$pval.perm$max.type$pval,
      # gSeg_MDT3$pval.perm$generalized$pval,
      # # MDT 5
      # gSeg_MDT5$pval.perm$ori$pval,
      # gSeg_MDT5$pval.perm$weighted$pval,
      # gSeg_MDT5$pval.perm$max.type$pval,
      # gSeg_MDT5$pval.perm$generalized$pval,
      # # MDT 7
      # gSeg_MDT7$pval.perm$ori$pval,
      # gSeg_MDT7$pval.perm$weighted$pval,
      # gSeg_MDT7$pval.perm$max.type$pval,
      # gSeg_MDT7$pval.perm$generalized$pval,
      # MDT 15
      gSeg_MDT15$pval.perm$ori$pval,
      gSeg_MDT15$pval.perm$weighted$pval,
      gSeg_MDT15$pval.perm$max.type$pval,
      gSeg_MDT15$pval.perm$generalized$pval,

      # # NNL 1
      # gSeg_NNL1$pval.perm$ori$pval,
      # gSeg_NNL1$pval.perm$weighted$pval,
      # gSeg_NNL1$pval.perm$max.type$pval,
      # gSeg_NNL1$pval.perm$generalized$pval,
      # # NNL 3
      # gSeg_NNL3$pval.perm$ori$pval,
      # gSeg_NNL3$pval.perm$weighted$pval,
      # gSeg_NNL3$pval.perm$max.type$pval,
      # gSeg_NNL3$pval.perm$generalized$pval,
      # # NNL 5
      # gSeg_NNL5$pval.perm$ori$pval,
      # gSeg_NNL5$pval.perm$weighted$pval,
      # gSeg_NNL5$pval.perm$max.type$pval,
      # gSeg_NNL5$pval.perm$generalized$pval,
      # # NNL 7
      # gSeg_NNL7$pval.perm$ori$pval,
      # gSeg_NNL7$pval.perm$weighted$pval,
      # gSeg_NNL7$pval.perm$max.type$pval,
      # gSeg_NNL7$pval.perm$generalized$pval,
      # NNL 15
      gSeg_NNL15$pval.perm$ori$pval,
      gSeg_NNL15$pval.perm$weighted$pval,
      gSeg_NNL15$pval.perm$max.type$pval,
      gSeg_NNL15$pval.perm$generalized$pval,

      # # MST 1
      # gSeg_MST1$pval.perm$ori$pval,
      # gSeg_MST1$pval.perm$weighted$pval,
      # gSeg_MST1$pval.perm$max.type$pval,
      # gSeg_MST1$pval.perm$generalized$pval,
      # # MST 3
      # gSeg_MST3$pval.perm$ori$pval,
      # gSeg_MST3$pval.perm$weighted$pval,
      # gSeg_MST3$pval.perm$max.type$pval,
      # gSeg_MST3$pval.perm$generalized$pval,
      # # MST 5
      # gSeg_MST5$pval.perm$ori$pval,
      # gSeg_MST5$pval.perm$weighted$pval,
      # gSeg_MST5$pval.perm$max.type$pval,
      # gSeg_MST5$pval.perm$generalized$pval,
      # # MST 7
      # gSeg_MST7$pval.perm$ori$pval,
      # gSeg_MST7$pval.perm$weighted$pval,
      # gSeg_MST7$pval.perm$max.type$pval,
      # gSeg_MST7$pval.perm$generalized$pval,
      # MST 15
      gSeg_MST15$pval.perm$ori$pval,
      gSeg_MST15$pval.perm$weighted$pval,
      gSeg_MST15$pval.perm$max.type$pval,
      gSeg_MST15$pval.perm$generalized$pval,

      ifelse(is.na(mean_change(data[[i]])),0,1)
    )
    sink()
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(means[l], as.numeric(colSums(results<=0.05,na.rm = T)/colSums(!is.na(results))))

  saveRDS(results_save,
          file=paste0(path,"\\results\\normal4skew_change_full.rds") )
  saveRDS(results_simplify,
          file=paste0(path,"\\results\\normal4skew_change_simple.rds") )
}


###############################
results <- data.frame(# 'MDT1ori'=1:trials, 'MDT1wei'=NA,
  # 'MDT1max'=NA, 'MDT1gen'=NA,
  # 'MDT3ori'=NA, 'MDT3wei'=NA,
  # 'MDT3max'=NA, 'MDT3gen'=NA,
  # 'MDT5ori'=NA, 'MDT5wei'=NA,
  # 'MDT5max'=NA, 'MDT5gen'=NA,
  # 'MDT7ori'=NA, 'MDT7wei'=NA,
  # 'MDT7max'=NA, 'MDT7gen'=NA,
  'MDT15ori'=1:trials, 'MDT15wei'=NA,
  'MDT15max'=NA, 'MDT15gen'=NA,

  # 'NNL1ori'=NA, 'NNL1wei'=NA,
  # 'NNL1max'=NA, 'NNL1gen'=NA,
  # 'NNL3ori'=NA, 'NNL3wei'=NA,
  # 'NNL3max'=NA, 'NNL3gen'=NA,
  # 'NNL5ori'=NA, 'NNL5wei'=NA,
  # 'NNL5max'=NA, 'NNL5gen'=NA,
  # 'NNL7ori'=NA, 'NNL7wei'=NA,
  # 'NNL7max'=NA, 'NNL7gen'=NA,
  'NNL15ori'=NA, 'NNL15wei'=NA,
  'NNL15max'=NA, 'NNL15gen'=NA,

  # 'MST1ori'=NA, 'MST1wei'=NA,
  # 'MST1max'=NA, 'MST1gen'=NA,
  # 'MST3ori'=NA, 'MST3wei'=NA,
  # 'MST3max'=NA, 'MST3gen'=NA,
  # 'MST5ori'=NA, 'MST5wei'=NA,
  # 'MST5max'=NA, 'MST5gen'=NA,
  # 'MST7ori'=NA, 'MST7wei'=NA,
  # 'MST7max'=NA, 'MST7gen'=NA,
  'MST15ori'=NA, 'MST15wei'=NA,
  'MST15max'=NA, 'MST15gen'=NA,

  'mean'=NA
)

results_save <- list()
means = c(0, 0.05,0.1,0.15,0.25,0.5,0.75,1)
results_simplify <- data.frame('mean'=means,
                               # 'MDT1ori'=NA, 'MDT1wei'=NA,
                               # 'MDT1max'=NA, 'MDT1gen'=NA,
                               # 'MDT3ori'=NA, 'MDT3wei'=NA,
                               # 'MDT3max'=NA, 'MDT3gen'=NA,
                               # 'MDT5ori'=NA, 'MDT5wei'=NA,
                               # 'MDT5max'=NA, 'MDT5gen'=NA,
                               # 'MDT7ori'=NA, 'MDT7wei'=NA,
                               # 'MDT7max'=NA, 'MDT7gen'=NA,
                               'MDT15ori'=NA, 'MDT15wei'=NA,
                               'MDT15max'=NA, 'MDT15gen'=NA,

                               # 'NNL1ori'=NA, 'NNL1wei'=NA,
                               # 'NNL1max'=NA, 'NNL1gen'=NA,
                               # 'NNL3ori'=NA, 'NNL3wei'=NA,
                               # 'NNL3max'=NA, 'NNL3gen'=NA,
                               # 'NNL5ori'=NA, 'NNL5wei'=NA,
                               # 'NNL5max'=NA, 'NNL5gen'=NA,
                               # 'NNL7ori'=NA, 'NNL7wei'=NA,
                               # 'NNL7max'=NA, 'NNL7gen'=NA,
                               'NNL15ori'=NA, 'NNL15wei'=NA,
                               'NNL15max'=NA, 'NNL15gen'=NA,

                               # 'MST1ori'=NA, 'MST1wei'=NA,
                               # 'MST1max'=NA, 'MST1gen'=NA,
                               # 'MST3ori'=NA, 'MST3wei'=NA,
                               # 'MST3max'=NA, 'MST3gen'=NA,
                               # 'MST5ori'=NA, 'MST5wei'=NA,
                               # 'MST5max'=NA, 'MST5gen'=NA,
                               # 'MST7ori'=NA, 'MST7wei'=NA,
                               # 'MST7max'=NA, 'MST7gen'=NA,
                               'MST15ori'=NA, 'MST15wei'=NA,
                               'MST15max'=NA, 'MST15gen'=NA,

                               'mean'=NA
)

data <- data_dist <- list()


for(l in 1:length(means)){

  cat(paste0('Mean, ',means[l],' (',trials,'): '))

  for(i in 1:trials){
    set.seed(123456 * (l+1) + i + 1)
    cat(paste0(i,', '))
    sink(nullfile())    # now suppresses

    ## Generate Data
    data[[i]] <- gen_FD_KL_Expansion(
      ns = c(n*break_location,n*(1-break_location)),
      eigsList = list(c(3,2,1,0.5),
                      c(3,2,1,0.5)),
      basesList = list(create.bspline.basis(nbasis=4, norder=4),
                       create.bspline.basis(nbasis=4, norder=4)),
      meansList = c(0,means[l]),
      distsArray = c('sn','sn'),
      skew = 10,
      evals = evalPts,
      kappasArray = 0,
      silent = T)
    data_dist[[i]] <- calculateDistanceMatrix(data[[i]], silent = T,
                                              errType='L2', dataUse = 'Orig')

    # MDT1 <- createMDT(data_dist[[i]],1)
    # MDT3 <- createMDT(data_dist[[i]],3)
    # MDT5 <- createMDT(data_dist[[i]],5)
    # MDT7 <- createMDT(data_dist[[i]],7)
    MDT15 <- createMDT(data_dist[[i]],15)
    # gSeg_MDT1 <- gseg1_new(n, MDT1, pval.perm = T, pval.appr = F)
    # gSeg_MDT3 <- gseg1_new(n, MDT3, pval.perm = T, pval.appr = F)
    # gSeg_MDT5 <- gseg1_new(n, MDT5, pval.perm = T, pval.appr = F)
    # gSeg_MDT7 <- gseg1_new(n, MDT7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_MDT15 <- list()
      gSeg_MDT15$pval.perm$ori$pval <- NA
      gSeg_MDT15$pval.perm$weighted$pval <- NA
      gSeg_MDT15$pval.perm$max.type$pval <- NA
      gSeg_MDT15$pval.perm$generalized$pval <- NA
    }


    # NNL1 <- nnl1(data_dist[[i]],1)
    # NNL3 <- nnl1(data_dist[[i]],3)
    # NNL5 <- nnl1(data_dist[[i]],5)
    # NNL7 <- nnl1(data_dist[[i]],7)
    NNL15 <- nnl1(data_dist[[i]],15)
    # gSeg_NNL1 <- gseg1_new(n, NNL1, pval.perm = T, pval.appr = F)
    # gSeg_NNL3 <- gseg1_new(n, NNL3, pval.perm = T, pval.appr = F)
    # gSeg_NNL5 <- gseg1_new(n, NNL5, pval.perm = T, pval.appr = F)
    # gSeg_NNL7 <- gseg1_new(n, NNL7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_NNL15 <- gseg1_new(n, NNL15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_NNL15 <- list()
      gSeg_NNL15$pval.perm$ori$pval <- NA
      gSeg_NNL15$pval.perm$weighted$pval <- NA
      gSeg_NNL15$pval.perm$max.type$pval <- NA
      gSeg_NNL15$pval.perm$generalized$pval <- NA
    }

    # MST1 <- mstree(as.dist(data_dist[[i]]), ngmax=1)
    # MST3 <- mstree(as.dist(data_dist[[i]]), ngmax=3)
    # MST5 <- mstree(as.dist(data_dist[[i]]), ngmax=5)
    # MST7 <- mstree(as.dist(data_dist[[i]]), ngmax=7)
    MST15 <- mstree(as.dist(data_dist[[i]]), ngmax=15)
    # gSeg_MST1 <- gseg1_new(n, MST1, pval.perm = T, pval.appr = F)
    # gSeg_MST3 <- gseg1_new(n, MST3, pval.perm = T, pval.appr = F)
    # gSeg_MST5 <- gseg1_new(n, MST5, pval.perm = T, pval.appr = F)
    # gSeg_MST7 <- gseg1_new(n, MST7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_MST15 <- gseg1_new(n, MST15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_MST15 <- list()
      gSeg_MST15$pval.perm$ori$pval <- NA
      gSeg_MST15$pval.perm$weighted$pval <- NA
      gSeg_MST15$pval.perm$max.type$pval <- NA
      gSeg_MST15$pval.perm$generalized$pval <- NA
    }

    results[i, ] <- c(
      # # MDT 1
      # gSeg_MDT1$pval.perm$ori$pval,
      # gSeg_MDT1$pval.perm$weighted$pval,
      # gSeg_MDT1$pval.perm$max.type$pval,
      # gSeg_MDT1$pval.perm$generalized$pval,
      # # MDT 3
      # gSeg_MDT3$pval.perm$ori$pval,
      # gSeg_MDT3$pval.perm$weighted$pval,
      # gSeg_MDT3$pval.perm$max.type$pval,
      # gSeg_MDT3$pval.perm$generalized$pval,
      # # MDT 5
      # gSeg_MDT5$pval.perm$ori$pval,
      # gSeg_MDT5$pval.perm$weighted$pval,
      # gSeg_MDT5$pval.perm$max.type$pval,
      # gSeg_MDT5$pval.perm$generalized$pval,
      # # MDT 7
      # gSeg_MDT7$pval.perm$ori$pval,
      # gSeg_MDT7$pval.perm$weighted$pval,
      # gSeg_MDT7$pval.perm$max.type$pval,
      # gSeg_MDT7$pval.perm$generalized$pval,
      # MDT 15
      gSeg_MDT15$pval.perm$ori$pval,
      gSeg_MDT15$pval.perm$weighted$pval,
      gSeg_MDT15$pval.perm$max.type$pval,
      gSeg_MDT15$pval.perm$generalized$pval,

      # # NNL 1
      # gSeg_NNL1$pval.perm$ori$pval,
      # gSeg_NNL1$pval.perm$weighted$pval,
      # gSeg_NNL1$pval.perm$max.type$pval,
      # gSeg_NNL1$pval.perm$generalized$pval,
      # # NNL 3
      # gSeg_NNL3$pval.perm$ori$pval,
      # gSeg_NNL3$pval.perm$weighted$pval,
      # gSeg_NNL3$pval.perm$max.type$pval,
      # gSeg_NNL3$pval.perm$generalized$pval,
      # # NNL 5
      # gSeg_NNL5$pval.perm$ori$pval,
      # gSeg_NNL5$pval.perm$weighted$pval,
      # gSeg_NNL5$pval.perm$max.type$pval,
      # gSeg_NNL5$pval.perm$generalized$pval,
      # # NNL 7
      # gSeg_NNL7$pval.perm$ori$pval,
      # gSeg_NNL7$pval.perm$weighted$pval,
      # gSeg_NNL7$pval.perm$max.type$pval,
      # gSeg_NNL7$pval.perm$generalized$pval,
      # NNL 15
      gSeg_NNL15$pval.perm$ori$pval,
      gSeg_NNL15$pval.perm$weighted$pval,
      gSeg_NNL15$pval.perm$max.type$pval,
      gSeg_NNL15$pval.perm$generalized$pval,

      # # MST 1
      # gSeg_MST1$pval.perm$ori$pval,
      # gSeg_MST1$pval.perm$weighted$pval,
      # gSeg_MST1$pval.perm$max.type$pval,
      # gSeg_MST1$pval.perm$generalized$pval,
      # # MST 3
      # gSeg_MST3$pval.perm$ori$pval,
      # gSeg_MST3$pval.perm$weighted$pval,
      # gSeg_MST3$pval.perm$max.type$pval,
      # gSeg_MST3$pval.perm$generalized$pval,
      # # MST 5
      # gSeg_MST5$pval.perm$ori$pval,
      # gSeg_MST5$pval.perm$weighted$pval,
      # gSeg_MST5$pval.perm$max.type$pval,
      # gSeg_MST5$pval.perm$generalized$pval,
      # # MST 7
      # gSeg_MST7$pval.perm$ori$pval,
      # gSeg_MST7$pval.perm$weighted$pval,
      # gSeg_MST7$pval.perm$max.type$pval,
      # gSeg_MST7$pval.perm$generalized$pval,
      # MST 15
      gSeg_MST15$pval.perm$ori$pval,
      gSeg_MST15$pval.perm$weighted$pval,
      gSeg_MST15$pval.perm$max.type$pval,
      gSeg_MST15$pval.perm$generalized$pval,

      ifelse(is.na(mean_change(data[[i]])),0,1)
    )
    sink()
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(means[l], as.numeric(colSums(results<=0.05,na.rm = T)/colSums(!is.na(results))))

  saveRDS(results_save,
          file=paste0(path,"\\results\\skew_strong_change_full1.rds") )
  saveRDS(results_simplify,
          file=paste0(path,"\\results\\skew_strong_change_simple1.rds") )
}


###############################
library(tidyverse)


plot_15tree_combine <- function(results_simplify){
  # Get all 15 trees
  results_simplify1 <- cbind(results_simplify[,1:2],
                             results_simplify[, stringr::str_detect(colnames(results_simplify),'15')])

  # Format
  data_plot <- results_simplify1 %>% #[,!stringr::str_detect(colnames(results_simplify1),'ori')] %>%
    pivot_longer(cols =MDT15ori:MST15gen)
  data_plot$graph <- substr(data_plot$name,1,3)
  data_plot$stat <- substr(data_plot$name,6,8)
  data_plot <- data_plot[data_plot$mean>=0,]
  data_plot$graph <- ifelse(data_plot$graph=="MDT",'MDP',data_plot$graph)

  ref_data_tmp <- results_simplify[,c(1:2,ncol(results_simplify))]
  ref_data <- rbind(
    cbind(ref_data_tmp,data.frame('graph'='NNL')),
    cbind(ref_data_tmp,data.frame('graph'='MDP')),
    cbind(ref_data_tmp,data.frame('graph'='MST')))
  ref_data <- ref_data[ref_data$mean>=0,]

  #cols <- rep(scales::hue_pal()(4)[c(1,2,4)],3)
  ref_data1 <- ref_data
  ref_data1$value <- ref_data1$`mean.1`
  ref_data1$`mean.1` <- NULL
  ref_data1$name <- 'mean'
  ref_data1$stat <- 'mean'

  data_plot <- rbind(data_plot,ref_data1)
  ggplot( data_plot ) +
    geom_hline(aes(yintercept=0.05), color='gray',
               linetype='dotted', linewidth=2) +
    facet_grid(cols=vars(graph),rows=vars(Type)) +
    geom_line(aes(x=mean,
                  y = value,
                  group=stat,
                  color=stat,
                  linetype=stat),
              linewidth=2) +
    geom_point(aes(x=mean,
                   y = value,
                   group=stat,
                   color=stat,
                   shape=stat),
               size=5) +
    scale_color_manual(
      values=c(scales::hue_pal()(4)[c(1,2,4,3)],'black'),
      breaks = c('gen','max','wei','ori','mean'),
      labels=c('Generalized','Maximum','Weighted','Original','ARS-18')
    ) +
    scale_shape_manual(
      values=c(16,15,17,18,NA),
      breaks = c('gen','max','wei','ori','mean'),
      labels=c('Generalized','Maximum','Weighted','Original','ARS-18')
    ) +
    scale_linetype_manual(
      values=c('twodash','dashed','longdash','dotdash','dotted'),
      #c(rep('solid',4),'dashed'),
      breaks = c('gen','max','wei','ori','mean'),
      labels=c('Generalized','Maximum','Weighted','Original','ARS-18')
    ) +
    theme_bw() +
    # theme(axis.title = element_text(size=28),
    #       axis.text = element_text(size=24),
    #       strip.text = element_text(size = 24),
    #       legend.position = 'bottom',
    #       legend.text = element_text(size=24),
    #       legend.title = element_blank())+
    theme(axis.title = element_text(size=30),
          axis.text = element_text(size=26),
          strip.text = element_text(size = 26),
          legend.position = 'bottom',
          legend.text = element_text(size=26),
          legend.title = element_blank())+
    xlab('Mean Change Size') +
    ylab('Null Rejection Rate') #+
}

results_simplify <-
  rbind(cbind(data.frame('Type'='Normal'),
              readRDS(paste0(path,"\\results\\normal4skew_change_simple.rds") )),
        cbind(data.frame('Type'='Weak Skew'),
              readRDS(paste0(path,"\\results\\skew_change_simple.rds") )),
        cbind(data.frame('Type'='Strong Skew'),
              readRDS(paste0(path,"\\results\\skew_strong_change_simple1.rds") ))
  )
results_simplify$mean.1 <- 1-results_simplify$mean.1
results_simplify$Type <- factor(results_simplify$Type,
                                levels = c('Normal','Weak Skew','Strong Skew'))

# Plot
tmp <- plot_15tree_combine(results_simplify)


# results_simplify$mean.1 <- 1-results_simplify$mean.1
# results_simplify1 <- results_simplify
# data_plot <- results_simplify1 %>% pivot_longer(cols =MDT15wei:mean.1)
#
# ggplot(na.omit(data_plot[data_plot$name!='mean.1',])) +
#   geom_hline(aes(yintercept=0.05), color='gray',linetype='dotted',linewidth=1.75) +
#   geom_line(aes(x=mean,
#                 y=value,
#                 color=name),linewidth=2) +
#   geom_point(aes(x=mean,
#                  y=value,
#                  color=name),size=4) +
#   #geom_vline(aes(xintercept=0)) +
#   geom_line(aes(x=mean,
#                 y=value), data=na.omit(data_plot[data_plot$name=='mean.1',]),
#             color='black', linewidth=1,linetype='dashed') +
#   guides(color="none") +
#   theme_bw() +
#   theme(axis.title = element_text(size=22),
#         axis.text = element_text(size=18))+
#   xlab('Mean Change Size') +
#   ylab('Null Rejection Rate') +
#   scale_color_manual(values=scales::hue_pal()(4)[c(1,2,4)])

png(paste0(path,"\\figures\\mean_power_skew.png"), width = 1800, height = 1000 )
print(tmp)
dev.off()
