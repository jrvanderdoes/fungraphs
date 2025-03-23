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

                      'mean'=NA#, 'cov'=NA
)

results_save <- list()
deps = c(0, 0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9)
results_simplify <- data.frame('dep'=deps,
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

                               'mean'=NA#, 'cov'=NA
)

data <- data_dep <- list()


for(l in 1:length(deps)){

  cat(paste0('dep, ',deps[l],' (',trials,'): '))

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
      meansList = c(0,0.5),
      distsArray = c('Normal'),
      kappasArray = deps[l],
      evals = evalPts,
      silent = T)
    data_dep[[i]] <- calculateDistanceMatrix(data[[i]], silent = T,
                                             errType='L2', dataUse = 'Orig')

    # MDT1 <- createMDT(data_dep[[i]],1)
    # MDT3 <- createMDT(data_dep[[i]],3)
    # MDT5 <- createMDT(data_dep[[i]],5)
    # MDT7 <- createMDT(data_dep[[i]],7)
    MDT15 <- createMDT(data_dep[[i]],15)
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


    # NNL1 <- nnl(data_dep[[i]],1)
    # NNL3 <- nnl(data_dep[[i]],3)
    # NNL5 <- nnl(data_dep[[i]],5)
    # NNL7 <- nnl(data_dep[[i]],7)
    NNL15 <- nnl1(data_dep[[i]],15)
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

    # MST1 <- mstree(as.dist(data_dep[[i]]), ngmax=1)
    # MST3 <- mstree(as.dist(data_dep[[i]]), ngmax=3)
    # MST5 <- mstree(as.dist(data_dep[[i]]), ngmax=5)
    # MST7 <- mstree(as.dist(data_dep[[i]]), ngmax=7)
    MST15 <- mstree(as.dist(data_dep[[i]]), ngmax=15)
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

      ifelse(is.na(mean_change(data[[i]])),0,1)#,
      # ifelse(is.na(cov_change(data[[i]])),0,1)
    )
    sink()
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(deps[l], as.numeric(colSums(results<=0.05,na.rm = T) /
                            colSums(!is.na(results))))

  saveRDS(results_save,
          file=paste0(path,"\\results\\dep_inc_mean_change_full.rds") )
  saveRDS(results_simplify,
          file=paste0(path,"\\results\\dep_inc_mean_change_simple.rds") )
}

##############################################################

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

  'mean'=NA#, 'cov'=NA
)

results_save <- list()
deps = c(0, 0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9)
results_simplify <- data.frame('dep'=deps,
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

                               'mean'=NA#, 'cov'=NA
)

data <- data_dep <- list()


for(l in 1:length(deps)){

  cat(paste0('dep, ',deps[l],' (',trials,'): '))

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
      meansList = c(0,0),
      distsArray = c('Normal'),
      kappasArray = deps[l],
      evals = evalPts,
      silent = T)
    data_dep[[i]] <- calculateDistanceMatrix(data[[i]], silent = T,
                                             errType='L2', dataUse = 'Orig')

    # MDT1 <- createMDT(data_dep[[i]],1)
    # MDT3 <- createMDT(data_dep[[i]],3)
    # MDT5 <- createMDT(data_dep[[i]],5)
    # MDT7 <- createMDT(data_dep[[i]],7)
    MDT15 <- createMDT(data_dep[[i]],15)
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


    # NNL1 <- nnl(data_dep[[i]],1)
    # NNL3 <- nnl(data_dep[[i]],3)
    # NNL5 <- nnl(data_dep[[i]],5)
    # NNL7 <- nnl(data_dep[[i]],7)
    NNL15 <- nnl1(data_dep[[i]],15)
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

    # MST1 <- mstree(as.dist(data_dep[[i]]), ngmax=1)
    # MST3 <- mstree(as.dist(data_dep[[i]]), ngmax=3)
    # MST5 <- mstree(as.dist(data_dep[[i]]), ngmax=5)
    # MST7 <- mstree(as.dist(data_dep[[i]]), ngmax=7)
    MST15 <- mstree(as.dist(data_dep[[i]]), ngmax=15)
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

      ifelse(is.na(mean_change(data[[i]])),0,1)#,
      # ifelse(is.na(cov_change(data[[i]])),0,1)
    )
    sink()
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(deps[l], as.numeric(colSums(results<=0.05,na.rm = T) /
                            colSums(!is.na(results))))

  saveRDS(results_save,
          file=paste0(path,"\\results\\dep_inc_null_full.rds") )
  saveRDS(results_simplify,
          file=paste0(path,"\\results\\dep_inc_null_simple.rds") )
}

##############################################################
library(tidyverse)
plot_15tree_size_adj <- function(data_full){

  # Get all 15 trees
  results_simplify <- data_full[stringr::str_detect(data_full$name,'15'),]

  # Format
  data_plot <- results_simplify #%>%
  #pivot_longer(cols =MDT15ori:MST15gen)
  data_plot$graph <- substr(data_plot$name,1,3)
  data_plot$stat <- str_replace(str_replace(data_plot$name,'[0-9]+',''),'(NNL)|(MST)|(MDT)','')
  data_plot$graph <- ifelse(data_plot$graph=="MDT",'MDP',data_plot$graph)

  ref_data_tmp <- data_full[data_full$name=='mean',c("name",'dep',"sizeAdjPower")]
  ref_data <- rbind(
    cbind(ref_data_tmp,data.frame('graph'='NNL')),
    cbind(ref_data_tmp,data.frame('graph'='MDP')),
    cbind(ref_data_tmp,data.frame('graph'='MST')))

  ref_data$value <- ref_data$mean
  ref_data$stat <- 'mean'
  ref_data$name <- 'mean'

  data_plot <- rbind(
    data_plot[,-c(3:5)],
    data.frame(ref_data)
  )


  ggplot( data_plot ) +
    geom_hline(aes(yintercept=0.05),
               color='gray',linetype='dotted',linewidth=2) +
    facet_grid(cols=vars(graph)) +
    geom_line(aes(x = dep,
                  y = sizeAdjPower,
                  group=stat,
                  color=stat,
                  linetype=stat),
              linewidth=2) +
    geom_point(aes(x=dep,
                   y = sizeAdjPower,
                   group=stat,
                   color=stat,
                   shape=stat),
               size=5) +
    # geom_line(aes(x=dep,
    #               y=sizeAdjPower),
    #           data=ref_data,
    #           color='black', linewidth=1,linetype='dotdash') +
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
      #values=c(rep('solid',4),'dashed'),
      values=c('twodash','dashed','longdash','dotdash','dotted'),
      breaks = c('gen','max','wei','ori','mean'),
      labels=c('Generalized','Maximum','Weighted','Original','ARS-18')
    ) +
    theme_bw() +
    xlim(c(0,1)) +
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
    xlab('Length of Sample') +
    xlab('Dependence Strength') +
    ylab('Null Rejection Rate')
}

# Plot
results_size <- readRDS(paste0(path,"\\results\\dep_inc_null_simple.rds"))
results_size$mean <- 1-as.numeric(results_size$mean)
data_size <- results_size %>% pivot_longer(cols = MDT15ori:mean)
data_size$value <- as.numeric(data_size$value)

results_power <- readRDS(paste0(path,"\\results\\dep_inc_mean_change_simple.rds"))
results_power$mean <- 1-as.numeric(results_power$mean)
data_power <- results_power %>% pivot_longer(cols = MDT15ori:mean)
data_power$value <- as.numeric(data_power$value)

## Full_data
data_power1 <- data_power
colnames(data_power1) <- c('dep','name','power')
data_size1 <- data_size
colnames(data_size1) <- c('dep','name','size')

data_full <- merge(data_power1,data_size1,by = c('name','dep'))
data_full$sizeTo5 <- data_full$size-0.05
data_full$sizeAdjPower <- data_full$power-data_full$sizeTo5

# Plot

tmp <- plot_15tree_size_adj(data_full)

png(paste0(path,"\\figures\\dep_size_adj_power.png"),width = 1600,height = 800)
print(tmp)
dev.off()
