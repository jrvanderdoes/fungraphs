Skewed Normal (weak) mean change (PLOT COMPLETE):
  ```{r}
trials <- 1000
evalPts <- seq(0,1,0.05)
n <- 50 #c(200, 500, 1000)
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
                      'MST15max'=NA, 'MST15gen'=NA,

                      'mean'=NA, 'cov'=NA
)

results_save <- list()
means = c(-1,-0.75,-0.5,-0.25,-0.15,-0.1,-0.05,0,
          0.05,0.1,0.15,0.25,0.5,0.75,1)
results_simplify <- data.frame('mean'=means,
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
                               'MST15max'=NA, 'MST15gen'=NA,

                               'mean'=NA, 'cov'=NA
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

    MDT1 <- createMDT(data_dist[[i]],1)
    MDT3 <- createMDT(data_dist[[i]],3)
    MDT5 <- createMDT(data_dist[[i]],5)
    MDT7 <- createMDT(data_dist[[i]],7)
    MDT15 <- createMDT(data_dist[[i]],15)
    gSeg_MDT1 <- gseg1_new(n, MDT1, pval.perm = T, pval.appr = F)
    gSeg_MDT3 <- gseg1_new(n, MDT3, pval.perm = T, pval.appr = F)
    gSeg_MDT5 <- gseg1_new(n, MDT5, pval.perm = T, pval.appr = F)
    gSeg_MDT7 <- gseg1_new(n, MDT7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_MDT15 <- list()
      gSeg_MDT15$pval.perm$ori$pval <- NA
      gSeg_MDT15$pval.perm$weighted$pval <- NA
      gSeg_MDT15$pval.perm$max.type$pval <- NA
      gSeg_MDT15$pval.perm$generalized$pval <- NA
    }


    NNL1 <- nnl(data_dist[[i]],1)
    NNL3 <- nnl(data_dist[[i]],3)
    NNL5 <- nnl(data_dist[[i]],5)
    NNL7 <- nnl(data_dist[[i]],7)
    NNL15 <- nnl(data_dist[[i]],15)
    gSeg_NNL1 <- gseg1_new(n, NNL1, pval.perm = T, pval.appr = F)
    gSeg_NNL3 <- gseg1_new(n, NNL3, pval.perm = T, pval.appr = F)
    gSeg_NNL5 <- gseg1_new(n, NNL5, pval.perm = T, pval.appr = F)
    gSeg_NNL7 <- gseg1_new(n, NNL7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_NNL15 <- gseg1_new(n, NNL15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_NNL15 <- list()
      gSeg_NNL15$pval.perm$ori$pval <- NA
      gSeg_NNL15$pval.perm$weighted$pval <- NA
      gSeg_NNL15$pval.perm$max.type$pval <- NA
      gSeg_NNL15$pval.perm$generalized$pval <- NA
    }

    MST1 <- mstree(as.dist(data_dist[[i]]), ngmax=1)
    MST3 <- mstree(as.dist(data_dist[[i]]), ngmax=3)
    MST5 <- mstree(as.dist(data_dist[[i]]), ngmax=5)
    MST7 <- mstree(as.dist(data_dist[[i]]), ngmax=7)
    MST15 <- mstree(as.dist(data_dist[[i]]), ngmax=15)
    gSeg_MST1 <- gseg1_new(n, MST1, pval.perm = T, pval.appr = F)
    gSeg_MST3 <- gseg1_new(n, MST3, pval.perm = T, pval.appr = F)
    gSeg_MST5 <- gseg1_new(n, MST5, pval.perm = T, pval.appr = F)
    gSeg_MST7 <- gseg1_new(n, MST7, pval.perm = T, pval.appr = F)
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
      gSeg_MST15$pval.perm$generalized$pval,

      ifelse(is.na(mean_change(data[[i]])),0,1),
      ifelse(is.na(cov_change(data[[i]])),0,1)
    )
    sink()
    saveRDS(results, paste0('C:/Users/jerem/OneDrive/Documents/',
                            'School/Waterloo/Research/',
                            'GraphCP/Data/Simulations/',
                            'tmp_t',trials,'_l',l,'_dsn2','.rds'))
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(means[l], as.numeric(colSums(results<=0.05,na.rm = T)/colSums(!is.na(results))))

  saveRDS(results_save, paste0('C:/Users/jerem/OneDrive/Documents/',
                               'School/Waterloo/Research/',
                               'GraphCP/Data/Simulations/power/',
                               'trials',trials,
                               '_ns',n,
                               '_bl',break_location,
                               '_dsn2',
                               '_results.rds'))
  saveRDS(results_simplify, paste0('C:/Users/jerem/OneDrive/Documents/',
                                   'School/Waterloo/Research/',
                                   'GraphCP/Data/Simulations/power/',
                                   'trials',trials,
                                   '_ns',n,
                                   '_bl',break_location,
                                   '_dsn2',
                                   '_results_simp.rds'))
}


# Plot

tmp <- plot_15tree_meanchange(results_simplify)

png("C:/Users/jerem/Downloads/updatedPlots/mean_power_sn2.png")
print(tmp)
dev.off()

```



```{r distChanges_old}
trials <- 1000
evalPts <- seq(0,1,0.05)
n <- 50 #c(200, 500, 1000)
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
                      'MST15max'=NA, 'MST15gen'=NA,

                      'mean'=NA, 'cov'=NA
)

results_save <- list()
dists = c('Exponential','t','Binomial')
results_simplify <- data.frame('dist'=dists,
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
                               'MST15max'=NA, 'MST15gen'=NA,

                               'mean'=NA, 'cov'=NA
)

data <- data_dist <- list()


for(l in 2:length(dists)){

  cat(paste0('Dist, ',dists[l],' (',trials,'): '))

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
      distsArray = c('Normal',dists[l]),
      evals = evalPts,
      kappasArray = 0,
      silent = T)
    data_dist[[i]] <- calculateDistanceMatrix(data[[i]], silent = T,
                                              errType='L2', dataUse = 'Orig')

    MDT1 <- createMDT(data_dist[[i]],1)
    MDT3 <- createMDT(data_dist[[i]],3)
    MDT5 <- createMDT(data_dist[[i]],5)
    MDT7 <- createMDT(data_dist[[i]],7)
    MDT15 <- createMDT(data_dist[[i]],15)
    gSeg_MDT1 <- gseg1_new(n, MDT1, pval.perm = T, pval.appr = F)
    gSeg_MDT3 <- gseg1_new(n, MDT3, pval.perm = T, pval.appr = F)
    gSeg_MDT5 <- gseg1_new(n, MDT5, pval.perm = T, pval.appr = F)
    gSeg_MDT7 <- gseg1_new(n, MDT7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_MDT15 <- list()
      gSeg_MDT15$pval.perm$ori$pval <- NA
      gSeg_MDT15$pval.perm$weighted$pval <- NA
      gSeg_MDT15$pval.perm$max.type$pval <- NA
      gSeg_MDT15$pval.perm$generalized$pval <- NA
    }


    NNL1 <- nnl(data_dist[[i]],1)
    NNL3 <- nnl(data_dist[[i]],3)
    NNL5 <- nnl(data_dist[[i]],5)
    NNL7 <- nnl(data_dist[[i]],7)
    NNL15 <- nnl(data_dist[[i]],15)
    gSeg_NNL1 <- gseg1_new(n, NNL1, pval.perm = T, pval.appr = F)
    gSeg_NNL3 <- gseg1_new(n, NNL3, pval.perm = T, pval.appr = F)
    gSeg_NNL5 <- gseg1_new(n, NNL5, pval.perm = T, pval.appr = F)
    gSeg_NNL7 <- gseg1_new(n, NNL7, pval.perm = T, pval.appr = F)
    if(n>30) {
      gSeg_NNL15 <- gseg1_new(n, NNL15, pval.perm = T, pval.appr = F)
    } else{
      gSeg_NNL15 <- list()
      gSeg_NNL15$pval.perm$ori$pval <- NA
      gSeg_NNL15$pval.perm$weighted$pval <- NA
      gSeg_NNL15$pval.perm$max.type$pval <- NA
      gSeg_NNL15$pval.perm$generalized$pval <- NA
    }

    MST1 <- mstree(as.dist(data_dist[[i]]), ngmax=1)
    MST3 <- mstree(as.dist(data_dist[[i]]), ngmax=3)
    MST5 <- mstree(as.dist(data_dist[[i]]), ngmax=5)
    MST7 <- mstree(as.dist(data_dist[[i]]), ngmax=7)
    MST15 <- mstree(as.dist(data_dist[[i]]), ngmax=15)
    gSeg_MST1 <- gseg1_new(n, MST1, pval.perm = T, pval.appr = F)
    gSeg_MST3 <- gseg1_new(n, MST3, pval.perm = T, pval.appr = F)
    gSeg_MST5 <- gseg1_new(n, MST5, pval.perm = T, pval.appr = F)
    gSeg_MST7 <- gseg1_new(n, MST7, pval.perm = T, pval.appr = F)
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
      gSeg_MST15$pval.perm$generalized$pval,

      ifelse(is.na(mean_change(data[[i]])),0,1),
      ifelse(is.na(cov_change(data[[i]])),0,1)
    )
    sink()
    saveRDS(results, paste0('C:/Users/jerem/OneDrive/Documents/',
                            'School/Waterloo/Research/',
                            'GraphCP/Data/Simulations/',
                            'tmp_t',trials,'_l',l,'_dists','.rds'))
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(dists[l], as.numeric(colSums(results<=0.05,na.rm = T)/colSums(!is.na(results))))

  saveRDS(results_save, paste0('C:/Users/jerem/OneDrive/Documents/',
                               'School/Waterloo/Research/',
                               'GraphCP/Data/Simulations/power/',
                               'trials',trials,
                               '_ns',n,
                               '_bl',break_location,
                               '_dists',
                               '_results.rds'))
  saveRDS(results_simplify, paste0('C:/Users/jerem/OneDrive/Documents/',
                                   'School/Waterloo/Research/',
                                   'GraphCP/Data/Simulations/power/',
                                   'trials',trials,
                                   '_ns',n,
                                   '_bl',break_location,
                                   '_dists',
                                   '_results_simp.rds'))
}


# Plot
results_simplify$mean <- 1-as.numeric(results_simplify$mean)
results_simplify$cov <- 1-as.numeric(results_simplify$cov)
data_plot <- results_simplify %>% pivot_longer(cols = MDT1ori:cov)
data_plot$value <- as.numeric(data_plot$value)

tmp <- ggplot(na.omit(data_plot)) +
  geom_line(aes(x=dist,
                y=value,
                color=name,group=name)) +
  geom_point(aes(x=dist,
                 y=value,
                 color=name,group=name)) +
  geom_hline(aes(yintercept=0.05)) +
  geom_vline(aes(xintercept=0)) +
  geom_line(aes(x=dist,
                y=value,group=name), data=na.omit(data_plot[data_plot$name=='mean',]),
            color='black', linewidth=1) +
  geom_line(aes(x=dist,
                y=value,group=name), data=na.omit(data_plot[data_plot$name=='cov',]),
            color='black', linewidth=1, linetype='dashed') +
  guides(color="none") +
  theme_bw() +
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=18))+
  xlab('Distribution Change') +
  ylab('Null Rejection Rate') +
  ylim(c(0,1))

png("C:/Users/jerem/Downloads/dist_power.png")
print(tmp)
dev.off()
```



## For 3 only
plot_15tree_meanchange <- function(results_simplify){
  results_simplify$mean.1 <- 1-results_simplify$mean.1
  results_simplify$cov <- 1-results_simplify$cov
  results_simplify1 <- results_simplify[,c(1,1+4*4+1:4, ncol(results_simplify)-1)]
  data_plot <- results_simplify1 %>% pivot_longer(cols =MDT15wei:mean.1)

  ggplot(na.omit(data_plot[data_plot$name!='mean.1',])) +
    geom_hline(aes(yintercept=0.05), color='gray',linetype='dotted',linewidth=1.75) +
    geom_line(aes(x=mean,
                  y=value,
                  color=name),linewidth=2) +
    geom_point(aes(x=mean,
                   y=value,
                   color=name),size=4) +
    #geom_vline(aes(xintercept=0)) +
    geom_line(aes(x=mean,
                  y=value), data=na.omit(data_plot[data_plot$name=='mean.1',]),
              color='black', linewidth=1,linetype='dashed') +
    guides(color="none") +
    theme_bw() +
    theme(axis.title = element_text(size=22),
          axis.text = element_text(size=18))+
    xlab('Mean Change Size') +
    ylab('Null Rejection Rate') +
    scale_color_manual(values=scales::hue_pal()(4)[c(1,2,4)])
}
