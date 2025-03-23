trials <- 1000
evalPts <- seq(0,1,0.05)
n <- 50
break_location <- 0.5


results_save <- list()
dists = c(10,5,4,3,2,1)
results_simplify <- data.frame('dist'=dists,
                               'MDT15ori'=NA, 'MDT15wei'=NA,
                               'MDT15max'=NA, 'MDT15gen'=NA,

                               'NNL15ori'=NA, 'NNL15wei'=NA,
                               'NNL15max'=NA, 'NNL15gen'=NA,

                               'MST15ori'=NA, 'MST15wei'=NA,
                               'MST15max'=NA, 'MST15gen'=NA,

                               'mean'=NA, 'cov'=NA
)

data <- data_dist <- list()


for(l in 1:length(dists)){

  cat(paste0('Dist, ',dists[l],' (',trials,'): '))
  results <- data.frame('MDT15ori'=1:trials, 'MDT15wei'=NA,
                        'MDT15max'=NA, 'MDT15gen'=NA,

                        'NNL15ori'=NA, 'NNL15wei'=NA,
                        'NNL15max'=NA, 'NNL15gen'=NA,

                        'MST15ori'=NA, 'MST15wei'=NA,
                        'MST15max'=NA, 'MST15gen'=NA,

                        'mean'=NA, 'cov'=NA
  )

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
      distsArray = c('Normal','t'),
      dof = dists[l],
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
    # if(n>30) {
    gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)
    # } else{
    #   gSeg_MDT15 <- list()
    #   gSeg_MDT15$pval.perm$ori$pval <- NA
    #   gSeg_MDT15$pval.perm$weighted$pval <- NA
    #   gSeg_MDT15$pval.perm$max.type$pval <- NA
    #   gSeg_MDT15$pval.perm$generalized$pval <- NA
    # }


    # NNL1 <- nnl(data_dist[[i]],1)
    # NNL3 <- nnl(data_dist[[i]],3)
    # NNL5 <- nnl(data_dist[[i]],5)
    # NNL7 <- nnl(data_dist[[i]],7)
    NNL15 <- nnl(data_dist[[i]],15)
    # gSeg_NNL1 <- gseg1_new(n, NNL1, pval.perm = T, pval.appr = F)
    # gSeg_NNL3 <- gseg1_new(n, NNL3, pval.perm = T, pval.appr = F)
    # gSeg_NNL5 <- gseg1_new(n, NNL5, pval.perm = T, pval.appr = F)
    # gSeg_NNL7 <- gseg1_new(n, NNL7, pval.perm = T, pval.appr = F)
    # if(n>30) {
    gSeg_NNL15 <- gseg1_new(n, NNL15, pval.perm = T, pval.appr = F)
    # } else{
    #   gSeg_NNL15 <- list()
    #   gSeg_NNL15$pval.perm$ori$pval <- NA
    #   gSeg_NNL15$pval.perm$weighted$pval <- NA
    #   gSeg_NNL15$pval.perm$max.type$pval <- NA
    #   gSeg_NNL15$pval.perm$generalized$pval <- NA
    # }

    # MST1 <- mstree(as.dist(data_dist[[i]]), ngmax=1)
    # MST3 <- mstree(as.dist(data_dist[[i]]), ngmax=3)
    # MST5 <- mstree(as.dist(data_dist[[i]]), ngmax=5)
    # MST7 <- mstree(as.dist(data_dist[[i]]), ngmax=7)
    MST15 <- mstree(as.dist(data_dist[[i]]), ngmax=15)
    # gSeg_MST1 <- gseg1_new(n, MST1, pval.perm = T, pval.appr = F)
    # gSeg_MST3 <- gseg1_new(n, MST3, pval.perm = T, pval.appr = F)
    # gSeg_MST5 <- gseg1_new(n, MST5, pval.perm = T, pval.appr = F)
    # gSeg_MST7 <- gseg1_new(n, MST7, pval.perm = T, pval.appr = F)
    # if(n>30) {
    gSeg_MST15 <- gseg1_new(n, MST15, pval.perm = T, pval.appr = F)
    # } else{
    #   gSeg_MST15 <- list()
    #   gSeg_MST15$pval.perm$ori$pval <- NA
    #   gSeg_MST15$pval.perm$weighted$pval <- NA
    #   gSeg_MST15$pval.perm$max.type$pval <- NA
    #   gSeg_MST15$pval.perm$generalized$pval <- NA
    # }

    results[i, ] <- c(
      # MDT 15
      gSeg_MDT15$pval.perm$ori$pval,
      gSeg_MDT15$pval.perm$weighted$pval,
      gSeg_MDT15$pval.perm$max.type$pval,
      gSeg_MDT15$pval.perm$generalized$pval,

      # NNL 15
      gSeg_NNL15$pval.perm$ori$pval,
      gSeg_NNL15$pval.perm$weighted$pval,
      gSeg_NNL15$pval.perm$max.type$pval,
      gSeg_NNL15$pval.perm$generalized$pval,

      # MST 15
      gSeg_MST15$pval.perm$ori$pval,
      gSeg_MST15$pval.perm$weighted$pval,
      gSeg_MST15$pval.perm$max.type$pval,
      gSeg_MST15$pval.perm$generalized$pval,

      ifelse(is.na(mean_change(data[[i]])),0,1),
      ifelse(is.na(cov_change(data[[i]])),0,1)
    )
    sink()
    # saveRDS(results, paste0('C:/Users/jerem/OneDrive/Documents/',
    #                         'School/Waterloo/Research/',
    #                         'GraphCP/Data/Simulations/',
    #                         'tmp_t',trials,'_l',l,'_dists_t','.rds'))
    saveRDS(results, paste0('C:/Users/jerem/Downloads/updatedPlots/sims/',
                            'tmp_t',trials,'_l',l,'_dists_t','.rds'))
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(dists[l], as.numeric(colSums(results<=0.05,na.rm = T)/colSums(!is.na(results))))

  saveRDS(results_save, paste0('C:/Users/jerem/Downloads/updatedPlots/sims/',
                               'trials',trials,
                               '_ns',n,
                               '_bl',break_location,
                               '_dists_t',
                               '_results.rds'))
  saveRDS(results_simplify, paste0('C:/Users/jerem/Downloads/updatedPlots/sims/',
                                   'trials',trials,
                                   '_ns',n,
                                   '_bl',break_location,
                                   '_dists_t',
                                   '_results_simp.rds'))
}
