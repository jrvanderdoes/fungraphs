trials <- 1000
evalPts <- seq(0,1,0.05)
n <- 100
break_locations <- c(0.05,0.1,0.2,0.25,0.5,0.75,0.8,0.9,0.95)

results <- data.frame(#'MDT1ori'=1:trials, 'MDT1wei'=NA,
                      #'MDT1max'=NA, 'MDT1gen'=NA,
                      #'MDT3ori'=NA, 'MDT3wei'=NA,
                      #'MDT3max'=NA, 'MDT3gen'=NA,
                      #'MDT5ori'=NA, 'MDT5wei'=NA,
                      #'MDT5max'=NA, 'MDT5gen'=NA,
                      #'MDT7ori'=NA, 'MDT7wei'=NA,
                      #'MDT7max'=NA, 'MDT7gen'=NA,
                      'MDT15ori'=1:trials, 'MDT15wei'=NA,
                      'MDT15max'=NA, 'MDT15gen'=NA,

                      #'NNL1ori'=NA, 'NNL1wei'=NA,
                      #'NNL1max'=NA, 'NNL1gen'=NA,
                      #'NNL3ori'=NA, 'NNL3wei'=NA,
                      #'NNL3max'=NA, 'NNL3gen'=NA,
                      #'NNL5ori'=NA, 'NNL5wei'=NA,
                      #'NNL5max'=NA, 'NNL5gen'=NA,
                      #'NNL7ori'=NA, 'NNL7wei'=NA,
                      #'NNL7max'=NA, 'NNL7gen'=NA,
                      'NNL15ori'=NA, 'NNL15wei'=NA,
                      'NNL15max'=NA, 'NNL15gen'=NA,

                      #'MST1ori'=NA, 'MST1wei'=NA,
                      #'MST1max'=NA, 'MST1gen'=NA,
                      #'MST3ori'=NA, 'MST3wei'=NA,
                      #'MST3max'=NA, 'MST3gen'=NA,
                      #'MST5ori'=NA, 'MST5wei'=NA,
                      #'MST5max'=NA, 'MST5gen'=NA,
                      #'MST7ori'=NA, 'MST7wei'=NA,
                      #'MST7max'=NA, 'MST7gen'=NA,
                      'MST15ori'=NA, 'MST15wei'=NA,
                      'MST15max'=NA, 'MST15gen'=NA,

                      'mean'=NA, 'cov'=NA
)

results_save <- list()
results_simplify <- data.frame('location'=break_locations,
                               #'MDT1ori'=NA, 'MDT1wei'=NA,
                               #'MDT1max'=NA, 'MDT1gen'=NA,
                               #'MDT3ori'=NA, 'MDT3wei'=NA,
                               #'MDT3max'=NA, 'MDT3gen'=NA,
                               #'MDT5ori'=NA, 'MDT5wei'=NA,
                               #'MDT5max'=NA, 'MDT5gen'=NA,
                               #'MDT7ori'=NA, 'MDT7wei'=NA,
                               #'MDT7max'=NA, 'MDT7gen'=NA,
                               'MDT15ori'=NA, 'MDT15wei'=NA,
                               'MDT15max'=NA, 'MDT15gen'=NA,

                               #'NNL1ori'=NA, 'NNL1wei'=NA,
                               #'NNL1max'=NA, 'NNL1gen'=NA,
                               #'NNL3ori'=NA, 'NNL3wei'=NA,
                               #'NNL3max'=NA, 'NNL3gen'=NA,
                               #'NNL5ori'=NA, 'NNL5wei'=NA,
                               #'NNL5max'=NA, 'NNL5gen'=NA,
                               #'NNL7ori'=NA, 'NNL7wei'=NA,
                               #'NNL7max'=NA, 'NNL7gen'=NA,
                               'NNL15ori'=NA, 'NNL15wei'=NA,
                               'NNL15max'=NA, 'NNL15gen'=NA,

                               #'MST1ori'=NA, 'MST1wei'=NA,
                               #'MST1max'=NA, 'MST1gen'=NA,
                               #'MST3ori'=NA, 'MST3wei'=NA,
                               #'MST3max'=NA, 'MST3gen'=NA,
                               #'MST5ori'=NA, 'MST5wei'=NA,
                               #'MST5max'=NA, 'MST5gen'=NA,
                               #'MST7ori'=NA, 'MST7wei'=NA,
                               #'MST7max'=NA, 'MST7gen'=NA,
                               'MST15ori'=NA, 'MST15wei'=NA,
                               'MST15max'=NA, 'MST15gen'=NA,

                               'mean'=NA, 'cov'=NA
)

data <- data_dist <- list()


for(l in 1:length(break_locations)){

  cat(paste0('location, ',break_locations[l],' (',trials,'): '))

  for(i in 1:trials){
    set.seed(123456 * (l+1) + i + 1)
    cat(paste0(i,', '))
    sink(nullfile())    # now suppresses

    ## Generate Data
    data[[i]] <- gen_FD_KL_Expansion(
      ns = c(n*break_locations[l],n*(1-break_locations[l])),
      eigsList = list(c(3,2,1,0.5),
                      c(3,2,1,0.5)),
      basesList = list(create.bspline.basis(nbasis=4, norder=4),
                       create.bspline.basis(nbasis=4, norder=4)),
      meansList = c(0,0.25),
      distsArray = c('Normal','Normal'),
      evals = evalPts,
      kappasArray = 0,
      silent = T)
    data_dist[[i]] <- calculateDistanceMatrix(data[[i]], silent = T,
                                              errType='L2', dataUse = 'Orig')

    #MDT1 <- createMDT(data_dist[[i]],1)
    #MDT3 <- createMDT(data_dist[[i]],3)
    #MDT5 <- createMDT(data_dist[[i]],5)
    #MDT7 <- createMDT(data_dist[[i]],7)
    MDT15 <- createMDT(data_dist[[i]],15)
    #gSeg_MDT1 <- gseg1_new(n, MDT1, pval.perm = T, pval.appr = F)
    #gSeg_MDT3 <- gseg1_new(n, MDT3, pval.perm = T, pval.appr = F)
    #gSeg_MDT5 <- gseg1_new(n, MDT5, pval.perm = T, pval.appr = F)
    #gSeg_MDT7 <- gseg1_new(n, MDT7, pval.perm = T, pval.appr = F)
    gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)


    #NNL1 <- nnl(data_dist[[i]],1)
    #NNL3 <- nnl(data_dist[[i]],3)
    #NNL5 <- nnl(data_dist[[i]],5)
    #NNL7 <- nnl(data_dist[[i]],7)
    NNL15 <- nnl(data_dist[[i]],15)
    #gSeg_NNL1 <- gseg1_new(n, NNL1, pval.perm = T, pval.appr = F)
    #gSeg_NNL3 <- gseg1_new(n, NNL3, pval.perm = T, pval.appr = F)
    #gSeg_NNL5 <- gseg1_new(n, NNL5, pval.perm = T, pval.appr = F)
    #gSeg_NNL7 <- gseg1_new(n, NNL7, pval.perm = T, pval.appr = F)
    gSeg_NNL15 <- gseg1_new(n, NNL15, pval.perm = T, pval.appr = F)

    #MST1 <- mstree(as.dist(data_dist[[i]]), ngmax=1)
    #MST3 <- mstree(as.dist(data_dist[[i]]), ngmax=3)
    #MST5 <- mstree(as.dist(data_dist[[i]]), ngmax=5)
    #MST7 <- mstree(as.dist(data_dist[[i]]), ngmax=7)
    MST15 <- mstree(as.dist(data_dist[[i]]), ngmax=15)
    #gSeg_MST1 <- gseg1_new(n, MST1, pval.perm = T, pval.appr = F)
    #gSeg_MST3 <- gseg1_new(n, MST3, pval.perm = T, pval.appr = F)
    #gSeg_MST5 <- gseg1_new(n, MST5, pval.perm = T, pval.appr = F)
    #gSeg_MST7 <- gseg1_new(n, MST7, pval.perm = T, pval.appr = F)
    gSeg_MST15 <- gseg1_new(n, MST15, pval.perm = T, pval.appr = F)

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

      ifelse(is.na(mean_change(data[[i]])),0,1),
      tryCatch({ifelse(is.na(cov_change(data[[i]])),0,1)},
               error=function(e){ifelse(is.na(cov_change(data[[i]])),0,1)})
    )
    sink()
    saveRDS(results, paste0('C:/Users/jerem/OneDrive/Documents/',
                            'School/Waterloo/Research/',
                            'GraphCP/Data/Simulations/CL_',
                            'tmp_t',trials,'_loc',l,'_dNorm','.rds'))
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(break_locations[l], as.numeric(colSums(results<=0.05,na.rm = T)/colSums(!is.na(results))))

  saveRDS(results_save, paste0('C:/Users/jerem/Downloads/funkycells_temp/CL_',
                               'trials',trials,
                               '_ns',n,
                               '_bl',break_locations[l],
                               '_dNorm',
                               '_results.rds'))
  saveRDS(results_simplify, paste0('C:/Users/jerem/Downloads/funkycells_temp/CL_',
                                   'trials',trials,
                                   '_ns',n,
                                   '_bl',break_locations[l],
                                   '_dNorm',
                                   '_results_simp.rds'))
}
