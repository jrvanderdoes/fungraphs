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
eigs <- c(1, 2, 4, 6, 8)

results_simplify <-
  data.frame('eigs'=1:length(eigs),
             'MDT15ori'=NA, 'MDT15wei'=NA,
             'MDT15max'=NA, 'MDT15gen'=NA,

             'NNL15ori'=NA, 'NNL15wei'=NA,
             'NNL15max'=NA, 'NNL15gen'=NA,

             'MST15ori'=NA, 'MST15wei'=NA,
             'MST15max'=NA, 'MST15gen'=NA,

             'cov'=NA
  )


results_save <- list()
for(l in 1:length(eigs)){
  cat(paste0('Eigs, ',eigs[l],' (',trials,'): '))
  data <- data_L2 <-list()
  results <-
    data.frame(
      'MDT15ori'=1:trials, 'MDT15wei'=NA,
      'MDT15max'=NA, 'MDT15gen'=NA,

      'NNL15ori'=NA, 'NNL15wei'=NA,
      'NNL15max'=NA, 'NNL15gen'=NA,

      'MST15ori'=NA, 'MST15wei'=NA,
      'MST15max'=NA, 'MST15gen'=NA,

      'cov'=NA
    )

  for(i in 1:trials){
    set.seed(123456 * (l+1) + i + 1)
    cat(paste0(i,', '))
    sink(nullfile())    # now suppresses

    ## Generate Data
    data <- gen_FD_KL_Expansion(
      ns = c(n*break_location,n*(1-break_location)),
      eigsList = list(c(3,2,1,0.5),
                      eigs[l]*c(3,2,1,0.5)),
      basesList = list(create.bspline.basis(nbasis=4, norder=4),
                       create.bspline.basis(nbasis=4, norder=4)),
      meansList = c(0,0),
      distsArray = c('Normal'),
      kappasArray = 0,
      evals = evalPts,
      silent = T)
    data_L2 <- calculateDistanceMatrix(data, silent = T,
                                       errType='L2', dataUse = 'Orig')

    dat <- data_L2

    MDT15 <- createMDT(dat,15)
    gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)

    NNL15 <- nnl1(dat,15)
    gSeg_NNL15 <- gseg1_new(n, NNL15, pval.perm = T, pval.appr = F)

    MST15 <- mstree(as.dist(dat), ngmax=15)
    gSeg_MST15 <- gseg1_new(n, MST15, pval.perm = T, pval.appr = F)

    cov_sim <- cov_change(data,kappa=0)

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

      ifelse(is.na(cov_sim),0,1)
    )
    sink()
  }

  results$cov <- !results$cov
  results_simplify[l,] <-
    c(eigs[l], as.numeric(colSums(results<=0.05,na.rm = T) /
                            colSums(!is.na(results))))
  results_save[[l]]  <- results

  saveRDS(results_save,
          file=paste0(path,"\\results\\cov_change_full.rds") )
  saveRDS(results_simplify,
          file=paste0(path,"\\results\\cov_change_simple.rds") )
}


#############################################
results_simplify <- readRDS(paste0(path,'\\results\\cov_change_simple.rds'))
library(tidyverse)

# Plot
plot_15tree_covchange <- function(results_simplify){
  data_plot1 <- results_simplify[,-ncol(results_simplify)] %>%
    pivot_longer(cols =MDT15ori:MST15gen)

  data_plot1$graph <- substr(data_plot1$name,1,3)
  data_plot1$stat <- substr(data_plot1$name,6,8)
  data_plot1$graph <- ifelse(data_plot1$graph=="MDT",'MDP',data_plot1$graph)

  # REFERENCE DATA
  ref_data_tmp <- results_simplify[,c(1,ncol(results_simplify))]
  ref_data <- rbind(
    cbind(ref_data_tmp,data.frame('graph'='NNL')),
    cbind(ref_data_tmp,data.frame('graph'='MDP')),
    cbind(ref_data_tmp,data.frame('graph'='MST')))

  ref_data$value <- ref_data$cov
  ref_data$stat <- 'cov'
  ref_data$name <- 'cov'

  # COMBINE AND PLOT
  data_plot <- rbind(
    data_plot1,
    ref_data[,-2])

  ggplot( data_plot ) +
    geom_hline(aes(yintercept=0.05), color='gray',
               linetype='dotted',linewidth=2) +
    facet_grid(cols=vars(graph),scales = 'free_y') +
    geom_line(aes(x=eigs,
                  y = value,
                  group=stat,
                  color=stat,
                  linetype=stat),
              linewidth=2) +
    geom_point(aes(x=eigs,
                   y = value,
                   group=stat,
                   color=stat,
                   shape=stat),
               size=5) +
    scale_color_manual(
      values=c(scales::hue_pal()(4)[c(1,2,4,3)],'black'),
      breaks = c('gen','max','wei','ori','cov'),
      labels=c('Generalized','Maximum','Weighted','Original','HRZ-22')
    ) +
    scale_shape_manual(
      values=c(16,15,17,18,NA),
      breaks = c('gen','max','wei','ori','cov'),
      labels=c('Generalized','Maximum','Weighted','Original','HRZ-22')
    ) +
    scale_linetype_manual(
      values=c('twodash','dashed','longdash','dotdash','dotted'),
      #c(rep('solid',4),'dashed'),
      breaks = c('gen','max','wei','ori','cov'),
      labels=c('Generalized','Maximum','Weighted','Original','HRZ-22')
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
    xlab('Length of Sample') +
    xlab('Eigenvalue Change Size') +
    ylab('Null Rejection Rate')
}

tmp <- plot_15tree_covchange(results_simplify)

png(paste0(path,"\\figures\\cov_change.png"),
    width = 1600,height = 800)
print(tmp)
dev.off()
