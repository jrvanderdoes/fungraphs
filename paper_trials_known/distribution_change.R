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


results_save <- list()
dists = c(5,3,1,0.5,0.1)
results_simplify <- data.frame('dist'=dists,
                               'MDT15ori'=NA, 'MDT15wei'=NA,
                               'MDT15max'=NA, 'MDT15gen'=NA,

                               'NNL15ori'=NA, 'NNL15wei'=NA,
                               'NNL15max'=NA, 'NNL15gen'=NA,

                               'MST15ori'=NA, 'MST15wei'=NA,
                               'MST15max'=NA, 'MST15gen'=NA,

                               'ce'=NA
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

                        'ce'=NA
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
      distsArray = c('Normal','gamma'),
      shape = c(0,dists[l]),
      evals = evalPts,
      kappasArray = 0,
      silent = T)
    data_dist[[i]] <- calculateDistanceMatrix(data[[i]], silent = T,
                                              errType='L2', dataUse = 'Orig')

    MDT15 <- createMDT(data_dist[[i]],15)
    gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)

    # Has issues with non-unique
    NNL15 <- nnl1(data_dist[[i]],15)
    gSeg_NNL15 <- gseg1_new(n, NNL15, pval.perm = T, pval.appr = F)

    MST15 <- mstree(as.dist(data_dist[[i]]), ngmax=15)
    gSeg_MST15 <- gseg1_new(n, MST15, pval.perm = T, pval.appr = F)

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

      NA
    )
    sink()
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(dists[l], as.numeric(colSums(results<=0.05,na.rm = T)/colSums(!is.na(results))))

  saveRDS(data,
          file=paste0(path,"\\results\\distribution_change_data_",l,".rds") )
  saveRDS(results_save,
          file=paste0(path,"\\results\\distribution_change_full.rds") )
  saveRDS(results_simplify,
          file=paste0(path,"\\results\\distribution_change_simple.rds") )
}

for(l in 1:length(dists)){
  results_simplify$ce[l] <-
    readRDS(paste0(path,"\\results\\distribution_change_simple_ce",l,"_s.rds"))$ce[l]
}
saveRDS(results_simplify,
        file=paste0(path,"\\results\\distribution_change_simple_comb.rds") )

### Plot
library(tidyverse)

plot_15tree_dist <- function(results_simplify){
  # Get all 15 trees
  results_simplify1 <- cbind(results_simplify[,1:2],
                             results_simplify[,stringr::str_detect(colnames(results_simplify),'15')])

  # Format
  data_plot <- results_simplify1 %>% #[,!stringr::str_detect(colnames(results_simplify1),'ori')] %>%
    pivot_longer(cols =MDT15ori:MST15gen)
  data_plot$graph <- substr(data_plot$name,1,3)
  data_plot$stat <- substr(data_plot$name,6,8)
  data_plot$graph <- ifelse(data_plot$graph=="MDT",'MDP',data_plot$graph)

  ref_data_tmp <- results_simplify[,c(1,ncol(results_simplify))]
  ref_data <- rbind(
    cbind(ref_data_tmp,data.frame('graph'='NNL')),
    cbind(ref_data_tmp,data.frame('graph'='MDP')),
    cbind(ref_data_tmp,data.frame('graph'='MST')))
  ref_data$name <- 'hrv25'
  ref_data$stat <- 'hrv25'
  ref_data$value <- ref_data$ce
  ref_data$ce <- NULL

  data_plot <- rbind(data_plot,ref_data)

  # cols <- rep(scales::hue_pal()(4)[c(1,2,4)],3)

  ggplot( data_plot ) +
    geom_hline(aes(yintercept=0.05), color='gray',linetype='dotted',linewidth=2) +
    facet_grid(rows=vars(graph)) +
    geom_line(aes(x=dist,
                  y = value,
                  group=stat,
                  color=stat,
                  linetype=stat),
              linewidth=2) +
    geom_point(aes(x=dist,
                   y = value,
                   group=stat,
                   color=stat,
                   shape=stat),
               size=5) +
    # geom_line(aes(x=dist,
    #               y=ce,
    #               group=stat,
    #               color=stat),
    #           data=ref_data,linewidth=1,linetype='dotdash') +
    scale_color_manual(
      values=c(scales::hue_pal()(4)[c(1,2,4,3)],'black'),
      breaks = c('gen','max','wei','ori','hrv25'),
      labels=c('Generalized','Maximum','Weighted','Original','HRV-25')
    ) +
    scale_shape_manual(
      values=c(16,15,17,18,NA),
      breaks = c('gen','max','wei','ori','hrv25'),
      labels=c('Generalized','Maximum','Weighted','Original','HRV-25')
    ) +
    scale_linetype_manual(
      values=c('twodash','dashed','longdash','dotdash','dotted'),
      #c(rep('solid',4),'dashed'),
      breaks = c('gen','max','wei','ori','hrv25'),
      labels=c('Generalized','Maximum','Weighted','Original','HRV-25')
    ) +
    theme_bw() +
    ylim(c(0,1)) +
    # theme(axis.title = element_text(size=42),
    #       axis.text = element_text(size=38),
    #       strip.text = element_text(size = 24),
    #       legend.position = 'bottom',
    #       legend.text = element_text(size=38),
    #       legend.title = element_blank()) +
    theme(axis.title = element_text(size=30),
          axis.text = element_text(size=26),
          strip.text = element_text(size = 26),
          legend.position = 'bottom',
          legend.text = element_text(size=26),
          legend.title = element_blank())+
    xlab('Length of Sample') +
    xlab('Shape Parameter') +
    ylab('Null Rejection Rate') +
    scale_x_reverse(breaks = 0:10)
}

results_simplify <- readRDS(
  file=paste0(path,"\\results\\distribution_change_simple_comb.rds") )
# Plot
tmp <- plot_15tree_dist(results_simplify)


png(paste0(path,"\\figures\\dist_power.png")
    ,width = 800,height = 1600)
print(tmp)
dev.off()
