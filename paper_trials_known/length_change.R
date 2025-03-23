library(devtools)
load_all()
library(fda)
library(ade4)
library(gSeg)
library(ggplot2)

path <- 'C://Users//jerem//OneDrive//Documents//School//Waterloo//Research//RPackages//fungraphs//paper_trials_known'

trials <- 1000
evalPts <- seq(0,1,0.05)
ns <- c(50, 100, 200)
break_location <- 0.5

results <- data.frame('MDT15ori'=1:trials, 'MDT15wei'=NA,
                      'MDT15max'=NA, 'MDT15gen'=NA,

                      'NNL15ori'=NA, 'NNL15wei'=NA,
                      'NNL15max'=NA, 'NNL15gen'=NA,

                      'MST15ori'=NA, 'MST15wei'=NA,
                      'MST15max'=NA, 'MST15gen'=NA,

                      'mean'=NA
)

results_save <- list()
results_simplify <- data.frame('Length'=rep(NA,length(ns)),
                               'MDT15ori'=NA, 'MDT15wei'=NA,
                               'MDT15max'=NA, 'MDT15gen'=NA,

                               'NNL15ori'=NA, 'NNL15wei'=NA,
                               'NNL15max'=NA, 'NNL15gen'=NA,

                               'MST15ori'=NA, 'MST15wei'=NA,
                               'MST15max'=NA, 'MST15gen'=NA,

                               'mean'=NA
)

for(l in 1:length(ns)){

  cat(paste0('Length, ',ns[l],' (',trials,'): '))

  for(i in 1:trials){
    set.seed(123456 * (l+1) + i)
    cat(paste0(i,', '))
    sink(nullfile())    # now suppresses

    ## Generate Data
    data_cp <- gen_FD_KL_Expansion(
      ns = c(ns[l]*break_location,ns[l]*(1-break_location)),
      eigsList = list(c(3,2,1,0.5),
                      c(3,2,1,0.5)),
      basesList = list(create.bspline.basis(nbasis=4, norder=4),
                       create.bspline.basis(nbasis=4, norder=4)),
      meansList = c(0,0.5),
      distsArray = c('Normal','gamma'),
      evals = evalPts,
      kappasArray = 0,
      silent = T)
    data_dist <- calculateDistanceMatrix(data_cp, silent = T,
                                         errType='L2', dataUse = 'Orig')

    MDT15 <- createMDT(data_dist,15)
    gSeg_MDT15 <- gseg1_new(ns[l], MDT15, pval.perm = T, pval.appr = F)

    NNL15 <- nnl1(data_dist,15)
    gSeg_NNL15 <- gseg1_new(ns[l], NNL15, pval.perm = T, pval.appr = F)

    MST15 <- mstree(as.dist(data_dist), ngmax=15)
    gSeg_MST15 <- gseg1_new(ns[l], MST15, pval.perm = T, pval.appr = F)

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

      ifelse(is.na(mean_change(data_cp)),0,1)
    )
    sink()
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(ns[l], as.numeric(colSums(results>0.05,na.rm = T)/colSums(!is.na(results))))

  saveRDS(results_save,
          file=paste0(path,"\\results\\length_change_full.rds") )
  saveRDS(results_simplify,
          file=paste0(path,"\\results\\length_change_simple.rds") )
}





##############################
results_simplify$mean <- 1-results_simplify$mean
results_simplify$cov <- 1-results_simplify$cov
# Get all 15 trees
results_simplify1 <- cbind(results_simplify[1],
                           results_simplify[,stringr::str_detect(colnames(results_simplify),'15')])

# Format
data_plot <- results_simplify1 %>% #[,!stringr::str_detect(colnames(results_simplify1),'ori')] %>%
  pivot_longer(cols =MDT15ori:MST15gen)
data_plot$graph <- substr(data_plot$name,1,3)
data_plot$stat <- substr(data_plot$name,6,8)
data_plot$graph <- ifelse(data_plot$graph=="MDT",'MDP',data_plot$graph)

ref_data_tmp <- results_simplify[,c(1,ncol(results_simplify)-1)]
ref_data <- rbind(
  cbind(ref_data_tmp,data.frame('graph'='NNL')),
  cbind(ref_data_tmp,data.frame('graph'='MDP')),
  cbind(ref_data_tmp,data.frame('graph'='MST')))

#cols <- rep(scales::hue_pal()(4)[c(1,2,4)],3)

ggplot( data_plot ) +
  geom_hline(aes(yintercept=0.05),
             color='gray',linetype='solid',linewidth=1.75) +
  facet_grid(cols=vars(graph)) +
  geom_line(aes(x=length,
                y = value,
                group=stat,
                color=stat,
                linetype=stat),
            linewidth=2) +
  geom_point(aes(x=length,
                 y = value,
                 group=stat,
                 color=stat,
                 shape=stat),
             size=5) +
  geom_line(aes(x=length,
                y=mean),
            data=ref_data,
            color='black', linewidth=1,linetype='dotdash') +
  scale_color_manual(
    values=c(scales::hue_pal()(4)[c(1,2,4,3)]),
    breaks = c('gen','max','wei','ori'),
    labels=c('Generalized','Maximum','Weighted','Original')
  ) +
  scale_shape_manual(
    values=c(16,15,17,18),
    breaks = c('gen','max','wei','ori'),
    labels=c('Generalized','Maximum','Weighted','Original')
  ) +
  scale_linetype_manual(
    values=c('twodash','dashed','longdash','dotdash'),
    #c(rep('solid',4),'dashed'),
    breaks = c('gen','max','wei','ori'),
    labels=c('Generalized','Maximum','Weighted','Original')
  ) +
  theme_bw() +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        strip.text = element_text(size = 24),
        legend.position = 'bottom',
        legend.text = element_text(size=24),
        legend.title = element_blank()) +
  xlab('Length of Sample') +
  ylab('Null Rejection Rate')
