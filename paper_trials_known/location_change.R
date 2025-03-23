library(devtools)
load_all()
library(fda)
library(ade4)
library(gSeg)
library(ggplot2)

path <- 'C://Users//jerem//OneDrive//Documents//School//Waterloo//Research//RPackages//fungraphs//paper_trials_known'


trials <- 1000
evalPts <- seq(0,1,0.05)
n <- 100
break_locations <- c(0.05,0.1,0.2,0.25,0.5,0.75,0.8,0.9,0.95)

###########################

results <- data.frame('MDT15ori'=1:trials, 'MDT15wei'=NA,
                      'MDT15max'=NA, 'MDT15gen'=NA,

                      'NNL15ori'=NA, 'NNL15wei'=NA,
                      'NNL15max'=NA, 'NNL15gen'=NA,

                      'MST15ori'=NA, 'MST15wei'=NA,
                      'MST15max'=NA, 'MST15gen'=NA,

                      'mean'=NA
)

results_save <- list()
results_simplify <- data.frame('location'=break_locations,
                               'MDT15ori'=NA, 'MDT15wei'=NA,
                               'MDT15max'=NA, 'MDT15gen'=NA,

                               'NNL15ori'=NA, 'NNL15wei'=NA,
                               'NNL15max'=NA, 'NNL15gen'=NA,

                               'MST15ori'=NA, 'MST15wei'=NA,
                               'MST15max'=NA, 'MST15gen'=NA,

                               'mean'=NA
)


data <- data_dist <- list()

mean_value <- 0.25
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
      meansList = c(0, mean_value),
      distsArray = c('Normal','Normal'),
      evals = evalPts,
      kappasArray = 0,
      silent = T)
    data_dist[[i]] <- calculateDistanceMatrix(data[[i]], silent = T,
                                              errType='L2', dataUse = 'Orig')

    MDT15 <- createMDT(data_dist[[i]],15)
    gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)

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

      ifelse(is.na(mean_change(data[[i]])),0,1)
    )
    sink()
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(break_locations[l], as.numeric(colSums(results<=0.05,na.rm = T)/colSums(!is.na(results))))

  saveRDS(results_save,
          file =paste0(path,"\\results\\location_sim_mean",mean_value,"_full.rds") )
  saveRDS(results_simplify,
          file =paste0(path,"\\results\\location_sim_mean",mean_value,"_simple.rds") )
}

#########################

results <- data.frame('MDT15ori'=1:trials, 'MDT15wei'=NA,
                      'MDT15max'=NA, 'MDT15gen'=NA,

                      'NNL15ori'=NA, 'NNL15wei'=NA,
                      'NNL15max'=NA, 'NNL15gen'=NA,

                      'MST15ori'=NA, 'MST15wei'=NA,
                      'MST15max'=NA, 'MST15gen'=NA,

                      'mean'=NA
)

results_save <- list()
results_simplify <- data.frame('location'=break_locations,
                               'MDT15ori'=NA, 'MDT15wei'=NA,
                               'MDT15max'=NA, 'MDT15gen'=NA,

                               'NNL15ori'=NA, 'NNL15wei'=NA,
                               'NNL15max'=NA, 'NNL15gen'=NA,

                               'MST15ori'=NA, 'MST15wei'=NA,
                               'MST15max'=NA, 'MST15gen'=NA,

                               'mean'=NA
)

data <- data_dist <- list()

mean_value <- 0.5
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
      meansList = c(0, mean_value),
      distsArray = c('Normal','Normal'),
      evals = evalPts,
      kappasArray = 0,
      silent = T)
    data_dist[[i]] <- calculateDistanceMatrix(data[[i]], silent = T,
                                              errType='L2', dataUse = 'Orig')

    MDT15 <- createMDT(data_dist[[i]],15)
    gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)

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

      ifelse(is.na(mean_change(data[[i]])),0,1)
    )
    sink()
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(break_locations[l], as.numeric(colSums(results<=0.05,na.rm = T)/colSums(!is.na(results))))

  saveRDS(results_save,
          file =paste0(path,"\\results\\location_sim_mean",mean_value,"_full.rds") )
  saveRDS(results_simplify,
          file =paste0(path,"\\results\\location_sim_mean",mean_value,"_simple.rds") )
}



# # Plot
#
# library(ggplot2)
# library(tidyverse)
#
# results_simplify$mean <- 1-results_simplify$mean
# results_simplify$cov <- 1-results_simplify$cov
# data_plot <- results_simplify %>% pivot_longer(cols = MDT1ori:cov)
#
# tmp <- ggplot(na.omit(data_plot)) +
#   geom_line(aes(x=location,
#                 y=value,
#                 color=name)) +
#   geom_point(aes(x=location,
#                  y=value,
#                  color=name)) +
#   geom_hline(aes(yintercept=0.05)) +
#   #geom_vline(aes(xintercept=0.5)) +
#   geom_line(aes(x=location,
#                 y=value), data=na.omit(data_plot[data_plot$name=='mean',]),
#             color='black', linewidth=1) +
#   geom_line(aes(x=location,
#                 y=value), data=na.omit(data_plot[data_plot$name=='cov',]),
#             color='black', linewidth=1, linetype='dashed') +
#   guides(color="none") +
#   theme_bw() +
#   theme(axis.title = element_text(size=22),
#         axis.text = element_text(size=18))+
#   xlab('Change Location') +
#   ylab('Null Rejection Rate')
#
# png("C:/Users/jerem/Downloads/location_power_n.png")
# print(tmp)
# dev.off()


############################

# break_locations <- c(0.05,0.1,0.2,0.25,0.5,0.75,0.8,0.9,0.95)
# n <- 100
# trials <- 1000

results_simplify5 <- readRDS(
  file =paste0(path,"\\results\\location_sim_mean",0.5,"_simple.rds") )
results_simplify25 <- readRDS(
  file =paste0(path,"\\results\\location_sim_mean",0.25,"_simple.rds") )


# Plot
library(tidyverse)
plot_15tree_loc <- function(results_simplify5,results_simplify25){
  results_simplify5$mean <- 1-results_simplify5$mean
  results_simplify25$mean <- 1-results_simplify25$mean
  # Get all 15 trees
  results_simplify5_1 <- cbind(results_simplify5[1],
                               results_simplify5[,stringr::str_detect(colnames(results_simplify5),'15')],
                               results_simplify5[,c(ncol(results_simplify5),ncol(results_simplify5)-1)])

  # Format
  data_plot5 <- results_simplify5_1 %>%
    pivot_longer(cols =MDT15ori:MST15gen)
  data_plot25 <- results_simplify25 %>%
    pivot_longer(cols =MDT15ori:MST15gen)

  data_plot5$graph <- substr(data_plot5$name,1,3)
  data_plot5$stat <- substr(data_plot5$name,6,8)
  data_plot5$graph <- ifelse(data_plot5$graph=="MDT",'MDP',data_plot5$graph)
  data_plot25$graph <- substr(data_plot25$name,1,3)
  data_plot25$stat <- substr(data_plot25$name,6,8)
  data_plot25$graph <- ifelse(data_plot25$graph=="MDT",'MDP',data_plot25$graph)

  ref_data_tmp <- results_simplify5[,c(1,ncol(results_simplify5))]
  ref_data5 <- rbind(
    cbind(ref_data_tmp,data.frame('graph'='NNL')),
    cbind(ref_data_tmp,data.frame('graph'='MDP')),
    cbind(ref_data_tmp,data.frame('graph'='MST')))
  ref_data_tmp <- results_simplify25[,c(1,ncol(results_simplify25))]
  ref_data25 <- rbind(
    cbind(ref_data_tmp,data.frame('graph'='NNL')),
    cbind(ref_data_tmp,data.frame('graph'='MDP')),
    cbind(ref_data_tmp,data.frame('graph'='MST')))

  ref_data5$value <- ref_data5$mean
  ref_data5$stat <- 'mean'
  ref_data5$name <- NULL#'mean'
  ref_data25$value <- ref_data25$mean
  ref_data25$stat <- 'mean'
  ref_data25$name <- NULL#'mean'

  data_plot <- rbind(
    data.frame('change'=0.5,data_plot5[,-c(2:3)]),
    data.frame('change'=0.25,data_plot25[,-c(2:3)]),
    data.frame(ref_data25[,-2],'change'=0.25),
    data.frame(ref_data5[,-2],'change'=0.5)
  )


  #cols <- rep(scales::hue_pal()(4)[c(1,2,3,4)],3)

  ggplot( data_plot ) +
    geom_hline(aes(yintercept=0.05), color='gray',
               linetype='dotted',linewidth=2) +
    facet_grid(cols=vars(graph),rows=vars(change),scales = 'free_y') +
    geom_line(aes(x=location,
                  y = value,
                  group=stat,
                  color=stat,
                  linetype=stat),
              linewidth=2) +
    geom_point(aes(x=location,
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
      #values=c(rep('solid',4),'dashed'),
      values=c('twodash','dashed','longdash','dotdash','dotted'),
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
    xlab('Length of Sample') +
    xlab('Change Location') +
    ylab('Null Rejection Rate')
}

tmp <- plot_15tree_loc(results_simplify5,results_simplify25)

png(paste0(path,"\\figures\\location_power_n2.png"),
    width = 1600,height = 800)
print(tmp)
dev.off()

