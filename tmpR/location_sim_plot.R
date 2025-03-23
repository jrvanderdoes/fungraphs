break_locations <- c(0.05,0.1,0.2,0.25,0.5,0.75,0.8,0.9,0.95)
n <- 100
trials <- 1000

results_simplify5 <- readRDS( paste0('C:/Users/jerem/OneDrive/Documents/',
                                    'School/Waterloo/Research/',
                                    'GraphCP/Data/Simulations/power/CL_',
                                    'trials',trials,
                                    '_ns',n,
                                    '_bl',break_locations[length(break_locations)],
                                    '_dNorm',
                                    '_results_simp.rds'))
results_simplify25 <- readRDS('C:\\Users\\jerem\\Downloads\\updatedPlots\\sims\\simp025.rds')


# Plot
plot_15tree_loc <- function(results_simplify5,results_simplify25){
  results_simplify5$mean <- 1-results_simplify5$mean
  results_simplify5$cov <- 1-results_simplify5$cov
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

  ref_data_tmp <- results_simplify5[,c(1,ncol(results_simplify5)-1)]
  ref_data5 <- rbind(
    cbind(ref_data_tmp,data.frame('graph'='NNL')),
    cbind(ref_data_tmp,data.frame('graph'='MDP')),
    cbind(ref_data_tmp,data.frame('graph'='MST')))
  ref_data_tmp <- results_simplify25[,c(1,ncol(results_simplify25)-1)]
  ref_data25 <- rbind(
    cbind(ref_data_tmp,data.frame('graph'='NNL')),
    cbind(ref_data_tmp,data.frame('graph'='MDP')),
    cbind(ref_data_tmp,data.frame('graph'='MST')))

  ref_data5$value <- ref_data5$mean
  ref_data5$stat <- 'mean'
  ref_data5$name <- 'mean'
  ref_data25$value <- ref_data25$mean
  ref_data25$stat <- 'mean'
  ref_data25$name <- 'mean'

  data_plot <- rbind(
    data.frame('change'=0.5,data_plot5[,-c(2:3)]),
    data.frame('change'=0.25,data_plot25[,-c(2:3)]),
    data.frame(ref_data25[,-2],'change'=0.25),
    data.frame(ref_data5[,-2],'change'=0.5)
  )


  #cols <- rep(scales::hue_pal()(4)[c(1,2,3,4)],3)

  ggplot( data_plot ) +
    geom_hline(aes(yintercept=0.05), color='gray',
               linetype='solid',linewidth=1.75) +
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
    theme(axis.title = element_text(size=28),
          axis.text = element_text(size=24),
          strip.text = element_text(size = 24),
          legend.position = 'bottom',
          legend.text = element_text(size=24),
          legend.title = element_blank())+
    xlab('Change Location') +
    ylab('Null Rejection Rate')
}

tmp <- plot_15tree_loc(results_simplify5,results_simplify25)

png("C:/Users/jerem/Downloads/location_power_n2.png",width = 1600,height = 800)
print(tmp)
dev.off()
