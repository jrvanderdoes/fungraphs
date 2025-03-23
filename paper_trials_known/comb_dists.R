results <- readRDS(file=paste0(path,"\\results\\distribution_change_simple_comb.rds") )
results_f <- readRDS(file=paste0(path,"\\results\\distribution_change_simple_new_comb.rds"))

results_comb <- rbind(cbind(data.frame(type='B-spline'),results),
                      cbind(data.frame(type='Fourier'),results_f))

plot_15tree_dist <- function(results_comb){
  # Get all 15 trees
  results_comb1 <- cbind(results_comb[,1:2],
                             results_comb[,stringr::str_detect(colnames(results_comb),'15')])

  # Format
  data_plot <- results_comb1 %>%
    pivot_longer(cols =MDT15ori:MST15gen)
  data_plot$graph <- substr(data_plot$name,1,3)
  data_plot$stat <- substr(data_plot$name,6,8)
  data_plot$graph <- ifelse(data_plot$graph=="MDT",'MDP',data_plot$graph)

  ref_data_tmp <- results_comb[,c(1,2,ncol(results_comb))]
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
    facet_grid(cols=vars(graph), rows=vars(type)) +
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
    scale_x_continuous(trans = c("log10", "reverse"))
  #   scale_x_reverse(breaks = 0:10)
  # scale_x_log10() +
}



png(paste0(path,"\\figures\\dist_power_comb.png")
    ,width = 1600,height = 800)
plot_15tree_dist(results_comb)
dev.off()
