
plot_15tree_dep <- function(results_simplify){
  results_simplify$mean <- 1-results_simplify$mean
  results_simplify$cov <- 1-results_simplify$cov
  # Get all 15 trees
  results_simplify1 <- cbind(results_simplify[1],
                             results_simplify[,stringr::str_detect(colnames(results_simplify),'15')])

  # Format
  data_plot <- results_simplify1[,!stringr::str_detect(colnames(results_simplify1),'ori')] %>%
    pivot_longer(cols =MDT15wei:MST15gen)
  data_plot$graph <- substr(data_plot$name,1,3)
  data_plot$stat <- substr(data_plot$name,6,8)
  data_plot$graph <- ifelse(data_plot$graph=="MDT",'MDP',data_plot$graph)

  ref_data_tmp <- results_simplify[,c(1,ncol(results_simplify)-1)]
  ref_data <- rbind(
    cbind(ref_data_tmp,data.frame('graph'='NNL')),
    cbind(ref_data_tmp,data.frame('graph'='MDP')),
    cbind(ref_data_tmp,data.frame('graph'='MST')))

  cols <- rep(scales::hue_pal()(4)[c(1,2,4)],3)

  ggplot( data_plot ) +
    geom_hline(aes(yintercept=0.05), color='gray',linetype='dotted',linewidth=1.75) +
    facet_grid(cols=vars(graph)) +
    geom_line(aes(x=dep,
                  y = value,
                  group=name,
                  color=name),
              linewidth=2) +
    geom_point(aes(x=dep,
                   y = value,
                   group=name,
                   color=name),
               size=4) +
    geom_line(aes(x=dep,
                  y=mean),
              data=ref_data,
              color='black', linewidth=1,linetype='dotdash') +
    scale_color_manual( values=cols ) +
    guides(color="none") +
    theme_bw() +
    theme(axis.title = element_text(size=28),
          axis.text = element_text(size=24),
          strip.text = element_text(size = 24))+
    xlab('Dependence Change') +
    ylab('Null Rejection Rate')
}

##########################

scale_color_manual(
  values=scales::hue_pal()(4)[c(1,2,4,3)],
  breaks = c('gen','max','wei','ori'),
  labels=c('Generalized','Maximum','Weighted','Original')
) +
  scale_shape_manual(
    values=c(16,15,17,18),
    breaks = c('gen','max','wei','ori'),
    labels=c('Generalized','Maximum','Weighted','Original')
  ) +

  ,
shape=stat),
size=5) +

  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        strip.text = element_text(size = 24),
        legend.position = 'bottom',
        legend.text = element_text(size=24),
        legend.title = element_blank())+
  xlab('Change Location') +


  ##################
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



  ggplot( data_plot ) +
    geom_hline(aes(yintercept=0.05), color='gray',linetype='dotted',linewidth=1.75) +
    facet_grid(cols=vars(graph)) +
    geom_line(aes(x = dep,
                  y = sizeAdjPower,
                  group=stat,
                  color=stat),
              linewidth=2) +
    geom_point(aes(x=dep,
                   y = sizeAdjPower,
                   group=stat,
                   color=stat,
                   shape=stat),
               size=5) +
    geom_line(aes(x=dep,
                  y=sizeAdjPower),
              data=ref_data,
              color='black', linewidth=1,linetype='dotdash') +
    scale_color_manual(
      values=scales::hue_pal()(4)[c(1,2,4,3)],
      breaks = c('gen','max','wei','ori'),
      labels=c('Generalized','Maximum','Weighted','Original')
    ) +
    scale_shape_manual(
      values=c(16,15,17,18),
      breaks = c('gen','max','wei','ori'),
      labels=c('Generalized','Maximum','Weighted','Original')
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

################################

ggplot( data_plot ) +
  geom_hline(aes(yintercept=0.05), color='gray',linetype='dotted',linewidth=1.75) +
  facet_grid(cols=vars(graph),rows=vars(Type)) +
  geom_line(aes(x=mean,
                y = value,
                group=stat,
                color=stat),
            linewidth=2) +
  geom_point(aes(x=mean,
                 y = value,
                 group=stat,
                 color=stat,
                 shape=stat),
             data_plot[data_plot$name!='mean.1',],
             size=5) +
  geom_line(aes(x=mean,
                y=mean.1),
            data=ref_data,
            color='black', linewidth=1,linetype='dotdash') +
  scale_color_manual(
    values=scales::hue_pal()(4)[c(1,2,4)],
    breaks = c('gen','max','wei'),
    labels=c('Generalized','Maximum','Weighted')
  ) +
  scale_shape_manual(
    values=c(16,15,17,18),
    breaks = c('gen','max','wei','ori'),
    labels=c('Generalized','Maximum','Weighted','Original')
  ) +
  theme_bw() +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        strip.text = element_text(size = 24),
        legend.position = 'bottom',
        legend.text = element_text(size=24),
        legend.title = element_blank())+
  xlab('Mean Change Size') +
  ylab('Null Rejection Rate')

############################
