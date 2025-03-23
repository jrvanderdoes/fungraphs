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
means = c(0,0.05,0.1,0.15,0.25,0.5,0.75,1)
dofs <- c(1,3,10)

for(d in 1:length(dofs)){

  results <- data.frame(# 'MDT1ori'=1:trials, 'MDT1wei'=NA,
                        # 'MDT1max'=NA, 'MDT1gen'=NA,
                        # 'MDT3ori'=NA, 'MDT3wei'=NA,
                        # 'MDT3max'=NA, 'MDT3gen'=NA,
                        # 'MDT5ori'=NA, 'MDT5wei'=NA,
                        # 'MDT5max'=NA, 'MDT5gen'=NA,
                        # 'MDT7ori'=NA, 'MDT7wei'=NA,
                        # 'MDT7max'=NA, 'MDT7gen'=NA,
                        'MDT15ori'=1:trials, 'MDT15wei'=NA,
                        'MDT15max'=NA, 'MDT15gen'=NA,

                        # 'NNL1ori'=NA, 'NNL1wei'=NA,
                        # 'NNL1max'=NA, 'NNL1gen'=NA,
                        # 'NNL3ori'=NA, 'NNL3wei'=NA,
                        # 'NNL3max'=NA, 'NNL3gen'=NA,
                        # 'NNL5ori'=NA, 'NNL5wei'=NA,
                        # 'NNL5max'=NA, 'NNL5gen'=NA,
                        # 'NNL7ori'=NA, 'NNL7wei'=NA,
                        # 'NNL7max'=NA, 'NNL7gen'=NA,
                        'NNL15ori'=NA, 'NNL15wei'=NA,
                        'NNL15max'=NA, 'NNL15gen'=NA,

                        # 'MST1ori'=NA, 'MST1wei'=NA,
                        # 'MST1max'=NA, 'MST1gen'=NA,
                        # 'MST3ori'=NA, 'MST3wei'=NA,
                        # 'MST3max'=NA, 'MST3gen'=NA,
                        # 'MST5ori'=NA, 'MST5wei'=NA,
                        # 'MST5max'=NA, 'MST5gen'=NA,
                        # 'MST7ori'=NA, 'MST7wei'=NA,
                        # 'MST7max'=NA, 'MST7gen'=NA,
                        'MST15ori'=NA, 'MST15wei'=NA,
                        'MST15max'=NA, 'MST15gen'=NA,

                        'mean'=NA
  )

  results_save <- list()
  results_simplify <- data.frame('mean'=means,
                                 # 'MDT1ori'=NA, 'MDT1wei'=NA,
                                 # 'MDT1max'=NA, 'MDT1gen'=NA,
                                 # 'MDT3ori'=NA, 'MDT3wei'=NA,
                                 # 'MDT3max'=NA, 'MDT3gen'=NA,
                                 # 'MDT5ori'=NA, 'MDT5wei'=NA,
                                 # 'MDT5max'=NA, 'MDT5gen'=NA,
                                 # 'MDT7ori'=NA, 'MDT7wei'=NA,
                                 # 'MDT7max'=NA, 'MDT7gen'=NA,
                                 'MDT15ori'=NA, 'MDT15wei'=NA,
                                 'MDT15max'=NA, 'MDT15gen'=NA,

                                 # 'NNL1ori'=NA, 'NNL1wei'=NA,
                                 # 'NNL1max'=NA, 'NNL1gen'=NA,
                                 # 'NNL3ori'=NA, 'NNL3wei'=NA,
                                 # 'NNL3max'=NA, 'NNL3gen'=NA,
                                 # 'NNL5ori'=NA, 'NNL5wei'=NA,
                                 # 'NNL5max'=NA, 'NNL5gen'=NA,
                                 # 'NNL7ori'=NA, 'NNL7wei'=NA,
                                 # 'NNL7max'=NA, 'NNL7gen'=NA,
                                 'NNL15ori'=NA, 'NNL15wei'=NA,
                                 'NNL15max'=NA, 'NNL15gen'=NA,

                                 # 'MST1ori'=NA, 'MST1wei'=NA,
                                 # 'MST1max'=NA, 'MST1gen'=NA,
                                 # 'MST3ori'=NA, 'MST3wei'=NA,
                                 # 'MST3max'=NA, 'MST3gen'=NA,
                                 # 'MST5ori'=NA, 'MST5wei'=NA,
                                 # 'MST5max'=NA, 'MST5gen'=NA,
                                 # 'MST7ori'=NA, 'MST7wei'=NA,
                                 # 'MST7max'=NA, 'MST7gen'=NA,
                                 'MST15ori'=NA, 'MST15wei'=NA,
                                 'MST15max'=NA, 'MST15gen'=NA,

                                 'mean'=NA
  )

  for(l in 1:length(means)){

    cat(paste0('Mean, ',means[l],' (',trials,'): '))

    for(i in 1:trials){
      set.seed(123456 * (l+1) + i + 1)
      cat(paste0(i,', '))
      sink(nullfile())    # now suppresses

      ## Generate Data
      data <- gen_FD_KL_Expansion(
        ns = c(n*break_location,n*(1-break_location)),
        eigsList = list(c(3,2,1,0.5),
                        c(3,2,1,0.5)),
        basesList = list(create.bspline.basis(nbasis=4, norder=4),
                         create.bspline.basis(nbasis=4, norder=4)),
        meansList = c(0,means[l]),
        distsArray = c('t','t'),
        dof = c(dofs[d],dofs[d]),
        evals = evalPts,
        kappasArray = 0,
        silent = T)
      data_dist <- calculateDistanceMatrix(data, silent = T,
                                                errType='L2', dataUse = 'Orig')

      # MDT1 <- createMDT(data_dist[[i]],1)
      # MDT3 <- createMDT(data_dist[[i]],3)
      # MDT5 <- createMDT(data_dist[[i]],5)
      # MDT7 <- createMDT(data_dist[[i]],7)
      MDT15 <- createMDT(data_dist,15)
      # gSeg_MDT1 <- gseg1_new(n, MDT1, pval.perm = T, pval.appr = F)
      # gSeg_MDT3 <- gseg1_new(n, MDT3, pval.perm = T, pval.appr = F)
      # gSeg_MDT5 <- gseg1_new(n, MDT5, pval.perm = T, pval.appr = F)
      # gSeg_MDT7 <- gseg1_new(n, MDT7, pval.perm = T, pval.appr = F)
      if(n>30) {
        gSeg_MDT15 <- gseg1_new(n, MDT15, pval.perm = T, pval.appr = F)
      } else{
        gSeg_MDT15 <- list()
        gSeg_MDT15$pval.perm$ori$pval <- NA
        gSeg_MDT15$pval.perm$weighted$pval <- NA
        gSeg_MDT15$pval.perm$max.type$pval <- NA
        gSeg_MDT15$pval.perm$generalized$pval <- NA
      }


      # NNL1 <- nnl1(data_dist[[i]],1)
      # NNL3 <- nnl1(data_dist[[i]],3)
      # NNL5 <- nnl1(data_dist[[i]],5)
      # NNL7 <- nnl1(data_dist[[i]],7)
      NNL15 <- nnl1(data_dist,15)
      # gSeg_NNL1 <- gseg1_new(n, NNL1, pval.perm = T, pval.appr = F)
      # gSeg_NNL3 <- gseg1_new(n, NNL3, pval.perm = T, pval.appr = F)
      # gSeg_NNL5 <- gseg1_new(n, NNL5, pval.perm = T, pval.appr = F)
      # gSeg_NNL7 <- gseg1_new(n, NNL7, pval.perm = T, pval.appr = F)
      if(n>30) {
        gSeg_NNL15 <- gseg1_new(n, NNL15, pval.perm = T, pval.appr = F)
      } else{
        gSeg_NNL15 <- list()
        gSeg_NNL15$pval.perm$ori$pval <- NA
        gSeg_NNL15$pval.perm$weighted$pval <- NA
        gSeg_NNL15$pval.perm$max.type$pval <- NA
        gSeg_NNL15$pval.perm$generalized$pval <- NA
      }

      # MST1 <- mstree(as.dist(data_dist[[i]]), ngmax=1)
      # MST3 <- mstree(as.dist(data_dist[[i]]), ngmax=3)
      # MST5 <- mstree(as.dist(data_dist[[i]]), ngmax=5)
      # MST7 <- mstree(as.dist(data_dist[[i]]), ngmax=7)
      MST15 <- mstree(as.dist(data_dist), ngmax=15)
      # gSeg_MST1 <- gseg1_new(n, MST1, pval.perm = T, pval.appr = F)
      # gSeg_MST3 <- gseg1_new(n, MST3, pval.perm = T, pval.appr = F)
      # gSeg_MST5 <- gseg1_new(n, MST5, pval.perm = T, pval.appr = F)
      # gSeg_MST7 <- gseg1_new(n, MST7, pval.perm = T, pval.appr = F)
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

        ifelse(is.na(mean_change(data)),0,1)#,
        # ifelse(is.na(cov_change(data[[i]])),0,1)
      )
      sink()
    }

    results_save[[l]] <- results
    results_simplify[l,] <-
      c(means[l], as.numeric(colSums(results<=0.05,na.rm = T)/colSums(!is.na(results))))

    saveRDS(results_save,
            file=paste0(path,"\\results\\t",dofs[d],"_change_full.rds") )
    saveRDS(results_simplify,
            file=paste0(path,"\\results\\t",dofs[d],"_change_simple.rds") )
  }
}

########################
library(tidyverse)

data_simple <-
  rbind(cbind(data.frame('dof'=dofs[1]),
              readRDS(paste0(path,"\\results\\t",dofs[1],"_change_simple.rds")) ),
        cbind(data.frame('dof'=dofs[2]),
              readRDS(paste0(path,"\\results\\t",dofs[2],"_change_simple.rds")) ),
        cbind(data.frame('dof'=dofs[3]),
              readRDS(paste0(path,"\\results\\t",dofs[3],"_change_simple.rds")) )
  )
data_simple$mean.1 <- 1-data_simple$mean.1

# Format
data_plot <- data_simple[,-ncol(data_simple)] %>%
  pivot_longer(cols =MDT15ori:MST15gen)
data_plot$graph <- substr(data_plot$name,1,3)
data_plot$stat <- str_replace(str_replace(data_plot$name,'[0-9]+',''),'(NNL)|(MST)|(MDT)','')
data_plot$graph <- ifelse(data_plot$graph=="MDT",'MDP',data_plot$graph)

ref_data_tmp <- data_simple[,c(1,2,ncol(data_simple))]
ref_data <- rbind(
  cbind(ref_data_tmp,data.frame('graph'='NNL')),
  cbind(ref_data_tmp,data.frame('graph'='MDP')),
  cbind(ref_data_tmp,data.frame('graph'='MST')))

ref_data$value <- ref_data$mean.1
ref_data$mean.1 <- NULL
ref_data$stat <- 'mean'
ref_data$name <- 'mean'

data_plot <- rbind(
  data_plot,
  data.frame(ref_data)
)
data_plot$dof <- paste0('DOF ',data_plot$dof)
# data_plot$dof <- factor(data_plot$dof,levels = c('DOF 1','DOF 3','DOF 10'))
data_plot$dof <- factor(data_plot$dof,levels = c('DOF 10','DOF 3','DOF 1'))


plt <- ggplot(na.omit(data_plot[data_plot$name!='mean',])) +
  facet_grid(cols=vars(graph),rows=vars(dof)) +
  geom_hline(aes(yintercept=0.05), color='gray',linetype='dotted',linewidth=2) +
  geom_line(aes(x=mean,
                y=value,
                color=stat,
                group=stat,
                linetype=stat ),linewidth=2) +
  geom_point(aes(x=mean,
                 y=value,
                 color=stat,
                 group=stat ),size=4) +
  #geom_vline(aes(xintercept=0)) +
  geom_line(aes(x=mean,
                y=value,
                col=stat,
                group=stat,
                linetype=stat ), data=na.omit(data_plot[data_plot$name=='mean',]),
            linewidth=2) +
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
  xlim(c(0,1)) +
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
  xlab('Mean Change Size') +
  ylab('Null Rejection Rate')

png(paste0(path,"/figures/mean_power_ts.png"),
    width = 1600,height = 800)
print(plt)
dev.off()

