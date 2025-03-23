library(devtools)
load_all()
library(fda)
library(ade4)
library(gSeg)
library(ggplot2)

path <- 'C://Users//jerem//OneDrive//Documents//School//Waterloo//Research//RPackages//fungraphs//paper_trials_known'

graph_change_o <- function(data){
  tmp <- gseg1(n, tree, pval.perm = T, B=100, pval.appr = F)

  ifelse(tmp$pval.perm$ori$pval<0.05,tmp$scanZ$ori$tauhat,NA)
}
graph_change_w <- function(data){
  tmp <- gseg1(n, tree, pval.perm = T, B=100, pval.appr = F)

  ifelse(tmp$pval.perm$wei$pval<0.05,tmp$scanZ$wei$tauhat,NA)
}
graph_change_m <- function(data){
  tmp <- gseg1(n, tree, pval.perm = T, B=100, pval.appr = F)

  ifelse(tmp$pval.perm$ori$pval<0.05,tmp$scanZ$ori$tauhat,NA)
}
graph_change_g <- function(data){
  tmp <- gseg1(n, tree, pval.perm = T, B=100, pval.appr = F)

  ifelse(tmp$pval.perm$ori$pval<0.05,tmp$scanZ$ori$tauhat,NA)
}

#########
bs_func_basic <- function(data, method, addAmt=0, silent=F, ...){
  potential_cp <- method(data, ...)

  # No Change Point Detected
  if(is.na(potential_cp)) return()

  # Display progress
  if(!silent)
    cat(paste0('ChangePoint Detected (',1+addAmt,'-' ,addAmt+ncol(data),' at ',
               addAmt+potential_cp,'): Segment Data and Re-Search\n'))

  return(c(
    bs_func_basic(data=data[,1:potential_cp],
            method=method,
            addAmt = addAmt,
            silent=silent,
            ...),
    potential_cp + addAmt,
    bs_func_basic(data=data[,(potential_cp+1):ncol(data)],
            method=method,
            addAmt = addAmt+potential_cp,
            silent=silent,
            ...)
  ))
}


ver_func_basic <- function(data, cps, method, ...){
  if(length(cps)==0){
    gc <- method(data,...)
    if(is.na(gc)) return()

    return(gc)
  }


  cps_new <- c()
  cps_ver <- c(0,cps,ncol(data))

  for(i in 2:(length(cps_ver)-1)){
    dat <- data[,(cps_ver[i-1]+1):cps_ver[i+1]]

    gc <- method(dat,...)
    if(!is.na(gc))
      cps_new <- c(cps_new,gc+cps_ver[i-1])
  }

  cps_new[!is.na(cps_new)]
}

############
graph_method_gen <- function(data,n,nTrees,treeType){
  if(treeType=='MDT'){
    tree <- createMDT(data,nTrees)
  }else if(treeType=='NNL'){
    tree <- nnl1(data,nTrees)
  }else if(treeType=='MST'){
    tree <- mstree(as.dist(data), ngmax=nTrees)
  }else{
    stop('treeType bad')
  }

  graph <- tryCatch({
    gseg1_new(n, tree,
                       pval.perm = T, pval.appr = F)
  },error=function(e){
    graph <- list()
    graph$pval.perm$max.type$pval <- 1
    graph
  })

  ifelse(graph$pval.perm$max.type$pval<0.05,
         graph$scanZ$max.type$tauhat,
         NA)
}

bs_func <- function(data, method, addAmt=0, silent=F, ...){
  data_dist <- calculateDistanceMatrix(
    data, silent = T, errType='L2', dataUse = 'Orig')

  potential_cp <- method(data=data_dist,n=ncol(data), ...)

  # No Change Point Detected
  if(is.na(potential_cp)) return()

  # Display progress
  if(!silent)
    cat(paste0('ChangePoint Detected (',1+addAmt,'-' ,addAmt+ncol(data),' at ',
               addAmt+potential_cp,'): Segment Data and Re-Search\n'))

  return(c(
    bs_func(data=data[,1:potential_cp],
            method=method,
            addAmt = addAmt,
            silent=silent,
            ...),
    potential_cp + addAmt,
    bs_func(data=data[,(potential_cp+1):ncol(data)],
            method=method,
            addAmt = addAmt+potential_cp,
            silent=silent,
            ...)
  ))
}

ver_func <- function(data, cps, method, ...){
  if(length(cps)==0){
    data_dist <- calculateDistanceMatrix(
      data=data, silent = T,
      errType='L2', dataUse = 'Orig')

    gc <- method(data=data_dist,n=ncol(data), ...)
    if(is.na(gc)) return()

    return(gc)
  }
  cps_new <- c()
  cps_ver <- c(0,cps,ncol(data))

  for(i in 2:(length(cps_ver)-1)){
    dat <- data[,(cps_ver[i-1]+1):cps_ver[i+1]]


    data_dist <- calculateDistanceMatrix(
      data=dat, silent = T,
      errType='L2', dataUse = 'Orig')

    gc <- method(data=data_dist,n=ncol(dat), ...)
    if(!is.na(gc))
      cps_new <- c(cps_new,gc+cps_ver[i-1])
  }

  cps_new[!is.na(cps_new)]
}

evalPts <- seq(0,1,0.05)
nSims <- 100
data <- list()
mst5_pts <- nnl5_pts <- mdt5_pts <-
  mean_pts <- cov_pts <- c()

change_locs <- c(30,70,100,50)
meanValues<- c(0,0.75,0.75, 0.75)
distValues <- c('Normal','Normal',
                'Exponential','Exponential')
eigMults <- c(1,1,1,3)

data_list <- list()
for(i in 1:nSims){
  set.seed(12345678 * i)
  cat(paste0(i,', '))
  sink(nullfile())    # now suppresses

  ## Generate Data
  data <- gen_FD_KL_Expansion(
    ns = change_locs,
    eigsList = list(eigMults[1]*c(3,2,1,0.5),
                    eigMults[2]*c(3,2,1,0.5),
                    eigMults[3]*c(3,2,1,0.5),
                    eigMults[4]*c(3,2,1,0.5)),
    basesList = list(create.bspline.basis(nbasis=4, norder=4),
                     create.bspline.basis(nbasis=4, norder=4),
                     create.bspline.basis(nbasis=4, norder=4),
                     create.bspline.basis(nbasis=4, norder=4)),
    meansList = meanValues,
    distsArray = distValues,
    evals = evalPts,
    kappasArray = 0,
    silent = T)
  data_list[[i]] <- data

  # MST5
  tmp <- bs_func(data = data,method=graph_method_gen,
                 nTrees=5,treeType='MST')
  pts <- ver_func(data = data,
                  cps = tmp,
                  method=graph_method_gen,
                  nTrees=5,
                  treeType='MST')
  mst5_pts <- c(mst5_pts,pts)

  # NNL
  tmp <- bs_func(data = data,method=graph_method_gen,
                 nTrees=5,treeType='NNL')
  pts <- ver_func(data = data,
                  cps = tmp,
                  method=graph_method_gen,
                  nTrees=5,
                  treeType='NNL')
  nnl5_pts <- c(nnl5_pts,pts)

  # MDT
  tmp <- bs_func(data = data,method=graph_method_gen,
                 nTrees=5,treeType='MDT')
  pts <- ver_func(data = data,
                  cps = tmp,
                  method=graph_method_gen,
                  nTrees=5,
                  treeType='MDT')
  mdt5_pts <- c(mdt5_pts,pts)

  # # Mean
  # tmp <- bs_func_basic(data = data,method=mean_change)
  # pts <- ver_func_basic(data = data,
  #                       method=mean_change,
  #                       cps = tmp)
  # mean_pts <- c(mean_pts,pts)
  #
  # # Cov
  # tmp <- bs_func_basic(data = data,method=cov_change)
  # pts <- ver_func_basic(data = data,
  #                       method=cov_change,
  #                       cps = tmp)
  # cov_pts <- c(cov_pts,pts)

  sink()

  data_list[[i]] <- data
}
saveRDS(data_list,file =paste0(path,"\\data\\mult_change.rds"))
saveRDS(mst5_pts,file = paste0(path,"\\results\\mult_change_mst.rds"))
saveRDS(mdt5_pts,file = paste0(path,"\\results\\mult_change_mdt.rds"))
saveRDS(nnl5_pts,file = paste0(path,"\\results\\mult_change_nnl.rds"))
# ce_pts <- readRDS(paste0(path,"\\results\\mult_change_ce.rds"))
# saveRDS(ce_pts,file = paste0(path,"\\results\\mult_change_ce.rds"))

#opts<- options()
#options(warn=2)
##################################
change_locs <- c(30,70,100,50)
meanValues<- c(0,0.75,0.75, 0.75)
distValues <- c('Normal','Normal',
                'Exponential','Exponential')
eigMults <- c(1,1,1,3)

### Hist
data_plot <- data.frame(type='MDT-5',
                        pts=readRDS(paste0(path,"\\results\\mult_change_mdt.rds")))
data_plot <- rbind(data_plot,
                   data.frame(type='NNL-5',
                              pts=readRDS(paste0(path,"\\results\\mult_change_nnl.rds"))))
data_plot <- rbind(data_plot,
                   data.frame(type='MST-5',
                              pts=readRDS(paste0(path,"\\results\\mult_change_mst.rds"))))
data_plot <- rbind(data_plot,
                   data.frame(type='HRV-25',
                              pts=readRDS(paste0(path,"\\results\\mult_change_ce_s.rds"))))
# data_plot <- rbind(data_plot,
#                    data.frame(type='MEAN',
#                               pts=mean_pts))
# data_plot <- rbind(data_plot,
#                    data.frame(type='COV',
#                               pts=cov_pts))
data_plot[data_plot$type=='MDT-5',"type"] <- 'MDP-5'
data_plot$type <- factor(data_plot$type,
                         levels=c('MDP-5','NNL-5','MST-5',
                                  'HRV-25'))
# saveRDS(data_plot,'C:/Users/jerem/OneDrive/Documents/School/Waterloo/Research/GraphCP/Data/Simulations/multiple_change_simulation_plot.rds')
detected_plot <-
  ggplot(data_plot,aes(x=pts)) +
  geom_vline(aes(xintercept=cumsum(change_locs)[1]),
             color='red', linetype='dashed') +
  geom_vline(aes(xintercept=cumsum(change_locs)[2]),
             color='red', linetype='dashed') +
  geom_vline(aes(xintercept=cumsum(change_locs)[3]),
             color='red', linetype='dashed') +
  # geom_vline(aes(xintercept=cumsum(change_locs)[4]),
  #            color='red', linetype='dashed') +
  facet_grid( rows=vars(type)) +
  geom_histogram(bins=200) +
  theme_bw() +
  # theme(legend.title = element_blank(),
  #       #plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
  #       axis.title = ggplot2::element_text(size=20),
  #       axis.text = ggplot2::element_text(size=14),
  #       legend.text = element_text(size=12),
  #       strip.text = element_text(size = 14)) +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        strip.text = element_text(size = 14),
        legend.position = 'bottom',
        legend.text = element_text(size=14),
        legend.title = element_blank())+
  xlab('Length of Sample') +
  xlab('Location of Detected Change') +
  ylab('Observation Count') +
  ylim(c(0,105)) +
  coord_cartesian(xlim=c(0,sum(change_locs)))#+
  #xlim(c(0,sum(change_locs)))

# ### Overview
# change_locs <- c(30,70,50,30,20)
# meanValues<- c(0,0.75,0.75, 0.75,0)
# distValues <- c('Normal','Normal',
#                 'Exponential','Exponential','Normal')
# eigMults <- c(1,1,1,3,1)


means <- rep(meanValues,times=change_locs)
vars <- rep(eigMults/4,times=change_locs)
dists <- rep(distValues,times=change_locs)
groups <- rep(1:length(change_locs),times=change_locs)

xs <- 1:sum(change_locs)
true_plot <- ggplot(data=data.frame(type=rep('Model',length(means))),
                    mapping = aes(x=xs)) +
  facet_grid( rows=vars(type)) +
  geom_line(aes(y=means,
                group=groups,
                linetype=dists,
                color=dists),
            linewidth=1.5) +
  geom_ribbon(aes(ymin=means-vars,
                  ymax=means+vars,
                  fill=dists,
                  group=groups),
              alpha=0.25) +
  theme_bw() +
  # theme(legend.title = element_blank(),
  #       #plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
  #       axis.title = ggplot2::element_text(size=20),
  #       axis.text = ggplot2::element_text(size=14),
  #       legend.text = element_text(size=12),
  #       strip.text = element_text(size = 14),
  #       legend.position = 'bottom') +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        strip.text = element_text(size = 14),
        legend.position = 'bottom',
        legend.text = element_text(size=14),
        legend.title = element_blank())+
  xlab('Length of Sample') +
  ylab(NULL)  +
  coord_cartesian(xlim=c(0,sum(change_locs))) +
  #xlim(0,sum(change_locs)) +
  xlab(NULL)#'Location in Data')

result <- patchwork::wrap_plots(true_plot,detected_plot) +
  patchwork::plot_layout(nrow = 2, guides = "collect",
                         heights = c(1,5)) +
  patchwork::plot_annotation(
    title = NULL,#'Multiple Change Point',
    # subtitle = subtitle,
    theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=22),
                           plot.subtitle = ggplot2::element_text(hjust = 0.5))
  ) &
  theme(legend.position='top')
  #ggplot2::theme(legend.position='none')

# png('C:/Users/jerem/Downloads/MultipleChange.png',
png(paste0(path,'//figures//MultipleChange.png'),
    width = 800,height = 800)
    #width = 1.5*480,height = 480)
print(result)
dev.off()

