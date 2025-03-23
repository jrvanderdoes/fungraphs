
path <- 'C://Users//jerem//OneDrive//Documents//School//Waterloo//Research//RPackages//fungraphs//paper_trials_known'


results_save <- list()
dists = c(5,3,1,0.5,0.1)
results_simplify <- data.frame('dist'=dists,
                               'ce'=NA
)

for(l in 1:length(dists)){
  data <- readRDS(file=paste0(path,"\\results\\distribution_change_data_",l,".rds") )

  results <- rep(NA,length(data))
  for(i in 1:length(data)){
    set.seed(1234 * i)
    cat(paste0(i,', '))

    tmp <- change(data[[i]],method='characteristic', statistic='Tn', critical='simulation',
                  type='single')
    results[i] <- tmp$pvalue
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(dists[l], as.numeric(sum(results<=0.05,na.rm = T)/sum(!is.na(results))))

  saveRDS(results_save,
          file=paste0(path,"\\results\\distribution_change_full_ce.rds") )
  saveRDS(results_simplify,
          file=paste0(path,"\\results\\distribution_change_simple_ce.rds") )
}

##########################

path <- 'C://Users//jerem//OneDrive//Documents//School//Waterloo//Research//RPackages//fungraphs//paper_trials_known'


results_save <- list()
dists = c(5,3,1,0.5,0.1)
results_simplify <- data.frame('dist'=dists,
                               'ce'=NA
)

for(l in 1:length(dists)){
  data <- readRDS(file=paste0(path,"\\results\\distribution_change_data_",l,"_new.rds") )

  results <- rep(NA,length(data))
  for(i in 1:length(data)){
    set.seed(1234 * i)
    cat(paste0(i,', '))

    tmp <- change(data[[i]],method='characteristic', statistic='Tn',
                  critical='simulation', type='single')
    results[i] <- tmp$pvalue
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(dists[l], as.numeric(sum(results<=0.05,na.rm = T)/sum(!is.na(results))))

  saveRDS(results_save,
          file=paste0(path,"\\results\\distribution_change_full_ce_new.rds") )
  saveRDS(results_simplify,
          file=paste0(path,"\\results\\distribution_change_simple_ce_new.rds") )
}


##################################

path <- 'C://Users//jerem//OneDrive//Documents//School//Waterloo//Research//RPackages//fungraphs//paper_trials_known'


results_save <- list()
dists = c(5,3,1,0.5,0.1)
results_simplify <- data.frame('dist'=dists,
                               'ce'=NA
)

for(l in 1:length(dists)){
  data <- readRDS(file=paste0(path,"\\results\\distribution_change_data_",l,".rds") )

  results <- rep(NA,length(data))
  for(i in 1:length(data)){
    set.seed(1234 * i)
    cat(paste0(i,', '))

    tmp <- change(data[[i]],method='characteristic', statistic='Tn',
                  critical='simulation', type='single')
    results[i] <- tmp$pvalue
  }

  results_save[[l]] <- results
  results_simplify[l,] <-
    c(dists[l], as.numeric(sum(results<=0.05,na.rm = T)/sum(!is.na(results))))

  saveRDS(results_save,
          file=paste0(path,"\\results\\distribution_change_full_ce",l,"_s.rds") )
  saveRDS(results_simplify,
          file=paste0(path,"\\results\\distribution_change_simple_ce",l,"_s.rds") )
}
##################################

## ELECTRICITY
errs <- projection_model(electricity,model = 'ets',TVE = 0.99,n.ahead = 0,
                         check.cp = FALSE, frequency = 7)
er <- errs$errors$data
colnames(er)<- errs$errors$labels
saveRDS(er,'C:\\Users\\jerem\\Downloads\\elec_errs.rds')
tmp <- electricity$data
colnames(tmp) <- electricity$labels
saveRDS(tmp,'C:\\Users\\jerem\\Downloads\\electricity.rds')
