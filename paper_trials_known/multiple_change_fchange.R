
path <- 'C://Users//jerem//OneDrive//Documents//School//Waterloo//Research//RPackages//fungraphs//paper_trials_known'
data <- readRDS(file=paste0(path,"\\data\\mult_change.rds"))

# ce_pts <- c()
# for(i in 1:length(data)){
#   set.seed(12345678 * i)
#   cat(paste0(i,', '))
#
#   tmp <- change(data[[i]],method='characteristic', statistic='Tn', critical='welch',
#                 type='segmentation')
#   ce_pts <- c(ce_pts, tmp$location)
# }
# saveRDS(ce_pts, "C:\\Users\\jerem\\Downloads\\graph_updates\\mult_change_ce.rds")

ce_pts <- c()
for(i in 1:length(data)){
  set.seed(12345678 * i)
  cat(paste0(i,', '))

  tmp <- change(data[[i]],method='characteristic', statistic='Tn', critical='simulation',
                type='segmentation')
  ce_pts <- c(ce_pts, tmp$location)
}
saveRDS(ce_pts, paste0(path,"\\results\\mult_change_ce_s.rds"))
