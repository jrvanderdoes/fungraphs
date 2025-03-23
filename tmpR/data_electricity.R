## https://www.jepx.jp/en/electricpower/market-data/spot/

tmp0 <- read.csv("C:\\Users\\jerem\\Downloads\\spot_summary_2020.csv",
                skip = 1,header = F)
tmp <- read.csv("C:\\Users\\jerem\\Downloads\\spot_summary_2021.csv",
                skip = 1,header = F)
tmp1 <- read.csv("C:\\Users\\jerem\\Downloads\\spot_summary_2022.csv",
                skip = 1,header = F)
tmp2 <- read.csv("C:\\Users\\jerem\\Downloads\\spot_summary_2023.csv",
                skip = 1,header = F)

j_elec <- rbind(tmp0,tmp,tmp1,tmp2)
j_elec <- j_elec[,c(1,2,6)]
colnames(j_elec) <- c('date','halfhr','price')
j_elec <- j_elec[j_elec$date >=as.Date('2023-01-01') &
                   j_elec$date <= as.Date('2023-12-31'),]

j_elec_plot <- j_elec %>%
  pivot_wider(names_from = date,id_cols = halfhr,values_from = price)
j_elec_plot <- j_elec_plot[,-1]
j_elec_plot <- as.data.frame(j_elec_plot)

data_dist <- calculateDistanceMatrix(
  j_elec_plot, silent = T, errType='L2', dataUse = 'Orig')

mm_bs1 <- bs_func(data_dist,graph_change)
mm_bs1

plot_fd(j_elec_plot,mm_bs1)
