This document analyzes Twitter Stock Data.
```{readData r}
data_path <- 'C:/Users/jerem/OneDrive/Documents/School/Waterloo/Research/GraphCP/Data/Examples/'

tmp <- read.csv2(paste0(data_path,"TWTR.txt"),sep=',')
tmp <- tmp[,c("X.DATE.","X.TIME.","X.CLOSE.")]
colnames(tmp) <- c('date','time','close')
tmp$time <- as.character(tmp$time)
tmp$time <- ifelse(nchar(tmp$time)<=5,
                   paste0("0",tmp$time),
                   tmp$time)
tmp[["dateTime"]] <- as.Date(paste(tmp$date,tmp$time),"%Y%m%d %H%M%S")

twtr_fd <-
  tmp %>%
  select('date','time','close') %>%
  pivot_wider(id_cols = 'time', id_expand = TRUE,
              names_from = 'date', values_from = 'close') %>%
  filter(!dplyr::row_number() %in% c(391)) %>%
  mutate(across(`20190102`:`20201231`, as.numeric)) %>%
  select(-c(time)) %>%
  as.data.frame() %>%
  linear_imputatation()
twtr_fd_cidr <- compute_cidr(twtr_fd)
```
