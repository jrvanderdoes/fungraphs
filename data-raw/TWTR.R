## code to prepare `TWTR` dataset goes here

tmp <- read.csv2(paste0(getwd(),'/data-raw/',"TWTR.txt"),sep=',')
tmp <- tmp[,c("X.DATE.","X.TIME.","X.CLOSE.")]
colnames(tmp) <- c('date','time','close')
tmp$time <- as.character(tmp$time)
tmp$time <- ifelse(nchar(tmp$time)<=5,
                   paste0("0",tmp$time),
                   tmp$time)
tmp[["dateTime"]] <- as.Date(paste(tmp$date,tmp$time),"%Y%m%d %H%M%S")

TWTR <-
  tmp %>%
  select('date','time','close') %>%
  pivot_wider(id_cols = 'time', id_expand = TRUE,
              names_from = 'date', values_from = 'close') %>%
  filter(!dplyr::row_number() %in% c(391)) %>%
  mutate(across(`20190102`:`20201231`, as.numeric)) %>%
  select(-c(time)) %>%
  as.data.frame() %>%
  linear_imputatation()

usethis::use_data(TWTR, overwrite = TRUE)
