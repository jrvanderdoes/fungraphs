---
  title: "Weather"
output: html_document
---

  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
compute_cidr <- function(dat){
  dat_cidr <- dat
  for(i in 1:nrow(dat_cidr)){
    dat_cidr[i,] <- 100*(log(dat[i,]) - log(dat[1,]))
  }
  dat_cidr
}

data_path <- 'C:/Users/jerem/OneDrive/Documents/School/Waterloo/Research/GraphCP/Data/'
```

This document analyzes Twitter Stock Data.
```{readData r}

tmp <- read.csv2(paste0(data_path,"Basil_Clean.txt"),sep=',')
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

```{computeCPs r}

mm_bs1 <- bs_func(twtr_fd_cidr,mean_change)
cv_bs1 <- bs_func(twtr_fd_cidr,cov_change)
gc_bs1 <- bs_func(twtr_fd_cidr,graph_change)
```

```{r plot_data}
plot_fd(twtr_fd_cidr, CPs = mm_bs1)
plot_fd(twtr_fd_cidr, CPs = cv_bs1)
plot_fd(twtr_fd_cidr, CPs = gc_bs1)
```
