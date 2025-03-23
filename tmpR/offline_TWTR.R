---
  title: "twitter"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{twitter}
%\VignetteEngine{knitr::rmarkdown}
%\VignetteEncoding{UTF-8}
---

  ```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(fungraphs)
```


```{r}
compute_cidr <- function(dat){
  dat_cidr <- dat
  for(i in 1:nrow(dat_cidr)){
    dat_cidr[i,] <- 100*(log(dat[i,]) - log(dat[1,]))
  }
  dat_cidr
}

data_path <- 'C:/Users/jerem/OneDrive/Documents/School/Waterloo/Research/GraphCP/Data/Examples/'
```

This document analyzes Twitter Stock Data.
```{readData r}
TWTR=readRDS('C:\\Users\\j53vande\\Downloads\\TWTR.rds')
tmp <- TWTR

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

########################################
graph_change <- function(data, treeType='MST', kTrees=15,
                         error='L2', dataUse='Orig'){

  data_dist <- calculateDistanceMatrix(
    data=data, silent = T,
    errType=error, dataUse = dataUse)

  if(treeType=='MST'){
    tree <- mstree(as.dist(data_dist), ngmax=kTrees)
  }else{
    stop('Sorry no other trees')
  }

  gc <- tryCatch({
    n = ncol(data)
    tmp <- gseg1(n, tree,
                 n0=max(0.05*n,20), n1=min(0.95*n,n-20),
                 pval.perm = T, B=100, pval.appr = F)
    ifelse(tmp$pval.perm$ori$pval<0.05,tmp$scanZ$ori$tauhat,NA)
  },error= function(x){NA})

  gc
}

bs_func <- function(data, method, addAmt=0, silent=F){
  if(ncol(data)<=30) return()

  potential_cp <- method(data)

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
            silent=silent),
    potential_cp + addAmt,
    bs_func(data=data[,(potential_cp+1):ncol(data)],
            method=method,
            addAmt = addAmt+potential_cp,
            silent=silent)
  ))
}

ver_func <- function(data, cps){
  cps_new <- c()
  cps_ver <- c(0,cps,ncol(data))

  for(i in 2:length(cps_ver)){
    dat <- data[,(cps_ver[i-1]+1):cps_ver[i]]


    data_dist <- calculateDistanceMatrix(
      data=dat, silent = T,
      errType='L2', dataUse = 'Orig')
    tree <- mstree(as.dist(data_dist), ngmax=5)

    gc <- tryCatch({
      n = ncol(dat)
      tmp <- gseg1(n, tree,
                   n0=max(0.05*n,20), n1=min(0.95*n,n-20),
                   pval.perm = T, B=100, pval.appr = F)
      ifelse(tmp$pval.perm$ori$pval<0.05,tmp$scanZ$ori$tauhat,NA)
    })
    cps_new <- c(cps,gc)

  }

  cps_new[!is.na(cps_new)]
}

set.seed(12345)
gc <- bs_func(data=twtr_fd_cidr, method=graph_change)
gc_ver <- ver_func(twtr_fd_cidr, gc)

# plot_fd(TWTR_final,gc_ver)

```

```{computeCPs r}

# mm_bs1 <- bs_func(twtr_fd_cidr,mean_change)
# cv_bs1 <- bs_func(twtr_fd_cidr,cov_change)
# gc_bs1 <- bs_func(twtr_fd_cidr,graph_change)
```

```{r plot_data}
# plot_fd(twtr_fd_cidr, CPs = mm_bs1)
# plot_fd(twtr_fd_cidr, CPs = cv_bs1)
# plot_fd(twtr_fd_cidr, CPs = gc_bs1)
```

```{r}

plot_dates <- function(data, curve_points = 1:nrow(data),
                       plot_title = NULL,
                       val_axis_title = '',
                       res_axis_title = '',
                       FD_axis_title = '',
                       FDReps = as.Date(colnames(data),'%Y%m%d'),
                       eye = list(x = -0.35, y = -1.5, z = 0.5),
                       aspectratio = list(x=0.75, y=0.75, z=0.5),
                       showticklabels = TRUE) {
  number <- length(data[1, ])
  valRange <- c(min(data), max(data))

  plotData <- data.frame(
    "resolution" = curve_points,
    "FDRep" = FDReps[1],
    "Value" = data[, 1]
  )

  for (i in 2:number) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = FDReps[i],
        "Value" = as.factor(data[, i])
      )
    )
  }

  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  tmpColors <- RColorBrewer::brewer.pal(11, "Spectral")
  tmpColors[6] <- "yellow"

  magrittr::`%>%`(
    magrittr::`%>%`(
      magrittr::`%>%`(
        plotly::plot_ly(plotData,
                        x = ~FDRep, y = ~resolution, z = ~Value,
                        type = "scatter3d", mode = "lines",
                        color = ~ as.factor(FDRep),
                        colors = tmpColors
        ),
        plotly::layout(
          scene = list(
            yaxis = list(
              title = res_axis_title,
              titlefont = list(size = 20),
              tickfont = list(size = 14),
              showticklabels = showticklabels
            ),
            xaxis = list(
              title = FD_axis_title,
              titlefont = list(size = 20),
              tickfont = list(size = 14),
              showticklabels = showticklabels
            ),
            zaxis = list(
              title = val_axis_title,
              titlefont = list(size = 20),
              tickfont = list(size = 14),
              showticklabels = showticklabels
            )
          )
        )
      ),
      plotly::layout(title = plot_title, scene = scene)
    ),
    plotly::layout(showlegend = FALSE)
  )

  # plotly::plot_ly(plotData,
  #                x = ~as.factor(FDRep), y = ~resolution, z = ~Value,
  #                type = 'scatter3d', mode = 'lines',
  #                color = ~as.factor(FDRep),
  #                colors = tmpColors) %>%
  #  #colors = c("red", "yellow", "blue")) %>%
  #  #colors='Spectral') %>%
  #  plotly::layout(
  #    scene = list(
  #      yaxis = list(title = "resolution"),
  #      xaxis = list(title = "FD Sims"),
  #      zaxis = list(title = "Value")
  #    )) %>%
  #  plotly::layout(title = plot_title, scene = scene) %>%
  #  plotly::layout(showlegend = FALSE)
  ## line = list(width = 4, color = ~as.factor(FDRep),
  ##            colorscale = list(c(0,'#BA52ED'), c(1,'#FCB040'))))
}

tmp <- plot_dates(twtr_fd_cidr)
plotly::save_image(tmp,
                   "C:/Users/j53vande/Downloads/twtr_cidr.png")
```

```{r}

plot_cps <- function(data, CPs, curve_points = 1:nrow(data),
                     plot_title = NULL,
                     val_axis_title = '',
                     res_axis_title = '',
                     FD_axis_title = '',
                     FDReps = as.Date(colnames(data),'%Y%m%d'),
                     eye = list(x = -0.35, y = -1.5, z = 0.5),
                     aspectratio = list(x=0.75, y=0.75, z=0.5),
                     showticklabels = TRUE) {
  number <- length(data[1, ])
  valRange <- c(min(data), max(data))

  plotData <- data.frame(
    "resolution" = NA,
    "FDRep" = NA,
    "Color" = NA,
    "Value" = NA
  )[-1, ]

  # Color Group to first CP
  for (j in 1:min(CPs)) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = j,
        "Color" = 1,
        "Value" = data[, j]
      )
    )
  }
  # Color Group from last CP
  for (j in (max(CPs) + 1):number) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = j,
        "Color" = length(CPs) + 1,
        "Value" = data[, j]
      )
    )
  }

  # Color Additional Groups
  if (length(CPs) > 1) {
    for (i in 2:length(CPs)) {
      for (j in (CPs[i - 1] + 1):CPs[i]) {
        plotData <- rbind(
          plotData,
          data.frame(
            "resolution" = curve_points,
            "FDRep" = j,
            "Color" = i,
            "Value" = data[, j]
          )
        )
      }
    }
  }

  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  # Get Colors
  tmpColors <- RColorBrewer::brewer.pal(min(9, max(3, length(CPs) + 1)), "Set1")
  if (length(CPs) > 9) {
    tmpColors <- rep(tmpColors, ceiling(c(length(CPs) + 1) / 9))[1:(length(CPs) + 1)]
  }

  tmpLines <- rep(c('solid','dashed'),length.out=length(CPs)+1)
  names(tmpLines) <- unique(levels(as.factor(as.factor(plotData$Color))))

  plotData$FDRep1 <- as.Date(colnames(data)[plotData$FDRep],'%Y%m%d')

  magrittr::`%>%`(
    magrittr::`%>%`(
      magrittr::`%>%`(
        plotly::plot_ly(plotData,
                        x = ~ as.factor(FDRep1), y = ~resolution, z = ~Value,
                        type = "scatter3d", mode = "lines",
                        split = ~ as.factor(FDRep),
                        color = ~ as.factor(Color),
                        linetype = ~ as.factor(Color),
                        colors = tmpColors,
                        linetypes = tmpLines
        ),
        plotly::layout(
          scene = list(
            yaxis = list(
              title = res_axis_title,
              titlefont = list(size = 20),
              tickfont = list(size = 14),
              showticklabels = showticklabels
            ),
            xaxis = list(
              title = FD_axis_title,
              titlefont = list(size = 20),
              tickfont = list(size = 14),
              showticklabels = showticklabels
            ),
            zaxis = list(
              title = val_axis_title,
              titlefont = list(size = 20),
              tickfont = list(size = 14),
              showticklabels = showticklabels
            )
          )
        )
      ),
      plotly::layout(title = plot_title, scene = scene)
    ),
    plotly::layout(showlegend = FALSE)
  )
}


tmp <- plot_cps(data = twtr_fd_cidr,CPs = gc_ver)
plotly::save_image(tmp,
                   "C:/Users/j53vande/Downloads/twtr_cidr_segmented.png")

```
