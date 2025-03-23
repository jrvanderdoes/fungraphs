library(tidyverse)
library(devtools);load_all()

path <- 'C://Users//jerem//OneDrive//Documents//School//Waterloo//Research//RPackages//fungraphs//paper_trials_known'

peds <- data.frame()

for(i in c(2019,2020,2021)){
  #for(i in c(2019,2020)){
  #for(i in c(2020,2021)){

  Jan <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         'January_',i,".csv"), na.strings="na")
  Feb <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "February_",i,".csv"), na.strings="na")
  Mar <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "March_",i,".csv"), na.strings="na")
  Apr <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "April_",i,".csv"), na.strings="na")
  May <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "May_",i,".csv"), na.strings="na")
  Jun <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "June_",i,".csv"), na.strings="na")
  Jul <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "July_",i,".csv"), na.strings="na")
  Aug <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "August_",i,".csv"), na.strings="na")
  Sep <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "September_",i,".csv"), na.strings="na")
  Oct <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "October_",i,".csv"), na.strings="na")
  Nov <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "November_",i,".csv"), na.strings="na")
  Dec <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "December_",i,".csv"), na.strings="na")
  if(i==2020){
    Nov$Date <- paste0(Nov$Date,'20')
    Dec$Date <- paste0(Dec$Date,'20')
  }
  if(i==2021){
    Jul$Date <- stringr::str_replace_all(
      as.character(as.Date(Jul$Date,format='%d/%m/%y')),'-','/')
    Jul$Date <- paste0(substr(Jul$Date,9,10),'/07/2021')
  }

  peds <- plyr::rbind.fill(peds,Jan,Feb,Mar,Apr,May,Jun,
                           Jul,Aug,Sep,Oct,Nov,Dec)
}
dropCols <- c()
for(i in 3:ncol(peds)){
  peds[,i] <- suppressWarnings(as.numeric(peds[,i]))
  peds[,i]<-ifelse(peds[,i]==-1,NA,peds[,i])
  if(sum(is.na(peds[,i]))>500)
    dropCols <- c(dropCols,i)
}
peds <- peds[,-dropCols]

# peds$total <- 0
# for(i in 1:nrow(peds)){
#   peds$total[i] <- sum(peds[i,-c(1:2,ncol(peds))],na.rm = T)
# }
peds$Date <- as.Date(peds$Date,'%d/%m/%Y')
peds$YearMonthDay <- paste0(year(peds$Date),'_',month(peds$Date),'_',day(peds$Date))
peds$MonthDayHour <- paste0(month(peds$Date),'_',day(peds$Date),'_',peds$Hour)
#peds <- peds[,c(ncol(peds),2,ncol(peds)-1)]

peds1 <- peds[(wday(peds$Date) %in% c(1,7)),]
peds_daily <- peds1 %>% pivot_wider(id_cols = Hour,
                                    names_from = YearMonthDay,
                                    values_from = Elizabeth.St.Lonsdale.St..South.)
peds_daily <- as.data.frame(peds_daily[,-1])
peds_daily <- linear_imputatation(peds_daily,use.prev.curve = T)


plot_dates <- function(data, curve_points = 1:nrow(data),
                       plot_title = NULL,
                       val_axis_title = '',
                       res_axis_title = '',
                       FD_axis_title = '',
                       FDReps = as.Date(colnames(data),'%Y_%m_%d'),
                       eye = list(x = -0.5, y = -1.5, z = 0.5),
                       aspectratio = list(x=1, y=0.75, z=0.5),
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
              showticklabels = FALSE
            ),
            xaxis = list(
              title = FD_axis_title,
              titlefont = list(size = 20),
              tickfont = list(size = 14),
              showticklabels = TRUE
            ),
            zaxis = list(
              title = val_axis_title,
              titlefont = list(size = 20),
              tickfont = list(size = 14),
              showticklabels = FALSE
            )
          )
        )
      ),
      plotly::layout(title = plot_title, scene = scene)
    ),
    plotly::layout(showlegend = FALSE)
  )

}

tmp <- plot_dates(peds_daily)
plotly::save_image(tmp,
        file=paste0(path,"/figures/pedestrian_counts.png"))

################################



graph_change <- function(data, treeType='MST', kTrees=15,
                         error='L2', dataUse='Orig', pvalues=FALSE){

  data_dist <- calculateDistanceMatrix(
    data=data, silent = T,
    errType=error, dataUse = dataUse)

  if(treeType=='MST'){
    tree <- mstree(as.dist(data_dist), ngmax=kTrees)
  }else{
    stop('Sorry no other trees')
  }

  gc <- tryCatch({
    sink(nullfile())
    tmp <- gseg1(ncol(data), tree, pval.perm = T, B=100, pval.appr = F)
    sink()

    if(pvalues) cat(tmp$pval.perm$max.type$pval,'\n')

    ifelse(tmp$pval.perm$max.type$pval<0.05, tmp$scanZ$max.type$tauhat, NA)
  },error= function(x){
    sink()
    NA
  })

  gc
}

bs_func <- function(data,method, addAmt=0, silent=F){

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
  cps_ver <- na.omit(c(0,cps,ncol(data)))

  for(i in 2:(length(cps_ver)-1)){
    gc <- graph_change( data[, (cps_ver[i-1]+1):cps_ver[i+1]],pvalues=TRUE )


    cps_new <- c(cps_new,cps_ver[i-1]+gc)

  }

  cps_new[!is.na(cps_new)]
}


set.seed(12345)
res <- bs_func(peds_daily,graph_change)
res1 <- ver_func(peds_daily,res)
colnames(peds_daily)[res1]

plot_cps <- function(data, CPs, curve_points = 1:nrow(data),
                     plot_title = NULL,
                     val_axis_title = '',
                     res_axis_title = '',
                     FD_axis_title = '',
                     FDReps = as.Date(colnames(data),'%Y_%m_%d'),
                     eye = list(x = -0.5, y = -1.5, z = 0.5),
                     aspectratio = list(x=1, y=0.75, z=0.5),
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

  plotData$FDRep1 <- as.Date(colnames(data)[plotData$FDRep],'%Y_%m_%d')

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
              showticklabels = FALSE
            ),
            xaxis = list(
              title = FD_axis_title,
              titlefont = list(size = 20),
              tickfont = list(size = 14),
              showticklabels = TRUE
            ),
            zaxis = list(
              title = val_axis_title,
              titlefont = list(size = 20),
              tickfont = list(size = 14),
              showticklabels = FALSE
            )
          )
        )
      ),
      plotly::layout(title = plot_title, scene = scene)
    ),
    plotly::layout(showlegend = FALSE)
  )
}


tmp <- plot_cps(data = peds_daily,CPs = res1)
plotly::save_image(tmp,
                   file=paste0(path,"/figures/pedestrian_counts_segments.png"))
