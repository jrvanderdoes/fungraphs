############################

linear_imputatation <- function(data, evalPts = 1:nrow(data),
                                use.prev.curve = FALSE) {
  # Look at all FDs
  for (i in 1:ncol(data)) {
    # If missing value
    if (sum(is.na(data[, i])) > 0) {
      # Set values if none in column and want to use others
      #   Assume filled last and will try to use future, but NA is possible
      if (sum(is.na(data[, i])) == nrow(data) &&
          use.prev.curve) {
        if (i > 1) data[1, i] <- data[nrow(data), i - 1]
        if (i < ncol(data)) data[nrow(data), i] <- data[1, i + 1]
      }

      # Set starting
      prevInfo <- c(data[1, i], 1)
      nextInfo <- c(data[nrow(data), i], nrow(data))

      # Check first
      if (is.na(prevInfo[1])) {
        for (j in 2:nrow(data)) {
          if (!is.na(data[j, i])) {
            prevInfo <- c(data[j, i], j)
            for (k in 1:j) {
              data[k, i] <- prevInfo[1]
            }
            break
          }
        }
      }
      # Check last
      if (is.na(nextInfo[1])) {
        for (j in (nrow(data) - 1):1) {
          if (!is.na(data[j, i])) {
            nextInfo <- c(data[j, i], j)
            for (k in nrow(data):j) {
              data[k, i] <- nextInfo[1]
            }
            break
          }
        }
      }

      # Fix Middle
      st <- prevInfo[2]
      en <- nextInfo[2]
      fill <- FALSE

      ## FIX HERE
      # if(!is.numeric(st) && !is.numeric(en)){
      #   warning(paste0('Error: Column ',j,' has no data in it.',
      #                  ' It is entirely dropped from data'))
      # } else if()

      for (j in (st + 1):en) {
        # Check if value needs to be interpolated
        if (is.na(data[j, i]) && !fill) {
          prevInfo <- c(data[j - 1, i], j - 1)
          fill <- TRUE
        }
        # Interpolate values if possible
        if (!is.na(data[j, i]) && fill) {
          nextInfo <- c(data[j, i], j)
          fill <- FALSE
          for (k in (prevInfo[2] + 1):(nextInfo[2] - 1)) {
            x1 <- evalPts[prevInfo[2]]
            x2 <- evalPts[k]
            x3 <- evalPts[nextInfo[2]]
            y1 <- prevInfo[1]
            y3 <- nextInfo[1]


            data[k, i] <- (x2 - x1) * (y3 - y1) / (x3 - x1) + y1
          }
        }
      }
    }
  }

  data
}


############################
#' Plot functional data
#'
#' \code{plot_fd} plots functional data either in an fd object or evaluated at certain points
#'
#' @param data A data.frame of evaluated fd objects (the columns being lines and
#'     the rows the evaluated points)
#' @param CPs (Optional) Vectors of numeric values indicating the location of
#'     change points. This will color each section differently. Default vallue
#'     is NULL.
#' @param curve_points (Optional) An vector containing the points at which the
#'     data was evaluated. Default is 1:nrow(data)
#' @param plot_title (Optional) String to title the plot. Default value is NULL.
#' @param val_axis_title (Optional) String to title the axis for values of
#'     observation. Default value is 'Value'.
#' @param res_axis_title (Optional) String to title the axis for the resolution
#'     range axis (observations in an FD object). Default value is 'resolution'.
#' @param FD_axis_title (Optional) String to title the axis for FD observations.
#'     Default value is 'Observation'.
#' @param FDReps (Optional) Vector of values for FD observation names. Default
#'     value is 1:ncol(data).
#' @param eye (Optional) List with certain parameters to determine the view of
#'     the resulting image. The default value is list(x = -1.5, y = -1.5, z = 1.5).
#' @param aspectratio (Optional) List with certain parameters to determine the
#'     image size / viewing parameters. The default value is
#'     list(x=1, y=1, z=1), which is used in the non-interactive case. If
#'     interactive is TRUE, then the default is c(2.5,.75,1) and all lists are
#'     converted to this value so enter a vector if you wish to modify it for
#'     the non-interactive case.
#' @param showticklabels (Optional) Boolean to indicate if the tick labels should
#'     be given in the image. Default value is TRUE.
#' @param interactive (Optional) Boolean value to indicate if the plot should
#'  be interactive or not. Recommended to be TRUE when the number of observations
#'  is high. Default is FALSE
#'
#' @return A plot for the data. It is an interactive plotly plot if
#'  interactive is FALSE (default) and a lattice plot if TRUE.
#' @export
#'
#' @examples
#' plot_fd(data = electricity[, 1:10])
#' plot_fd(data = electricity[, 1:50], CPs = c(25))
#' plot_fd(
#'   data = electricity, CPs = c(50, 150, 220, 300),
#'   interactive = FALSE, showticklabels = FALSE
#' )
plot_fd <- function(data, CPs = NULL, curve_points = 1:nrow(data), plot_title = NULL,
                    val_axis_title = "Value", res_axis_title = "resolution",
                    FD_axis_title = "Observations", FDReps = 1:ncol(data),
                    eye = list(x = -1.5, y = -1.5, z = 1.5),
                    aspectratio = list(x = 1, y = 1, z = 1),
                    showticklabels = TRUE, interactive = TRUE) {
  if (!interactive) {
    fdPlot <- .plot_evalfd_highdim(
      data = data, curve_points = curve_points,
      CPs = CPs, plot_title = plot_title,
      val_axis_title = val_axis_title,
      res_axis_title = res_axis_title,
      FD_axis_title = FD_axis_title,
      aspectratio = aspectratio,
      showticklabels = showticklabels
    )
  } else if (!is.null(CPs) && length(stats::na.omit(CPs)) > 0) {
    fdPlot <- .plot_evalfd_3dlines_cps(
      data = data, curve_points = curve_points,
      CPs = CPs[order(CPs)], plot_title = plot_title,
      val_axis_title = val_axis_title,
      res_axis_title = res_axis_title,
      FD_axis_title = FD_axis_title,
      eye = eye, aspectratio = aspectratio,
      showticklabels = showticklabels
    ) # ,FDReps)
  } else {
    fdPlot <- .plot_evalfd_3dlines(
      data = data, curve_points = curve_points,
      plot_title = plot_title,
      val_axis_title = val_axis_title,
      res_axis_title = res_axis_title,
      FD_axis_title = FD_axis_title, FDReps = FDReps,
      eye = eye, aspectratio = aspectratio,
      showticklabels = showticklabels
    )
  }

  fdPlot
}


#' Plot Functional Data Without Change Points
#'
#' This (internal) function to plot the function data, with no coloring based on
#'     change points.
#'
#' @inheritParams plot_fd
#'
#' @return A plotly plot
#'
#' @noRd
.plot_evalfd_3dlines <- function(data, curve_points, plot_title = NULL,
                                 val_axis_title = "Value",
                                 res_axis_title = "Resolution",
                                 FD_axis_title = "Observation",
                                 FDReps = 1:ncol(data),
                                 eye = list(x = -1.5, y = -1.5, z = 1.5),
                                 aspectratio = list(x = 1, y = 1, z = 1),
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
              showticklabels = showticklabels
            ),
            xaxis = list(
              title = FD_axis_title,
              showticklabels = showticklabels
            ),
            zaxis = list(
              title = val_axis_title,
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


#' Plot Functional Data With Change Points
#'
#' This (internal) function to plot the function data, with coloring based on
#'     change points.
#'
#' @inheritParams plot_fd
#' @param CPs Vectors of numeric values indicating the location of change points.
#'     This will color each section differently.
#'
#' @return A plotly plot
#'
#' @noRd
.plot_evalfd_3dlines_cps <- function(data, curve_points, CPs,
                                     plot_title = NULL,
                                     val_axis_title = "Value",
                                     res_axis_title = "resolution",
                                     FD_axis_title = "FD Sims",
                                     eye = list(x = -1.5, y = -1.5, z = 1.5),
                                     aspectratio = list(x = 1, y = 1, z = 1),
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

  magrittr::`%>%`(
    magrittr::`%>%`(
      magrittr::`%>%`(
        plotly::plot_ly(plotData,
                        x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                        type = "scatter3d", mode = "lines",
                        split = ~ as.factor(FDRep),
                        color = ~ as.factor(Color),
                        colors = tmpColors
        ),
        plotly::layout(
          scene = list(
            yaxis = list(
              title = res_axis_title,
              showticklabels = showticklabels
            ),
            xaxis = list(
              title = FD_axis_title,
              showticklabels = showticklabels
            ),
            zaxis = list(
              title = val_axis_title,
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


#' Plot With Surface for Speed
#'
#' @inheritParams plot_fd
#' @param aspectratio (Optional) List with certain parameters to determine the
#'     image size parameters. The default value is c(2.5,.75,1). Also any lists
#'     are converted to this (as it is likely from default call `in plot_fd()`).
#'
#' @return A lattice plot
#'
#' @noRd
.plot_evalfd_highdim <- function(data, curve_points, CPs = NULL,
                                 plot_title = NULL, val_axis_title = NULL,
                                 res_axis_title = NULL,
                                 FD_axis_title = NULL,
                                 aspectratio = c(2.5, .75, 1),
                                 showticklabels = FALSE) {
  if (!requireNamespace("lattice", quietly = TRUE)) {
    stop(
      "Package \"lattice\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # This resets ratio, since default was likely used
  if (is.list(aspectratio)) aspectratio <- c(2.5, 0.75, 1)

  number <- length(data[1, ])
  valRange <- c(
    floor(min(data)),
    ceiling(max(data))
  )

  plotData <- data.frame(
    "resolution" = curve_points,
    "FDRep" = 1,
    "Value" = data[, 1]
  )

  for (i in 2:number) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = i,
        "Value" = data[, i]
      )
    )
  }

  ## Setup up Colors
  #   Rainbow for no CPs, colored for CPs
  plotData[["color"]] <- rep(1:ncol(data), each = nrow(data))
  if (!is.null(CPs)) {
    tmp_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(CPs) + 1)), "Set1")
    if (length(CPs) > 9) {
      tmp_colors <- rep(tmp_colors, ceiling(c(length(CPs) + 1) / 9))[1:(length(CPs) + 1)]
    }

    CPs <- unique(c(1, CPs, ncol(data)))
    colors_plot <- rep(tmp_colors[1], ncol(data))
    for (i in 2:(length(CPs) - 1)) {
      colors_plot[CPs[i]:CPs[i + 1]] <- tmp_colors[i]
    }
  } else {
    colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
    colors_plot[6] <- "yellow"
    colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(data))
    plotData[["color"]] <- rep(colors_plot,
                               each = nrow(data)
    )
  }

  ## Set up Tick Labels
  z_range <- round(range(plotData$Value), -2)
  if (showticklabels) {
    scale_info <- list(
      col = "black", arrows = FALSE, cex = 0.75, cex.title = 1.5,
      x = list(
        at = seq(min(curve_points),
                 max(curve_points),
                 length.out = 5
        ),
        # labels=.specify_decimal(rev(seq(max(curve_points),
        #                                 min(curve_points),
        #                                 length.out=4)),2)),
        labels = NULL
      ),
      y = list(
        at = -seq(1, number, length.out = 8),
        labels = .specify_decimal(
          seq(1, number, length.out = 8), 0
        )
      ),
      z = list(
        at = seq(z_range[1], z_range[2],
                 length.out = 5
        ),
        labels = .specify_decimal(
          seq(z_range[1], z_range[2], length.out = 5), 0
        )
      )
    )
  } else {
    scale_info <- list(
      col = "black", arrows = FALSE, cex = 0.75, cex.title = 1.5,
      x = list(
        at = seq(min(curve_points),
                 max(curve_points),
                 length.out = 5
        ),
        labels = NULL
      ),
      y = list(
        at = -seq(1, number, length.out = 8),
        labels = NULL
      ),
      z = list(
        at = seq(z_range[1], z_range[2],
                 length.out = 5
        ),
        labels = NULL
      )
    )
  }

  color <- NULL # Check note removal
  # lattice:::wireframe
  lattice::cloud(
    x = Value ~ (resolution) * (-FDRep),
    data = plotData,
    type = "l", groups = color,
    par.box = c(col = "transparent"),
    par.settings =
      list(
        axis.line = list(col = "transparent"),
        superpose.line = list(col = colors_plot)
      ),
    # screen=list(z = 90, x = -75,y=-45),
    # trellis.par.set(list(axis.text=list(cex=2)),
    #                "axis.line",list(col=NA)),
    zlim = valRange,
    aspect = aspectratio,
    drape = TRUE, colorkey = FALSE,
    scales = scale_info,
    xlab = list(res_axis_title, rot = 30),
    ylab = list(paste0("\n", FD_axis_title), rot = -30),
    zlab = list(val_axis_title, rot = 90, just = 0.75),
    main = plot_title
  )
}


############################
library(tidyverse)
library(devtools);load_all()

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

#plot_fd(peds_daily)

##################################################

graph_change <- function(data, treeType='MST', kTrees=5,
                         error='L2', dataUse='Orig'){

  n = ncol(data)
  n0 = round(max(0.05*n,14))
  n1 = round(min(0.95*n,n-14))
  if(n0>=n1) return(NA)

  data_dist <- calculateDistanceMatrix(
    data=data, silent = T,
    errType=error, dataUse = dataUse)

  if(treeType=='MST'){
    tree <- mstree(as.dist(data_dist), ngmax=kTrees)
  }else{
    stop('Sorry no other trees')
  }

  gc <- tryCatch({
    tmp <- gseg1(n, tree,
                 n0=n0, n1=n1,
                 pval.perm = T, B=300, pval.appr = F)
    ifelse(tmp$pval.perm$generalized$pval<0.05,
           tmp$scanZ$generalized$tauhat,
           NA)
  },
  error = function(x){
    NA
  })

  gc
}

bs_func <- function(data, addAmt=0, silent=F){
  potential_cp <- graph_change(data)

  # No Change Point Detected
  if(is.na(potential_cp)) return()

  # Display progress
  if(!silent)
    cat(paste0('ChangePoint Detected (',1+addAmt,'-' ,addAmt+ncol(data),' at ',
               addAmt+potential_cp,'): Segment Data and Re-Search\n'))

  return(c(
    bs_func(data=data[,1:potential_cp],
            addAmt = addAmt,
            silent=silent),
    potential_cp + addAmt,
    bs_func(data=data[,(potential_cp+1):ncol(data)],
            addAmt = addAmt+potential_cp,
            silent=silent)
  ))
}

ver_func <- function(data, cps){
  cps_new <- c()
  cps_ver <- na.omit(c(0,cps,ncol(data)))

  for(i in 2:length(cps_ver)){
    gc <- graph_change( data[, (cps_ver[i-1]+1):cps_ver[i]] )


    cps_new <- c(cps_new,cps_ver[i-1]+gc)

  }

  cps_new[!is.na(cps_new)]
}

# set.seed(12345)
# gc <- bs_func(data=peds_daily)
# gc_ver <- ver_func(data=peds_daily, cps=gc)
# # gc_ver1 <- ver_func(data=peds_daily, cps=gc_ver)
# # gc_ver2 <- ver_func(data=peds_daily, cps=gc_ver1)
#
# plot_fd(peds_daily, gc_ver)


#######################################################


electricity_fd <- fda::Data2fd(1:24, as.matrix(peds_daily))
# sum(fda::eval.fd(1:24,electricity_fd)-electricity)
#   eval.fd: eval.basis(1:24,electricity_fd$basis) %*% electricity_fd$coefs

# Play with more components
nPCs <- 5
electricity_fpca <- fda::pca.fd(electricity_fd, nharm = nPCs)
electricity_fpca_comp <- electricity_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for (i in 1:nPCs) {
  ts_dat[[i]] <- ts(electricity_fpca_comp[, i], freq = 7)
  # comps[[i]] <- forecast::ets(ts_dat[[i]])
  # comps_resids[[i]] <- resid(comps[[i]])
  comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps_resids[[i]] <- resid(comps[[i]])
}

electricity_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
#   Want: 24 x 365
#     forecast: 365 x 3
#    coefs: 26  x 3
#      coefs %*% comps: 26 x 365
#    Eval: 24 x 26
#       eval %*% orig: 24 x 365
orig_coefs <- electricity_fpca$harmonics$coefs %*% t(electricity_fpca_forecast)
eval_fd_vals <- fda::eval.basis(1:24, electricity_fd$basis) %*% orig_coefs


set.seed(12345)
gc <- bs_func(data=eval_fd_vals)
# gc_ver <- ver_func(data=eval_fd_vals, cps=gc)
# gc; gc_ver;

plot_fd(peds_daily, gc)
