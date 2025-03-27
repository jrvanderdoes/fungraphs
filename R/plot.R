#' Plot functional data
#'
#' Function to display functional data. See vastly improved visualization in the
#'  package fChange.
#'
#' @param data Numeric data.frame with evaled points on rows and fd objects in columns
#' @param curve_points Evaluation points
#' @param plot_title,val_axis_title,res_axis_title,FD_axis_title Titles for use
#'  in the plot
#' @param FDReps Point to label individual curves
#' @param eye Angle to view plot
#' @param aspectratio Relative size of each dimension
#' @param showticklabels Boolean if ticks should be shown on plot
#'
#' @return A plotly object.
#' @export
plot_fd <- function(data, curve_points = 1:nrow(data),
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
