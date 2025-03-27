#' FD Generation - KL Expansion
#'
#' Generation of FD using an autoregressive karhunen-loeve expansion. This can
#'  include change points in any combination of the following:
#'  \itemize{
#'    \item Mean
#'    \item Distribution
#'    \item Eigenvalue(s)
#'    \item Eigenvector(s)
#'  }
#'  In this sense, the function creates m 'groups' of discretely observed
#'  functions with similar properties. See updated
#'  version in the package fChange (\url{https://github.com/jrvanderdoes/fChange}).
#'
#' @param ns Vector of Numerics. Each value in N is the number of observations
#'  for a given group.
#' @param eigsList Vector of eigenvalues. Length 1 or m.
#' @param basesList A list of bases (eigenfunctions), length m.
#' @param meansList A vector of means, length 1 or N.
#' @param distsArray A vector of distributions, length 1 or m.
#' @param evals A vector of points indicating the points to evaluate the
#'     functions on.
#' @param kappasArray Numeric \[0,1\] indicating strength of VAR(1) process.
#' @param burnin A numeric value indicating the number of burnin trials.
#' @param silent A Boolean that toggles running output.
#' @param dof Numeric for degrees of freedom with t-distribution.
#' @param shape Numeric for shape of gamma distribution.
#' @param ... Additional parameters to send in. Unused.
#'
#' @return List with (1) data (N-by-m) and (2) previous errors.
#' @export
#'
#' @examples
#' data_KL <- gen_FD_KL_Expansion(ns = c(20,20),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(0,0),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
gen_FD_KL_Expansion <- function(ns,
                                eigsList,
                                basesList,
                                meansList,
                                distsArray,
                                evals,
                                kappasArray = c(0),
                                burnin = 100,
                                silent = F,
                                dof=NULL,
                                shape=NULL,
                                ...){
  # ns is a vector with length m for the number of data runs until next CP
  # - i.e. c(10,10,10) has 10 length TS then CP followed by 10 and another CP
  # eigsList is a list of vectors giving the eigenvalues for each distribution
  # - Length of list be 1 or m
  # basesList is a list of bases (eigf) for each distribution
  # - Length should be 1 or m
  # - Define basis on (0,1)
  # meansList is a list of means for each distribution
  # - Length should be 1 or m
  # distsArray is a vector of distributions to each run
  # - Length should be 1 or m
  # evals is the vector of points to evaluate
  # kappasArray is a vector of Kappas for VAR
  # - Length should be 1 or m

  ## Functions
  checkLength <- function(dataList, name, m){
    if(length(dataList) == 1){
      dataList <- list(rep(dataList[[1]], m))
    } else if(length(dataList) != m){
      stop(paste(name,'is length',length(dataList),'not 1 or',m,'\n'))
    }

    dataList
  }

  getPsiList <- function(D, m,eigsList,kappasArray){
    psi <- list()
    normsSD <- stats::rnorm(max(D), mean=0, sd=1)

    for(i in 1:m){
      groupSD <- normsSD[1:D[i]] * sqrt(eigsList[[i]])
      psi0 <- groupSD %*% t(groupSD)
      psi0 <- psi0 / sqrt(sum(psi0^2)) ## TODO:: Check this
      psi[[i]] <- kappasArray[i] * psi0
    }

    psi
  }

  ## Verification
  m <- length(ns)

  eigsList <- checkLength(eigsList, 'eigsList', m)
  basesList <- checkLength(basesList, 'basesList', m)
  meansList <- checkLength(meansList, 'meansList', m)
  distsArray <- unlist(checkLength(distsArray, 'distsArray', m))
  kappasArray <- unlist(checkLength(kappasArray, 'kappaArray', m))
  if(!is.null(dof)) dof <- unlist(checkLength(dof, 'dof', m))
  if(!is.null(shape)) shape <- unlist(checkLength(shape, 'shape', m))

  # Run Code
  data <- data.frame(matrix(NA, ncol = sum(ns), nrow = length(evals)))
  addIdx <- 0

  # Setup psi
  Ds <- 1:m
  for(i in 1:m){
    Ds[i] <- length(eigsList[[i]])
  }
  psi <- getPsiList(Ds, m,eigsList,kappasArray)

  # Burnin for VAR
  peps <- data.frame(matrix(0,ncol=length(evals),nrow=Ds[1]))
  for(j in 1:burnin){
    waste <- generateData_KL_Expansion(
      eigs = eigsList[[1]],
      basis = basesList[[1]],
      means = meansList[[1]],
      dist = distsArray[1],
      evals = evals,
      peps = peps,
      psi = psi[[1]],
      dof=dof[1],
      shape=shape[1],
      ...)

    peps <- waste[[2]]
  }

  for(i in 1:m){
    if(!silent)
      cat(paste0('Running setup ', i, '/',m,'\n'))

    for(j in 1:ns[i]){

      # If Num of Eigs increases or decreases (only at CPs)
      psiDim1 <- dim(psi[[i]])[1]
      pepDim1 <- dim(peps)[1]
      if(psiDim1 != pepDim1){
        if(psiDim1>pepDim1){
          # Bind row of 0s to the bottom if didn't have a value previously
          for(k in 1:(psiDim1-pepDim1)){
            peps <- rbind(peps,0)
          }
        }else if(psiDim1<pepDim1){
          # Drop Rows if not needed
          for(k in 1:(pepDim1-psiDim1)){
            peps <- peps[-nrow(peps),]
          }
        }
      }

      result <- generateData_KL_Expansion(
        eigs = eigsList[[i]],
        basis = basesList[[i]],
        means = meansList[[i]],
        dist = distsArray[i],
        evals = evals,
        peps = peps,
        psi = psi[[i]],
        dof=dof[i],
        shape=shape[i],
        ...)

      data[,addIdx + j] <- result[[1]]
      peps <- result[[2]]
    }
    addIdx <- addIdx + j
  }

  data
}


#' @noRd
#' @keywords internal
generateData_KL_Expansion <- function(eigs, basis, means, dist,
                                      evals, peps, psi, dof=NULL, shape=NULL,
                                      ...){

  ## Functions
  generateXi <- function(dist, sd,dof=3,skew=2,shape=1,rate=1){
    ## This function give centered distributions with eig^2 var

    xi <- 0

    if(dist == 'Normal'){

      xi <- stats::rnorm(1,mean=0, sd=sd)

    }else if(dist == 'sn'){
      if(!requireNamespace("fGarch", quietly = TRUE))
        stop('Need "fGarch" package for skewed normal errors.')
      xi <- fGarch::rsnorm(1,sd = sd,xi = skew)

    }else if(dist == 'Binomial'){

      if(sd==0)
        return(0)

      mean <- 10 * sd^2 # arbitrary, must exceed var
      p <- 1 - sd^2/mean
      size <- round(mean/p)

      xi <- stats::rbinom(n=1,size=size,p=p) - mean

    }else if(dist == 'Exponential'){

      xi <- stats::rexp(1,rate = 1/sd) - sd

    }else if(dist == 't'){
      xi <- stats::rt(1, dof)
      if(dof>2)
        xi <- xi * sqrt(sd^2 * (dof - 2) / dof)
    }else if (dist == "gamma") {
      xi <- ( stats::rgamma(1,shape = shape, rate = rate) - shape/rate ) /
        sqrt( shape * (1/rate)^2 ) * sd
    }else if(dist == 'Laplace'){
      if(!requireNamespace("jmuOutlier", quietly = TRUE))
        stop('Need "jmuOutlier" package for laplace errors.')
      xi <- jmuOutlier::rlaplace(1,0,sd)
    }else{
      stop(paste('Sorry, dist',dist,'not implemented yet'))
    }

    xi
  }

  ## Code

  # Setup
  n <- length(evals)
  D <- length(eigs)

  X <- rep(0, n)
  Zeta <- data.frame(matrix(NA,ncol=n,nrow=D)) # Matrix with col as time, row as dimension
  eps <- Zeta

  # Verify
  if(length(means)==1){
    means = rep(means,n)
  }else if(length(means)!=n){
    stop(paste('Length of means is',length(means),'not 1 or',n))
  }

  # Generate Noise
  xi <- rep(NA,D)
  for(j in 1:D){
    xi[j] <- generateXi(dist=dist, sd=sqrt(eigs[j]), dof=dof, shape=shape, ...)
  }

  # Generate Data
  for(t in 1:n){
    for(j in 1:D){
      Zeta[j,t] <- xi[j] * fda::eval.basis(evals[t], basis)[j]
    }

    eps[,t] <- Zeta[,t] + psi %*% peps[,t]
    X[t] <- means[t] + sum(eps[,t])
  }

  list(X, eps)
}



# #' Plot - FD Object Plot
# #'
# #' @noRd
# #' @keywords internal
# plot_fd_3dsurf <- function(fdobj, titleVal){
#
#   singleRange <- fdobj$basis$rangeval
#   singleRangeSeq <- seq(singleRange[1],singleRange[2],length.out=100)
#   number <- length(fdobj$fdnames$reps)
#
#   fd_eval <- eval.fd(singleRangeSeq,
#                      data_fd)
#   valRange <- c(min(fd_eval),max(fd_eval))
#
#   plotData <- data.frame('evalRange'=singleRangeSeq,
#                          'FDRep'=1,
#                          'Value'=fd_eval[,1])
#
#   for(i in 2:number){
#     plotData <- rbind(plotData,
#                       data.frame('evalRange'=singleRangeSeq,
#                                  'FDRep'=i,
#                                  'Value'=fd_eval[,i])
#     )
#   }
#
#   #trellis.par.set("axis.line",list(col='black'))
#   trellis.par.set("axis.line",list(col=NA))
#   wireframe(x=Value ~ (evalRange)*(-FDRep),
#             data= plotData,
#             #trellis.par.set(list(axis.text=list(cex=2)),
#             #                "axis.line",list(col=NA)),
#             zlim=valRange,
#             aspect=c(3,.75,1),
#             drape=TRUE,colorkey = FALSE,
#             scales = list(arrows=FALSE,cex= 0.75,
#                           cex.title=1.5,
#                           x = list(at=seq(singleRange[1],
#                                           singleRange[2],
#                                           0.2),
#                                    labels=rev(seq(singleRange[2],
#                                                   singleRange[1],
#                                                   -0.2))),
#                           y = list(at=-seq(1, number, 9),
#                                    labels=seq(1, number, 9))),
#             xlab="Eval Range",ylab="\nFD Reps",zlab="Value",
#             main=titleVal)
# }
#
#
# #' Title
# #'
# #' @noRd
# #' @keywords internal
# plot_fd_3dsurf_plotly <- function(fdobj, titleVal){
#
#   singleRange <- fdobj$basis$rangeval
#   singleRangeSeq <- seq(singleRange[1],singleRange[2],length.out=100)
#   number <- length(fdobj$fdnames$reps)
#
#   fd_eval <- eval.fd(singleRangeSeq,
#                      data_fd)
#   valRange <- c(min(fd_eval),max(fd_eval))
#
#   plotData <- data.frame('evalRange'=singleRangeSeq,
#                          'FDRep'=1,
#                          'Value'=fd_eval[,1])
#
#   for(i in 2:number){
#     plotData <- rbind(plotData,
#                      data.frame('evalRange'=singleRangeSeq,
#                                  'FDRep'=i,
#                                  'Value'=fd_eval[,i])
#     )
#   }
#
#   plot_ly() %>%
#     add_trace(data = plotData,
#               x=plotData$FDRep,
#               y=plotData$evalRange,
#               z=plotData$Value, type="mesh3d" )
#
#   tmpData <-
#    as.matrix(plotData %>%
#                 pivot_wider(., names_from ='evalRange',values_from= 'Value')
#     )[,-1]
#
#   scene = list(camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))
#
#   plot_ly(x = 1:number,
#          y = singleRangeSeq,
#           z = t(tmpData)) %>%
#     add_surface() %>%
#     layout(
#       scene = list(
#         yaxis = list(title = "EvalRange"),
#         xaxis = list(title = "FD Sims"),
#         zaxis = list(title = "Value")
#       )) %>%
#     layout(title = titleVal, scene = scene)
# }
#
#
# #' Title
# #'
# #' @param fdobj
# #' @param titleVal
# #'
# #' @return
# #' @export
# #'
# #' @examples
# plot_fd_3dlines_plotly <- function(fdobj, titleVal){
#
#   singleRange <- fdobj$basis$rangeval
#   singleRangeSeq <- seq(singleRange[1],singleRange[2],length.out=100)
#   number <- length(fdobj$fdnames$reps)
#
#   fd_eval <- eval.fd(singleRangeSeq,
#                      data_fd)
#   valRange <- c(min(fd_eval),max(fd_eval))
#
#   plotData <- data.frame('evalRange'=singleRangeSeq,
#                          'FDRep'=1,
#                          'Value'=fd_eval[,1])
#
#   for(i in 2:number){
#     plotData <- rbind(plotData,
#                       data.frame('evalRange'=singleRangeSeq,
#                                  'FDRep'=i,
#                                  'Value'=fd_eval[,i])
#     )
#   }
#
#   scene = list(camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))
#
#   tmpColors <- brewer.pal(11,"Spectral")
#   tmpColors[6] <- 'yellow'
#
#   plot_ly(plotData,
#           x = ~as.factor(FDRep), y = ~evalRange, z = ~Value,
#           type = 'scatter3d', mode = 'lines',
#           color = ~as.factor(FDRep),
#           colors = tmpColors) %>%
#     #colors = c("red", "yellow", "blue")) %>%
#     #colors='Spectral') %>%
#     layout(
#       scene = list(
#         yaxis = list(title = "EvalRange"),
#         xaxis = list(title = "FD Sims"),
#         zaxis = list(title = "Value")
#       )) %>%
#     layout(title = titleVal, scene = scene) %>%
#     layout(showlegend = FALSE)
#   #line = list(width = 4, color = ~as.factor(FDRep),
#   #            colorscale = list(c(0,'#BA52ED'), c(1,'#FCB040'))))
# }
#
#
# #' Title
# #'
# #' @noRd
# #' @keywords internal
# plot_evalfd_3dsurf <- function(fd_eval, singleRangeSeq,
#                                titleVal=NULL){
#
#   number <- length(fd_eval[1,])
#   valRange <- c(floor(min(fd_eval)),
#                 ceiling(max(fd_eval)))
#
#   plotData <- data.frame('evalRange'=singleRangeSeq,
#                          'FDRep'=1,
#                          'Value'=fd_eval[,1])
#
#   for(i in 2:number){
#     plotData <- rbind(plotData,
#                       data.frame('evalRange'=singleRangeSeq,
#                                  'FDRep'=i,
#                                  'Value'=fd_eval[,i])
#     )
#   }
#
#   #trellis.par.set("axis.line",list(col='black'))
#   trellis.par.set("axis.line",list(col=NA))
#   wireframe(x=Value ~ (evalRange)*(-FDRep),
#             data= plotData,
#             #trellis.par.set(list(axis.text=list(cex=2)),
#             #                "axis.line",list(col=NA)),
#             zlim=valRange,
#             aspect=c(3,.75,1),
#             drape=TRUE,colorkey = FALSE,
#             scales = list(arrows=FALSE,cex= 0.75,
#                           cex.title=1.5,
#                           x = list(at=seq(min(singleRangeSeq),
#                                           max(singleRangeSeq),
#                                           0.2),
#                                    labels=rev(seq(max(singleRangeSeq),
#                                                   min(singleRangeSeq),
#                                                   -0.2))),
#                           y = list(at=-seq(1, number, 9),
#                                    labels=seq(1, number, 9))),
#             xlab="Eval Range",ylab="\nFD Reps",zlab="Value",
#             main=titleVal)
# }
#
#
# #' Title
# #'
# #' @param fd_eval
# #' @param singleRangeSeq
# #' @param titleVal
# #'
# #' @return
# #' @export
# #'
# #' @examples
# plot_evalfd_3dsurf_plotly <- function(fd_eval, singleRangeSeq,
#                                       titleVal=NULL){
#
#   number <- length(fd_eval[1,])
#   valRange <- c(min(fd_eval),max(fd_eval))
#
#   plotData <- data.frame('evalRange'=singleRangeSeq,
#                          'FDRep'=1,
#                          'Value'=fd_eval[,1])
#
#   for(i in 2:number){
#     plotData <- rbind(plotData,
#                       data.frame('evalRange'=singleRangeSeq,
#                                  'FDRep'=i,
#                                  'Value'=fd_eval[,i])
#     )
#   }
#
#   plot_ly() %>%
#     add_trace(data = plotData,
#               x=plotData$FDRep,
#               y=plotData$evalRange,
#               z=plotData$Value, type="mesh3d" )
#
#   tmpData <-
#     as.matrix(plotData %>%
#                 pivot_wider(., names_from ='evalRange',values_from= 'Value')
#     )[,-1]
#
#   scene = list(camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))
#
#   plot_ly(x = 1:number,
#           y = singleRangeSeq,
#           z = t(tmpData)) %>%
#     add_surface() %>%
#     layout(
#       scene = list(
#         yaxis = list(title = "EvalRange"),
#         xaxis = list(title = "FD Sims"),
#         zaxis = list(title = "Value")
#       )) %>%
#     layout(title = titleVal, scene = scene)
# }
#
#
# #' Title
# #'
# #' @param fd_eval
# #' @param singleRangeSeq
# #' @param titleVal
# #'
# #' @return
# #' @export
# #'
# #' @examples
# plot_evalfd_3dlines_plotly <- function(fd_eval, singleRangeSeq,
#                                        titleVal=NULL){
#
#   number <- length(fd_eval[1,])
#   valRange <- c(min(fd_eval),max(fd_eval))
#
#   plotData <- data.frame('evalRange'=singleRangeSeq,
#                          'FDRep'=1,
#                          'Value'=fd_eval[,1])
#
#   for(i in 2:number){
#     plotData <- rbind(plotData,
#                       data.frame('evalRange'=singleRangeSeq,
#                                  'FDRep'=i,
#                                  'Value'=fd_eval[,i])
#     )
#   }
#
#   scene = list(camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))
#
#   tmpColors <- brewer.pal(11,"Spectral")
#   tmpColors[6] <- 'yellow'
#
#   plot_ly(plotData,
#           x = ~as.factor(FDRep), y = ~evalRange, z = ~Value,
#           type = 'scatter3d', mode = 'lines',
#           color = ~as.factor(FDRep),
#           colors = tmpColors) %>%
#     #colors = c("red", "yellow", "blue")) %>%
#     #colors='Spectral') %>%
#     layout(
#       scene = list(
#         yaxis = list(title = "EvalRange"),
#         xaxis = list(title = "FD Sims"),
#         zaxis = list(title = "Value")
#       )) %>%
#     layout(title = titleVal, scene = scene) %>%
#     layout(showlegend = FALSE)
#   #line = list(width = 4, color = ~as.factor(FDRep),
#   #            colorscale = list(c(0,'#BA52ED'), c(1,'#FCB040'))))
# }
