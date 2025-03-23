rlaplace = function(n,mu,sigma){
  U = runif(n,0,1)
  #This will give negative value half of the time
  sign = ifelse(rbinom(n,1,.5)>.5,1,-1)
  y = mu + sign*sigma/sqrt(2)*log(1-U)
  y
}


data <- gen_FD_KL_Expansion(
  ns = c(500,500),
  eigsList = list(c(5,3,1,0.5),
                  c(5,3,1,0.5)),
                  #c(2,1.75,1.5,1.25)/2),
  basesList = list(create.bspline.basis(nbasis=4, norder=4),
                   create.bspline.basis(nbasis=4, norder=4)),
  meansList = c(0,0),
  distsArray = c('Normal','Laplace'),
  evals = seq(0,1,length.out=20),
  kappasArray = 0,
  silent = T)
#plot_fd(data)
getPts(data, method = graph_method_gen_giveRawData,
       nTrees=5, treeType='MDT')
