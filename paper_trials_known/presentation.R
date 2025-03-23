library(ggraph)
library(igraph)
library(devtools)
library(ade4)

load_all()

path <- 'C:/Users/jerem/Downloads/'

set.seed(1234)
points <- t(MASS::mvrnorm(20,mu = rep(0,2),diag(1,2,2)))
data_dist <- calculateDistanceMatrix(
  points, silent = T, errType='L2', dataUse = 'Orig')

nnl_tree2 <- nnl1(data_dist,2)

edge_cols <- apply(nnl_tree2,MARGIN = 1,FUN = function(x){
  if(min(x)<=10) return('red')

  'black'
})

graph <- tidygraph::as_tbl_graph(nnl_tree2)
edge_colors <- data.frame(
  edge_id = seq_len(ecount(graph)),
  color = edge_cols#c(rep("black", 20), rep("red", ecount(graph) - 20))
)

png(paste0(path,"graph_random.png"),width = 800,height = 800)
ggraph(graph,layout=t(points)) +
  geom_node_point(size=10) +
  geom_edge_link(aes(color = edge_colors$color,
                     linetype =ifelse(edge_colors$edge_id <= 20,
                                      'dashed', 'solid')),width=2) +
  scale_edge_color_identity() +
  theme_graph(background = 'white') +
  theme(legend.position="none")
dev.off()


set.seed(1234)
points[,1:10] <- points[,1:10] +5
data_dist <- calculateDistanceMatrix(
  points, silent = T, errType='L2', dataUse = 'Orig')

nnl_tree2 <- nnl1(data_dist,2)

edge_cols <- apply(nnl_tree2,MARGIN = 1,FUN = function(x){
  if(min(x)<=10) return('red')

  'black'
})

graph <- tidygraph::as_tbl_graph(nnl_tree2)
edge_colors <- data.frame(
  edge_id = seq_len(ecount(graph)),
  color = edge_cols#c(rep("black", 20), rep("red", ecount(graph) - 20))
)

png(paste0(path,"graph_mean.png"),width = 800,height = 800)
ggraph(graph,layout=t(points)) +
  geom_node_point(size=10) +
  geom_edge_link(aes(color = edge_colors$color,
                     linetype =ifelse(edge_colors$edge_id <= 20,
                                      'dashed', 'solid')),width=2) +
  scale_edge_color_identity() +
  theme_graph(background = 'white') +
  theme(legend.position="none")
dev.off()
