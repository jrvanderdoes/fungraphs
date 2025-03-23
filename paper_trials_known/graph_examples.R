library(ggraph)
library(devtools)
library(ade4)

load_all()

path <- 'C://Users//jerem//OneDrive//Documents//School//Waterloo//Research//RPackages//fungraphs//paper_trials_known'

set.seed(1234)
points <- t(MASS::mvrnorm(20,mu = rep(0,2),diag(1,2,2)))
data_dist <- calculateDistanceMatrix(
  points, silent = T, errType='L2', dataUse = 'Orig')

mdt_tree <- createMDT(data_dist,1)
# nnl_tree <- RANN::nn2(t(points),t(points),k=2)$nn.idx#nnl(data_dist,1)
nnl_tree <- nnl1(data_dist,1)
# nnl_tree <- nnl1(data_dist,1)
mst_tree <- mstree(as.dist(data_dist), ngmax=1)

# MDT
graph <- tidygraph::as_tbl_graph(mdt_tree)

png(paste0(path,"\\figures\\mdp.png"),width = 800,height = 800)
ggraph(graph,layout=t(points)) +
  geom_node_point(size=10) +
  geom_edge_link(width=2) +
  theme_graph(background = 'white')
dev.off()

# MST
graph <- tidygraph::as_tbl_graph(neig2mat(mst_tree))
png(paste0(path,"\\figures\\mst.png"),width = 800,height = 800)
ggraph(graph,layout=t(points)) +
  geom_node_point(size=10) +
  geom_edge_link(width=2) +
  theme_graph(background = 'white')
dev.off()


# NNL
graph <- tidygraph::as_tbl_graph(nnl_tree)
png(paste0(path,"\\figures\\nnl.png"),width = 800,height = 800)
ggraph(graph,layout=t(points)) +
  geom_node_point(size=10) +
  geom_edge_link(width=2) +
  theme_graph(background = 'white')
dev.off()

##################
library(igraph)

set.seed(1234)
nnl_tree2 <- nnl1(data_dist,2)

graph <- tidygraph::as_tbl_graph(nnl_tree2)
edge_colors <- data.frame(
  edge_id = seq_len(ecount(graph)),
  color = c(rep("black", 20), rep("red", ecount(graph) - 20))
)

png(paste0(path,"\\figures\\nnl_2.png"),width = 800,height = 800)
ggraph(graph,layout=t(points)) +
  geom_node_point(size=10) +
  geom_edge_link(aes(color = edge_colors$color,
                     linetype =ifelse(edge_colors$edge_id <= 20,
                                      'dashed', 'solid')),width=2) +
  scale_edge_color_identity() +
  theme_graph(background = 'white') +
  theme(legend.position="none")
dev.off()






# library(igraph)
# library(ggraph)
# library(tidygraph)
