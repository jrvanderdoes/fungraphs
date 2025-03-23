#' @noRd
#' @keywords internal
nnl1 = function(distance, K) {
  n <- nrow(distance)
  num_neigs <- 1

  distance = as.matrix(distance)
  maxdis = 1e5*max(distance)
  nodes = dim(distance)[1]
  tempE = vector("list",K)
  E = matrix(0,0,2)
  for (i in 1:K) {
    tempE[[i]] <- matrix(NA,ncol=1+num_neigs,nrow=n)
    for(j in 1:n){
      tempE[[i]][j,] <- c(j,
                          FastKNN::k.nearest.neighbors(j, distance, k = num_neigs)
                          )
    }
    # tempE[[i]] = nnlink1(distance)
    E = rbind(E,tempE[[i]])
    E = unique(E)
    if (i == K) break
    for (e in 1:dim(tempE[[i]])[1]) {
      e1 = tempE[[i]][e,1]
      e2 = tempE[[i]][e,2]
      distance[e1,e2] = distance[e2,e1] = maxdis
    }
  }
  return(E)
}



#' @noRd
#' @keywords internal
nnlink1 = function(distance){
  nodes = dim(distance)[1]
  temp = nnlink_Com(distance)
  Com = temp$Com
  adj = temp$adj
  edgenum = temp$edgenum
  Components = length(Com)
  while (1){
    if(Components ==1){
      E = matrix(0,edgenum,2)
      e = 1
      for (i in 1:nodes){
        if (length(adj[[i]])>0){
          for (j in 1:length(adj[[i]])){
            E[e,1] = i
            E[e,2] = adj[[i]][j]
            # Remove NA and ones with distance to self
            #   Should already be gone, but temp fix for now
            if(is.na(E[e,1]) || is.na(E[e,2]) || E[e,1] == E[e,2]){
              stop('sadas')
            }

            adj[[E[e,2]]] = setdiff(adj[[E[e,2]]],i)
            e = e+1
          }
        }
      }
      break
    } else{
      newdist = ID_edges_candidate = matrix(0,Components,Components)
      edges_candidate = vector("list",Components*(Components-1)/2)
      edge_com = 1
      for (i in 1:(Components-1)){
        for (j in (i+1):Components){
          tt = getComdist(Com[[i]],Com[[j]],distance)
          newdist[i,j] = newdist[j,i] = tt$mindis
          edges_candidate[[edge_com]] = tt$minID
          ID_edges_candidate[i,j] = ID_edges_candidate[j,i] = edge_com
          edge_com = edge_com + 1
        }
      }
      temp2 = nnlink_Com(newdist)
      adj2 = temp2$adj
      for (i in 1:Components){
        if (length(adj2[[i]])>0){
          for (j in 1:length(adj2[[i]])){
            e1 = i
            e2 = adj2[[i]][j]
            id_com = ID_edges_candidate[e1,e2]
            addedge = edges_candidate[[id_com]]
            for (adde in 1:dim(addedge)[1]){
              adj[[addedge[adde,1]]] = c(adj[[addedge[adde,1]]],addedge[adde,2])
              adj[[addedge[adde,2]]] = c(adj[[addedge[adde,2]]],addedge[adde,1])
            }
            adj2[[e2]] = setdiff(adj2[[e2]],i)
          }
        }
      }
      edgenum = 0
      for (i in 1:nodes){
        adj[[i]] = sort(unique(adj[[i]]))
        edgenum = edgenum + length(adj[[i]])
      }
      edgenum = edgenum/2
      Com2  = temp2$Com
      Components = length(Com2)
      newCom = vector("list",Components)
      for (i in 1:Components){
        for (j in 1:length(Com2[[i]])){
          newCom[[i]] = c(newCom[[i]], Com[[Com2[[i]][j]]])
        }
        newCom[[i]] = sort(newCom[[i]])
      }
      Com = newCom
    }
  }

  return(E)
}


#' @noRd
#' @keywords internal
getComdist = function(g1,g2,distance){
  tempdis = distance[g1,g2]
  mindis = min(tempdis)
  minid = which(tempdis==mindis, arr.ind=TRUE)
  minID = cbind(g1[minid[,1]],g2[minid[,2]])

  return(list(mindis=mindis, minID=minID))
}


#' @noRd
#' @keywords internal
nnlink_Com = function(distance){
  nodes = dim(distance)[1]
  adj = vector("list",nodes)
  Com = vector("list",nodes)
  Components = 0

  for (i in 1:nodes) adj[[i]] = rep(0,0)
  for (i in 1:nodes){
    edgeid = which(distance[i,] == min(distance[i,-i]))
    adj[[i]] = c(adj[[i]],edgeid)

    for (j in 1:length(edgeid)){
      adj[[edgeid[j]]] = c(adj[[edgeid[j]]],i)
    }
  }
  edgenum = 0
  for (i in 1:nodes){
    adj[[i]] = sort(unique(adj[[i]]))
    edgenum = edgenum + length(adj[[i]])
  }
  edgenum = edgenum/2
  visited = rep(0,nodes)
  for (i in 1:nodes){
    if (visited[i] == 0){
      visited = dfs(i,visited,adj)
      Components = Components + 1
      Com[[Components]] = which(visited==1)
    }
  }

  if (Components > 1){
    for (i in Components:2){
      Com[[i]] = sort(setdiff(Com[[i]],Com[[i-1]]))
    }
    Com[[1]] = sort(Com[[1]])
  }else{
    Com[[1]] = sort(Com[[1]])
  }
  Com = Com[1:Components]

  return(list(Com=Com, adj=adj, edgenum=edgenum))
}


#' @noRd
#' @keywords internal
dfs = function(s,visited,adj){
  visited[s] = 1
  for (i in 1:length(adj[[s]])){
    if (visited[adj[[s]][i]] == 0){
      visited = dfs(adj[[s]][i],visited,adj)
    }
  }
  return(visited)
}

