require(igraph)
require(Matrix)
require(plyr)
require(MASS)
require(data.table)
Rcpp::sourceCpp('./Rcpp_PriorityQueue_Ball_Grow.cpp')

genMeshCoord3D <- function(dimVec,offset) {
out <- matrix(data=0,nrow=dimVec[1]*dimVec[2]*dimVec[3],ncol=3)
ii=1;
for (zz in 1:dimVec[3]) {
for (yy in 1:dimVec[2]) {
for (xx in 1:dimVec[1]) {
out[ii,]=c(xx,yy,zz)+runif(3,-offset,offset);
ii=ii+1;
}
}
}
return(out)
}

genMesh2D <- function(dimVec,cOffset=0.0,periodic=FALSE,eWeight=1,vWeight=0,
                      directed=FALSE) {
  #dimVec - the dimension sizes of the graph, must be array of 2 integers
  #cOffset - random offset for coordinates (helpful for periodic or multi edge)
  #periodic - true / false - whether mesh is periodic
  #eWeight - weights of edges. Can be scalar or 1D array of values
  #           be sure array dimension matches number of edges
  #vWeight - weights of vertices. Can be scalar or array of values
  myGraph <- graph.lattice(dimVec,circular=periodic)
  myGraph$coordinates <- genMeshCoord2D(dimVec,cOffset)
  V(myGraph)$id <- seq(vcount(myGraph))
  V(myGraph)$weight <- vWeight
  if (directed) {
    myGraph<-as.directed(myGraph)
    myEdgeList <- E(myGraph)
    E(myGraph)$id <- seq(ecount(myGraph))
    E(myGraph)$weight <- eWeight
  } else {
    myEdgeList <- E(myGraph)
    E(myGraph)$id <- seq(ecount(myGraph))
    E(myGraph)$weight <- eWeight    
  }
  return(myGraph)
}

genMesh3D <- function(dimVec,cOffset=0.0,periodic=FALSE,eWeight=1,vWeight=0) {
  #dimVec - the dimension sizes of the graph, must be array of 2 integers
  #cOffset - random offset for coordinates (helpful for periodic or multi edge)
  #periodic - true / false - whether mesh is periodic
  #eWeight - weights of edges. Can be scalar or 1D array of values
  #           be sure array dimension matches number of edges
  #vWeight - weights of vertices. Can be scalar or array of values
  myGraph <- graph.lattice(dimVec,circular=periodic)
  myGraph$coordinates <- genMeshCoord3D(dimVec,cOffset)
  myEdgeList <- E(myGraph)
  V(myGraph)$id <- seq(vcount(myGraph))
  E(myGraph)$id <- seq(ecount(myGraph))
  V(myGraph)$weight <- vWeight
  E(myGraph)$weight <- eWeight
  return(myGraph)
}

genMesh2DLap <- function(dimVec,periodic=FALSE,eWeight=-1,vWeight=0,verbose=FALSE) {
  xm = dimVec[1]; ym=dimVec[2]
  coordToLinMap = array(data=0,dim=dimVec)
  linToCoordMap = array(data=0,dim=c(xm*ym,2))
  yiL=1:ym
  for (xi in 1:xm) {
    coordToLinMap[xi,]=xi+xm*((1:ym)-1)
    xiL=as(array(data=xi,dim=ym),"numeric");
    linToCoordMap[xi+xm*((1:ym)-1),] = t(rbind(xiL,yiL))
  }
  if (periodic) {
    tList = as(array(data=0,dim=4*xm*ym),"numeric")
    sList = as(array(data=0,dim=4*xm*ym),"numeric")
    for (li in 1:(xm*ym)) {
      xi=linToCoordMap[li,1];yi=linToCoordMap[li,2];
      xim=((xi-2)%%xm)+1;xip=(xi%%xm)+1
      yim=((yi-2)%%ym)+1;yip=(yi%%ym)+1
      tList[4*(li-1)+(1:4)] = c(coordToLinMap[xim,yi],coordToLinMap[xip,yi],
                                      coordToLinMap[xi,yim],coordToLinMap[xi,yip])
      sList[4*(li-1)+(1:4)] = c(li,li,li,li)
    }
  } else {
    tList=c();sList=c()
    for (yi in 1:(ym-1)) {
      for (xi in 1:(xm-1)) {
        li = coordToLinMap[xi,yi]
        lxip = coordToLinMap[xi+1,yi]
        lyip = coordToLinMap[xi,yi+1]
        sLTemp = c(li,li)
        tLTemp = c(lxip,lyip)
        tList = append(tList,c(tLTemp,sLTemp))
        sList = append(sList,c(sLTemp,tLTemp))
      }
    }
    for (xi in 1:(xm-1)) {
      li = coordToLinMap[xi,ym]
      lxip=coordToLinMap[xi+1,ym]
      tLTemp=c(lxip,li)
      sLTemp=c(li,lxip)
      tList=append(tList,tLTemp)
      sList=append(sList,sLTemp)
    }
    for (yi in 1:(ym-1)) {
      li=coordToLinMap[xm,yi]
      lyip=coordToLinMap[xm,yi+1]
      tLTemp=c(lyip,li)
      sLTemp=c(li,lyip)
      tList=append(tList,tLTemp)
      sList=append(sList,sLTemp)
    }
  }
  if (verbose) print(rbind(tList,sList))
  outmat=sparseMatrix(i=tList,j=sList,x=eWeight)
  for (ii in 1:dim(outmat)[1]) {
    outmat[ii,ii] = vWeight-sum(outmat[ii,])
  }
  return(outmat)
}

genMeshCoord2D <- function(dimVec,offset) {
out <- matrix(data=0,nrow=dimVec[1]*dimVec[2],ncol=2)
ii=1;
for (yy in 1:dimVec[2]) {
for (xx in 1:dimVec[1]) {
out[ii,]=c(xx,yy)+runif(2,-offset,offset);
ii=ii+1;
}
}
return(out)
}

genMeshCoord3D <- function(dimVec,offset) {
  out <- matrix(data=0,nrow=dimVec[1]*dimVec[2]*dimVec[3],ncol=3)
  ii=1;
  for(zz in 1:dimVec[3]){
    for (yy in 1:dimVec[2]) {
      for (xx in 1:dimVec[1]) {
        out[ii,]=c(xx,yy,zz)+runif(3,-offset,offset);
        ii=ii+1;
      }
    }
  }
  return(out)
}

getSubMeshCoords <- function(coords,idVec) {
  subCoords <- coords[idVec,]
  return(subCoords)
}

plotSubgraph <- function(g,s,vSize=10,lSize=.625,arrowSize=.4,highlight="#DC143C",
                         vLabel=c(),eWidth=1) {
  myPlotGraph <- g
  sgVmap <- array(data=0,dim=length(V(s))); sgEmap <- array(data=0,dim=length(E(s)))
  sVids<-V(s)$id; sEids<-E(s)$id
  gVids<-V(g)$id; gEids<-E(g)$id
  for (ii in 1:length(V(s))) {
    tempind = match(sVids[ii],gVids,nomatch=0)
    if (tempind !=0) {sgVmap[ii]=tempind} else {
      stop("Vertex IDs of g and s do not match at ",ii,"!")}
  }
  for (ii in 1:length(E(s))) {
    tempind = match(sEids[ii],gEids,nomatch=0)
    if (tempind !=0) {sgEmap[ii]=tempind} else {
      stop("Edge IDs of g and s do not match at ",ii,"!")}
  }
  V(myPlotGraph)$color <- "#00FFFF" #cyan
  V(myPlotGraph)[sgVmap]$color <- highlight
  E(myPlotGraph)$color <- "gray"
  E(myPlotGraph)[sgEmap]$color <- highlight #crimson
  if (length(vLabel) > 0){
    plot.igraph(myPlotGraph,layout=myPlotGraph$coordinates[gVids,],
                vertex.size=vSize,vertex.label.cex=lSize,
                edge.arrow.size=arrowSize,vertex.label=vLabel,
                edge.width=eWidth)
  }else{
    plot.igraph(myPlotGraph,layout=myPlotGraph$coordinates[gVids,],
              vertex.size=vSize,vertex.label.cex=lSize,
              edge.arrow.size=arrowSize,edge.width=eWidth)}
}

extractTree <- function(g,roots) {
r <- graph.bfs(g,root=roots,neimode='all',order=TRUE,father=TRUE)
h <- graph(rbind(r$order,r$father[r$order])[,-1],directed=FALSE)
return(h)
}

findBall <- function(g,r,d,weighted=FALSE,directed=FALSE,verbose=FALSE,distances=FALSE,parents=FALSE) {
#input: 
#	graph g - an igraph graph object
#	root node r - the integer index for the root node of the ball
#	distance d - the maximum distance allowed (i.e. the ball radius)
#output: an integer vector containing the indices of all nodes in the ball
  
#if the graph is unweighted, we use a simple breadth first search based
#approach. If not, we will use the Rcpp growBall algorithm
  if (!weighted) {
	  ball <- c() #initiallize the ball to the root node
	  myBFS <- graph.bfs(g,root=r,dist=TRUE)
	  for (i in 1:length(myBFS$dist)){
		  if (myBFS$dist[i] <= d) {
			  ball = append(ball,i)
		  }
	  }
	  return(ball)
  } else {
    #since the graph is weighted, findBall will just act as a
    #wrapper for the growBall algorithm.
    #if the graph g is not a directed graph, it needs to be made into one.
    if (!directed) {
      gd = as.directed(g);
    } else {
      gd = g;
    }
    gvl = V(gd);
    gel = get.edgelist(gd);
    gelSources = gel[,1]; gelTargets=gel[,2];
    gew = E(gd)$weight;
    ball = growBall(vl=gvl,elSources=gelSources,elTargets=gelTargets,
                    ew=gew,rv=r,r=d,verbose=verbose,distances=distances,
                    parents=parents);
    return(ball)
  }
}

findCone <- function(g,s,t,d) {
  cone <- c()
  mySourceBFS <- graph.bfs(g,root=s,dist=TRUE)
  myTargetBFS <- graph.bfs(g,root=t,dist=TRUE)
  myDistList <- abs(mySourceBFS$dist[t]+ myTargetBFS$dist-mySourceBFS$dist )
  for (i in 1:length(myDistList)){
    if (myDistList[i] <= d) {
      cone = append(cone,i)
    }
  }
  return(cone)
}

findPetal <- function(g,s,t,r,verbose=FALSE,distances=FALSE,parents=FALSE) {
  gDir <- as.directed(g,mode="mutual")
  gDirElist <- get.edgelist(gDir)
  gDirEweights <- E(gDir)$weight
  gDistMap <- growBall(vl=V(g),elSources=gDirElist[,1],elTargets=gDirElist[,2],
                       ew=gDirEweights,rv=s,r=9999,distances=TRUE,parents=TRUE)
  gDistMap[1,]=gDistMap[1,]+1
  gDist = gDistMap[1,]
  gDist[gDistMap[1,]] = gDistMap[2,]
  gPar=gDistMap[1,]; gPar[gDistMap[1,]] = gDistMap[3,]+1
  path=c(t)
  tempv=t
  while (tempv != s){
    tempv=gPar[tempv]
    path=append(path,tempv)
  }
  if(verbose) print("indices of edge-sources matching members of path-edge source vertices: ")
  pSourceInds=which(gDirElist[,1] %in% path[2:length(path)])
  if(verbose) {
    print(pSourceInds); print("The corresponding vertex indices are: ")
    print(gDirElist[pSourceInds,1])
    print("indices of edge-targets along path originating from sources along path:")
  }
  pTargetInds=pSourceInds[which(gDirElist[pSourceInds,2] %in% path[1:(length(path)-1)])]
  if(verbose) { 
    print(pTargetInds); print("the corresponding vertex indices are: ");
    print(gDirElist[pTargetInds,2]);
    print("path (target --> source): "); print(path); print("gDist: ");
    print(gDist); print("gDistMap"); print(gDistMap);
  }
  cEweights = gDirEweights - (gDist[gDirElist[,2]]-gDist[gDirElist[,1]])
  #gDirEweights[gDistMap[1,]]=(gDirEweights[gDistMap[1,]]^-1 + gDistMap[2,])^-1
  if(verbose) { print("cEweights: "); print(cEweights) }
  if (floor(r) > 0) {
    cEweights[pTargetInds[floor(r)]] = cEweights[pTargetInds[floor(r)]] / 2
  }
  petalMap = growBall(vl=V(g),elSources=gDirElist[,1],elTargets=gDirElist[,2],
                  ew=cEweights,rv=t,r=r/2,distances=distances,parents=parents)
  petalMap[1,] = petalMap[1,]+1
  if(!(distances || parents)) {
    return(petalMap[1,])
  } else {
    return(petalMap)
  }
}

extractGraphData <- function(G,giveEdata=TRUE,giveBallMap=TRUE,giveDistMap=TRUE,
                             giveParMap=TRUE,rv=V(G)[floor(length(V(G))/2)]) {
  outlist=list()
  if (is.null(E(G)$weight)) {
    warning("No edge weights defined, using weight of 1")
    E(G)$weight<-1
  }
  gEdata <- t(rbind(t(get.edgelist(G)),E(G)$weight))
  gBallMap <- growBall(V(G),gEdata[,1],gEdata[,2],gEdata[,3],rv=rv,r=999999,
                       distances=giveDistMap,parents=giveParMap)
  if (giveDistMap) {
    gDist <- gBallMap[1,]; gDist[gBallMap[1,]+1] <- gBallMap[2,]
    outlist[["dist"]] = gDist
    if (giveParMap) {
      gPar <- gBallMap[3,]; gPar[gBallMap[1,]+1] <- gBallMap[3,]+1
      outlist[["par"]] <- gPar
    }
  } else {
    if (giveParMap) {
      gPar <- gBallMap[1,]; gPar[gBallMap[1,]+1] <- gBallMap[2,]+1
      outlist[["par"]] <- gPar
    }
  }
  if (giveBallMap) {
    outlist[["ball"]] <- gBallMap
  }
  if (giveEdata) {
    outlist[["Edata"]] <- gEdata
  }
  return(outlist)
}
  
#generate a mapping from l2 to l1. l2 must be a subset of l1.
#If l1 contains duplicates, the mapping may not function properly
#this should be slightly faster than the unordered version
#for long lists
orderedSublistMap <- function(l1,l2) {
  l1Len = length(l1); l2Len = length(l2)
  if (l2Len > l1Len) stop("ERROR length of l2 greater than length of l1")
  ii=1
  lmap= array(data=0,dim=l2Len)
  tempind = match(l2[ii],l1,nomatch=0);
  if (tempind == 0) stop("no match found first item! l2 is not sublist of l1")
  while (tempind != 0 && tempind < l1Len && ii < l2Len) {
    lmap[ii]=tempind; ii=ii+1
    temp2=match(l2[ii],l1[(tempind+1):l1Len],nomatch=0)
    if (temp2 != 0) {tempind = tempind+temp2} else {
      tempind=0;
      warning("no match found for ",ii,"th item, l2 may not be a sublist of l1!")
    }
  }
  if (ii <= l2Len && tempind < l1Len && tempind != 0) {
    lmap[ii] = tempind 
  }
  if (ii < l2Len) warning("lmap is short!")
  return(lmap)
}

unorderedSublistMap <- function(l1,l2,wanings=TRUE) {
  l1Len=length(l1);l2Len=length(l2)
  if (l2Len > l1Len) stop("length of l2 is greater than length of l1")
  lmap=array(data=0,dim=l2Len)
  for (ii in 1:l2Len) {
    lmap[ii] = match(l2[ii],l1,nomatch=0)
  }
  if (match(0,lmap,nomatch=0) != 0) {
    warning("unmatched entries in lmap!")
  }
  return(lmap)
}

constructPath <- function(x0,t,vPar,maxLen=1000) {
  path=c(t)
  tempv=t
  it=0
  while (tempv != x0 && it < maxLen) {
    tempv=vPar[tempv]
    path=append(path,tempv)
    it=it+1
  }
  return(path)
}

getPathEdgeInds <- function(path,Edata,Par) {
  tempInds=which(Edata[,1] %in% path[1:(length(path)-1)])
  tempInds=tempInds[which(Edata[tempInds,2] %in% path[2:length(path)])]
  pxEinds=c()
  for (ii in 1:length(tempInds)) {
    if (Edata[tempInds[ii],2] == Par[xEdata[tempInds[ii],1]]) {
      pxEinds=append(pxEinds,tempInds[ii])
    }
  }
  pEorder=array(data=0,dim=length(pxEinds))
  for (ii in 1:length(pxEinds)) {
    pEorder[ii] = match(path[ii],Edata[pxEinds,1],nomatch=0)
    if (pEorder[ii] == 0) stop("error matching pathing indices")
  }
  pxEinds = pxEinds[pEorder]
  return(pxEinds)
}

mapPathVerticesToSubgraph <- function(xpath,gVids,sVids) {
  indList = which(gVids[xpath] %in% sVids)
  if (max(indList[2:length(indList)]-indList[1:(length(indList)-1)]) > 1 ) {
    warning("the path in y is not contiguous!")
  }
  return(match(gVids[xpath[indList]],sVids))
}

mapPathEdgesToSubgraph <- function(pgEinds,gEids,sEids) {
  matchList <- match(gEids[pgEinds],sEids,nomatch=0)
  zeroIds<-which(matchList %in% 0)
  if (length(zeroIds) > 0){
    if (max(zeroIds[2:length(zeroIds)]-zeroIds[1:(length(zeroIds)-1)])>1) {
      warning("Y seems to intersect multiple non-contiguos segements of path in X!") 
    }
  }
  matchList <- matchList[which(!(matchList %in% 0))]
  return(matchList)
}

mapEdgeWeightListToSubgraph <- function(xEweights,xEids,yEids) {
  yEweights <- array(data=0,dim=length(yEids))
  count=0
  matchList<-match(xEids,yEids,nomatch=0)
  xIdMatchLocs <- which(!(matchList %in% 0))
  yIdMatchLocs <- matchList[xIdMatchLocs]
  count=length(yEids)-length(yIdMatchLocs)
  if (count > 0) warning(count," edges in Y did not get mapped")
  yEweights <- array(data=0,dim=length(yEids))
  yEweights[yIdMatchLocs] <- xEweights[xIdMatchLocs]
  #for (ii in 1:length(xEids)) {
  #  tempind = match(xEids[ii],yEids,nomatch=0)
  #  if (tempind > 0) {
  #    yEweights[tempind] = xEweights[ii]
  #    count=count+1
  #  }
  #}
  return(yEweights)
}

calcSubgraphSurfaceWeight <- function(G,S,gEweights) {
  sVids <- V(S)$id; sEids <- E(S)$id
  gVids <- V(G)$id; gEids <- E(G)$id
  sComp <- induced.subgraph(G,which(!(gVids %in% sVids)))
  sCompEids <- E(sComp)$id
  return(sum(gEweights[which(!(gEids[which(!(gEids %in% sEids))] %in% sCompEids))]))
}

#convert edge data for an edge weight only directed graph into matrix representation
#of the corresponding vertex weight only undirected graph
#Edata should be a m by 3 matrix c1=source vertices, c2=target vertices, c3=weights
edgeDataToLap <- function(Edata,Vweights=c()) {
  outmat = sparseMatrix(i=Edata[,1],j=Edata[,2],x=-Edata[,3])
  vWeighted = FALSE
  if (!is.null(Vweights) && length(Vweights) == dim(outmat)[1]) vWeighted = TRUE
  for (i in 1:(dim(outmat)[1])) {
    outmat[i,i] = -sum(outmat[i,])
    if (vWeighted) outmat[i,i] = outmat[i,i] + Vweights[i]
  }
  return(outmat)
}

#Creates a petal in sugraph Y of X with target node t, source node x0
#and max radius R for reference see the article by Abraham and Ofer in:
#STOC 2012: "Using Petal-Decompositions to Build a Low Stretch Spanning Tree"
#Returns the resulting petal Wr and the corresponding node along the
#path from x0 to t, pr.
#The petal is grown iteratively using an Rcpp function 'expandBall'
#
#Y must be a subgraph of X created using induced.subgraph, and X must
#have id values mapped to its edges and vertices to allow proper mapping
#of distance metrics between X and Y.
createPetal <- function(X,Y,t,x0,R,
                        xDistMap=c(),xDist=c(),xPar=c(),
                        xcEweights=c(), path=c(), pxEinds=c(),
                        showTimings=FALSE) {
  ptm <- proc.time()
  if (!is.directed(X) || !is.directed(Y)) {
    stop("X and Y must both be directed graphs!")
  }
  L = ceiling(log(log(length(V(Y)))))
  if (showTimings) print(paste("L=",L,"time: "))
  if (showTimings) print(ptm-proc.time())
  X_vlength = length(V(X))
  #Edge and vertex ids will be used to map between X and Y
  if (showTimings) print("building id arrays")
  xEids <- E(X)$id 
  yEids <- E(Y)$id 
  xVids <- V(X)$id
  yVids <- V(Y)$id
  yEdata <- get.edgelist(Y)
  #make sure that X is directed then generate edge data for X
  if (!is.directed(X)) {
    stop("X is not a directed graph, it will be converted to one!")
    #xEdata <- get.edgelist(as.directed(X))
    #xEdata = t(rbind(t(xEdata),E(as.directed(X))$weight))
  }
  else {
    xEdata <- t(rbind(t(get.edgelist(X)),E(X)$weight))
  }
  #build distance map for X using the growBall routine if data was not provided
  if (showTimings) print("generating X distance maps, time: ")
  if (showTimings) print(ptm-proc.time())
  #build xDistMap and its corresponding vertex data lists if not provided
  if (length(xDistMap) < 3*X_vlength) {
    xDistMap <- growBall(vl=V(X),elSources=xEdata[,1],elTargets=xEdata[,2],
                       ew=xEdata[,3],rv=x0,r=9999,distances=TRUE,parents=TRUE)
    #adjust for off by 1 error if needed
    if (min(xDistMap[1,])==0) {
      xDistMap[1,] = xDistMap[1,]+1
    }
    if (min(xDistMap[3,]) == 0) {
      xPar <- xDistMap[1,]; xPar[xDistMap[1,]] <- xDistMap[3,]+1
    } else {
      xPar <- xDistMap[1,]; xPar[xDistMap[1,]] <- xDistMap[3,]
    }
    xDist <- xDistMap[1,];xDist[xDistMap[1,]] <- xDistMap[2,]
  }
  #If xDistMap was given, check to see if the remapped parent and distance
  #lists were also provided or build them if not
  if (length(xDist) < X_vlength) {
    #fix off by 1 error in xDistMap[1,] if needed
    if (min(xDistMap[1,])==0) {
      xDistMap[1,] = xDistMap[1,]+1
    }
    xDist <- xDistMap[1,];xDist[xDistMap[1,]] <- xDistMap[2,]
  }
  if (length(xPar) < X_vlength) {
    #fix off by 1 error in xDistMap[1,] if needed
    if (min(xDistMap[1,])==0) {
      xDistMap[1,] = xDistMap[1,]+1
    }
    #check for off by 1 error in xDistMap[3,] and act accordingly
    if (min(xDistMap[3,]==0)) {
      xPar <- xDistMap[1,]; xPar[xDistMap[1,]] <- xDistMap[3,]+1      
    } else {
      xPar <- xDistMap[1,]; xPar[xDistMap[1,]] <- xDistMap[3,]
    }
  }
  #compute and store edge volume of X for later use
  net_X_ew = sum(xDist) 
  #use distance map to construct cone metric edge weights if not provided
  if (length(xcEweights) < length(xEids)) {
    xcEweights = xEdata[,3] - (xDist[xEdata[,2]]-xDist[xEdata[,1]]) 
  }
  #build the path from x0 to t in X if not provided
  if (length(path)<=1){
    if (showTimings) print("generating path in X time: ")
    if (showTimings) print(ptm-proc.time())
    path = constructPath(x0,t,xPar)
  }

  #build edge list for path in X if not provided
  if (length(pxEinds) < 1) {
    pxEinds = getPathEdgeInds(path,xEdata,xPar)    
  }
  #map path edges in X to path edges in Y
  if (showTimings) print("mapping path Edges to Y time: ")
  if (showTimings) print(ptm-proc.time())
  pyEinds <- mapPathEdgesToSubgraph(pxEinds,xEids,yEids)
  #map path vertices in X to path vertices in Y
  if (showTimings) print("mapping path Vertices to Y time: ")
  if (showTimings) print(ptm-proc.time())
  ypath = mapPathVerticesToSubgraph(path,xVids,yVids)
  #yEdata <- get.edgelist(Y) #appears unnecessary
  #map conemetric of x0 in X onto subgraph Y
  if (showTimings) print("mapping path conemetric to Y time: ")
  if (showTimings) print(ptm-proc.time())
  ycEweights <- mapEdgeWeightListToSubgraph(xcEweights,xEids,yEids)
  #construct path distance list
  if (showTimings) print("calculating path distances in Y time: ")
  if (showTimings) print(ptm-proc.time())
  ypath_d = array(data=0,dim=length(ypath))
  for (ii in 1:length(pyEinds)) {
    ypath_d[ii+1] = ypath_d[ii]+ycEweights[pyEinds[ii]]
  }
  #begin building reference petal in Y
  if (showTimings) print("growing first petal time: ")
  if (showTimings) print(ptm-proc.time())
  xEsum = sum(xEdata[,3])
  wrad = R/2
  Wnew <- growBall(V(Y),yEdata[,1],yEdata[,2],ycEweights,
                   rv=ypath[1] ,r=wrad, verbose=FALSE, 
                   distances=TRUE,parents=TRUE)
  ii=1; foundLimit=FALSE
  m = length(E(X))/2
  pi=1 #move along path starting at target node
  while(ii <= L && !foundLimit) {
    if (showTimings) print(paste("starting iteration ",ii,"time: "))
    if (showTimings) print(ptm-proc.time())
    if (showTimings) print(paste("petal radius:",wrad,"path vertex: ",pi))
    Wa <- Wnew
    ppV=ypath[pi] #edge index of the pi'th point in the path
    ppE=pyEinds[pi] #edge index leading from pi-1 to pi
    ycEweights[ppE]=ycEweights[ppE]/2 
    wrad = (1 + pi / L)*R/2 #radius of petal to create
    Wnew <- expandBall(V(Y),t(rbind(t(yEdata[,1:2]),ycEweights)),rv=ypath[1],
                       ballData=Wa,r=wrad, pv=ppV)
    wEsum = sum(E(induced.subgraph(Y,Wnew[1,]+1))$weight)
    if (wEsum <= 2*xEsum / 2^(log(m)^(1-pi/L))) foundLimit=TRUE
    ii=ii+1
    pi=pi+1
  }
  if (showTimings) print(paste("finding Wa required",ii,"iterations ","path vertex #",pi-1,"time: "))
  if (showTimings) print(ptm-proc.time())
  alpha=(1+(pi-1)/L)*R/2 #end of step 2 (FINALLY!)
  wrad = alpha
  Wr = Wa
  chi = wEsum / sum(E(induced.subgraph(Y,Wa[1,]+1))$weight)
  #need to do steps 3 through 5#
  WrG <- induced.subgraph(Y,Wr[1,]+1)
  WrDEsum = calcSubgraphSurfaceWeight(X,WrG,xEdata[,3])
  WrEsum = sum(E(WrG)$weight)
  if (showTimings) print(paste("chi: ",chi,"WrSurfaceArea: ",WrDEsum,"WrVolume: ",WrEsum))
  if (WrDEsum >= (WrEsum * 8*L*log(chi)/R)) {
    Wr = Wnew
    WrG <- induced.subgraph(Y,Wr[1,]+1)
    WrDEsum = calcSubgraphSurfaceWeight(X,WrG,xEdata[,3])
    WrEsum = sum(E(WrG)$weight)
    while ((WrDEsum >= (WrEsum * 8*L*log(chi)/R)) && wrad <=R) {
      pi = pi+1
      ii=ii+1
      ppV=ypath[pi] #edge index of the pi'th point in the path
      ppE=pyEinds[pi] #edge index leading from pi-1 to pi
      ycEweights[ppE]=ycEweights[ppE]/2 
      wrad = (1 + pi / L)*R/2 #radius of petal to create
      Wr <- expandBall(V(Y),t(rbind(t(yEdata[,1:2]),ycEweights)),rv=ypath[1],
                         ballData=Wr,r=wrad, pv=ppV)
      WrG <- induced.subgraph(Y,Wr[1,]+1)
      WrDEsum = calcSubgraphSurfaceWeight(X,WrG,xEdata[,3])
      WrEsum = sum(E(WrG)$weight)
      if (showTimings) print(paste("finished",ii,"iterations","time: "))
      if (showTimings) print(ptm-proc.time())
      if (showTimings) print(paste("petal radius:",wrad,"path vertex: #",pi))
    }
    if (showTimings) print(paste("Wr was Wa. Found in ",ii,"iterations"))
    if (showTimings) print(paste("chi: ",chi,"WrSurfaceArea: ",WrDEsum,"WrVolume: ",WrEsum))
    return(list(Wr,yVids[match(yVids[ypath[pi]],xVids)]))
  } else {
    if (showTimings) print(paste("Wr was Wa. Found in ",ii,"iterations","path vertex #",pi-1,"time: "))
    if (showTimings) print(ptm-proc.time())
    if (showTimings) print(paste("chi: ",chi,"WrSurfaceArea: ",WrDEsum,"WrVolume: ",WrEsum))
    return(list(Wr,yVids[match(yVids[ypath[pi-1]],xVids)])) 
  }
}

findNeighbors <- function(v,Edata) { #find list of vertices connected to v
  return(Edata[which(Edata[,1] %in% v),2])
}

mostDistantNeighbor <- function(v,Edata,Vdist) {
  adjVinds = findNeighbors(v,Edata)
  return(adjVinds[which.max(Vdist[adjVinds])])
}

petalDecomposition <- function(X,x0,t) {
  decompList = c(); petalList = c(); portalList = c(); targetList = c()
  xEdata = t(rbind(t(get.edgelist(X)),E(X)$weight))
  xDistMap = growBall(V(X),xEdata[,1],xEdata[,2],xEdata[,3],x0,r=999999,
                      distances=TRUE,parents=TRUE)
  xDist<-xDistMap[1,];xDist[xDistMap[1,]+1]<-xDistMap[2,]
  xPar <- xDistMap[1,];xPar[xDistMap[1,]+1]<-xDistMap[3,]+1
  radx0 = max(xDist); r0 = radx0/2
  #Build First petal
  t1 = t
  ##find taret for first petal, should be closest vertex
  ##which is a distance of at least r0 from x0
  if (xDist[t1] < r0) { #step 1, check length of path to t
    while (xDist[t1] < r0) {
      t1 <- mostDistantNeighbor(t1,xEdata,xDist)
    }
  }
  else {
    foundt = FALSE
    while ((xDist[t1] > r0) && !foundt) {
      ttemp = xPar[t1]
      if (xDist[ttemp] < r0) {
        foundt=TRUE
      } else {
        t1 = ttemp
      }
    }
  }
  #add t1 to the targetList and build the first petal
  targetList[[1]] = t1
  WrL = createPetal(X,X,t,x0,r0/2,xDistMap,xDist,xPar)
  return(WrL)
  
}

eVec<-function(i,n){
  outvec = array(data=0,dim=n)
  outvec[i]=1
  return(outvec)
}
genRPmat <- function(pvec,nrow,ncol,sparse=TRUE) { #row perm matrix
  if (sparse) {
    return(genSRPmat(pvec,nrow,ncol))
  } else {
    return(genRowPmat(pvec,nrow,ncol))
  }
}
genRowPmat <- function(pVec,nrow,ncol) {#non-sparse format row perm. matrix
  if (length(pVec) < nrow) stop("permutation vector is short!")
  if (min(pVec) < 1) stop("negative permutations indices not supported")
  if (max(pVec) > ncol) stop("some values in permutation vector are out of bounds")
  pMat = array(data=0,dim=(c(nrow,ncol)))
  for (i in 1:nrow) {pMat[i,]=eVec(pVec[i],ncol)}
  return(pMat)
}
genSRPmat <- function(pVec,nrow,ncol) {#sparse row permuation matrix
  if (length(pVec) < nrow) stop("permutation vector is short!")
  if (min(pVec) < 1) stop("negative permutations indices not supported")
  if (max(pVec) > ncol) stop("some values in permutation vector are out of bounds")
  return(as(sparseMatrix(i=1:nrow,j=pVec,dims=c(nrow,ncol)),"dgCMatrix"))
}

rowNZcountVec <- function(mat) {
  nrow=dim(mat)[1]
  rzcVec <- as(array(data=0,dim=nrow),"numeric")
  for (i in 1:nrow) {
    rzcVec[i] = nnzero(mat[i,])
  }
  return(rzcVec)
}

greedGraphElimination <- function(V,E,w) {
  
}

treeDist_u_to_v <- function(vDist,vPar,u,v,root,upath=c(),vpath=c(),
                            giveupath=FALSE) {
  upath = constructPath(root,u,vPar)
  vpath = constructPath(root,v,vPar)
  uvMatchList=match(upath,vpath,nomatch=0)
  uvlink=vpath[uvMatchList[which(uvMatchList>0)[1]]]
  return(vDist[u]+vDist[v]-2*vDist[uvlink])
}
totalTreeStretch <- function(G,Tree,dist,par,gEdata,root) {
  stretchIdList = which(!(E(as.undirected(G))%in%E(as.undirected(Tree))))
  sum=0
  for (edge in stretchIdList) {
    sum=sum+treeDist_u_to_v(dist,par,gEdata[edge,1],gEdata[edge,2],root)
  }
  return(sum)
}
treeStretchList <- function(G,Tree,dist,par,gEdata,root) {
  Gdir<-as.directed(G);Treedir<-as.directed(Tree)
  stretchIdList = which(!(E(Gdir)%in%E(Treedir)))
  stretchLenList = as(array(data=0,dim=length(stretchIdList)),"numeric")
  for (ii in 1:length(stretchIdList)) {
    edge=stretchIdList[ii]
    stretchLenList[ii] = treeDist_u_to_v(dist,par,
                            gEdata[edge,1],gEdata[edge,2],root)
  }
  return(rbind(stretchIdList,stretchLenList))
}

genEidMap <- function(Edata,Eids) {
  sparseMatrix(i=Edata[,1],j=Edata[,2],x=Eids)
}
distMapToTree <- function(order,par,EIdMap=c(),biDirected=FALSE) {
  if (biDirected) {
    gTreeEdata <- rbind(c(order[-1],(par[order])[-1]),
                        c((par[order])[-1],order[-1]))
    gTree <- graph(gTreeEdata,directed=TRUE)
    if (nnzero(EIdMap)>=length(E(gTree))) {
      print("mapping Eids")
      E(gTree)$id <- EIdMap[t(gTreeEdata)]
    }
  } else {
    gTreeEdata <- rbind(order,par[order])[,-1]
    gTree <- graph(gTreeEdata,directed=FALSE)
    if (nnzero(EIdMap)>=length(E(gTree))) {
      E(gTree)$id <- EIdMap[t(gTreeEdata)]
    }
  }
  return(gTree)
}
getSubgraphEdgeIds <- function(G,S) {
  sElist <- get.edgelist(S)
  sEmatchList <- array(data=0,dim=length(E(S))*2)
  sEmatchList[2*(1:length(E(S)))-1]<-sElist[,1]
  sEmatchList[2*(1:length(E(S)))]<-sElist[,2]
  return(get.edge.ids(G,sEmatchList))
}
genBallTree <- function(G,Gdata=list(),rv,directed=NULL) {
  if (length(Gdata) < 1) {
    Gdata <- extractGraphData(G,rv=rv)
  }
  gTree<-distMapToTree(Gdata$ball[1,]+1,Gdata$par)
  V(gTree)$id <- V(G)$id
  if (!is.null(directed)) {
    if (directed == TRUE) {
      gTree <- as.directed(gTree)
    } else {
      gTree <- as.undirected(gTree)
    }
  }
  E(gTree)$id <- getSubgraphEdgeIds(G,gTree)
  return(gTree)
}

cLind <- function(i,j,jmax) {
  if (i>j) stop("i must be less than j!")
  if (j<=jmax) {
    return((j-1)*j/2+(i%%j)+1)    
  } else {
    
  }
}
partialChol<-function(A,jmax) {
  
}

lapToGraphData <- function(A,ztol=1e-15,vWeight=FALSE,verbose=FALSE) {
  if (dim(A)[1] != dim(A)[2]) stop("Only square matrices are supported")
  if (dim(A)[1] < 2) stop("minimum matrix size is 2x2")
  if (verbose) print("parsing row 1")
  sIds = which(abs(A[2:dim(A)[1],1]) > ztol)+1
  if (length(sIds > 0)){
    tIds = as(array(data=1,dim=length(sIds)),"numeric")
    weights = A[1,sIds] 
  }
  for (k in 2:(dim(A)[1]-1)) {
    if (verbose) print(paste("parsing row ",k,sep=""))
    rIds = c(which(abs(A[k,1:k-1]) > ztol),which(abs(A[k,(k+1):dim(A)[1]]) > ztol)+k)
    if (length(rIds > 0)) {
      sIds = append(sIds,rIds)
      weights= append(weights,A[k,rIds])
      tIds = append(tIds,as(array(data=k,dim=length(rIds)),"numeric")) 
    }
  }
  if (verbose) print(paste("parsing row ",k,sep=""))
  rIds = which(abs(A[dim(A)[1],1:(dim(A)[1]-1)]) > ztol)
  if (length(rIds > 0)) {
    sIds = append(sIds,rIds)
    weights= append(weights,A[dim(A)[1],rIds])
    tIds = append(tIds,as(array(data=dim(A)[1],dim=length(rIds)),"numeric")) 
  }
  return(rbind(rbind(sIds,tIds),abs(weights)))
}
gLapToTreeLap <- function(A,ztol=1e-15,root,offEdgeList=FALSE,
                          vWeight=FALSE,verbose=FALSE) {
  if (dim(A)[1] != dim(A)[2]) stop("Only square matrices are supported")
  if (dim(A)[1] < 2) stop("minimum matrix size is 2x2")
  aEdata = lapToGraphData(A,ztol)
  aDistMap = growBall(c(1:dim(A)[1]),aEdata[1,],aEdata[2,],abs(aEdata[3,]),rv=root,r=999999,
                      distances=TRUE,parents=TRUE)

  iList = c((1+aDistMap[1,][-1]),(1+aDistMap[3,][-1]))
  jList = c((1+aDistMap[3,][-1]),(1+aDistMap[1,][-1]))
  xList = as(array(data=0,dim=length(iList)),"numeric")
  for (i in 1:length(iList)) {
    xList[i] = A[iList[i],jList[i]]
  }
  outmat=as(sparseMatrix(i=iList,j=jList,x=xList,dims=dim(A)),"dgCMatrix")
  for (j in 1:(dim(outmat)[2]-1)) {
    outmat[j,j] = -sum(outmat[which(abs(outmat[1:(j-1),j]) > ztol ),j])
    outmat[j,j] = outmat[j,j] - sum(outmat[j+which(abs(outmat[(j+1):(dim(outmat)[2]),j]) > ztol ),j])
  }
  j=dim(outmat)[2];outmat[j,j] = -sum(outmat[which(abs(outmat[1:(j-1),j]) > ztol ),j])
  if (offEdgeList) {
    outlist=list()
    outlist[["mat"]] = outmat
    Apar <- 
    for (ii in 1:dim(A)[1]) {
      
    }
  } else {
    return(outmat)
  }
}

getOffTreeEdgeIds <- function(GEids,GtreeEids) {
  return(which(!(GEids%in%GtreeEids)))
}
getTreeStretchData <- function(G,Gtree,GtreeDist) {
  
}
sampleOffTreeEdges <- function(A,Atree,AEdat,offEdgeData,sampleSize) {
  
}

cgSol <- function(A,b,x=c(),tol=1e-6,imax=10000,verbose=FALSE,
                  cWarn=TRUE) {
  l=0
  b = array(data=as(b,"numeric"),dim=c(length(b),1))
  if (length(x)<1) {
    x=array(data=0,dim=c(length(b),1))
    r = b - A%*%x
  } else {
    x = array(data=as(x,"numeric"),dim=(c(length(x),1)) )
  }
  r0= b - A%*%array(data=0,dim=(c(length(x),1)))
  r= b - A%*%x
  if (sum(abs(b)) == 0) {
    return(b)
  } else {
    norm0 = sum(abs(r0))
  }
  p=r
  norm = c(sum(abs(r)))
  if (norm0 == 0 || norm[1] == 0) {
    return(r0)
    stop()
  }
  it=1
  converged=FALSE
  while (!converged && it <= imax) {
    Ap = A%*%p
    alpha = t(r)%*%p / (t(p)%*%Ap)
    #print(paste("iteration: ",it," dim alpha: ",dim(alpha)," dim p: ",dim(p)))
    #print(paste("alpha: ",alpha))
    #print(paste("p: ",p))
    x = x + t(alpha * t(p))
    r0=r
    r = r - t(alpha*t(Ap))
    beta = t(r)%*%r / (t(r0)%*%r0)
    p = r + t(beta*t(p))
    it = it +1
    norm=append(norm,sum(abs(r)))
    if (norm[it] / norm0 <= tol) converged = TRUE
  }
  if (!converged && cWarn) warning("Did not converge!")
  if (verbose) print(paste("Required ",it," iterations",sep=""))
  #if (verbose) print(norm/norm0)
  return(x)
}

tccgSol <- function(A,Atree,b,x=c(),tol=1e-6,imax=1000,verbose=FALSE,f=1,
                    MP=c(),pAtreeR=c(),itimes=0) {
  
  #Mchol <- chol(Atree + 1.0e-12 * diag(min(dim(Atree))))
  #Minv <- chol2inv(Mchol)
  #create permutation and cholesky decomposition of the tree if bad or not provided
  if (verbose && itimes > 0) {
    ttemp <- proc.time()
  }
  if ( is.null(MP) || is.null(pAtreeR)) {
    MP = as(genSRPmat(order(rowNZcountVec(Atree)),dim(Atree)[1],dim(Atree)[2]),"dgCMatrix")
    pAtree = MP%*%(Atree+1e-10*diag(dim(A)[1]))%*%t(MP)
    pAtreeR = chol(pAtree)     
  } else {
    if (sum(as(c(dim(MP)!=dim(A),dim(pAtreeR)!=dim(A)),"integer")) != 0) {
      MP = as(genSRPmat(order(rowNZcountVec(Atree)),dim(Atree)[1],dim(Atree)[2]),"dgCMatrix")
      pAtree = MP%*%(Atree+1e-10*diag(dim(A)[1]))%*%t(MP)
      pAtreeR = chol(pAtree) 
    }    
  }
  
  #Minv = diag(min(dim(Minv)))
  
  b = as(array(data=b,dim=c(length(b),1)),"numeric")
  if (length(x)<1) {
    x=as(array(data=0,dim=c(length(b),1)),"numeric")
  } else {
    x = as(array(data=x,dim(c(length(x),1))),"numeric")
  }
  
  r = b - A %*% x
  #z = Minv %*% r
  z = f*t(MP)%*%backsolve(pAtreeR,forwardsolve(t(pAtreeR),MP%*%r))
  p=z
  
  norm0 = sum(abs(r))
  norm= c()
  it=1
  converged=FALSE
  if (verbose && itimes > 0) {
    print(paste("setup time: ",(proc.time()-ttemp)[3],sep=""))
  }
  while (!converged && it <= imax) {
    if (verbose && it < itimes){
      ttemp<-proc.time()
    }
    alpha = t(r)%*%z / (t(p)%*%A%*%p)
    x = x + t(alpha*t(p))
    r0 = r
    r= r - t(alpha*t(A%*%p))
    norm = append(norm,sum(abs(r)))
    if (norm[it]/norm0 < tol) {
      converged = TRUE
    } else {
      if (verbose && it < itimes) {
        ttemp=proc.time()
      }
      z0=z
      z = f*t(MP)%*%backsolve(pAtreeR,forwardsolve(t(pAtreeR),MP%*%r))
      beta = t(z)%*%(r-r0) / (t(z0)%*%r0)
      p = z + t(beta * t(p))
      if(verbose && it < itimes) {
        print(paste("iteration ",it," time: ",(proc.time()-ttemp)[3],sep=""))
      }
      it = it +1
    }
  }
  if (!converged) warning("did not converge!")
  if (verbose) print(paste("iterations: ",it))
  #if (verbose) print(norm/norm0)
  return(x)
}

fastTreeDist <- function(vert1,vert2,root,par,distMap) {
  path1 = findTreePath(vert1-1,root-1,par-1,length(par))+1
  path2 = findTreePath(vert2-1,root-1,par-1,length(par))+1
  matchVert=path1[which(!(match(path1,path2,nomatch=0)%in%0))[1]]
  return(distMap[vert1]+distMap[vert2]-2*distMap[matchVert])
}

buildGraphConditioners <- function(G,rv,Gtree=NULL,Gdata=list(),
                                   Gtreedata=list()) {
  if (is.null(Gdata$par) || is.null(Gdata.ball)
      || is.null(Gdata.dist) || is.null(Gdata.Edata) ||
        !(is.directed(G))) {
    if (!is.directed(G)) {
      G<-as.directed(G)
      E(G)$id <- 1:length(E(G))
    }
    Gdata <- extractGraphData(G,rv=rv) 
  }
  if (is.null(Gtree)) {
    #Special tree building routines can go here for now, just
    #using the search / ball grow tree method
    Gtree = genBallTree(G,Gdata,rv,directed=TRUE)
    E(Gtree)$weight <- E(G)$weight[E(Gtree)$id]
  }
  if (is.null(E(Gtree)$weight)) {
    E(Gtree)$weight <- E(G)$weight[E(Gtree)$id]
  }
  GtreeEdata <- t(rbind(t(get.edgelist(Gtree)),E(Gtree)$weight))
  if (is.null(Gtreedata$par) || is.null(Gtreedata.ball)
      || is.null(Gtreedata.dist) || is.null(Gtreedata.Edata) ) {
    Gtreedata <- extractGraphData(as.directed(Gtree),rv=rv) 
  }
  offTreeEdgeInds = which(!(E(G)$id %in% E(Gtree)$id))
  offTreeDists = mapply(fastTreeDist, Gdata$Edata[offTreeEdgeInds,1], 
                        Gdata$Edata[offTreeEdgeInds,2],
                            MoreArgs=list(root=rv,par=Gtreedata$par,
                                               distMap=Gdata$dist))
  probNorm = sum(E(G)$weight[offTreeEdgeInds] / offTreeDists)
  offEdgeProbs = E(G)$weight[offTreeEdgeInds] / offTreeDists / probNorm
  indMat = sparseMatrix(i=Gdata$Edata[,1],j=Gdata$Edata[,2],x=1:length(E(G)))
  Cinds = mapply(function(ii) indMat[Gdata$Edata[ii,2],Gdata$Edata[ii,1]],
                 offTreeEdgeInds)
  offTreeData = as.data.table(list(Einds=offTreeEdgeInds,
                        Cinds=Cinds,
                        Dinds=Gdata$Edata[offTreeEdgeInds,1],
                        probs = offEdgeProbs,
                       counts=0*offTreeEdgeInds))
  offDiagData = as.data.table(list(Dinds=unique(Gdata$Edata[offTreeEdgeInds,1]),
                                    counts=0*unique(Gdata$Edata[offTreeEdgeInds,1])))
  setkey(offTreeData,Einds)
  setkey(offDiagData,Dinds)
  diagInds = t(rbind(offDiagData$Dinds,offDiagData$Dinds))
  offInds = Gdata$Edata[offTreeData$Einds,1:2]
  Plist=list()
  Plist[[1]] <- edgeDataToLap(GtreeEdata,V(Gtree)$weight)
  treeLap <- edgeDataToLap(Gtreedata$Edata,V(Gtree)$weight)
  tempPbase <- edgeDataToLap(Gdata$Edata,V(G)$weight)
  tempPbase[Gdata$Edata[offTreeEdgeInds,1:2]] = 0
  tempPbase[Gtreedata$Edata[,1:2]] = treeLap[Gtreedata$Edata[,1:2]]
  diag(tempPbase) <- diag(treeLap)
  tempP=tempPbase
  logN = floor(log(length(V(G))))
  eM = length(E(G))
  Mdiff = length(offTreeDists)
  if (logN > 1) {
    #for (ii in 1:(logN^2)) {
      tempP = tempPbase * logN
      addCounts = count(sample(offTreeData$Einds,size = eM/logN^2,replace=TRUE,
                        prob=offTreeData$probs))
      offTreeData[Einds %in% addCounts$x]$counts = offTreeData[Einds %in% addCounts$x]$counts +
                                            addCounts$freq
      offTreeData$counts = offTreeData$counts + 
              offTreeData$counts[match(offTreeData$Cinds,offTreeData$Einds)]
      diagCounts = aggregate(offTreeData$counts,by=list(Dinds=offTreeData$Dinds),sum)
      offDiagData[Dinds %in% diagCounts$Dinds]$counts = offDiagData[Dinds %in% diagCounts$Dinds]$counts +
          diagCounts$x
      tempP[offInds] = tempP[offInds] - offTreeData$counts * 
          Gdata$Edata[offTreeData$Einds,3]
      tempP[diagInds] = (tempP[diagInds] +
          offDiagData$counts * Gdata$Edata[offDiagData$Dinds,3])
    #} 
  }
  return(list(mat=tempP,tlap=treeLap,offInds=offInds,
              diagInds=diagInds,offData=offTreeData,diagData=offDiagData,
              Esample=addCounts))
}

spannerListPvec <- function(matList) {
  len = length(matList)
  size = dim(matList[[1]])[1]
  cVecList <- lapply(matList,function(mat) rowNZcountVec(mat))
  mInd <- 1
  pVec<-order(cVecList[[len]])
  matchVec <- rep(1,len)
  for (i in (len-1):1) {
    tempVec <- cVecList[[i]][pVec]
    tempVec <- tempVec[mInd:size]
    toVec<-order(tempVec)
    pVec[mInd:size] <- pVec[mInd:size][toVec]
    mVal <- match(TRUE,(4 <= tempVec[toVec]),nomatch=0)
    if (mVal == 0) {
      matchVec[i]<-len
      mInd<-len
    } else {
      mInd<-mInd+mVal-1
      matchVec[i]<-mInd
    }
  }
  return(list(pVec=pVec,mInds=matchVec))
}
decomposeSpanner <- function(mat) {
  rnzVec <- rowNZcountVec(mat)
  permVec <- order(rnzVec)
  dVal <- match(4,rnzVec[permVec],nomatch=0)-1
  if (dVal > 1) {
    MP <- as(genSRPmat(permVec,dim(mat)[1],dim(mat)[2]),"dgCMatrix")
    pmat <- MP%*%mat%*%t(MP)
    lList <- matLL(pmat,dVal)
    return(list(P=MP,L=lList$L,C=lList$C,d=dVal))
  }
}
decomposeSpannerList <- function(matList) {
  len <- length(matList)
  
}

graphDataToConditioners <- function(Gdata) {
  Alap <- edgeDataToLap(Gdata$Edata)
  treeEdata <- rbind(c((Gdata$ball[1,]+1)[-1],(Gdata$ball[3,]+1)[-1]),
                     c((Gdata$ball[3,]+1)[-1],(Gdata$ball[1,]+1)[-1]) )
  treeEdata <- t(rbind(treeEdata,-Alap[t(treeEdata)]))
  Atree <- edgeDataToLap(treeEdata)
  offEdgeData <- subset(subset(summary(Alap-Atree),i!=j),x!=0)
  Adata <- summary(Alap)
  Adata$id <- 1:nrow(Adata)
  EidMat <- sparseMatrix(i=Adata$i,j=Adata$j,x=Adata$id)
  offEdgeData$id <- mapply(function(i,j) EidMat[i,j],offEdgeData$i,offEdgeData$j)
  offEdgeData$tdists <- mapply(fastTreeDist, offEdgeData$i, offEdgeData$j,
                        MoreArgs=list(root=Gdata$ball[1,1]+1,par=Gdata$par,
                                      distMap=Gdata$dist))
  probNorm <- sum(-1 * offEdgeData$x / offEdgeData$tdists)
  offEdgeData$prob <- 1 / offEdgeData$tdists / probNorm
  loffEdges<-subset(offEdgeData,i>j)
  #sampling portion
  matList <- list();
  matN <- length(diag(Alap))
  LogN <- floor(log(matN))
  baseMat<- Alap * 0 + Atree
  Adata[,5]<- 0
  sfac=LogN^2
  matList[[1]] <- baseMat * sfac
  for (ii in 1:(2*LogN)){ #need to rework this to only add to empty off tree edges
    idSample <- count(sample(x=loffEdges$id,prob=loffEdges$prob,replace=TRUE,
                           size=LogN^2))
    Adata[,ii+5] <- Adata[,ii+4]
    Adata[idSample$x,ii+5] <- (abs(  (abs(Adata[idSample$x,ii+5]) >0 ) + 
                                     (abs(idSample$freq) >0 ) ) >0 ) * 
                              Adata$x[idSample$x]
    cids <- EidMat[t(rbind(Adata[idSample$x,"j"],Adata[idSample$x,"i"]))]
    Adata[cids,ii+5] <- Adata[idSample$x,ii+5]
    dcounts <-  count(data.frame(i=Adata$i[c(idSample$x,cids)],
                                c=Adata[c(idSample$x,cids),ii+5]
                                ),vars="i",wt_var="c")
    Adata[EidMat[t(rbind(dcounts$i,dcounts$i))],ii+5] <- -1*dcounts$freq
    tempMat <- baseMat * sfac
    tempMat[t(rbind(Adata$i,Adata$j))] =  tempMat[t(rbind(Adata$i,Adata$j))] +
                                          Adata[,ii+5]
    matList[[ii+1]] = tempMat
  }
  matList[[2*LogN+1]] = baseMat*(sfac - 1) + Alap
  pList <- list()
  lList <- list()
  cList <- list()
  if (TRUE){
    MPdata <- spannerListPvec(matList)
    MP <- genSRPmat(MPdata$pVec,matN,matN)
    lList <- list()
    cList <- list()
    for (ii in 1:(length(matList)-1)) {
      #print(MPdata$mInds[[ii]]-1)
      lcTemp <- matLL(MP%*%matList[[ii]]%*%t(MP),MPdata$mInds[ii]-1)
      lList[[ii]] <- lcTemp$L
      cList[[ii]] <- lcTemp$C
    }
  }
  liList <- lapply(lList,function(mat) as(solve(mat),"dgCMatrix"))
  return(list(offEdgeData=offEdgeData,Adata=Adata,
              #counts=data.frame(i=ivec,j=jvec,c=cvec),dcounts=diagCounts,
              matList=matList,idMat=EidMat,MPdata=MPdata,MP=MP,
              lList=lList,cList=cList,liList=liList))
}

#need to add automated testing 

#clpData is list with
# L, C, MP s.t. conditioner = t(MP).t(L).C.L.MP
# ind : the length of non unity entries in the diagonal of C
clpccgSol <- function(A,clpData,b,x=c(),tol=1e-6,imax=1000,
                      verbose=FALSE,itimes=1) {
  
  #Mchol <- chol(Atree + 1.0e-12 * diag(min(dim(Atree))))
  #Minv <- chol2inv(Mchol)
  #create permutation and cholesky decomposition of the tree if bad or not provided
  if (verbose && itimes > 0) {
    ttemp <- proc.time()
  }

N <- dim(A)[1]
MP <- clpData$MP
L <- clpData$L
C <- clpData$C
Linv <- solve(clpData$L)
Linv <- (abs(Linv) > tol)*Linv
ind <- clpData$ind

cits <- ceiling(2*log(N-ind))

#Minv = diag(min(dim(Minv)))

b = as(array(data=b,dim=c(length(b),1)),"numeric")
if (length(x)<1) {
  x=as(array(data=0,dim=c(length(b),1)),"numeric")
} else {
  x = as(array(data=x,dim(c(length(x),1))),"numeric")
}

r = b - A %*% x
#z = Minv %*% r
z = Linv%*%(MP%*%b)
z[ind:N] <- cgSol(C,z[ind:N],imax=cits,cWarn=FALSE)
z = 
p=z

norm0 = sum(abs(r))
norm= c()
it=1
converged=FALSE
if (verbose && itimes > 0) {
  print(paste("setup time: ",(proc.time()-ttemp)[3],sep=""))
}
while (!converged && it <= imax) {
  if (verbose && it < itimes){
    ttemp<-proc.time()
  }
  alpha = t(r)%*%z / (t(p)%*%A%*%p)
  x = x + t(alpha*t(p))
  r0 = r
  r= r - t(alpha*t(A%*%p))
  norm = append(norm,sum(abs(r)))
  if (norm[it]/norm0 < tol) {
    converged = TRUE
  } else {
    if (verbose && it < itimes) {
      ttemp=proc.time()
    }
    z0=z
    z = f*t(MP)%*%backsolve(pAtreeR,forwardsolve(t(pAtreeR),MP%*%r))
    beta = t(z)%*%(r-r0) / (t(z0)%*%r0)
    p = z + t(beta * t(p))
    if(verbose && it < itimes) {
      print(paste("iteration ",it," time: ",(proc.time()-ttemp)[3],sep=""))
    }
    it = it +1
  }
}
if (!converged) warning("did not converge!")
if (verbose) print(paste("iterations: ",it))
#if (verbose) print(norm/norm0)
return(x)
}

hccgSol <- function(A,Atree,b,x=c(),tol=1e-6,imax=1000,verbose=FALSE,f=1,
                    MPlist=list(),pAtreeRlist=list(),level=0) {
      
  #build hierachy of preconditioners if not given
  if (length(MPlist)==0 || length(pAtreeRlist==0)) {
    #create permutation and cholesky decomposition of the tree
    MP = as(genSRPmat(order(rowNZcountVec(Atree)),dim(Atree)[1],dim(Atree)[2]),"dgCMatrix")
    pAtree = MP%*%(Atree+1e-10*diag(dim(A)[1]))%*%t(MP)
    pAtreeR = chol(pAtree)     
    
    #create hierachy of preconditioners MPlist and pAtreeRlist with tree as level 1
    MPlist = list(); pAtreeRlist = list();
    MPlist[[1]] = MP; pAtreeRlist[[1]] = pAtreeR 
  }
  
  b = as(array(data=b,dim=c(length(b),1)),"numeric")
  if (length(x)<1) {
    x=as(array(data=0,dim=c(length(b),1)),"numeric")
  } else {
    x = as(array(data=x,dim(c(length(x),1))),"numeric")
  }
  
  r = b - A %*% x
  #z = Minv %*% r
  z = f*t(MP)%*%backsolve(pAtreeR,forwardsolve(t(pAtreeR),MP%*%r))
  p=z
  
  norm0 = sum(r^2)
  norm= c()
  it=0
  converged=FALSE
  while (!converged && it <= imax) {
    alpha = t(r)%*%z / (t(p)%*%A%*%p)
    x = x + t(alpha*t(p))
    r0 = r
    r= r - t(alpha*t(A%*%p))
    norm = append(norm,sum(r^2))
    if (norm[it+1]/norm0 < tol) {
      converged = TRUE
    } else {
      z0=z
      z = f*t(MP)%*%backsolve(pAtreeR,forwardsolve(t(pAtreeR),MP%*%r))
      beta = t(z)%*%(r-r0) / (t(z0)%*%r0)
      p = z + t(beta * t(p))
      it = it + 1
    }
  }
  if (!converged) warning("did not converge!")
  if (verbose) print(paste("iterations: ",it))
  #if (verbose) print(norm/norm0)
  return(x)
}

#assumes a symmetric diagonally dominant matrix
#this implies that we can bound the eigenvalues by 0 at minimum
#and by twice the largest diagonal entry at maximum
tcChebySol <- function(A,Atree,b,x0=c(),tol=1.0e-10,imax=1000,verbose=FALSE,
                       MP=c(),pAtreeR=c(),f=1) {
  lmin=0;lmax=2*max(diag(A))
  d = (lmax + lmin)/2
  c = (lmax - lmin)/2
  
  if (dim(A)[1] != dim(A)[2]) stop("Only square matrices supported!")
  if (dim(Atree)[1] != dim(Atree)[2]) stop("Only square matrices supported!")
  
  if ((length(MP) != length(A)) || (length(pAtreeR) != length(A))) {
    #print("building tree R'R")
    MP = as(genSRPmat(order(rowNZcountVec(Atree)),dim(Atree)[1],dim(Atree)[2]),"dgCMatrix")
    pAtree = MP%*%(Atree+1e-10*diag(dim(A)[1]))%*%t(MP)
    pAtreeR = chol(pAtree) 
  }
  #print("loading x and b")
  it=0
  b = as(array(data=b,dim=c(length(b),1)),"numeric")
  if (length(x0)!=length(bvec)) {
    x=as(array(data=0,dim=c(length(b),1)),"numeric")
  } else {
    #print(paste("length of x0:",length(x0)))
    x = as(array(data=x0,dim=c(length(x0),1)),"numeric")
  }
  r= b - A%*%as(array(data=0,dim=c(length(b),1)),"numeric")
  if (sum(abs(b)) == 0) {
    return(b)
  } else {
    norm0 = sum(abs(r))
  }
  it=1
  r=b-A%*%x
  z = f*t(MP)%*%backsolve(pAtreeR,forwardsolve(t(pAtreeR),MP%*%r))
  p=z
  alpha = 2 / d
  x = x + t(alpha*t(p))
  r = b - A%*%x
  norm=c(sum(abs(r)))
  converged=FALSE
  if (norm[it] / norm0 <= tol) converged=TRUE
  it = it+1
  #print("starting main loop")
  if (imax >= 2) {
    while (it <= imax && !converged) {
      z = f*t(MP)%*%backsolve(pAtreeR,forwardsolve(t(pAtreeR),MP%*%r))
      beta = (c*alpha/2)^2
      alpha = 1 / (d-beta)
      p=z+t(beta*t(p))
      x=x+t(alpha*t(p))
      r = b - A%*%x
      norm=append(norm,sum(abs(r)))
      #print(paste(it,norm[it]/norm0))
      if ( (norm[it] / norm0) <= tol) converged=TRUE
      it=it+1
    }
  }
  if (!converged) warning("Did not converge!")
  if (verbose) print(t(norm/norm0))
  print(paste(it," iterations used",sep=""))
  return(x)
}

scChebySol <- function(A,Atree,b,x0=c(),tol=1.0e-10,imax=1000,verbose=FALSE) {
  lmin=0;lmax=2*max(diag(A))
  d = (lmax + lmin)/2
  c = (lmax - lmin)/2
  
  if (dim(A)[1] != dim(A)[2]) stop("Only square matrices supported!")
  if (dim(Atree)[1] != dim(A)[2]) stop("Only square matrices supported!")
  
  print("setting up tree and spine")
  shA = A+Atree
  
  MP = as(genSRPmat(order(rowNZcountVec(Atree)),dim(Atree)[1],dim(Atree)[2]),"dgCMatrix")
  pAtree = MP%*%(2*Atree+(1e-2)*tol*diag(dim(A)[1]))%*%t(MP)
  pAtreeR = chol(pAtree)
  
  it=0
  b = as(array(data=b,dim=c(length(b),1)),"numeric")
  if (length(x0)!=length(bvec)) {
    x=as(array(data=0,dim=c(length(b),1)),"numeric")
  } else {
    x = as(array(data=x0,dim=c(length(x0),1)),"numeric")
  }
  r= b - A%*%as(array(data=0,dim=c(length(b),1)),"numeric")
  if (sum(abs(b)) == 0) {
    return(b)
  } else {
    norm0 = sum(abs(r))
  }
  it=1
  r=b-A%*%x
  print("starting first backsolve"); flush.console()
  z = 2*tcChebySol(shA,Atree*2,b=r,x,MP=MP,pAtreeR=pAtreeR)
  p=z
  alpha = 2 / d
  x = x + t(alpha*t(p))
  r = b - A%*%x
  norm=c(sum(abs(r)))
  converged=FALSE
  print("starting solver loop")
  if (norm[it] / norm0 <= tol) converged=TRUE
  it = it+1
  if (imax >= 2) {
    while (it <= imax && !converged) {
      z = 2*tcChebySol(shA,2*Atree,b=r,x,MP=MP,pAtreeR=pAtreeR,tol=tol/norm[it-1])
      beta = (c*alpha/2)^2
      alpha = 1 / (d-beta)
      p=z+t(beta*t(p))
      x=x+t(alpha*t(p))
      r = b - A%*%x
      norm=append(norm,sum(abs(r)))
      print(paste(it,norm[it]/norm0))
      if ( (norm[it] / norm0) <= tol) converged=TRUE
      it=it+1
    }
  }
  if (!converged) warning("Did not converge!")
  if (verbose) print(t(norm/norm0))
  return(x)
}

cmcgSol <- function(Cmat,Ccond,b,x=c(),tol=1e-6,imax=1000,verbose=FALSE,itimes=1) {
  
  b = as(array(data=b,dim=c(length(b),1)),"numeric")
  if (length(x)<1) {
    x=as(array(data=0,dim=c(length(b),1)),"numeric")
  } else {
    x = as(array(data=x,dim(c(length(x),1))),"numeric")
  }
  
  r = b - Cmat %*% x
  #z = Minv %*% r
  z = cgSol(Ccond,r,x,tol,imax,verbose=FALSE)
  p=z
  
  ttemp<-proc.time()
  
  norm0 = sum(abs(r))
  norm= c()
  it=1
  converged=FALSE
  if (verbose && itimes > 0) {
    print(paste("setup time: ",(proc.time()-ttemp)[3],sep=""))
  }
  while (!converged && it <= imax) {
    if (verbose && it < itimes){
      ttemp<-proc.time()
    }
    alpha = t(r)%*%z / (t(p)%*%Cmat%*%p)
    x = x + t(alpha*t(p))
    r0 = r
    r= r - t(alpha*t(Cmat%*%p))
    norm = append(norm,sum(abs(r)))
    if (norm[it]/norm0 < tol) {
      converged = TRUE
    } else {
      if (verbose && it < itimes) {
        ttemp=proc.time()
      }
      z0=z
      z = cgSol(Ccond,r,p,tol,imax,verbose)
      beta = t(z)%*%(r-r0) / (t(z0)%*%r0)
      p = z + t(beta * t(p))
      if(verbose && it < itimes) {
        print(paste("iteration ",it," time: ",(proc.time()-ttemp)[3],sep=""))
      }
      it = it +1
    }
  }
  if (!converged) warning("did not converge!")
  if (verbose) print(paste("iterations: ",it))
  #if (verbose) print(norm/norm0)
  return(x)
  
}

clpcgSol <- function(A,cmat,limat,pmat,cval=1,b,x=c(),tol=1e-6,imax=1000,verbose=FALSE,itimes=1,
                     cWarn=TRUE) {
  b =array(data=as(b,"numeric") ,dim=c(length(b),1))
  if (length(x)<1) {
    x=array(data=0,dim=c(length(b),1))
  } else {
    x =array(data=as(x,"numeric"),dim(c(length(x),1)))
  }
  nr = dim(A)[1]
  A=A+diag(nr)*tol/log(nr)
  cits = ceiling((nr-cval)+1)
  #print(cits)
  r = b - A %*% x
  #z = Minv %*% r

  z = (limat)%*%(pmat)%*%r
  z[cval:nr] = cgSol(cmat[cval:nr,cval:nr]+diag(nr-cval+1)*tol/log(nr-cval+1),
                     b=array(z[cval:nr],dim=c(nr-cval+1,1)),
                     x=array(z[cval:nr],dim=c(nr-cval+1,1)),
                     tol,imax=cits,verbose=FALSE,cWarn=FALSE)
  z = t(pmat)%*%t(limat)%*%z
  p=z
  
  ttemp<-proc.time()
  norm0 = sum(abs(r))
  if (norm0 == 0) {
    return(r)
    stop()
  }
  norm= c(norm0)
  ccount=1
  it=1
  converged=FALSE
  if (verbose && itimes > 0) {
    #print(paste("setup time: ",(proc.time()-ttemp)[3],sep=""))
  }
  while (!converged && it <= imax) {
    if (verbose && it < itimes){
      ttemp<-proc.time()
    }
    it = it +1
    alpha = as((t(r)%*%z / (t(p)%*%A%*%p)),"numeric")
    x = x + t(alpha*t(p))
    r0 = r
    r= r - t(alpha*t(A%*%p))
    norm = append(norm,sum(abs(r)))
    if (is.nan(norm[it]) ) {
      print(t(z))
      stop("ERROR!, norm undefined!")
    }
    if (norm[it]/norm0 < tol) {
      converged = TRUE
      print(paste("converge in",it-1,"iterations"))
    } else {
      if (verbose && it < itimes) {
        ttemp=proc.time()
      }
      z0=z
      if (norm[it]/norm0 > norm[it-1]/norm0) {
        cits = max(1,cits-2)        
        #ccount=1
      }
      else {
        cits = max(ceiling(((nr-cval+1) - cits)/2)+cits,
                   (nr-cval+1))
        #ccount=ccount+1
      }
      z = (limat)%*%(pmat)%*%r
      z[cval:nr] = cgSol(cmat[cval:nr,cval:nr]+diag(nr-cval+1)*tol/log(nr-cval+1),
                         b=array(z[cval:nr],dim=c(nr-cval+1,1)),
                         x=array(z[cval:nr],dim=c(nr-cval+1,1)),
                         tol,imax=cits,verbose=FALSE,cWarn=FALSE)
      z = t(pmat)%*%t(limat)%*%z
      
      beta = t(z)%*%(r-r0) / (t(z0)%*%r0)
      p = z + t(beta * t(p))

      if(verbose && it < itimes) {
        print(paste("iteration ",it-1," time: ",(proc.time()-ttemp)[3],sep=""))
      }

    }
  }
  if (!converged && cWarn) warning("did not converge!")
  if (verbose) print(paste("iterations: ",it,sep=""))
  #if (verbose) print(norm/norm0)
  return(x)
  
}

rccgSol <- function(A,b,x=c(),tol=1e-6,imax=1000,verbose=FALSE,f=1,
                    cDataList) {
  
  #cDataList: list of L and associated C matrices.
  #   If the list is null, then standard CG is run on the
  #   input matrix.
  #   If C is not a simple diagonal matrix then cutVal
  #   should also be present to indicate where to cut off
  #   the elimination backsolve and start next level
  #   May also contain P (a permutation matrix) which will
  #   be applied if present.
  
  if (is.null(cDataList)) {
    #run cg
  } else {
    if (is.null(cDataList[[length(cDataList)]]$P)) { 
      
    } else {
      #Apply permutations
      
    }
  }
  
  b = as(array(data=b,dim=c(length(b),1)),"numeric")
  if (length(x)<1) {
    x=as(array(data=0,dim=c(length(b),1)),"numeric")
  } else {
    x = as(array(data=x,dim(c(length(x),1))),"numeric")
  }
  
  r = b - A %*% x
  #z = Minv %*% r
  z = f*t(MP)%*%backsolve(pAtreeR,forwardsolve(t(pAtreeR),MP%*%r))
  p=z
  
  norm0 = sum(r^2)
  norm= c()
  it=0
  converged=FALSE
  while (!converged && it <= imax) {
    alpha = t(r)%*%z / (t(p)%*%A%*%p)
    x = x + t(alpha*t(p))
    r0 = r
    r= r - t(alpha*t(A%*%p))
    norm = append(norm,sum(r^2))
    if (norm[it+1]/norm0 < tol) {
      converged = TRUE
    } else {
      z0=z
      z = f*t(MP)%*%backsolve(pAtreeR,forwardsolve(t(pAtreeR),MP%*%r))
      beta = t(z)%*%(r-r0) / (t(z0)%*%r0)
      p = z + t(beta * t(p))
      it = it + 1
    }
  }
  if (!converged) warning("did not converge!")
  if (verbose) print(paste("iterations: ",it))
  #if (verbose) print(norm/norm0)
  return(x)
}