#GMG solver
require(MASS)
require(Matrix)
require(plyr)

build2DIndexGrids <- function(fineDims) {
  levels = ceiling(log(min(fineDims)) / (log(2)) - 2)
  lDims = fineDims
  GeoList = list()
  VecList = list()
  DimList = list()
  if (levels > 0 ) {
    for (level in 1:levels) {
      DimList[[level]]=lDims
      mapGeoToVec = array(1:(lDims[1]*lDims[2]),dim=lDims)
      mapVecToGeo = array(0,dim=c(lDims[1]*lDims[2],2))
      for(ii in 1:lDims[1]) {
        for (jj in 1:lDims[2]) {
          mapVecToGeo[mapGeoToVec[ii,jj],] = c(ii,jj)
        }
      }
      GeoList[[level]] = mapGeoToVec; VecList[[level]] = mapVecToGeo
      lDims = floor(lDims / 2) 
    }
  } else {
    stop("Error: Grid dimesnions are too small!")
  }
  return(list(geoList=GeoList,vecList=VecList,dimList=DimList))
}
build3DIndexGrids <- function(fineDims) {
  levels = ceiling(log(min(fineDims)) / (log(2)) - 2)
  lDims = fineDims
  GeoList = list()
  VecList = list()
  DimList = list()
  if (levels > 1 ) {
    for (level in 1:levels) {
      DimList[[level]]=lDims
      mapGeoToVec = array(1:(lDims[1]*lDims[2]*lDims[3]),dim=lDims)
      mapVecToGeo = array(0,dim=c(lDims[1]*lDims[2]*lDims[3],3))
      for(ii in 1:lDims[1]) {
        for (jj in 1:lDims[2]) {
          for (kk in 1:lDims[3]) {
            mapVecToGeo[mapGeoToVec[ii,jj,kk],] = c(ii,jj,kk)
          }
        }
      }
      GeoList[[level]] = mapGeoToVec; VecList[[level]] = mapVecToGeo
      lDims = floor(lDims / 2) 
    }
  } else {
    stop("Error: Grid dimesnions are too small!")
  }
  return(list(geoList=GeoList,vecList=VecList,dimList=DimList))
}

build2DInterpolate <- function(vecList,geoList,dimList,periodic=FALSE) {
  #This method uses simple linear interpolation. No operator based interpolation
  #variable name notes: i and j will refer to matrix indices, x and y
  #                     will refer to geometric indices.
  #                     lgeo and lvec are fine grid geometric to vector
  #                     index and vector to geometric index mapping arrays
  #                     lpgeo and lpvec are the equivalent coarse grid
  #                     mappings.
  levels = length(vecList)
  if (levels < 2) {
    stop("Error: Must have at least 2 levels")
  }
  lgeo = geoList[[level]]; lpgeo = geoList[[level+1]]
  lvec = vecList[[level]]; lpvec = vecList[[level+1]]
  ldim = dimList[[level]]; lpdim = dimList[[level+1]]
  lgvol = dim(lvec)[1]; lpgvol = dim(lpvec)[1]
  #interp will ultimately become a dgCMatrix class object
  #during construction we will want to store it in terms of the
  #row indices, column indices, and data values needed to construct it
  #its dimensions will ultimately be length(lvec) by length(lpvec)
  interpi = c();interpj = c(); interpx = c() #hold i,j, and x of interp matrix
  #Handle Type1 points first. These are where coarse and fine grid coincide
  #Performed using injection, ie LH=Lh
  interpi = c(interpi,c(1:dim(lpvec)[1]))
  interpj= c(interpj,lgeo[lpvec[interpi,]*2])
  interpx=c(interpx,rep(1,dim(lpvec)[1]))
  #Handle Type2 points next. These lie on grid lines between coarse grids
  #Performed by summing over adjoining coarse grids
  #First instance where periodicity comes into play will branch code here
  #if (periodic) { 
    #to be handled later.
  #} else {
    #for 2D, there will be two line directions, X and Y.
    #we start with terms for fine grids lying upon grid lines
    #start with x i+1/2 terms
    #xCoarse, yCoarse, xFine, and yFine are geometric grid indices
    xCoarse = c(1:floor((ldim[1]-1)/2));yCoarse=c(1:floor(ldim[2]/2))
    interpi = c(interpi,as(array(lpgeo[xCoarse,yCoarse],
                       dim=c(1,length(xCoarse)*length(yCoarse))),"numeric"))
    xFine = 2*xCoarse+1;yFine=2*yCoarse
    interpj = c(interpj,as(array(lgeo[xFine,yFine],
                       dim=c(1,length(xFine)*length(yFine))),"numeric"))
    interpx = c(interpx,rep(.5,length(xFine)*length(yFine)))
    #next x i-1/2 terms
    xCoarse = c(1:floor((ldim[1]-1)/2))
    interpi = c(interpi,as(array(lpgeo[xCoarse,yCoarse],
                                 dim=c(1,length(xCoarse)*length(yCoarse))),"numeric"))
    xFine = 2*xCoarse-1
    interpj = c(interpj,as(array(lgeo[xFine,yFine],
                                 dim=c(1,length(xFine)*length(yFine))),"numeric"))
    interpx = c(interpx,rep(.5,length(xFine)*length(yFine)))
    #now do y j+1/2 terms
    xCoarse = c(1:floor(ldim[1]/2));yCoarse=c(1:floor((ldim[2]-1)/2))
    interpi = c(interpi,as(array(lpgeo[xCoarse,yCoarse],
                                 dim=c(1,length(xCoarse)*length(yCoarse))),"numeric"))
    xFine = 2*xCoarse; yFine=2*yCoarse+1
    interpj = c(interpj,as(array(lgeo[xFine,yFine],
                                 dim=c(1,length(xFine)*length(yFine))),"numeric"))
    interpx = c(interpx,rep(.5,length(xFine)*length(yFine)))
    #now do y j-1/2 terms
    yCoarse = c(1:floor((ldim[1])/2))
    interpi = c(interpi,as(array(lpgeo[xCoarse,yCoarse],
                                 dim=c(1,length(xCoarse)*length(yCoarse))),"numeric"))
    yFine = 2*yCoarse-1
    interpj = c(interpj,as(array(lgeo[xFine,yFine],
                                 dim=c(1,length(xFine)*length(yFine))),"numeric"))
    interpx = c(interpx,rep(.5,length(xFine)*length(yFine)))
    #For 2D the last piece is fine grids not lying on any grid lines. These would be
    #type III and type IV points in 3D. They require summing over all adjacent coarse grids
    #start with x i+1 y j+1
    xCoarse = c(1:floor((ldim[1]-1)/2)); ycoarse = c(1:floor((ldim[2]-1)/2))
    interpi = c(interpi,as(array(lpgeo[xCoarse,yCoarse],
                                 dim=c(1,length(xCoarse)*length(yCoarse))),"numeric"))
    xFine = 2*xCoarse+1; yFine = 2*yCoarse+1
    interpj = c(interpj,as(array(lgeo[xFine,yFine],
                                 dim=c(1,length(xFine)*length(yFine))),"numeric"))
    interpx = c(interpx,rep(.5,length(xFine)*length(yFine)))
    #x i-1, y j+1
    xCoarse = c(1:floor(ldim[1]/2))
    interpi = c(interpi,as(array(lpgeo[xCoarse,yCoarse],
                                 dim=c(1,length(xCoarse)*length(yCoarse))),"numeric"))
    xFine = 2*xCoarse - 1
    interpj = c(interpj,as(array(lgeo[xFine,yFine],
                                 dim=c(1,length(xFine)*length(yFine))),"numeric"))
    interpx = c(interpx,rep(.5,length(xFine)*length(yFine)))
    #x i-1, y j-1
    yCoarse = c(1:floor(ldim[2]/2))
    interpi = c(interpi,as(array(lpgeo[xCoarse,yCoarse],
                                 dim=c(1,length(xCoarse)*length(yCoarse))),"numeric"))
    yFine = 2*yCoarse - 1
    interpj = c(interpj,as(array(lgeo[xFine,yFine],
                                 dim=c(1,length(xFine)*length(yFine))),"numeric"))
    interpx = c(interpx,rep(.5,length(xFine)*length(yFine)))
    #x i+1, y j-1
    xCoarse = c(1:floor((ldim[1]-1)/2))
    interpi = c(interpi,as(array(lpgeo[xCoarse,yCoarse],
                                 dim=c(1,length(xCoarse)*length(yCoarse))),"numeric"))
    xFine = 2*xCoarse + 1
    interpj = c(interpj,as(array(lgeo[xFine,yFine],
                                 dim=c(1,length(xFine)*length(yFine))),"numeric"))
    interpx = c(interpx,rep(.5,length(xFine)*length(yFine)))
  #}
  
}

gridIndGen2D <- function(Iinds,Jinds) {
  ivec <- rep(Iinds,length(Jinds))
  jvec <- mapply(function(val) rep(val,length(Iinds)),Jinds)
  return(rbind(ivec,jvec))
}

oneLevel2Dinterp <- function(lvec,lgeo,ldim,lpvec,lpgeo,lpdim) {
  #Non periodic only for now
  #if (periodic) {} else {}
  lgvol = dim(lvec)[1]; lpgvol = dim(lpvec)[1]
  #interp will ultimately become a dgCMatrix class object during construction we
  #will want to store it in terms of the row indices, column indices, and data
  #values needed to construct it its dimensions will ultimately be length(lvec)
  #by length(lpvel)
  
  #number of coarse grid points along x and y including upper edges
  xna = floor((ldim[1])/2); yna = floor((ldim[2])/2)
  #number of coarse grid points along x and y excluding upper edges
  xnm = floor((ldim[1]-1)/2); ynm = floor((ldim[2]-1)/2)
  
  ptotal = 4*xna*yna+2*xnm*yna+2*xna*ynm+xnm*ynm #total points in interp matrix
  pcount = 0 #counter to track number points added so far at each step
  padd = 0 #counter to track number of points to be added in each step
  
  interpi = rep(0,ptotal);interpj = rep(0,ptotal); interpx = rep(0,ptotal) #hold i,j, and x of interp matrix
  #Handle Type1 points first. These are where coarse and fine grid coincide
  #Performed using injection, ie LH=Lh
  #Deal with fine grids at 
  #x i, y j; x i-1, y i-1; x i-1, y j; x i, y j-1
  xCoarse = 1:xna; yCoarse = 1:yna
  padd = xna*yna
  #x i, y j fine grids
  xFine = 2*xCoarse; yFine = 2*yCoarse
  interpi[(pcount+1):(pcount+padd)] = as(array(lpgeo[xCoarse,yCoarse],
                                               dim=c(1,padd)),"numeric")
  interpj[(pcount+1):(pcount+padd)] = as(array(lgeo[xFine,yFine],
                                               dim=c(1,padd)),"numeric")
  interpx[(pcount+1):(pcount+padd)] = rep(1,padd)
  pcount = pcount+padd
  #x i-1, y j-1 fine grids
  xFine = 2*xCoarse-1; yFine = 2*yCoarse-1
  interpi[(pcount+1):(pcount+padd)] = as(array(lpgeo[xCoarse,yCoarse],
                                           dim=c(1,padd)),"numeric")
  interpj[(pcount+1):(pcount+padd)] = as(array(lgeo[xFine,yFine],
                               dim=c(1,padd)),"numeric")
  interpx[(pcount+1):(pcount+padd)] = rep(.25,padd)
  pcount = pcount+padd
  #x i-1, y j fine grids
  xFine = 2*xCoarse-1; yFine = 2*yCoarse
  interpi[(pcount+1):(pcount+padd)] = as(array(lpgeo[xCoarse,yCoarse],
                                           dim=c(1,padd)),"numeric") 
  interpj[(pcount+1):(pcount+padd)] = as(array(lgeo[xFine,yFine],
                                           dim=c(1,padd)),"numeric")
  interpx[(pcount+1):(pcount+padd)] = rep(.5,padd)
  pcount = pcount+padd
  #x i, y j-1 fine grids
  xFine = 2*xCoarse; yFine = 2*yCoarse-1
  interpi[(pcount+1):(pcount+padd)] = as(array(lpgeo[xCoarse,yCoarse],
                                           dim=c(1,padd)),"numeric") 
  interpj[(pcount+1):(pcount+padd)] = as(array(lgeo[xFine,yFine],
                                           dim=c(1,padd)),"numeric")
  interpx[(pcount+1):(pcount+padd)] = rep(.5,padd)
  pcount = pcount+padd
  #Now handle fine grids for x i+1, y j-1; x i+1, y j
  xCoarse = 1:xnm; yCoarse = 1:yna
  padd = xnm*yna
  #x i+1, y j-1
  xFine = 2*xCoarse+1; yFine = 2*yCoarse-1
  interpi[(pcount+1):(pcount+padd)] = as(array(lpgeo[xCoarse,yCoarse],
                                           dim=c(1,padd)),"numeric")
  interpj[(pcount+1):(pcount+padd)] = as(array(lgeo[xFine,yFine],
                                           dim=c(1,padd)),"numeric")
  interpx[(pcount+1):(pcount+padd)] = rep(.25,padd)
  pcount = pcount+padd
  #x i+1, y j
  xFine = 2*xCoarse+1; yFine = 2*yCoarse
  interpi[(pcount+1):(pcount+padd)] = as(array(lpgeo[xCoarse,yCoarse],
                                           dim=c(1,padd)),"numeric")
  interpj[(pcount+1):(pcount+padd)] = as(array(lgeo[xFine,yFine],
                                           dim=c(1,padd)),"numeric")
  interpx[(pcount+1):(pcount+padd)] = rep(.5,padd)
  pcount = pcount+padd
  #Now handle fine grids for x i-1, y j+1; x i, y j+1
  xCoarse = 1:xna; yCoarse = 1:ynm
  padd = xna*ynm
  #x i-1, y j+1
  xFine = 2*xCoarse-1; yFine=2*yCoarse+1
  interpi[(pcount+1):(pcount+padd)] = as(array(lpgeo[xCoarse,yCoarse],
                                           dim=c(1,padd)),"numeric")
  interpj[(pcount+1):(pcount+padd)] = as(array(lgeo[xFine,yFine],
                                           dim=c(1,padd)),"numeric")
  interpx[(pcount+1):(pcount+padd)] = rep(.25,padd)
  pcount = pcount+padd
  #x i, y j+1
  xFine = 2*xCoarse; yFine = 2*yCoarse+1
  interpi[(pcount+1):(pcount+padd)] = as(array(lpgeo[xCoarse,yCoarse],
                                           dim=c(1,padd)),"numeric")
  interpj[(pcount+1):(pcount+padd)] = as(array(lgeo[xFine,yFine],
                                           dim=c(1,padd)),"numeric")
  interpx[(pcount+1):(pcount+padd)] = rep(.5,padd)
  pcount = pcount+padd
  #Now handle fine grids for x i+1, y j+1
  xCoarse = 1:xnm; yCoarse = 1:ynm
  #x i+1, y j+1
  xFine = 2*xCoarse+1; yFine = 2*yCoarse+1
  interpi[(pcount+1):(pcount+padd)] = as(array(lpgeo[xCoarse,yCoarse],
                                           dim=c(1,padd)),"numeric")
  interpj[(pcount+1):(pcount+padd)] = as(array(lgeo[xFine,yFine],
                                           dim=c(1,padd)),"numeric")
  interpx[(pcount+1):(pcount+padd)] = rep(.25,padd)
  pcount = pcount+padd
  
  #still need to handle ragged lower edges
  
  return(sparseMatrix(i=interpi,j=interpj,x=interpx,dims=c(lpgvol,lgvol)))
}

oneLevel3Dinterp <- function(lvec,lgeo,ldim,lpvec,lpgeo,lpdim,periodic=TRUE) {
  lgvol = ldim[1]*ldim[2]*ldim[3]; lpgvol=lpdim[1]*lpdim[2]*lpdim[3]
  #interpi - indices of coarse grid entries in the matrix
  #interpj - indices of fine grid entries in the matrix
  gridFrame = data.frame(id=c(1:lgvol),xi=lvec[,1],yi=lvec[,2],zi=lvec[,3],
                         pv=(lvec[,1]%%2 + lvec[,2]%%2 + lvec[,3]%%2),
                         ev=((lvec[,1]%%ldim[1]<2)+(lvec[,2]%%ldim[2]<2)+
                               (lvec[,3]%%ldim[3]<2)))
  for(di in c(-1,1)) {
    for (dj in c(-1,1)) {
      for (dk in c(-1,1)) {
        parity = sum(abs(di)+abs(dj)+abs(dk))
        factor = 1 / (2*parity)
        interpi = 
        interpj = 
        
      }
    }
  }
  
}

vcycleGMG <- function(Amat,bvec,gridDims,levels=1,scaleExp=1,tol=1e-6) {
  RmatList <- buildRestrict()
  ImatList <- buildInterpolate()
  SmatList <- buildRelax()
  
  bvec = Rmat
  
  
}