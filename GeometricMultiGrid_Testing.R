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

build2DInterpolate <- function(vecList,geoList,dimList,periodic=FALSE) {
  levels = length(vecList)
  if (levels < 2) {
    stop("Error: Must have at least 2 levels")
  }
  lgeo = geoList[[level]]; lpgeo = geoList[[level+1]]
  lvec = vecList[[level]]; lpvec = vecList[[level+1]]
  ldim = dimList[[level]]; lpdim = dimList[[level+1]]
  #interp will ultimately become a dgCMatrix class object
  #during construction we will want to store it in terms of the
  #row indices, column indices, and data values needed to construct it
  #its dimensions will ultimately be length(lvec) by length(lpvec)
  interpi = c();interpj = c(); interpx = c()
  #Handle Type1 points first. These are where coarse and fine grid coincide
  #Performed using injection, ie LH=Lh
  interpj = c(interpi,c(1:length(lpvec)))
  interpi= lgeo[lpvec[interpi]*2]
  interpx=interpx(rep(1,length(lpvec)))
  #Handle Type2 points next. These lie on grid lines between coarse grids
  #Performed by summing over adjoining coarse grids
  #First instance where periodicity comes into play will branch code here
  if (periodic) { 
    #to be handled later.
  } else {
    #for 2D, there will be two line directions, X and Y.
    #start with x i-1/2 terms
    iCoarse = c(1:floor((ldim[1]-1)/2));jCoarse=c(1:floor(ldim[2]/2))
    interpj = c(interpj,as(array(lpgeo[iCoarse,jCoarse],
                       dim=c(1,length(iCoarse)*length(jCoarse))),"numeric"))
    iFine = 2*iCoarse-1;jFine=2*jCoarse
    interpi = c(interpi,as(array(lgeo[iFine,jFine],
                       dim=c(1,length(iFine)*length(jFine))),"numeric"))

  }
  
}

gridIndGen2D <- function(Iinds,Jinds) {
  ivec <- rep(Iinds,length(Jinds))
  jvec <- mapply(function(val) rep(val,length(Iinds)),Jinds)
  return(rbind(ivec,jvec))
}

vcycleGMG <- function(Amat,bvec,gridDims,levels=1,scaleExp=1,tol=1e-6) {
  RmatList <- buildRestrict()
  ImatList <- buildInterpolate()
  SmatList <- buildRelax()
  
  bvec = Rmat
  
  
}