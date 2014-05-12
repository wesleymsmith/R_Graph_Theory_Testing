require(Matrix);
require(MASS);

stepLU <- function(A,i,n) {
  A[(i+1):n,i] = A[(i+1):n,i] / A[i,i]
  A[(i+1):n,(i+1):n] = A[(i+1):n,(i+1):n]-
    matrix(A[(i+1):n,i],nrow=n-i,ncol=1) %*%
    matrix(A[i,(i+1):n],nrow=1,ncol=n-i)
  return(A)
}
matLU <- function(A,ni=dim(A)[1]) {
  n=dim(A)[1]
  for(i in c(1:(ni-1))) {
    A=stepLU(A,i,n) 
  }
  return(A)
}

stepLL <- function(A,i,n) {
  A[i,i] = sqrt(A[i,i])
  val = 1 / A[i,i]
  A[(i+1):n,i] = val * A[(i+1):n,i]
  A[(i+1):n,(i+1):n] = A[(i+1):n,(i+1):n] -
    matrix(A[(i+1):n,i],nrow=n-i,ncol=1) %*%
    matrix(A[i,(i+1):n],nrow=1,ncol=n-i)
  return(A)
}
matLL <- function(A,ni=dim(A)[1]-1) {
  n = dim(A)[1]
  AA=A
  for ( i in c( 1:(ni) ) ) {
    A[i,i] = sqrt(AA[i,i])
    val = 1 / A[i,i]
    A[(i+1):n,i] = val * AA[(i+1):n,i]
    A[i,(i+1):n] = A[(i+1):n,i]
    A[(i+1):n,(i+1):n]= A[(i+1):n,(i+1):n] - 
      matrix(A[(i+1):n,i],nrow=n-i,ncol=1) %*% 
      matrix(A[i,(i+1):n],nrow=1,ncol=n-i)
  }
  Asubset = subset(summary(A),i>=j)
  Lsubset=subset(Asubset,j<=ni)
  Lmat = sparseMatrix(i=c(Lsubset$i,(ni+1):n),j=c(Lsubset$j,(ni+1):n),
                      x=c(Lsubset$x,rep(1,n-ni)))
  Cmat = sparseMatrix(i=c(1:ni),j=c(1:ni),x=rep(1,ni),dims=c(n,n))
  Cmat[(ni+1):n,(ni+1):n] = AA[(ni+1):n,(ni+1):n] - A[(ni+1):n,1:ni] %*% A[1:ni,(ni+1):n]
  return(list(L=Lmat,C=Cmat))
}