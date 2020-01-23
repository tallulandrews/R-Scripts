# From: https://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r

sparse.cor <- function(x){
  require("Matrix")
  n <- nrow(x)
  m <- ncol(x)
  ii <- unique(x@i)+1 # rows with a non-zero element

  Ex <- colMeans(x)
  nozero <- as.vector(x[ii,]) - rep(Ex,each=length(ii))        # colmeans

  covmat <- ( crossprod(matrix(nozero,ncol=m)) +
              crossprod(t(Ex))*(n-length(ii))
            )/(n-1)
  sdvec <- sqrt(diag(covmat))
  covmat/crossprod(t(sdvec))
}

sparse.cor2 <- function(x){
    n <- nrow(x)

    covmat <- (crossprod(x)-2*(colMeans(x) %o% colSums(x))
              +n*colMeans(x)%o%colMeans(x))/(n-1)

    sdvec <- sqrt(diag(covmat)) # standard deviations of columns
    covmat/crossprod(t(sdvec)) # correlation matrix
}

sparse.cor3 <- function(x){
    memory.limit(size=10000)
    n <- nrow(x)

    cMeans <- colMeans(x)
    cSums <- colSums(x)

    # Calculate the population covariance matrix.
    # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
    # The code is optimized to minize use of memory and expensive operations
    covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
    crossp <- as.matrix(crossprod(x))
    covmat <- covmat+crossp

    sdvec <- sqrt(diag(covmat)) # standard deviations of columns
    covmat/crossprod(t(sdvec)) # correlation matrix
}

sparse.cor4 <- function(x){
    n <- nrow(x)
    cMeans <- colMeans(x)
    covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
    sdvec <- sqrt(diag(covmat)) 
    cormat <- covmat/tcrossprod(sdvec)
    list(cov=covmat,cor=cormat)
}
