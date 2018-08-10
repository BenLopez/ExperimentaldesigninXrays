# Gcor # 

# Compute the Gaussian correlation function for two matrices with a different theta for each input.  theta must be a vector.

CF_GaussianCorrelationFunction <- function(X, Y, theta){
  if(length(theta)>1){
    Xnew <- X%*%diag(1/theta)
    Ynew <- Y%*%diag(1/theta)
  }else{
    Xnew <- matrix(X/theta, ncol=1)
    Ynew <- matrix(Y/theta, ncol=1)
  }
  gcmat <- as.matrix(pdist(Xnew, Ynew))^2
  gcmat <- exp(-gcmat)
  return(gcmat)}

# Gcorla # 

# Compute the Gaussian correlation function for one matrix with a different theta for each input.  theta must be a vector.  This just avoids the error that is returnede by
# pdist when the two matrices X and Y are the same.

CF_GaussianCorrelationFunctionLa <- function(X, theta){
  if(length(theta)>1){
    Xnew <- X%*%diag(1/theta)
  }else{
    Xnew <- matrix(X/theta, ncol=1)
  }
  gcmat <- as.matrix(dist(Xnew))^2
  gcmat <- exp(-gcmat)
  return(gcmat)}
