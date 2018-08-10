CF_ExponentialFamily <- function(x , xstar , l , p)
{
  # Expenential family correlation length family
  
  # Turn into matrices
x <- as.matrix(x)
xstar <- as.matrix(xstar)
l <- as.matrix(l^p[1])
p <- as.matrix(p)

sizex <- dim(x)
sizexstar <- dim(xstar)


distmatrix = matrix(0, sizex[1] , sizexstar[1] , sizex[2])

if(sizex[2] > 1){

  for (i in 1:sizex[2])
{
  distmatrix[  ,  , i] <- (( as.matrix(pdist( as.matrix((x[,i])) , as.matrix( (xstar[,i]) )  )) )^p[1])/l[i]
}

} else
{
  distmatrix <- (( as.matrix(pdist( as.matrix((x[,i])) , as.matrix( (xstar[,i]) )  )) )^p[1])/l[1]
}


if(sizex[2] > 1)
{
distmatrix < - apply(distmatrix , 3 , sum)
}

KXXstar <- exp( -0.5 * distmatrix )
return(KXXstar)
}


BayesLinearEmulatorGLSEstimates <- function(y , x , xstar , w , l , p , h )
{

  KXX <- CF_ExponentialFamily(x , x , l , p)
  KXstarX <-  CF_ExponentialFamily( xstar , x , l , p)
  KXstarXstar <-  CF_ExponentialFamily(xstar , xstar , l , p)
  H <- h(x)
  sizeh = dim(H)
  Hstar <- h(xstar)
  
  D <- KXX + w*diag(sizex[1])
  L = t(chol(D))
  
  alpha <- solve(L , H)
  Q <- t(chol(t(alpha)%*%alpha))
  # Calulate approximated adjusted expectation of Beta
  Betastar <- solve(t(Q) , solve(Q))%*%t(alpha)%*%solve( L , y)
  
  zeta <- solve(L , y - H%*%Betastar)
  
  # Calulate approximated adjusted expectation of sigma
  Sigmastar <- 1/(sizex[1] - sizeh[2]) * ( t(zeta)%*%zeta )
  
  omega <-solve(L, t(KXstarX))
  kappa <-  solve( Q , t( Hstar - t(omega)%*%alpha ))
  
  E_z_y <- Hstar%*%Betastar + KXstarX %*% solve(t(L) , zeta) 
  V_z_y <- Sigmastar[1]*(  KXstarXstar  - t(omega)%*%omega + t(kappa)%*%kappa )
  
  return(list(E_z_y , V_z_y))

  }


