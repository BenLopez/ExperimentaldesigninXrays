### X-Ray Design ###

library("tgp")
library("MASS")

# Import the data.
xrx <- read.table("X_Sam", header = FALSE, sep=",")[,c(1,3,2)]   # proportion pmma, thickness, KV
dimnames(xrx)[[2]] <- c("V1", "V2", "V3")
xry <- read.table("Y_Sam", header = FALSE, sep=",")[,1:6]

# Ranges
rangee <- apply(xrx, 2, range)

# Could Standardise
me <- (rangee[2,] - rangee[1,])/2
te <- (rangee[2,] + rangee[1,])/2

xrxhat <- t((t(xrx)-te)/me)

# Log transform the output.
lxry <- log(xry)

# Ratio of outputs i from 2-6 to 1.
fx <- matrix(0, nrow=296, ncol=6)
fx[,1] <- lxry[,1]
for(i in 2:6){
   fx[,i] <- lxry[,i] - lxry[,1]
}

# Not sure how to choose s2f, s2n or l.
s2n <- 0.0000001
l <- c(0.25, 0.5, 15)

# Model discrepancy and measurement error.
s2md <- 0.01*diag(6)
s2e <- 0.01*diag(6)

# First of all we need a model - let's start with the full quadratic model.
y <- rep(1, nrow(xrx))
FLmod <- lm(y ~ ., data=xrx)

# The active variable names as a vector.
actv <- names(FLmod$coefficients)[-1]

# Construct the formula for considering full quadratic models with interaction terms.
FO <- paste(actv, collapse="+")
SO <- paste( "(",FO,"):(",FO,")","+",paste("I(",actv,"^2)",collapse="+",sep=""))
form1 <- as.formula(paste("y~",SO))
FQmod <- lm(form1, data=xrx)

# Use a 2/3, 1/3 training/test set split to diagnose the emulator.
sammy <- sample(296)
dg <- MOGPEmulator(X=xrx[sammy[1:200],], Xstar=xrx[sammy[201:296],], Y=lxry[sammy[1:200],], l=l, s2n=s2n, model=FQmod)
mahala <- rep(0, 96)
for(j in 1:96){mahala[j] <- mahalanobis(lxry[sammy[200+j],], dg$m[j,], dg$vax[,,j])}
plot(mahala)
abline(h=qchisq(0.99, 6))

# Now we need a training set of points - we can add some extra points for t=0 p=a for a= 0.2, 0.4, 0.6, 0.8, 1.
th0 <- which(xrx[,2]==0)
ep0 <- matrix(0, nrow=40, ncol=3)
ep0[,1] <- c(rep(0.2, 8), rep(0.4, 8), rep(0.6, 8), rep(0.8, 8), rep(1, 8))
ep0[,3] <- rep(seq(from=50, to=120, by=10), 5)
# Add these extra points onto the bottom of the old matrix.
xtr <- rbind(xrx, ep0)
epy0 <- kronecker(rep(1, 5), as.matrix(fx[th0,]))
fxtr <- data.frame(rbind(fx, epy0))

# Now we need a test set of points.
ll <- c(0, 0, 50)   # lower limits of LH.
ul <- c(1, 15, 120)    # upper limits of LH.
xLH <- lhs(6000, matrix(c(ll, ul), nrow=3, byrow=FALSE))   # LH for test points
# Trim down the LH.
trim <- which(xLH[,2] <= 3 + 12*xLH[,1]^2)
xte <- data.frame(xLH[trim,])  # Test points.
dimnames(xte)[[2]] <- c("V1", "V2", "V3")

# Now we emulate.
fxu <- MOGPEmulator(X=xtr, Xstar=xte, Y=fxtr, l=l, s2n=s2n, model=FQmod)

# Pairs plot and density plot of outputs.
durt3 <- c("gray35", "white", "cyan", "cyan2", "cornflowerblue", "blue", "blue2", "blue3", "darkorchid4", "darkorchid3", "darkorchid1")
pairs(fxu$m[1:1000,], pch=16)
par(mfrow=c(7, 7))
for(i in 1:7){
   for(j in 1:7){
      if(i==j){
         par(mar=c(0.3, 0.3, 0.3+0.7, 0.3))
         hist(fxu$m[,i], main=dimnames(fxtr)[[2]][i], xlab="", ylab="", cex.main=1, xaxt="n", yaxt="n", col="blue")
      }
      else{
         par(mar=rep(0.3, 4))
         smo <- kde2d(x=fxu$m[isubx,i], y=fxu$m[isubx,j], n=25, lims=c(min(fxu$m[isubx,i]), max(fxu$m[isubx,i]), min(fxu$m[isubx,j]), max(fxu$m[isubx,j])))
         levels <- seq(0, max(smo$z)*1.01, len=length(durt3)) # levels chosen specific to each plot
         plot(0, xlim=c(min(fxu$m[isubx,1]), max(fxu$m[isubx,1])), ylim=c(min(fxu$m[isubx,2:7]), max(fxu$m[isubx,2:7])), ty="n", axes=F, frame.plot=T, xlab="", ylab="", xaxs="i", yaxs="i")
         .filled.contour(as.double(smo$x), as.double(smo$y), smo$z, as.double(levels), col=durt3)
   }
   }
}

# Let's have a go at design.
zsamx <- zsamp(Efx=fxu$m, nz=3, s2fx=fxu$vax, s2md, s2e, df=0, Emd=0)

# Multivariate Utility Calculator
umu <- UCmult(Efx=fxu$m, zsam=zsamx, zh=0, s2fx=fxu$vax, s2md=s2md, s2e=s2e, impl=3, Emd=0)
umuB <- UCmultB(Efx=fxu$m, zsam=zsamx, zh=0, s2fx=fxu$vax, s2md=s2md, s2e=s2e, impl=3, Emd=0)

# Experiment 1 tells us:
subse <- 1
umu <- UCmult(Efx=fxu$m[,subse], zsam=zsamx[,subse], zh=0, s2fx=fxu$vax[subse,subse,], s2md=s2md[subse], s2e=s2e[subse], impl=3, Emd=0)

# Best experiment given that we have experiment 1 anyway.
u2s <- matrix(0, nrow=10, ncol=10)
for(i in 2:10){
   for(j in 1:(i-1)){
      subse <- c(1,i+1,j+1)
      u2s[i,j] <- UCmult(Efx=fxu$m[,subse], zsam=zsamx[,subse], zh=0, s2fx=fxu$vax[subse,subse,], s2md=s2md[subse,subse], s2e=s2e[subse,subse], impl=3, Emd=0)$U
   }
}
for(i in 1:10){
   subse <- c(1,i+1)
u2s[i,i] <- UCmult(Efx=fxu$m[,subse], zsam=zsamx[,subse], zh=0, s2fx=fxu$vax[subse,subse,], s2md=s2md[subse,subse], s2e=s2e[subse,subse], impl=3, Emd=0)$U
}

# Emulator Correlation
fxucorr <- matrix(0, nrow=11, ncol=11)
for(i in 1:11){
   for(j in 1:11){
      fxucorr[i,j] <- fxu$SigmaHat[i,j]/sqrt(fxu$SigmaHat[i,i]*fxu$SigmaHat[j,j])
   }
}

# Variance Resolution Utility Example

ft <- rep(0,7)
for(i in 2:7){
   subse <- c(1,i)
   ft[i] <- Restest(xv = xte, 
                 Efx = fxu$m[,subse], 
                 zsam = zsamx[,subse], 
                 zh = 0, 
                 s2fx = fxu$vax[subse, subse,], 
                 s2md = s2md[subse, subse], 
                 s2e = s2e[subse, subse], 
                 impl = 3, 
                 Emd = 0, 
                 type = "trace", 
                 wei = c(4,1,0), 
                 svx = "FALSE")
}

hist(ft$vx, xlim=c(0,1), xaxs="i", freq=FALSE)










