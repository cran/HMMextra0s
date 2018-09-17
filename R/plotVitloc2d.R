# plotVitloc2d plots the locations of the observations with differet colors indicating different states.
# object is a list containing y (the estimated Viterbi path) and 
## v (the estimated probability of each time point being in each state)
# CI.level is a scalar or a vector the confidence level for the ellipse contour of each state 
# npoints is the number of points used in the ellipse. Default is 100
# cols: a vector defines the colors to be used for each state. If col=NA, then the default colors will be used. 
require(ellipse)
plotVitloc2d <- function(object, R, Z, HMMest, CI.level=0.95, npoints=100, cols=NA,
        cex.lab=1.5, cex.axis=1.5, cex=1, cex.text=2){
## plot tremor locations according to states
    y <- object$y
    colV <- y[Z==1]
    mu <- HMMest$mu
    sig <- HMMest$sig
    m <- nrow( mu )
    colVs <- 1:m
    if (any(is.na(cols))){
      colV[colV==1] <- gray(0.2)
      colVs[1] <- gray(0.2)
    }else{
      for (j in 1:max(y)){
        colV[colV==j] <- cols[j]
      }
      colVs <- cols
    }
    prange <- NULL
    for (inds in 1:m){
        ellipf <- function(CI.level){ ellipse(sig[,,inds],centre=mu[inds,],level=CI.level,npoints=npoints)}
        ellip <- apply(t(CI.level),2,ellipf)
        nellip <- ncol(ellip)
        for (nj in 1:nellip)
        prange <- cbind(prange,t(apply(matrix(ellip[,nj],nrow=100)[,c(2,1)],2,range)))
    }
    plot(R[Z==1,2],R[Z==1,1],col=colV,pch=y[Z==1],xlab="Longitude",ylab="Latitude",
        cex.lab=cex.lab,cex.axis=cex.axis,cex=cex,xlim=range(c(prange[1,],R[Z==1,2])),
        ylim=range(c(prange[2,],R[Z==1,1])))
    for (inds in 1:m){
        text(mu[inds,2],mu[inds,1],col=1,paste(inds),cex=cex.text)
        ellipf <- function(CI.level){ ellipse(sig[,,inds],centre=mu[inds,],level=CI.level,npoints=npoints)}
        ellip <- apply(t(CI.level),2,ellipf)
        nellip <- ncol(ellip)
        for (nj in 1:nellip)
          points(matrix(ellip[,nj],nrow=100)[,c(2,1)],col=colVs[inds],lwd=2, type="l")
    }
}
