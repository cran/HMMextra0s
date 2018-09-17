# plotVitpath2d plots the Viterbi path and the probability of each time point being in each state.
# object is a list containing y (the estimated Viterbi path) and 
## v (the estimated probability of each time point being in each state)
# varb: an integer indicates how many hours of data will be ploted on each page. The default is 8780.
# yearstart=2005, and yearend=2012 by default
## yearstart is the starting year of the data used
## yearend is the end year of the data used
# cols: a vector defines the colors to be used for each state. If col=NA, then the default colors will be used. 
plotVitpath2d <- function(object, R, Z, HMMest, len.dat=96432, varb=8780,
         yearstart=2005, yearend=2012, cols=NA, cex.lab=1.5, cex.axis=1.5){
## plot viterbi path
  mu <- HMMest$mu
  sig <- HMMest$sig
  m <- nrow( mu )
  nn <- nrow( R )
  y <- object$y
  v <- object$v
  Indlt <- 1:nn
  MUlat <- MUlon <- NULL
  for (il in 1:nn){
     MUlat[il] <- mu[y[il],1]
     MUlon[il] <- mu[y[il],2]
  }
  legdx <- legx <- NULL
  for (i in yearstart:yearend){
    for (j in 1:12){
     legdx <- append(legdx,(julian(as.Date(paste(i,j,"01",sep="-")))-julian(as.Date(paste(yearstart,"-01-01",sep=""))))*24+1)
     legx <- append(legx,paste(i,j,sep="-"))
    }
  }
  colVs <- 1:m
  if (any(!is.na(cols))){
     colVs <- cols
  }
  for (i in 1:ceiling(len.dat/varb)){
     par(mfrow=c(4,1), mar=c(5.1, 5.1, 1, 1))
     plot(Indlt[Z==1],R[Z==1,1],xlim=c((i-1)*varb,i*varb),xlab="Time (year-month)",ylab="Latitude",cex.lab=cex.lab,axes=FALSE)
     axis(1,legdx,legx,cex.axis=cex.axis)
     axis(2,cex.axis=cex.axis)
     box()
     lines(1:nn,MUlat,col="red")
     plot(Indlt[Z==1],R[Z==1,2],xlim=c((i-1)*varb,i*varb),xlab="Time (year-month)",ylab="Longitude",cex.lab=cex.lab,axes=FALSE)
     axis(1,legdx,legx,cex.axis=cex.axis)
     axis(2,cex.axis=cex.axis)
     box()
     lines(1:nn,MUlon,col="red")
     plot(1:nn,y,type="l",xlim=c((i-1)*varb,i*varb),xlab="Time (year-month)",ylab="Viterbi Path",cex.lab=cex.lab,axes=FALSE)
     axis(1,legdx,legx,cex.axis=cex.axis)
     axis(2,cex.axis=cex.axis)
     box()
     plot(1:nn,rep(1,nn), type="l",xlim=c((i-1)*varb,i*varb),ylim=c(-0.2,1),xlab="Time (year-month)",ylab="State Probability",cex.lab=cex.lab,axes=FALSE)
     axis(1,legdx,legx,cex.axis=cex.axis)
     axis(2,cex.axis=cex.axis)
     box()
     points(1:nn,v[,1],type="l")
     polygon(c(1:nn,seq(nn,1,-1)),c(rep(0,nn),rev(v[,1])))
     points(1:nn,v[,1]+v[,2],type="l")
     polygon(c(1:nn,seq(nn,1,-1)),c(v[,1],rev(v[,1]+v[,2])),col=colVs[2])
     if (m >= 4){
       for (ji in 3:(m-1)){
         points(1:nn,apply(v[,1:ji],1,sum),type="l")
         polygon(c(1:nn,seq(nn,1,-1)),c(apply(v[,1:(ji-1)],1,sum),rev(apply(v[,1:(ji)],1,sum))),col=colVs[ji])
       }
     }
     polygon(c(1:nn,seq(nn,1,-1)),c(apply(v[,1:(m-1)],1,sum),rep(1,nn)),col=colVs[m])
     readline("Press <enter> to continue")
   }
}
