\name{Kii}
\alias{Kii}
\title{
Tremor data in the Kii region in 2002 and 2003 for use in function hmm0norm2d
}
\description{
A data frame containing a subset (in years 2002 and 2003) of Kii tremor data used in Wang et al. (2018). The columns are named \code{"year"}, \code{"month"}, \code{"day"}, \code{"hour"}, \code{"lat"}, \code{"lon"}.

We provide some R code in the Examples below for how to convert this dataset into the variables \code{R} and \code{Z} used in the function \code{\link{hmm0norm2d}}. This dataset can be obtained directly from the Slow Earthquake Database http://www-solid.eps.s.u-tokyo.ac.jp/~sloweq/.

If you have your own way to convert the data into the variables \code{R} and \code{Z}, then you can go to the function \code{\link{hmm0norm2d}} directly. 
}
\usage{
data(Kii)
}
\format{
A data frame with 1112 rows, each row representing the hour in which tremor events occurred. It contains the following variables:
\describe{
\item{year, month, day, hour}{time of tremor occurrence.}
\item{lat}{latitude of the tremor event in that hour.}
\item{lon}{longitude of the tremor event in that hour.}
}}

\references{
%% ~put references to the literature/web site here ~
Wang, T., Zhuang, J., Buckby, J., Obara, K. and Tsuruoka, H. (2018) Identifying the recurrence patterns of non-volcanic tremors using a 2D hidden Markov model with extra zeros. Journal of Geophysical Research, \doi{10.1029/2017JB015360}.
Obara, K., Tanaka, S., Maeda, T., & Matsuzawa, T. (2010) Depth-dependent activity of non-volcanic tremor in southwest Japan, Geophysical Research Letters, 37, L13306. \doi{10.1029/2010GL043679}.
Maeda, T., & Obara. K. (2009) Spatio-temporal distribution of seismic energy radiation from low-frequency tremor in western Shikoku, Japan, J. Geophys. Res., 114, B00A09, \doi{10.1029/2008JB006043}. 
}
\examples{
data(Kii)
year <- Kii$year
month <- Kii$month
day <- Kii$day
hour <- Kii$hour
lat <- Kii$lat
lon <- Kii$lon

## Transform the time into days*100+hour. Can use other transformation.
## The purpose is to make sure that each hour of a day has a unique number.
xd <- NULL
for (i in 1:nrow(Kii))
  xd[i] <- julian(as.Date(paste(year[i],month[i],day[i],sep="-")))*100+hour[i]

## Create a unique number for each hour in the years 2002 and 2003
## This is to match with xd above, so that we can create the Z variable
## which is 0 for the hours without any tremor occurrence and
## 1 for the hours with tremor events.
a <- seq( julian(as.Date("2002-01-01")), julian(as.Date("2002-12-31")), 1 )*100
b <- seq( julian(as.Date("2003-01-01")), julian(as.Date("2003-12-31")), 1 )*100
aa <- rep(a,each=24)+rep(0:23,times=length(a))
bb <- rep(b,each=24)+rep(0:23,times=length(b))

## Combine all the tremor events which occurred 
## in the same hour to be one tremor cluster.
## Kii has maximum 4 events in the same hour 
## so we used the code below.
## One can adjust the code for regions with more 
## tremor events in the same hour.
## indt: actual time as in each hour
Time <- c(aa,bb)
lt <- length(Time)
indt <- 1:lt

Tim <- Lat <- Lon <- NULL
j <- 1
while (j <= nrow(Kii)-3){
  i <- j
  if (xd[i+3]==xd[i] & xd[i+2]==xd[i] & xd[i+1]==xd[i]){
    Tim <- append(Tim,xd[i])
    Lat <- append(Lat,mean(lat[i:(i+3)]))
    Lon <- append(Lon,mean(lon[i:(i+3)]))
    j <- i+4
  }else{
    if (xd[i+2]==xd[i] & xd[i+1]==xd[i]){
      Tim <- append(Tim,xd[i])
      Lat <- append(Lat,mean(lat[i:(i+2)]))
      Lon <- append(Lon,mean(lon[i:(i+2)]))
      j <- i+3
    }else{
      if (xd[i+1]==xd[i]){
        Tim <- append(Tim,xd[i])
        Lat <- append(Lat,mean(lat[i:(i+1)]))
        Lon <- append(Lon,mean(lon[i:(i+1)]))
        j <- i+2
      }else{
        Tim <- append(Tim,xd[i])
        Lat <- append(Lat,lat[i])
        Lon <- append(Lon,lon[i])
        j <- i+1
      }
    }
  }
}
Tim <- append(Tim,xd[(nrow(Kii)-1):nrow(Kii)])
Lat <- append(Lat,lat[(nrow(Kii)-1):nrow(Kii)])
Lon <- append(Lon,lon[(nrow(Kii)-1):nrow(Kii)])

## Create a data frame in which each hour is a point
## Those hours when there was no tremor, we set the 
## number of tremors as 0

data1 <- array(0,dim=c(lt,3))
Thour <- NULL
for (i in 1:length(Tim)){
  use <- Time==Tim[i]
  idtem <- (1:lt)[use]
  Thour <- append(Thour,idtem)
  data1[idtem,2] <- Lat[i]
  data1[idtem,3] <- Lon[i]
}
data1[,1] <- indt ## Every hour is one time point

###########################################################
###########   Data for time series analysis   #############  
###########################################################
lt <- length(indt)
Z <- rep(0,lt)
Z[Thour] <- 1
R <- data1[,2:3]


###########################################################
# Setting up initial values for analysing real-world data
## nk is the number of states for the fitted model
### In this example we use nk=3
###########################################################

LL <- -10^200 ## A very small value to compare with
## the log likelihood from the model

nk = 3

gamma <- array(NA,dim=c(nk,nk))
mu <- array(NA,dim=c(nk,2))
sig <- array(NA,dim=c(2,2,nk))
pie <- array(NA,dim=c(1,nk))

kk <- 1
N <- 2
while(kk<N)
{
  temp <- matrix(runif(nk*nk,0,1),ncol=nk)
  diag(temp) = diag(temp) + rpois(1,6) * apply(temp, 1, sum)
  temp <- temp * matrix(rep(1/apply(temp, 1, sum), ncol(temp)), ncol=ncol(temp), byrow=FALSE)
  gamma <- temp
  
  R1min <- min((R[,1])[R[,1]>=1e-6])
  R1max <- max((R[,1])[R[,1]>=1e-6])
  R2min <- min((R[,2])[R[,2]>=1e-6])
  R2max <- max((R[,2])[R[,2]>=1e-6])
  temp <- cbind(runif(nk,R1min,R1max),runif(nk,R2min,R2max))
  temp <- temp[order(temp[,2]),]
  mu <- temp
  
  sdR1 <- sd((R[,1])[R[,1]>=1e-6])
  sdR2 <- sd((R[,2])[R[,2]>=1e-6])
  for (j in 1:nk){
    temp <- matrix(runif(4,0.0001,max(sdR1,sdR2)), ncol=2)
    temp[1,2] <- temp[2,1] <- runif(1,-1,1)* sqrt(prod(diag(temp)))
    sig[, ,j] <- temp
  }
  
  pie <- matrix(sort(c(runif(1, 0, 0.01),runif(nk-1, 0, 1))), nrow = 1, byrow = TRUE )
  
  delta <- c(6,runif(nk-1, 0,1)) 
  delta <- delta/sum(delta)
  
  tryCatch({
    temp <- hmm0norm2d(R, Z, pie, gamma, mu, sig, delta)
    kk<-kk+1
    if( LL <= temp$LL){
      HMMest <- temp
      LL =HMMest$LL
      eval(parse(text=paste('HMM',kk,'est = HMMest',sep="")))
#      eval(parse(text=paste('save(HMM',kk,'est, file="HMM',kk,'est.image")',sep='')))
## Uncomment the line above if you would like to save the result into a .image file.
    }
  }, error=function(e){})
  print(kk)
}

}
\seealso{
\code{\link{hmm0norm2d}}
}
\keyword{datasets}



