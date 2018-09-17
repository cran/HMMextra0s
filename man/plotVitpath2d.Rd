\name{plotVitpath2d}
\alias{plotVitpath2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Viterbi Path of a Bivariate HMM With Extra Zeros
}
\description{
Plot the 2-D data, Viterbi path and the probability of each time point being in each state over time.
}
\usage{
plotVitpath2d(object, R, Z, HMMest, len.dat=96432, varb=8780,
  yearstart=2005, yearend=2012, cols=NA, cex.lab=1.5, cex.axis=1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{is a list containing \code{y} (the estimated Viterbi path) and \code{v} (the estimated probability of each time point being in each state). This object is returned from running \code{Viterbi.hmm0norm2d(R, Z, HMMest)}.}
  \item{R}{is the observed data. \code{R} is a \eqn{T * 2} matrix, where \eqn{T} is the number of observations.}
  \item{Z}{is the binary data with the value 1 indicating that an event was observed and 0 otherwise. \code{Z} is a vector of length \eqn{T}.}
  \item{HMMest}{is a list which contains pie, gamma, sig, mu, and delta (the bivariate HMM parameter estimates).}
  \item{len.dat}{is the length of the data, that is, the number of time points. Default is 96432.}
  \item{varb}{is an integer indicating the length of data that will be ploted on each page. The default is 8780.}
  \item{yearstart}{is the starting year of the data used. Default is 2005.}
  \item{yearend}{is the end year of the data used. Default is 2012.}
  \item{cols}{is a vector defines the colors to be used for different states. If col=NA, then the default colors will be used.}
  \item{cex.lab}{specifies the size of the axis label text.}
  \item{cex.axis}{specifies the size of the tick label numbers/text.}
}
\details{
The returned object has four panels. Top two panels: Observed latitudes and longitudes with the center \eqn{\hat{\mu}_i} of each state overlaid as the red lines; third panel: tracked most likely state sequence of the HMM; bottom panel: the estimated probability of the data being in each state, with blank representing the probability of being in the last state (typically the plot looks better if the last state represents the background state with the minimum proportion of tremor occurrence). Some example plots are in the supplementary file of the reference Wang et al. (2018).
}

\references{
%% ~put references to the literature/web site here ~
Wang, T., Zhuang, J., Buckby, J., Obara, K. and Tsuruoka, H. (2018) Identifying the recurrence patterns of non-volcanic tremors using a 2D 
hidden Markov model with extra zeros. Journal of Geophysical Research, \doi{10.1029/2017JB015360}.
}
\author{
Ting Wang and Jiancang Zhuang
}
\examples{
pie <- c(0.008,0.2,0.4)
gamma <- matrix(c(0.99,0.007,0.003,
                  0.02,0.97,0.01,
                  0.04,0.01,0.95),byrow=TRUE, nrow=3)
mu <- matrix(c(35.03,137.01,
               35.01,137.29,
               35.15,137.39),byrow=TRUE,nrow=3)
sig <- array(NA,dim=c(2,2,3))
sig[,,1] <- matrix(c(0.005, -0.001,
                   -0.001,0.01),byrow=TRUE,nrow=2)
sig[,,2] <- matrix(c(0.0007,-0.0002,
                    -0.0002,0.0006),byrow=TRUE,nrow=2)
sig[,,3] <- matrix(c(0.002,0.0018,
                     0.0018,0.003),byrow=TRUE,nrow=2)
delta <- c(1,0,0)
y <- sim.hmm0norm2d(mu,sig,pie,gamma,delta, nsim=5000)
R <- y$x
Z <- y$z
HMMEST <- hmm0norm2d(R, Z, pie, gamma, mu, sig, delta)
Viterbi3 <- Viterbi.hmm0norm2d(R,Z,HMMEST)
plotVitpath2d(Viterbi3, R, Z,HMMEST,len.dat=5000,varb=5000,yearstart=2005, yearend=2005)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot,classification,Viterbi}



