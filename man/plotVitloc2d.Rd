\name{plotVitloc2d}
\alias{plotVitloc2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Classified 2-D Data of a Bivariate HMM With Extra Zeros
}
\description{
Plot the classified 2-D data with different colours representing different hidden states (or different clusters) obtained from the Viterbi path and confidence contours.
}
\usage{
plotVitloc2d(object, R, Z, HMMest, CI.level=0.95, npoints=100, cols=NA, 
cex.lab=1.5, cex.axis=1.5, cex=1, cex.text=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{is a list containing \code{y} (the estimated Viterbi path) and \code{v} (the estimated probability of each time point being in each state). This object is returned from running \code{Viterbi.hmm0norm2d(R, Z, HMMest)}.}
  \item{R}{is the observed data. \code{R} is a \eqn{T * 2} matrix, where \eqn{T} is the number of observations.}
  \item{Z}{is the binary data with the value 1 indicating that an event was observed and 0 otherwise. \code{Z} is a vector of length \eqn{T}.}
  \item{HMMest}{is a list which contains pie, gamma, sig, mu, and delta (the bivariate HMM parameter estimates).}
  \item{CI.level}{is a scalar or a vector, the confidence level for the ellipse contour of each state. Default is 0.95.}
  \item{npoints}{is the number of points used in the ellipse. Default is 100.}
  \item{cols}{is a vector defines the colors to be used for different states. If col=NA, then the default colors will be used.}
  \item{cex.lab}{specifies the size of the axis label text.}
  \item{cex.axis}{specifies the size of the tick label numbers/text.}
  \item{cex}{specifies the size of the points.}
  \item{cex.text}{specifies the size of the text indicting the state number.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\references{
%% ~put references to the literature/web site here ~
Wang, T., Zhuang, J., Buckby, J., Obara, K. and Tsuruoka, H. (2018) Identifying the recurrence patterns of non-volcanic tremors using a 
2D hidden Markov model with extra zeros. Journal of Geophysical Research, \doi{10.1029/2017JB015360}.
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
plotVitloc2d(Viterbi3, R, Z,HMMEST)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot,2-D classification,Viterbi}

