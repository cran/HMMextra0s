\name{hmm0norm2d}
\alias{hmm0norm2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter Estimation of a bivariate HMM with Extra Zeros
}
\description{
Calculates the parameter estimates of an HMM with bivariate observations having extra zeros.
}
\usage{
hmm0norm2d(R, Z, pie, gamma, mu, sig, delta, tol=1e-6, print.level=1, fortran = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{is the observed data. \code{R} is a \eqn{T * 2} matrix, where \eqn{T} is the number of observations.}
  \item{Z}{is the binary data with the value 1 indicating that an event was observed and 0 otherwise. \code{Z} is a vector of length \eqn{T}.}
  \item{pie}{is a vector of length \eqn{m}, the \eqn{j}th element of which is the probability of \eqn{Z=1} when the process is in state \eqn{j}.}
  \item{gamma}{is the transition probability matrix (\eqn{m * m}) of the hidden Markov chain.}
  \item{mu}{is an \eqn{m * 2} matrix, the \eqn{j}th row of which is the mean of the bivariate (Gaussian) distribution of the observations in state \eqn{j}.}
  \item{sig}{is a \eqn{2 * 2 * m} array. The matrix \code{sig[,,j]} is the variance-covariance matrix of the bivariate (Gaussian) distribution of the observations in state \eqn{j}.}
  \item{delta}{is a vector of length \eqn{m}, the initial distribution vector of the Markov chain.}
  \item{tol}{is the tolerance for testing convergence of the iterative estimation process. Default is 1e-6. For initial test of model fit to your data, a larger tolerance (e.g., 1e-3) should be used to save time.}
  \item{print.level}{controls the amount of output being printed. Default is 1. If \code{print.level=1}, only the log likelihoods and the differences between the log likelihoods at each step of the iterative estimation process, and the final estimates are printed. If \code{print.level=2}, the log likelihoods, the differences between the log likelihoods, and the estimates at each step of the iterative estimation process are printed.}
  \item{fortran}{is logical, and determines whether Fortran code is used; default is \code{TRUE}.}
}
\details{
  Setting up initial values for the real world data can be challenging, especially when the model is large (the number of states is big). 
  In the example below, we include a simple way to set up initial values. If the model is large, the model fitting process should be repeated for
  many different initial values. In the example below, we set the number of initial values to be \eqn{N=2} for the ease of compilation. For real-world data 
  analysis, taking the 2D model for the tremor data in Wang et al. (2018) for example, we used at least \eqn{N=1000} initial values for the large models
  with more than 15 hidden states. 
}
\value{
  \item{pie}{is the estimated probability of \eqn{Z=1} when the process is in each state.}
  \item{mu}{is the estimated mean of the bivariate (Gaussian) distribution of the observations in each state.}
  \item{sig}{is the estimated variance-covariance matrix of the bivariate (Gaussian) distribution of the observations in each state.}
  \item{gamma}{is the estimated transition probability matrix of the hidden Markov chain.}
  \item{delta}{is the estimated initial distribution vector of the Markov chain.}
  \item{LL}{is the log likelihood.}
}
\references{
%% ~put references to the literature/web site here ~
Wang, T., Zhuang, J., Buckby, J., Obara, K. and Tsuruoka, H. (2018) Identifying the recurrence patterns of non-volcanic tremors using a 2D hidden Markov model with extra zeros. Journal of Geophysical Research, \doi{10.1029/2017JB015360}.
}
\author{
Ting Wang
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
pie <- c(0.002,0.2,0.4)
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
yn <- hmm0norm2d(R, Z, pie, gamma, mu, sig, delta)
yn

# Setting up initial values when analysing real-world data
## nk is the number of states for the fitted model
### In this example we use nk=3

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
      eval(parse(text=paste('save(HMM',kk,'est, file="HMM',kk,'est.image")',sep='')))
    }
   }, error=function(e){})
 print(kk)
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
\keyword{methods}
