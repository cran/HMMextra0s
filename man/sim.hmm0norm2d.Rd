\name{sim.hmm0norm2d}
\alias{sim.hmm0norm2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of a Bivariate HMM with Extra Zeros
}
\description{
Simulates the observed process and the associated binary variable of a bivariate HMM with extra zeros.
}
\usage{
sim.hmm0norm2d(mu, sig, pie, gamma, delta, nsim = 1, mc.hist = NULL, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pie}{is a vector of length \eqn{m}, the \eqn{j}th element of which is the probability of \eqn{Z=1} when the process is in state \eqn{j}.}
  \item{gamma}{is the transition probability matrix (\eqn{m * m}) of the hidden Markov chain.}
  \item{mu}{is an \eqn{m * 2} matrix, the \eqn{j}th row of which is the mean of the bivariate (Gaussian) distribution of the observations in state \eqn{j}.}
  \item{sig}{is a \eqn{2 * 2 * m} array. The matrix \code{sig[,,j]} is the variance-covariance matrix of the bivariate (Gaussian) distribution of the observations in state \eqn{j}.}
  \item{delta}{is a vector of length \eqn{m}, the initial distribution vector of the Markov chain.}
  \item{nsim}{is an integer, the number of observations to simulate.}
  \item{mc.hist}{is a vector containing the history of the hidden Markov chain. This is mainly used for forecasting. If we fit an HMM to the data, and obtained the Viterbi path for the data, we can let \code{mc.hist} equal to the Viterbi path and then forecast futrue steps by simulation.}
  \item{seed}{is the seed for simulation. Default \code{seed=NULL}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{x}{is the simulated observed process.}
  \item{z}{is the simulated binary data with the value 1 indicating that an event was observed and 0 otherwise.}
  \item{mcy}{is the simulated hidden Markov chain.}
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
## Simulating a sequence of data without using any history.
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


## Forecast future tremor occurrences and locations when tremor occurs.
R <- y$x
Z <- y$z
HMMEST <- hmm0norm2d(R, Z, pie, gamma, mu, sig, delta)
Viterbi3 <- Viterbi.hmm0norm2d(R,Z,HMMEST)
y <- sim.hmm0norm2d(mu,sig,pie,gamma,delta,nsim=2,mc.hist=Viterbi3$y)
# This only forecasts two steps forward when we use nsim=2. 
# One can increase nsim to get longer simulated forecasts.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{simulation}
