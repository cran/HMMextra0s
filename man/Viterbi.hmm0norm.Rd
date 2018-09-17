\name{Viterbi.hmm0norm}
\alias{Viterbi.hmm0norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Viterbi Path of a 1-D HMM with Extra Zeros
}
\description{
Finds the most probable sequence of hidden states of an observed process.
}
\usage{
Viterbi.hmm0norm(R, Z, HMMest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{is the observed data. \code{R} is a \eqn{T * 1} matrix, where \eqn{T} is the number of observations.}
  \item{Z}{is the binary data with the value 1 indicating that an event was observed and 0 otherwise. \code{Z} is a vector of length \eqn{T}.}
  \item{HMMest}{is a list which contains pie, gamma, sig, mu, and delta (the HMM parameter estimates).}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{y}{is the estimated Viterbi path.} 
  \item{v}{is the estimated probability of each time point being in each state.}
}
\references{
%% ~put references to the literature/web site here ~
Wang, T., Zhuang, J., Obara, K. and Tsuruoka, H. (2016) Hidden Markov Modeling of Sparse Time Series from Non-volcanic Tremor Observations. Journal of the Royal Statistical Society, Series C, Applied Statistics, 66, Part 4, 691-715.
}
\author{
Ting Wang
}
\examples{
pie <- c(0.002,0.2,0.4)
gamma <- matrix(c(0.99,0.007,0.003,
                  0.02,0.97,0.01,
                  0.04,0.01,0.95),byrow=TRUE, nrow=3)
mu <- matrix(c(0.3,0.7,0.2),nrow=1)
sig <- matrix(c(0.2,0.1,0.1),nrow=1)
delta <- c(1,0,0)
y <- sim.hmm0norm(mu,sig,pie,gamma,delta, nsim=5000)
R <- as.matrix(y$x,ncol=1)
Z <- y$z
HMMEST <- hmm0norm(R, Z, pie, gamma, mu, sig, delta)
Viterbi3 <- Viterbi.hmm0norm(R,Z,HMMEST)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{inversion}

