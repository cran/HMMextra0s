\name{cumdist.hmm0norm}
\alias{cumdist.hmm0norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cumulative distribution of an HMM with Extra Zeros
}
\description{
Calculates the cumulative distribution of an HMM with extra zeros.
}
\usage{
cumdist.hmm0norm(x,HMMest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{x is a value at which the cumulative distribution is evaluated.}
  \item{HMMest}{is a list which contains pie, gamma, sig, mu, and delta (the HMM parameter estimates).}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{prob}{is the calculated cumulative distribution.} 
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Wang, T., Zhuang, J., Obara, K. and Tsuruoka, H. (2016) Hidden Markov Modeling of Sparse Time Series from Non-volcanic Tremor Observations. Journal of the Royal Statistical Society, Series C, Applied Statistics, 66, Part 4, 691-715.
}
\author{
Ting Wang
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
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
xx <- seq(0,1,0.05)
cumdist <- apply(t(xx),2,cumdist.hmm0norm,HMMest=HMMEST)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ goodness-of-fit }

