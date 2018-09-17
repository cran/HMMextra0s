\name{hmm0norm}
\alias{hmm0norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter Estimation of an HMM with Extra Zeros
}
\description{
Calculates the parameter estimates of a 1-D HMM with observations having extra zeros.
}
\usage{
hmm0norm(R, Z, pie, gamma, mu, sig, delta, tol=1e-6, print.level=1, fortran = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{is the observed data. \code{R} is a \eqn{T * 1} matrix, where \eqn{T} is the number of observations.}
  \item{Z}{is the binary data with the value 1 indicating that an event was observed and 0 otherwise. \code{Z} is a vector of length \eqn{T}.}
  \item{pie}{is a vector of length \eqn{m}, the \eqn{j}th element of which is the probability of \eqn{Z=1} when the process is in state \eqn{j}.}
  \item{gamma}{is the transition probability matrix (\eqn{m * m}) of the hidden Markov chain.}
  \item{mu}{is a \eqn{1 * m} matrix, the \eqn{j}th element of which is the mean of the (Gaussian) distribution of the observations in state \eqn{j}.}
  \item{sig}{is a \eqn{1 * m} matrix, the \eqn{j}th element of which is the standard deviation of the (Gaussian) distribution of the observations in state \eqn{j}.}
  \item{delta}{is a vector of length \eqn{m}, the initial distribution vector of the Markov chain.}
  \item{tol}{is the tolerance for testing convergence of the iterative estimation process. The default tolerance is 1e-6. For initial test of model fit to your data, a larger tolerance (e.g., 1e-3) should be used to save time.}
  \item{print.level}{controls the amount of output being printed. Default is 1. If \code{print.level=1}, only the log likelihoods and the differences between the log likelihoods at each step of the iterative estimation process, and the final estimates are printed. If \code{print.level=2}, the log likelihoods, the differences between the log likelihoods, and the estimates at each step of the iterative estimation process are printed.}
  \item{fortran}{is logical, and determines whether Fortran code is used; default is \code{TRUE}.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  \item{pie}{is the estimated probability of \eqn{Z=1} when the process is in each state.}
  \item{mu}{is the estimated mean of the (Gaussian) distribution of the observations in each state.}
  \item{sig}{is the estimated standard deviation of the (Gaussian) distribution of the observations in each state.}
  \item{gamma}{is the estimated transition probability matrix of the hidden Markov chain.}
  \item{delta}{is the estimated initial distribution vector of the Markov chain.}
  \item{LL}{is the log likelihood.}
}
\references{
%% ~put references to the literature/web site here ~
Wang, T., Zhuang, J., Obara, K. and Tsuruoka, H. (2016) Hidden Markov Modeling of Sparse Time Series from Non-volcanic Tremor Observations. Journal of the Royal Statistical Society, Series C, Applied Statistics, 66, Part 4, 691-715.
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
mu <- matrix(c(0.3,0.7,0.2),nrow=1)
sig <- matrix(c(0.2,0.1,0.1),nrow=1)
delta <- c(1,0,0)
y <- sim.hmm0norm(mu,sig,pie,gamma,delta, nsim=5000)
R <- as.matrix(y$x,ncol=1)
Z <- y$z
yn <- hmm0norm(R, Z, pie, gamma, mu, sig, delta)
yn
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
\keyword{methods}
