\name{HMMextra0s-package}
\alias{HMMextra0s-package}
\alias{HMMextra0s}
\docType{package}
\title{
\packageTitle{HMMextra0s} Hidden Markov Models (HMMs) with Extra Zeros
}
\description{
The DESCRIPTION file:
\packageDESCRIPTION{HMMextra0s}
\packageIndices{HMMextra0s}
This package contains functions to estimate the parameters of the HMMs with extra zeros using \code{\link{hmm0norm}} (1-D HMM) 
and \code{\link{hmm0norm2d}} (2-D HMM), to calculate the cumulative distribution of the 1-D HMM using \code{\link{cumdist.hmm0norm}}, 
to estimate the Viterbi path using \code{\link{Viterbi.hmm0norm}} (1-D HMM) and \code{\link{Viterbi.hmm0norm2d}} (2-D HMM), to 
simulate this class of models using \code{\link{sim.hmm0norm}} (1-D HMM) and \code{\link{sim.hmm0norm2d}} (2-D HMM), to plot the 
classified 2-D data with different colours representing different hidden states using \code{\link{plotVitloc2d}}, and to plot the 
Viterbi path using \code{\link{plotVitloc2d}}.
}
\details{
This package is used to estimate the parameters, carry out simulations, and estimate the Viterbi path for 1-D and 2-D HMMs with extra 
zeros as defined in the two publications in the reference (also briefly defined below). It contains examples using simulated data for 
how to set up initial values for a data analysis and how to plot the results.

An HMM is a statistical model in which the observed process is dependent on an unobserved Markov chain. A Markov chain is a sequence 
of states which exhibits a short-memory property such that the current state of the chain is dependent only on the previous state in 
the case of a first-order Markov chain. Assume that the Markov chain has \eqn{m} states, where \eqn{m} can be estimated from the data. 
Let \eqn{S_t \in \{1,\cdots,m\}}{S_t in {1,...,m}} denote the state of the Markov chain at time \eqn{t}. The probability of a first-order Markov chain 
in state \eqn{j} at time \eqn{t} given the previous states is \eqn{P(S_t=j|S_{t-1},\cdots,S_{1})=P(S_t=j|S_{t-1})}{P(S_t=j|S_{t-1},...,S_1)=P(S_t=j|S_{t-1})}. 
These states are not observable. The observation \eqn{Y_t} at time \eqn{t} depends on the state \eqn{S_t} of the Markov chain. 

In this framework, we are interested in estimating the transition probability matrix \eqn{\Gamma=(\gamma_{ij})_{m\times m}}{Gamma=(gamma_{ij})} of the 
Markov chain that describes the migration pattern and the density function \eqn{f(y_t|S_t=i)} that gives the distribution feature of 
observations in state \eqn{i}, where \eqn{\gamma_{ij}=P(S_t=j|S_{t-1}=i)}{gamma_{ij}=P(S_t=j|S_{t-1}=i)}.

Let \eqn{Z_t} be a Bernoulli variable, with \eqn{Z_t=1} if an event is present at \eqn{t}, and \eqn{Z_t=0}, otherwise.
Let \eqn{\mathbf{X}_t}{X_t} be the response variable (e.g., location of the tremor cluster in 2D space) at time \eqn{t}. 
We set \eqn{P(Z_t=0|S_t=i)=1-p_i} and \eqn{P(Z_t=1|S_t=i)=p_i}. We assume that, given \eqn{Z_t=1} and \eqn{S_t=i}, 
\eqn{\mathbf{X}_t}{X_t} follows a univariate or bivariate normal distribution, e.g. for a bivariate normal,
\deqn{
f(\mathbf{x}_t|Z_t=1, S_t=i)=\frac{1}{2\pi |\bm{\Sigma}_i|^{1/2}}\exp\left(-\frac{1}{2}
(\mathbf{x}_t-\bm{\mu}_i)^T\bm{\Sigma}_i^{-1}(\mathbf{x}_t-\bm{\mu}_i)\right).
}{f(x_t | Z_t=1,S_t=i)=1/(2 pi|Sigma_i|^(1/2)) exp[-(x_t-mu_i)^T Sigma_i^(-1)(x_t-mu_i)/2].}
The joint probability density function of \eqn{Z_t} and \eqn{\mathbf{X}_t}{X_t} conditional on the system being in state \eqn{i} at time \eqn{t} is
\deqn{
f(\mathbf{x}_t,z_t | S_t=i)=(1-p_i)^{1-z_t}\left[p_i\frac{1}{2\pi |\bm{\Sigma}_i|^{1/2}}\exp\left(-\frac{1}{2}(\mathbf{x}_t-\bm{\mu}_i)^T\bm{\Sigma}_i^{-1}(\mathbf{x}_t-\bm{\mu}_i)\right)\right]^{z_t},
}{f(x_t,z_t=1|S_t=i)=(1-p_i)^(1-z_t)[p_i/(2 pi|Sigma_i|^(1/2)) exp[-(x_t-mu_i)^T Sigma_i^(-1)(x_t-mu_i)/2]]^(z_t).}
where \eqn{p_i}, \eqn{\bm{\mu}_i=E(\mathbf{X}_t|S_t=i,Z_t=1)}{mu_i=E(X_t|S_t=i,Z_t=1)} and \eqn{\bm{\Sigma}_i=Var(\mathbf{X}_t|S_t=i,Z_t=1)}{Sigma_i=Var(X_t|S_t=i,Z_t=1)} are parameters to be 
estimated.
}
\author{
\packageAuthor{HMMextra0s}

Maintainer: \packageMaintainer{HMMextra0s}
}
\references{
Wang, T., Zhuang, J., Obara, K. and Tsuruoka, H. (2016) Hidden Markov Modeling of Sparse Time Series from Non-volcanic Tremor Observations. Journal of the Royal Statistical Society, Series C, Applied Statistics, 66, Part 4, 691-715.

Wang, T., Zhuang, J., Buckby, J., Obara, K. and Tsuruoka, H. (2018) Identifying the recurrence patterns of non-volcanic tremors using a 2D hidden Markov model with extra zeros. Journal of Geophysical Research, \doi{10.1029/2017JB015360}.

Some of the functions in the package are based on those of the R package ``HiddenMarkov":

Harte, D. (2017) HiddenMarkov: Hidden Markov Models. R package version 1.8-11. Statistics Research Associates, Wellington. URL: \url{http://www.statsresearch.co.nz/dsh/sslib/}
}
\keyword{ package }



