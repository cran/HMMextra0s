##28/09/2015
###########################################################
##################   Viterbi algorithm   ##################
###########################################################
## HMMest is a list which contains pie, gamma, sig, mu, delta (the HMM parameter estimates)
## y is the estimated hidden states 
## v is the estimated probability of each time point being in each state
#R is the tremor location data.
###R is a nn*n matrix. If R is nn*1, should use 
###R <- as.matrix(R) to make R as a matrix.
###Z is the binary data indicating if tremors are observed
Viterbi.hmm0norm2d <- function(R, Z, HMMest){
  pie <- HMMest$pie
  mu <- HMMest$mu
  sig <- HMMest$sig
  n <- ncol( mu )
  m <- nrow( mu )
  nn <- nrow( R )
  pRS <- matrix( 1, nn, m )
  for (k in 1:m)
  {
    pRS[,k] <- ( pie[k] * dmvnorm(R, mean = mu[k,], sigma = sig[,,k]) )^Z * (1-pie[k])^(1-Z)
  }
  muh <- mu
  sigh <- sig
  gammah <- HMMest$gamma
  deltah <- HMMest$delta
  nu <- matrix(NA, nrow = nn, ncol = m)
  y <- rep(NA, nn)
  nu[1, ] <- log(deltah) + log(pRS[1,])
  logPi <- log(gammah)
  for (i in 2:nn) {
      matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
      nu[i, ] <- apply(matrixnu + logPi, 2, max) + log(pRS[i,])
  }
  if (any(nu[nn, ] == -Inf)) 
      stop("Problems With Underflow")
  y[nn] <- which.max(nu[nn, ])
  for (i in seq(nn - 1, 1, -1)) y[i] <- which.max(logPi[, y[i + 1]] + nu[i, ])
  logalpha <- matrix(as.double(rep(0, nn * m)), nrow = nn)
  lscale <- as.double(0)
  phi <- as.double(deltah)
  for (t in 1:nn)
  {
    if (t > 1) 
       phi <- phi %*% gammah
    phi <- phi * pRS[t,]
    sumphi <- sum(phi)
    phi <- phi/sumphi
    lscale <- lscale + log(sumphi)
    logalpha[t, ] <- log(phi) + lscale
  }
  LL <- lscale
#
##Scaled backward variable
  logbeta <- matrix(as.double(rep(0, nn * m)), nrow = nn)
  phi <- as.double(rep(1/m, m))
  lscale <- as.double(log(m))
  for (t in seq(nn - 1, 1, -1)) 
  {
     phi <- gammah %*% (pRS[t+1,] * phi)
     logbeta[t, ] <- log(phi) + lscale
     sumphi <- sum(phi)
     phi <- phi/sumphi
     lscale <- lscale + log(sumphi)
  }
#
##E-step
###Calculate v_t(j)
  v <- exp(logalpha + logbeta - LL)
## y is the estimated hidden states 
## v is the estimated probability of each time point being in each state
  return(list(y=y,v=v))
}


