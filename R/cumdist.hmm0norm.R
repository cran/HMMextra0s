## pai0 is the stationary distribution of gamma
## x is a vector of values at which the cumulative distribution is evaluated.
cumdist.hmm0norm <- function(x,HMMest){
   pie <- HMMest$pie
   mu <- HMMest$mu
   sig <- HMMest$sig
   gamma <- HMMest$gamma
   pitemp <- eigen( t(gamma) )$vectors
   pai0 <- Re(pitemp[,1]/sum(pitemp[,1]))
   m <- ncol(mu)
   prob <- 0
   for (k in 1:m)
   {
     prob <- prob + pai0[k]*(pie[k] * pnorm(x,mu[k],sig[k]) + (1-pie[k]))
   }	
   return(prob)
}
