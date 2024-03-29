sim.hmm0norm <-
function(mu,sig,pie,gamma,delta, nsim=1, seed=NULL){
    if (!is.null(seed)) set.seed(seed)
    m <- ncol(gamma)
    if (sum(delta)!=1) stop("Invalid delta")
    if (any(delta==1))
        initial <- (1:m)[as.logical(delta)]
    else
        initial <- sample(m, 1, prob=delta)
    xx <- rep(NA, nsim)
    xx[1] <- initial
	if (nsim>1){
      for (i in 2:nsim)
        xx[i] <- sample(x=1:m, size=1, prob=gamma[(xx[i-1]),])
    }
	mcy <- xx
    x <- z <- rep(NA, nsim)
    U <- runif(nsim)
    for (i in 1:nsim) {
      if (U[i]<pie[mcy[i]]){
        x[i] <- rnorm(1,mu[mcy[i]],sig[mcy[i]])
        z[i] <- 1
      }else{
        x[i] <- 0
        z[i] <- 0
      }
    }
    return(list(x=x,z=z,mcy=mcy))
}
