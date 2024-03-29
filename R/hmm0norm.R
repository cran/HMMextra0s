hmm0norm <-
function( R, Z, pie, gamma, mu, sig, delta, tol=1e-6, print.level=1, fortran = TRUE)
{
  n <- nrow( sig )
  m <- ncol( sig )
  nn <- nrow( R )
  count = 1
  repeat
  {
#
    pRS <- matrix( 1, nn, m )
    if (fortran!=TRUE){
        for (k in 1:m)
        {
          pRS[,k] <- (pie[k] * exp( - (R[,1] - mu[,k])^2 / (2 * sig[,k]^2) ) /
                       (sqrt(2 * pi) * sig[,k]) )*(Z) + (1-pie[k])*(1-Z)
        }
    } else {
        prsloop <- .Fortran("prsloop", m, nn, pie, R[,1], mu[1,], sig[1,], Z,
                            pRS, PACKAGE="HMMextra0s")
        pRS <- prsloop[[8]]
    }
#
##Scaled forward variable
    logalpha <- matrix(as.double(rep(0, nn * m)), nrow = nn)
    lscale <- as.double(0)
    phi <- as.double(delta)
    gamma <- matrix(as.double(gamma),nrow=m)
    pRS <- matrix(as.double(pRS),nrow=nn)
    if (fortran!=TRUE){
        #  loop1 using R code
      for (t in 1:nn)
      {
        if (t > 1) 
           phi <- phi %*% gamma
        phi <- phi * pRS[t,]
        sumphi <- sum(phi)
        phi <- phi/sumphi
        lscale <- lscale + log(sumphi)
        logalpha[t, ] <- log(phi) + lscale
      }
	  if (count > 1.5){
        LLn = lscale
        diffL = LLn - LL 
        print(format(LL,digits=12))
        print(format(LLn,digits=12))
        print(format(diffL,digits=12))
        if (LLn < LL) stop ('worse likelihood')
        if (diffL <= tol) break
      } 
      LL <- lscale
    }else{
        if (!is.double(gamma)) stop("gamma is not double precision")
        memory0 <- rep(as.double(0), m)
        loop1 <- .Fortran("loop1", m, nn, phi, pRS, gamma, logalpha,
                          lscale, memory0, PACKAGE="HMMextra0s")
        logalpha <- loop1[[6]]
      if (count > 1.5){
        LLn = loop1[[7]]
        diffL = LLn - LL 
        print(format(LL,digits=12))
        print(format(LLn,digits=12))
        print(format(diffL,digits=12))
        if (LLn < LL) stop ('worse likelihood')
        if (diffL <= tol) break
      } 
        LL <- loop1[[7]]
    }
#
##Scaled backward variable
    logbeta <- matrix(as.double(rep(0, nn * m)), nrow = nn)
    phi <- as.double(rep(1/m, m))
    lscale <- as.double(log(m))
    if (fortran!=TRUE){
        #  loop2 using R code
      for (t in seq(nn - 1, 1, -1)) 
      {
         phi <- gamma %*% (pRS[t+1,] * phi)
         logbeta[t, ] <- log(phi) + lscale
         sumphi <- sum(phi)
         phi <- phi/sumphi
         lscale <- lscale + log(sumphi)
      }
    } else{
        memory0 <- rep(as.double(0), m)
        loop2 <- .Fortran("loop2", m, nn, phi, pRS, gamma, logbeta,
                          lscale, memory0, PACKAGE="HMMextra0s")
        logbeta <- loop2[[6]]
    }
#
##E-step
    if (fortran!=TRUE) {
###Calculate v_t(j)
        v <- exp(logalpha + logbeta - LL)
###Calculate w_t(i,j) 
        w <- array( NA, c( nn-1, m, m ) )
        for (k in 1:m) {
            logprob <- log( pRS[-1, k] )
            logPi <- matrix(log(gamma[, k]), byrow = TRUE, nrow = nn - 
                1, ncol = m)
            logPbeta <- matrix(logprob + logbeta[-1, k],
                byrow = FALSE, nrow = nn - 1, ncol = m)
            w[, , k] <- logPi + logalpha[-nn, ] + logPbeta - LL
        }
        w <- exp(w)
    } else {
        v <- array(0.0, c( nn, m ))
        w <- array(0.0, c( nn-1, m, m ) )
# logalpha can contain -Inf values where phi=0; set NAOK=TRUE as the Fortran code
# will handle such cases safely
        estep <- .Fortran("estep", m, nn, logalpha, logbeta, LL,
                           pRS, gamma, v, w, NAOK=TRUE, PACKAGE="HMMextra0s")
        v <- estep[[8]]
        w <- estep[[9]]
    }

#
##M-step
###Estimate pie_{i}, mu_{ij} and sigma_{ij}
    hatpie <- rep(0, m)
    hatmu <- matrix( 0, n, m )
    hatsig <- matrix( 0, n, m )
    if (fortran!=TRUE) {
        for (j in 1:m) {
            hatpie[j] <- ( v[,j] %*% Z ) / sum( v[,j] )
            hatmu[,j] <- ( v[Z==1,j] %*% R[Z==1,] ) / sum( v[Z==1,j] )
            hatsig[,j] <- sqrt( ( v[Z==1,j] %*% ( t(t(R[Z==1,]) - hatmu[,j]) )^2 ) /
                                sum( v[Z==1,j] ) )
        }
    } else {
        mstep1d <- .Fortran("mstep1d", n, m, nn, v, Z, R,
                            hatpie, hatmu, hatsig,
                            PACKAGE="HMMextra0s")
        hatpie <- mstep1d[[7]]
        hatmu <- mstep1d[[8]]
        hatsig <- mstep1d[[9]]
    }

###Estimate gamma_{ij}
#
    hatgamma <- colSums( w ) / replicate(m, colSums( v ) - v[nn,])
#
###Estimate delta
    hatdelta <- v[1,]
#
    pie <- hatpie
    mu <- hatmu
    sig <- hatsig
    gamma <- hatgamma
    delta <- hatdelta
    count = count + 1    
#
    if (print.level==2){
      print(pie)
      print(mu)
      print(sig)
      print(gamma)
      print(delta)
	}
  }  
#  write(format(t(nvtj)),'vtj.txt',ncolumns=m)
  print(pie)
  print(mu)
  print(sig)
  print(gamma)
  print(delta)
  return(list(pie=pie,mu=mu, sig=sig, gamma=gamma, delta=delta, LL=LLn))
}
