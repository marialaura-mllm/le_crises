
#### Lee-Miller (2001) with simulation
LM <- function(E, LMX, Y, a, sex, t1, tF, nS){
    
    n1 <- length(t1)
    nF <- length(tF)
    
    Alpha <- LMX[,ncol(LMX)] #### Last year
    
    M.tilde <- LMX - Alpha
    
    svdM <- svd(M.tilde,nu=1,nv=1)
    
    Beta <- c(svdM$u)
    Kappa1 <- c(svdM$v)
    
    sum.Beta <- sum(Beta)
    Beta <- Beta/sum.Beta
    
    Kappa1 <- Kappa1 - mean(Kappa1)
    Kappa1 <- Kappa1*sum.Beta
    
    Kappa <- rep(NA,n1)
    
    e0.obs.real <- apply(exp(LMX),2,lifetable.e0,x=a,sex=sex)
    
    # Adjusting kt for the e0
    koptim <- function(par,Alpha,Beta,e0.obs,x,sex){
      Kappa <- par[1]
      lmx.lc <- Alpha+Beta*Kappa
      e0.est <- lifetable.e0(mx=exp(lmx.lc),x=a,sex=sex) #e0 adjustment
      diff.lc <- abs(e0.obs-e0.est)
      return(diff.lc)
    }
    
    for (i in 1:n1){
      KappaSecStep <- optimize(f=koptim, interval=c(-150,150), Alpha=Alpha, 
                               Beta=Beta, e0.obs=e0.obs.real[i], x=a, sex=sex)
      Kappa[i] <- KappaSecStep$minimum
    }
    
    Kts <- ts(c(Kappa), start = t1[1])
    modK <- Arima(Kts, order=c(0,1,0), include.drift=TRUE)
    
    
    # Simulation of kts with bootstrapping
    #SIMe0 <- c()
    
    SIMnmx <- c()
    for(s in 1:nS){
      kappa.sim <- simulate(modK, nsim=nF,future=TRUE, bootstrap=TRUE)
      kappa.matrix <- matrix(kappa.sim)
      
      OneF <- rep(1,nF)
      nmx <- Alpha%*%t(OneF) + Beta%*%t(kappa.matrix)
      
      ## Calculate le by horizon
      #e0 <- apply(exp(nmx),2,e0.mx,x=a,sex=sex)
      #SIMe0 <- rbind(SIMe0, e0) # Saving the e0
      
      SIMnmx <- rbind(SIMnmx, nmx) # Saving nmx
      colnames(SIMnmx) <- tF
    }
    
    return(SIMnmx)
    
}


# Example of application
#LMX_fore <- LM(E=E, LMX=LMX, Y=D, a=a, sex="M", t1=t1, tF=tF, nS=nS)













