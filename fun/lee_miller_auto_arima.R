
#### Lee-Miller (2001) with simulation
LM_AA <- function(E, LMX, Y, a, sex, t1, tF, nS){
    
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
    
    ## adjusting Kappa in the last period
    ## (it should be 0 for fitted e0 to match exactly the observed one)
    Kappa[n1] <- 0
    # plot(t1,Kappa)
    
    Kts <- ts(c(Kappa), start = t1[1])
    modK <- auto.arima(y=Kts)
    # plot(forecast(modK))
    predK <- forecast(modK, h = nF)
    # Simulation of kts with bootstrapping
    #SIMe0 <- c()
    
    SIMnmx <- c()
    E0 <- matrix(NA,nF,nS)
    for(s in 1:nS){
      kappa.sim <- simulate(modK, nsim=nF,future=TRUE, bootstrap=TRUE)
      kappa.matrix <- matrix(kappa.sim)
      
      OneF <- rep(1,nF)
      nmx <- Alpha%*%t(OneF) + Beta%*%t(kappa.matrix)
      
      ## life exp
      E0[,s] <- apply(exp(nmx),2,lifetable.e0,x=a,sex=sex)
      
      ## Calculate le by horizon
      #e0 <- apply(exp(nmx),2,e0.mx,x=a,sex=sex)
      #SIMe0 <- rbind(SIMe0, e0) # Saving the e0
      
      SIMnmx <- rbind(SIMnmx, nmx) # Saving nmx
      colnames(SIMnmx) <- tF
    }
    
    # matplot(tF,E0,t="l",lty=1)
    
    ## median e0
    e0 <- apply(E0,1,median)
    # plot(tF,e0)
    
    ## theoretical e0 (non simulation based)
    nmx <- Alpha%*%t(OneF) + Beta%*%t(predK$mean)
    e0.med <- apply(exp(nmx),2,lifetable.e0,x=a,sex=sex)
    
    ## output
    out <- list(SIMnmx=SIMnmx,e0=e0,e0.med=e0.med)
    
    return(out)
    
}


# Example of application
#LMX_fore <- LM(E=E, LMX=LMX, Y=D, a=a, sex="M", t1=t1, tF=tF, nS=nS)













