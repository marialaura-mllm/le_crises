
#### Lee-Miller (2001) with simulation
LM <- function(E, LMX, Y, a, sex, t1, tF, nS){
    
    n1 <- length(t1)
    nF <- length(tF)
    
    Alpha <- LMX[,ncol(LMX)] #### Last year
    # Alpha <- apply(LMX,1,mean)
    
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
    
    for (i in 1:n1){
      KappaSecStep <- optimize(f=koptim, interval=c(-150,150), Alpha=Alpha, 
                               Beta=Beta, e0.obs=e0.obs.real[i], x=a, sex=sex)
      Kappa[i] <- KappaSecStep$minimum
    }
    
    ## adjusting Kappa in the last period
    ## (it should be 0 for fitted e0 to match exactly the observed one)
    Kappa[n1] <- 0
    
    # sum(Kappa)
    Kts <- ts(c(Kappa), start = t1[1])
    modK <- Arima(Kts, order=c(0,1,0), include.drift=TRUE)
    
    
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
    
    ## median e0
    e0 <- apply(E0,1,median)
    
    ## output
    out <- list(SIMnmx=SIMnmx,Beta=Beta,Kappa=Kappa,e0=e0)
    
    return(out)
    
}

# Adjusting kt for the e0
koptim <- function(par,Alpha,Beta,e0.obs,x,sex){
  Kappa <- par[1]
  lmx.lc <- Alpha+Beta*Kappa
  e0.est <- lifetable.e0(mx=exp(lmx.lc),x=a,sex=sex) #e0 adjustment
  diff.lc <- abs(e0.obs-e0.est)
  return(diff.lc)
}


#### Li-Lee (2005) with simulation
## here for females and males
## with 1=Females, 2 = Males
## here modified as to avoid jump-off bias
LiLee <- function(LMX1, LMX2, a, t1, tF, nS){
  
  n1 <- length(t1)
  nF <- length(tF)
  m <- length(a)
  
  ## common group
  ## average of death rates
  LMX <- log((exp(LMX1)+exp(LMX2))/2)
  ## life expectancy for the whole group
  e0.group <- apply(exp(LMX),2,lifetable.e0,x=a,sex="Female")
  
  ## extracting common factor: fit LM to the whole group
  fitGroup <- LM(LMX=LMX,a=a,t1=t1,tF=tF,nS=nS,sex = "Female")
  Beta <- fitGroup$Beta
  Kappa <- fitGroup$Kappa
  
  ##  group specific alpha_x 
  # Alpha1 <- apply(LMX1,1,mean)
  # Alpha2 <- apply(LMX2,1,mean)
  Alpha1 <- LMX1[,ncol(LMX)]
  Alpha2 <- LMX2[,ncol(LMX)]
  
  ##  group specific beta_x kappa_t
  LMXres1 <- LMX1 - Alpha1 - Beta%*%t(Kappa)
  LMXres2 <- LMX2 - Alpha2 - Beta%*%t(Kappa)
  ## group 1
  LCsvd1 <- svd(LMXres1,nu=1,nv=1)
  Beta1 <- c(LCsvd1$u)
  Kappa1 <- c(LCsvd1$v)
  ## group 2
  LCsvd2 <- svd(LMXres2,nu=1,nv=1)
  Beta2 <- c(LCsvd2$u)
  Kappa2 <- c(LCsvd2$v)
  
  ## forecasting the common Kappa
  Kts <- ts(c(Kappa), start = t1[1])
  modK <- auto.arima(y=Kts)
  # predK <- forecast(modK,h=nF)
  # k.point.foreMEAN <- matrix(predK$mean, nrow=nF, ncol=1)
  ## forecasting population-specific Kappas
  ## group 1
  K1ts <- ts(c(Kappa1), start = t1[1])
  modK1 <- Arima(K1ts, order=c(1,1,0), include.drift=FALSE)
  # predK1 <- forecast(modK1,h=nF)
  # k1.point.foreMEAN <- matrix(predK1$mean, nrow=nF, ncol=1)
  ## group 2
  K2ts <- ts(c(Kappa2), start = t1[1])
  modK2 <- Arima(K2ts, order=c(1,1,0), include.drift=FALSE)
  # predK2 <- forecast(modK2,h=nF)
  # k2.point.foreMEAN <- matrix(predK2$mean, nrow=nF, ncol=1)
  
  SIMnmx1 <- SIMnmx2 <- c()
  E01 <- E02 <- matrix(NA,nF,nS)
  for(s in 1:nS){
    ## common kappa
    kappa.sim <- simulate(modK, nsim=nF,future=TRUE, bootstrap=TRUE)
    kappa.matrix <- matrix(kappa.sim)
    ## kappa1
    kappa1.sim <- simulate(modK1, nsim=nF,future=TRUE, bootstrap=TRUE)
    kappa1.matrix <- matrix(kappa1.sim)
    ## kappa2
    kappa2.sim <- simulate(modK2, nsim=nF,future=TRUE, bootstrap=TRUE)
    kappa2.matrix <- matrix(kappa2.sim)
    
    ## computing forecast rates for the two populations
    nmx1 <- matrix(Alpha1,m,nF) + Beta%*%t(kappa.matrix) + Beta1%*%t(kappa1.matrix)
    nmx2 <- matrix(Alpha2,m,nF) + Beta%*%t(kappa.matrix) + Beta2%*%t(kappa2.matrix)
    ## computing corresponding life expectancy
    e0.fore1 <- apply(exp(nmx1),2,lifetable.e0,x=a,sex="Female")
    e0.fore2 <- apply(exp(nmx2),2,lifetable.e0,x=a,sex="Male")
    
    ## Calculate le by horizon
    #e0 <- apply(exp(nmx),2,e0.mx,x=a,sex=sex)
    #SIMe0 <- rbind(SIMe0, e0) # Saving the e0
    
    SIMnmx1 <- rbind(SIMnmx1, nmx1) # Saving nmx1
    SIMnmx2 <- rbind(SIMnmx2, nmx2) # Saving nmx1
    colnames(SIMnmx1) <- colnames(SIMnmx2) <- tF
    
    E01[,s] <- e0.fore1
    E02[,s] <- e0.fore2
    
  }
  
  ## median e0
  e01 <- apply(E01,1,median)
  e02 <- apply(E02,1,median)
  

  
  ## output
  out <- list(SIMnmx1=SIMnmx1,SIMnmx2=SIMnmx2,E01=E01,E02=E02,
              e01=e01,e02=e02)
  
  return(out)
  
}












