
################################################################################
# Function created to calculate measures of life table
# This is a generic function that can be used for different age groups' sizes
# x is a vector for the age groups and mx is a vector for the death rates 
lifetable.e0 <- function(x, mx, sex="Male", ax=NULL){
  m <- length(x)
  n <- c(diff(x), NA)
  # Considering different values for ax considering the sex and first age groups
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="Female"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="Male"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  e0 <- ex[1]
  return(e0) # In this case it will return only the life expectancy at birth
}



################################################################################
# Customized function running a linear interpolation on the log death rates
lin_inter <- function(data){
  
  # Filter for the interpolation period+year before and after for interpolation
  data_interp <- data %>% 
    filter(year %in% 2019:(table_interp[country==countries[c],]$last_year+1))
  
  lmx <- data_interp$lmx
  n <- length(lmx)
  
  # Changing values of NA(interp) to linear interpolation obtained with extremes
  data$lmx[is.na(data$lmx)] <- 
    cumsum(c(lmx[1L], rep((lmx[n]-lmx[1L])/(n - 1L), n - 1L)))[2:(n-1)]
  
  data
  
}



################################################################################
# Function to change period lmx to cohort using Schmertmann (2024)
# Needs matrix with lmx (age in rows) by year (column) as input
# a = age groups

period_cohort <- function(lmx, a, first_cohort, last_cohort){
  
  lmx_cohort <- matrix(nrow=length(a), 
                       ncol=length(seq(first_cohort,last_cohort)))
  colnames(lmx_cohort) <- seq(first_cohort, last_cohort)
  
  for (i in 1:(ncol(lmx_cohort))){ #For each cohort
    
    # Each time excludes the year before to take the values on the diagonal
    aux_1 <- lmx[,i:ncol(lmx)]
    aux_2 <- lmx[,(i+1):ncol(lmx)]
    
    # Taking the mean between consecutive years as in Schmertmann (2024)
    lmx_cohort[,i] <- (diag(aux_1) + diag(aux_2))/2
  }
  
  lmx_cohort <- lmx_cohort %>% as.data.frame()
  return(lmx_cohort)
  
}



################################################################################













