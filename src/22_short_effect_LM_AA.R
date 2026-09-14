
#'Evaluating the Impact of Mortality Crises on Life Expectancy: 
#Period and Cohort Perspectives'
# Short-term effect - Lee-Miller Auto ARIMA

# This code is to be used to replicate the results of the paper "Evaluating the 
#Impact of Mortality Crises on Life Expectancy: Period and Cohort Perspectives". 
#Specifically, it can be used to obtain cohort life expectancy at birth 
#estimates considering a short-term disturbance using as example the COVID-19 
#pandemic.



#-------------------------------------------------------------------------------

# Cleaning the workspace
rm(list=ls(all=TRUE))

# Packages
require(tidyverse)
require(forecast)

# Functions
source("fun/functions.R")
# List of functions includes 'lifetable.e0' for calculating life expectancy at 
#birth and 'lin_inter' to obtain linear interpolation between two years.

# Lee-Carter functions
source("fun/lee_miller_auto_arima.R")
# Function to apply the Lee-Miller (2001) with auto arima.



#-------------------------------------------------------------------------------
# Reading data from WPP already cleaned and filtered by country

mx1dt_wpp <- 
  read.table("dat/01_wpp_cleaned.txt", header = T)

# mx1dt_wpp: Age-specific death rates for selected countries with WPP data.
#Years from 1950 to 2023.
# Variables: "country", "year", "age", "sex", "mx", "death" and "exp"



#-------------------------------------------------------------------------------
# Obtaining cohort life expectancy for all countries using the Lee-Miller AA
# Short-term disturbance (2020-2023)


summary_e0 <- tibble()
summary_lmx <- tibble()
period_e0 <- tibble()

filter_country <- mx1dt_wpp %>% distinct(country) %>% pull()

filter_sex <- c("Male", "Female")

# Loop for all countries
set.seed(1234)
for (c in 1:length(filter_country)) {
  cat("analysing",filter_country[c],'\n')
  a <- as.integer(seq(0,100,1))
  
  aux_all <- 
    mx1dt_wpp %>% 
    filter(country == filter_country[c])
  
  # Exposures
  aux_E <-
    aux_all %>% 
    select(-c(country, death, mx)) %>% 
    pivot_wider(names_from = year, values_from = exp) %>% 
    arrange(sex, age)
  
  # Deaths
  aux_D <-
    aux_all %>% 
    select(-c(country, exp, mx)) %>% 
    pivot_wider(names_from = year, values_from = death) %>% 
    arrange(sex, age)
  
  # Death rates
  aux_MX <-
    aux_all %>% 
    select(-c(country, exp, death)) %>% 
    pivot_wider(names_from = year, values_from = mx) %>% 
    arrange(sex, age)
  
  # Initial cohort
  aux_ini_cohort <- aux_all %>%  select(year) %>%  min
  
  
  # Loop by sex
  for (s in 1:length(filter_sex)) {
    
    #########################       1950-2019       ############################
    
    E <- as.matrix(aux_E %>% 
                     filter(sex == filter_sex[s]) %>% 
                     select(-c(age, sex, `2020`, `2021`, `2022`, `2023`)))
    D <- as.matrix(aux_D %>%
                     filter(sex == filter_sex[s]) %>% 
                     select(-c(age, sex, `2020`, `2021`, `2022`, `2023`)))
    LMX <- log(as.matrix(aux_MX %>% 
                           filter(sex == filter_sex[s]) %>% 
                           select(-c(age, sex,`2020`,`2021`,`2022`,`2023`))))
    
    e0.obs <- apply(exp(LMX),2,lifetable.e0,x=a,sex=filter_sex[s])
    
    t1 <- as.integer(seq(1950,2019,1)) # Baseline period
    tF <- as.integer(seq(2020,2101,1)) # Forecast horizon
    nS <- 100 # Number of simulations
    
    # Lee-Miller AA function from "lc_functions.R"
    LMX_fore <- LM_AA(E=E, LMX=LMX, Y=D, a=a, sex=filter_sex[s],
                      t1=t1, tF=tF, nS=nS)
    
    e0.fore <- LMX_fore$e0.med
    # plot(t1,e0.obs,ylim=range(e0.obs,e0.fore),xlim=range(t1,tF))
    # lines(tF,e0.fore)
    # lines(tF,LMX_fore$e0,col=2)
    
    ## saving period forecasts
    period_e0_obs <- 
      tibble(country = filter_country[c],
             sex = filter_sex[s],
             year = t1,
             e0 = e0.obs,LM=NA)
    period_e0_fore <- 
      tibble(country = filter_country[c],
             sex = filter_sex[s],
             year = tF,e0=NA,
             LM_AA = e0.fore)
    period_e0_temp <- period_e0_obs %>% 
      bind_rows(period_e0_fore)
    
    period_e0 <- period_e0 %>% 
      bind_rows(period_e0_temp)
    
    ## computing cohort 
    for (sim in 1:nS) {
      
      LMX_fore_sim <- LMX_fore$SIMnmx[((101*sim)-100):(101*sim),]
      
      # Death rates from 1950 to end of forecast 
      lmx_all <- cbind(LMX, LMX_fore_sim)
      
      ## Observed LMX 2020-2023
      LMX_obs <- log(as.matrix(aux_MX %>% 
                                 filter(sex == filter_sex[s]) %>% 
                                 select(c(`2020`, `2021`, `2022`, `2023`))))
      
      # Death rates from 1950 to end of forecast with observed 2020-23 
      lmx_all_obs <- cbind(LMX, LMX_obs, LMX_fore_sim[,5:82])
      
      #### COHORT ####
      # From Period to Cohort using function from "functions.R"
      lmx_cohort <- 
        period_cohort(lmx=lmx_all, a=a, 
                      first_cohort=aux_ini_cohort, last_cohort=2000)
      
      lmx_cohort_obs <- 
        period_cohort(lmx=lmx_all_obs, a=a, 
                      first_cohort=aux_ini_cohort, last_cohort=2000)
      
      
      e0_cohort <- apply(exp(lmx_cohort), 2, 
                         lifetable.e0, x=a, sex=filter_sex[s])
      
      e0_cohort_obs <- apply(exp(lmx_cohort_obs), 2, 
                             lifetable.e0, x=a, sex=filter_sex[s])
      
      # Saving results cohort e0
      summary_e0_sex <- 
        data.frame(country = filter_country[c],
                   sex = filter_sex[s],
                   simulation = paste0("sim_", sim),
                   cohort = seq(aux_ini_cohort, 2000),
                   e0 = e0_cohort,
                   e0_obs = e0_cohort_obs) %>% 
        as_tibble()
      
      
      # Saving results period lmx (with forecasts)
      summary_lmx_sex <- 
        cbind(country = filter_country[c],
              sex = filter_sex[s],
              simulation = paste0("sim_", sim),
              lmx_all %>% 
                as.data.frame() %>% 
                mutate(age = seq(0,100)) %>% 
                pivot_longer(-age, names_to = "year", values_to = "lmx") %>% 
                arrange(year, age)) %>% 
        left_join(lmx_all_obs %>% 
                    as.data.frame() %>% 
                    mutate(age = seq(0,100)) %>% 
                    pivot_longer(-age, names_to="year", values_to="lmx_obs") %>% 
                    arrange(year, age),
                  by = c("age", "year")) %>% 
        as_tibble()
      
      
      
      summary_e0 <- 
        summary_e0 %>%
        bind_rows(summary_e0_sex)
      
      summary_lmx <- 
        summary_lmx %>%
        bind_rows(summary_lmx_sex)
      
    }
    
  }
  
}

# summary_e0: Cohort life expectancy at birth by birth cohort and sex for  
#selected countries. Cohorts vary from 1950 until 2000. Results 
#were obtained forecasting with the Lee-Miller AA using data until 2019.

# summary_lmx: Period log age-specific death rates for selected countries.
#Years vary from 1950 until 2100. Results were obtained forecasting 
#with the Lee-Miller AA using data until 2019.




#-------------------------------------------------------------------------------
# Saving results for Short-term disturbance (2020-2023)

# Cohort life expectancy estimates
write.table(
  summary_e0 %>% 
    rename(e0_short = e0_obs,
           e0_short_baseline = e0),
  row.names = F,
  "out/22_e0_short_LM_AA.txt")

# Death rates
write.table(
  summary_lmx %>% 
    rename(lmx_short = lmx_obs,
           lmx_short_baseline = lmx),
  row.names = F,
  "out/22_lmx_short_LM_AA.txt")


# period e0 obs & forecast
save(period_e0,file="out/22_period_e0_LM_AA.Rdata")


