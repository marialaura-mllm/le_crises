
#'Evaluating the Impact of Mortality Crises on Life Expectancy: 
#Period and Cohort Perspectives'
# Short-term effect - Li-Lee modified different number of simulations

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
require(ggpubr)

# Functions
source("fun/functions.R")
# List of functions includes 'lifetable.e0' for calculating life expectancy at 
#birth and 'lin_inter' to obtain linear interpolation between two years.

# Lee-Carter functions
source("fun/li_lee_mod.R")
# Function to apply the Li-Lee (2005).



#-------------------------------------------------------------------------------
# Reading data from WPP already cleaned and filtered by country

mx1dt_wpp <- 
  read.table("dat/01_wpp_cleaned.txt", header = T)

# mx1dt_wpp: Age-specific death rates for selected countries with WPP data.
#Years from 1950 to 2023.
# Variables: "country", "year", "age", "sex", "mx", "death" and "exp"



#-------------------------------------------------------------------------------
# Obtaining cohort life expectancy for all countries using the Li-Lee modified
# Short-term disturbance (2020-2023)


summary_e0 <- tibble()
summary_lmx <- tibble()
period_e0 <- tibble()

e0_period <- tibble()
e0_cohort <- tibble()

# Filters for loops
filter_country <- mx1dt_wpp %>% distinct(country) %>% pull()
filter_sim <- c(100, 1000, 10000)

set.seed(1234)
# Loop for different number of simulations
for (n_sim in 1:length(filter_sim)) {
  cat("# of simulations:",filter_sim[n_sim],'\n')

  # Loop for all countries
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
    
    ## time periods
    t1 <- as.integer(seq(1950,2019,1)) # Baseline period
    tF <- as.integer(seq(2020,2101,1)) # Forecast horizon
    nS <- filter_sim[n_sim] # Number of simulations
    
    ## data for females
    LMX1 <- log(as.matrix(aux_MX %>% 
                          filter(sex == "Female") %>% 
                          select(-c(age, sex, `2020`, `2021`, `2022`, `2023`))))
    e01.obs <- apply(exp(LMX1),2,lifetable.e0,x=a,sex="Female")
    ## data for males
    LMX2 <- log(as.matrix(aux_MX %>% 
                         filter(sex == "Male") %>% 
                         select(-c(age, sex, `2020`, `2021`, `2022`, `2023`))))
    e02.obs <- apply(exp(LMX2),2,lifetable.e0,x=a,sex="Male")
    
    
    ## fitting Lee-Miller (comparison purposes)
    ## pop 1
    fitLM <- LM(LMX=LMX1,a=a,t1=t1,tF=tF,nS=nS,sex = "Female")
    e0.fore1.lm <- fitLM$e0
    ## pop 2
    fitLM <- LM(LMX=LMX2,a=a,t1=t1,tF=tF,nS=nS,sex = "Male")
    e0.fore2.lm <- fitLM$e0
    
    ## fitting Li-Lee
    fit.Li.Lee <- LiLee(LMX1 = LMX1,LMX2 = LMX2,a = a,
                        t1 = t1,tF = tF,nS = nS)
    e0.fore1.ll <- fit.Li.Lee$e01
    e0.fore2.ll <- fit.Li.Lee$e02
    
    
    ## saving central period forecasts 
    period_e0_obs_fem <- 
      tibble(country = filter_country[c],
             sex = "Female",
             year = t1,
             e0 = e01.obs,LM=NA,LL=NA)
    period_e0_fore_fem <- 
      tibble(country = filter_country[c],
             sex = "Female",
             year = tF,e0=NA,
             LM = e0.fore1.lm,
             LL=e0.fore1.ll)
    period_e0_fem <- period_e0_obs_fem %>% 
      bind_rows(period_e0_fore_fem)
    
    period_e0_obs_mal <- 
      tibble(country = filter_country[c],
             sex = "Male",
             year = t1,
             e0 = e02.obs,LM=NA,LL=NA)
    period_e0_fore_mal <- 
      tibble(country = filter_country[c],
             sex = "Male",
             year = tF,e0=NA,
             LM = e0.fore2.lm,
             LL=e0.fore2.ll)
    period_e0_mal <- period_e0_obs_mal %>% 
      bind_rows(period_e0_fore_mal)
    
    period_e0_temp <- period_e0_fem %>% 
      bind_rows(period_e0_mal)
    
    period_e0 <- period_e0 %>% 
      bind_rows(period_e0_temp)

    
    # Simulating
    for (sim in 1:nS) {
  
      LMX1_fore_sim <- fit.Li.Lee$SIMnmx1[((101*sim)-100):(101*sim),]
      LMX2_fore_sim <- fit.Li.Lee$SIMnmx2[((101*sim)-100):(101*sim),]
      
      # Death rates from 1950 to end of forecast 
      lmx1_all <- cbind(LMX1, LMX1_fore_sim)
      lmx2_all <- cbind(LMX2, LMX2_fore_sim)
      
      ## Observed LMX 2020-2023
      LMX1_obs <- log(as.matrix(aux_MX %>% 
                                 filter(sex == "Female") %>% 
                                 select(c(`2020`, `2021`, `2022`, `2023`))))
      LMX2_obs <- log(as.matrix(aux_MX %>% 
                                  filter(sex == "Male") %>% 
                                  select(c(`2020`, `2021`, `2022`, `2023`))))
      
      
      # Death rates from 1950 to end of forecast with observed 2020-23 
      lmx1_all_obs <- cbind(LMX1, LMX1_obs, LMX1_fore_sim[,5:82])
      lmx2_all_obs <- cbind(LMX2, LMX2_obs, LMX2_fore_sim[,5:82])
      
      #### COHORT ####
      # From Period to Cohort using function from "functions.R"
      
      ## for females
      lmx1_cohort <- 
        period_cohort(lmx=lmx1_all, a=a, 
                      first_cohort=aux_ini_cohort, last_cohort=2000)
      
      lmx1_cohort_obs <- 
        period_cohort(lmx=lmx1_all_obs, a=a, 
                      first_cohort=aux_ini_cohort, last_cohort=2000)
      
      
      e01_cohort <- apply(exp(lmx1_cohort), 2, 
                         lifetable.e0, x=a, sex="Female")
      
      e01_cohort_obs <- apply(exp(lmx1_cohort_obs), 2, 
                             lifetable.e0, x=a, sex="Female")
      
      # Saving results cohort e0
      summary_e0_fem <- 
        data.frame(country = filter_country[c],
                   sex = "Female",
                   simulation = paste0("sim_", sim),
                   cohort = seq(aux_ini_cohort, 2000),
                   e0 = e01_cohort,
                   e0_obs = e01_cohort_obs) %>% 
        as_tibble()
      
      
      # Saving results period lmx (with forecasts)
      summary_lmx_fem <- 
        cbind(country = filter_country[c],
              sex = "Female",
              simulation = paste0("sim_", sim),
              lmx1_all %>% 
                as.data.frame() %>% 
                mutate(age = seq(0,100)) %>% 
                pivot_longer(-age, names_to = "year", values_to = "lmx") %>% 
                arrange(year, age)) %>% 
        left_join(lmx1_all_obs %>% 
                    as.data.frame() %>% 
                    mutate(age = seq(0,100)) %>% 
                    pivot_longer(-age, names_to="year", values_to="lmx_obs") %>% 
                    arrange(year, age),
                  by = c("age", "year")) %>% 
        as_tibble() %>% 
        filter(year >=2020 & year <= 2023) %>% 
        # Life Expectancy calculations
        group_by(year, sex, country, simulation) %>% 
        mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx)),
               e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_obs)),
               change_e0_period = (e0_observed-e0_fore2019)*12)
      
      ## for males
      lmx2_cohort <- 
        period_cohort(lmx=lmx2_all, a=a, 
                      first_cohort=aux_ini_cohort, last_cohort=2000)
      
      lmx2_cohort_obs <- 
        period_cohort(lmx=lmx2_all_obs, a=a, 
                      first_cohort=aux_ini_cohort, last_cohort=2000)
      
      
      e02_cohort <- apply(exp(lmx2_cohort), 2, 
                          lifetable.e0, x=a, sex="Male")
      
      e02_cohort_obs <- apply(exp(lmx2_cohort_obs), 2, 
                              lifetable.e0, x=a, sex="Male")
      
      # Saving results cohort e0
      summary_e0_mal <- 
        data.frame(country = filter_country[c],
                   sex = "Male",
                   simulation = paste0("sim_", sim),
                   cohort = seq(aux_ini_cohort, 2000),
                   e0 = e02_cohort,
                   e0_obs = e02_cohort_obs) %>% 
        as_tibble()
      
      
      # Saving results period lmx (with forecasts)
      summary_lmx_mal <- 
        cbind(country = filter_country[c],
              sex = "Male",
              simulation = paste0("sim_", sim),
              lmx2_all %>% 
                as.data.frame() %>% 
                mutate(age = seq(0,100)) %>% 
                pivot_longer(-age, names_to = "year", values_to = "lmx") %>% 
                arrange(year, age)) %>% 
        left_join(lmx2_all_obs %>% 
                    as.data.frame() %>% 
                    mutate(age = seq(0,100)) %>% 
                    pivot_longer(-age, names_to="year", values_to="lmx_obs") %>% 
                    arrange(year, age),
                  by = c("age", "year")) %>% 
        as_tibble() %>% 
        filter(year >=2020 & year <= 2023) %>% 
        # Life Expectancy calculations
        group_by(year, sex, country, simulation) %>% 
        mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx)),
               e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_obs)),
               change_e0_period = (e0_observed-e0_fore2019)*12)
      
      
      ## combine
      summary_e0 <- 
        summary_e0 %>%
        bind_rows(summary_e0_fem) %>% 
        bind_rows(summary_e0_mal) 
      
      summary_lmx <- 
        summary_lmx %>%
        bind_rows(summary_lmx_fem) %>% 
        bind_rows(summary_lmx_mal)
      
    }
    
  }


e0_period <- 
  e0_period %>% 
  bind_rows(summary_lmx %>% 
              group_by(year, sex, country) %>% 
              summarise(change_p_med = median(change_e0_period, na.rm=T),
                        change_p_low = quantile(change_e0_period,
                                                prob=(1-0.95)/2, na.rm=T),
                        change_p_upp = quantile(change_e0_period,
                                                prob=1-(1-0.95)/2, na.rm=T)) %>% 
              mutate(n_sim = filter_sim[n_sim])) %>% 
  group_by(sex, country, n_sim) %>% 
  # Most substantial impact
  arrange(change_p_med) %>% 
  slice_head(n=1)

e0_cohort <-
  e0_cohort %>% 
  bind_rows(summary_e0 %>% 
              mutate(change_e0_cohort = (e0_obs-e0)*12) %>% 
              group_by(cohort, sex, country) %>% 
              summarise(change_c_med = median(change_e0_cohort, na.rm=T),
                        change_c_low = quantile(change_e0_cohort,
                                                prob=(1-0.95)/2, na.rm=T),
                        change_c_upp = quantile(change_e0_cohort,
                                                prob=1-(1-0.95)/2, na.rm=T)) %>% 
              mutate(n_sim = filter_sim[n_sim])) %>% 
  group_by(sex, country, n_sim) %>% 
  # Most substantial impact
  arrange(change_c_med) %>% 
  slice_head(n=1)


}

# summary_e0: Cohort life expectancy at birth by birth cohort and sex for  
#selected countries. Cohorts vary from 1950 until 2000. Results 
#were obtained forecasting with Lee-Miller using data until 2019.

# summary_lmx: Period log age-specific death rates for selected countries.
#Years vary from 1950 until 2100. Results were obtained forecasting 
#with Lee-Miller using data until 2019.




#-------------------------------------------------------------------------------
# Saving results for simulations


write.table(
  e0_cohort,
  row.names = F,
  "out/06_e0_cohort_sim.txt")

write.table(
  e0_period,
  row.names = F,
  "out/06_e0_period_sim.txt")




