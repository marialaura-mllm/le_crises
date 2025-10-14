
#'How period life expectancy can distort our interpretation of mortality crises'
# Short-term effect with HMD data

# This code is to be used to replicate Figure S3 of the paper "How period life
#expectancy can distort our interpretation of mortality crises". Specifically, 
#it can be used to obtain cohort life expectancy at birth estimates considering 
#a short-term disturbance using the example of the COVID-19 pandemic and data 
#from HMD.



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
source("fun/lee_miller.R")
# Function to apply the Lee-Miller (2001).




#-------------------------------------------------------------------------------
# Reading data from WPP and HMD already cleaned and filtered by country

mx1dt_hmd <- 
  read.table("dat/02_wpp_hmd_cleaned.txt", header = T)

# mx1dt_hmd: Age-specific death rates for selected countries combining HMD and WPP.
#Years vary from beginning of HMD until 2100. For "death" and "exp", data starts
#in 1950, as this information comes from WPP.
# Variables: "country", "year", "age", "sex", "mx", "death" and "exp"




#-------------------------------------------------------------------------------
# Obtaining cohort life expectancy for all countries using Lee-Miller forecast
# Short-term disturbance (2020-2023)


summary_e0 <- tibble()
summary_lmx <- tibble()

filter_country <- mx1dt_hmd %>% distinct(country) %>% pull()

filter_sex <- c("Male", "Female")

# Loop for all countries
set.seed(1234)
for (c in 1:length(filter_country)) {
  
  a <- as.integer(seq(0,100,1))
  
  aux_all <- 
    mx1dt_hmd %>% 
    filter(country == filter_country[c])
  
  # Exposures
  aux_E <-
    aux_all %>% 
    select(-c(country, death, mx)) %>% 
    filter(year >= 1950) %>% 
    pivot_wider(names_from = year, values_from = exp) %>% 
    arrange(sex, age)
  
  # Deaths
  aux_D <-
    aux_all %>% 
    select(-c(country, exp, mx)) %>% 
    filter(year >= 1950) %>% 
    pivot_wider(names_from = year, values_from = death) %>% 
    arrange(sex, age)
  
  # Death rates
  aux_MX <-
    aux_all %>% 
    select(-c(country, exp, death)) %>% 
    filter(year >= 1950) %>% 
    pivot_wider(names_from = year, values_from = mx) %>% 
    arrange(sex, age)
  
  # Initial cohort
  aux_ini_cohort <- aux_all %>%  select(year) %>%  min
  
  
  # Loop by sex
  for (s in 1:length(filter_sex)) {
    
    #########################       1950-2019       ############################
    
    # LMX with HMD
    LMX_HMD <-
      aux_all %>% 
      select(-c(country, exp, death)) %>% 
      filter(year <= 2019) %>% 
      pivot_wider(names_from = year, values_from = mx) %>% 
      arrange(sex, age) %>%  
      filter(sex == filter_sex[s]) %>% 
      select(-c(age, sex)) %>% 
      as.matrix() %>% 
      log()
    
    
    E <- as.matrix(aux_E %>% 
                     filter(sex == filter_sex[s]) %>% 
                     select(-c(age, sex, `2020`, `2021`, `2022`, `2023`)))
    D <- as.matrix(aux_D %>%
                     filter(sex == filter_sex[s]) %>% 
                     select(-c(age, sex, `2020`, `2021`, `2022`, `2023`)))
    LMX <- log(as.matrix(aux_MX %>% 
                           filter(sex == filter_sex[s]) %>% 
                           select(-c(age, sex, `2020`, `2021`, `2022`, `2023`))))
    
    t1 <- as.integer(seq(1950,2019,1)) # Baseline period
    tF <- as.integer(seq(2020,2111,1)) # Forecast horizon
    nS <- 100 # Number of simulations
    
    # Lee-Miller function from "lc_functions.R"
    LMX_fore <- LM(E=E, LMX=LMX, Y=D, a=a, sex=filter_sex[s], t1=t1, tF=tF, nS=nS)
    
    for (sim in 1:nS) {
      
      LMX_fore_sim <- LMX_fore[((101*sim)-100):(101*sim),]
      
      # Death rates from beginning HMD to end of forecast 
      lmx_all <- cbind(LMX_HMD, LMX_fore_sim)
      
      ## Observed LMX 2020-2023
      LMX_obs <- log(as.matrix(aux_MX %>% 
                                 filter(sex == filter_sex[s]) %>% 
                                 select(c(`2020`, `2021`, `2022`, `2023`))))
      
      # Death rates from beginning HMD to end of forecast with observed 2020-23 
      lmx_all_obs <- cbind(LMX_HMD, LMX_obs, LMX_fore_sim[,5:92])
      
      #### COHORT ####
      # From Period to Cohort using function from "functions.R"
      lmx_cohort <- 
        period_cohort(lmx=lmx_all, a=a, 
                      first_cohort=aux_ini_cohort, last_cohort=2010)
      
      lmx_cohort_obs <- 
        period_cohort(lmx=lmx_all_obs, a=a, 
                      first_cohort=aux_ini_cohort, last_cohort=2010)
      
      
      e0_cohort <- apply(exp(lmx_cohort), 2, 
                         lifetable.e0, x=a, sex=filter_sex[s])
      
      e0_cohort_obs <- apply(exp(lmx_cohort_obs), 2, 
                             lifetable.e0, x=a, sex=filter_sex[s])
      
      # Saving results cohort e0
      summary_e0_sex <- 
        data.frame(country = filter_country[c],
                   sex = filter_sex[s],
                   simulation = paste0("sim_", sim),
                   cohort = seq(aux_ini_cohort, 2010),
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
#selected countries. Cohorts vary from beginning of HMD until 2000. Results 
#were obtained forecasting with Lee-Miller using data until 2019.

# summary_lmx: Period log age-specific death rates for selected countries.
#Years vary from beginning of HMD until 2100. Results were obtained forecasting 
#with Lee-Miller using data until 2019.




#-------------------------------------------------------------------------------
# Saving results for Short-term disturbance (2020-2023) with HMD data

# Cohort life expectancy estimates
write.table(
  summary_e0 %>% 
    rename(e0_short = e0_obs,
           e0_short_baseline = e0),
  row.names = F,
  "out/e0_short_hmd.txt")

# Death rates
write.table(
  summary_lmx %>% 
    rename(lmx_short = lmx_obs,
           lmx_short_baseline = lmx),
  row.names = F,
  "out/lmx_short_hmd.txt")





