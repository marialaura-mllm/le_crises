
#'Evaluating the Impact of Mortality Crises on Life Expectancy: 
#Period and Cohort Perspectives'
# Long-term effect - Lee-Miller auto ARIMA

# This code is to be used to replicate the results of the paper "Evaluating the 
#Impact of Mortality Crises on Life Expectancy: Period and Cohort Perspectives". 
#Specifically, it can be used to obtain cohort life expectancy at birth 
#estimates considering a long-term disturbance using as example the COVID-19 
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
# Long-term disturbance


summary_e0_long <- tibble()
summary_lmx_long <- tibble()

filter_country <- mx1dt_wpp %>% distinct(country) %>% pull()

filter_sex <- c("Male", "Female")

# Loop for all countries
set.seed(1234)
for (c in 1:length(filter_country)) {
  
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
    
    #########################       1950-2023       ############################
    
    E <- as.matrix(aux_E %>% 
                     filter(sex == filter_sex[s]) %>% 
                     select(-c(age, sex)))
    D <- as.matrix(aux_D %>%
                     filter(sex == filter_sex[s]) %>% 
                     select(-c(age, sex)))
    LMX <- log(as.matrix(aux_MX %>% 
                           filter(sex == filter_sex[s]) %>% 
                           select(-c(age, sex))))
    
    t1 <- as.integer(seq(1950,2023,1)) # Baseline period
    tF <- as.integer(seq(2024,2101,1)) # Forecast horizon
    nS <- 100 # Number of simulations
    
    # Lee-Miller function from "lc_functions.R"
    LMX_fore <- LM_AA(E=E, LMX=LMX, Y=D, a=a, sex=filter_sex[s],
                   t1=t1, tF=tF, nS=nS)
    
    for (sim in 1:nS) {
      
      LMX_fore_sim <- LMX_fore$SIMnmx[((101*sim)-100):(101*sim),]
      
      # Death rates from beginning HMD to end of forecast 
      lmx_all <- cbind(LMX, LMX_fore_sim)
      
      
      #### COHORT ####
      # From Period to Cohort using function from "functions.R"
      lmx_cohort <- 
        period_cohort(lmx=lmx_all, a=a, 
                      first_cohort=aux_ini_cohort, last_cohort=2000)
      
      e0_cohort <- apply(exp(lmx_cohort), 2, 
                         lifetable.e0, x=a, sex=filter_sex[s])
      
      # Saving results cohort e0
      summary_e0_sex <- 
        data.frame(country = filter_country[c],
                   sex = filter_sex[s],
                   simulation = paste0("sim_", sim),
                   cohort = seq(aux_ini_cohort, 2000),
                   e0 = e0_cohort) %>% 
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
        as_tibble()
      
      
      
      summary_e0_long <- 
        summary_e0_long %>%
        bind_rows(summary_e0_sex)
      
      summary_lmx_long <- 
        summary_lmx_long %>%
        bind_rows(summary_lmx_sex)
      
    }
    
  }
  
}

# summary_e0_long: Cohort life expectancy at birth by birth cohort and sex 
#for selected countries. Cohorts vary from 1950 to 2000. Results 
#were obtained forecasting with the Lee-Miller AA using data until 2023.

# summary_lmx_long: Period log age-specific death rates for selected 
#countries. Years vary from 1950 to 2100. Results were obtained 
#forecasting with the Lee-Miller AA using data until 2023.





#-------------------------------------------------------------------------------
# Saving results for Long-term disturbance

# Getting results from short to include information on baseline mortality
summary_e0_short <- 
  read.table("out/22_e0_short_LM_AA.txt", header = T)


# Cohort life expectancy estimates
write.table(
  # Baseline from Short-term
  summary_e0_short %>% 
    select(country, sex, simulation_19=simulation, cohort, 
           fore_19=e0_short_baseline) %>% 
    # Results from Long-term
    left_join(summary_e0_long %>% 
                select(country, sex, simulation_23=simulation, cohort, 
                       fore_23=e0),
              by = c("country", "sex", "cohort"), relationship= "many-to-many"),
  row.names = F,
  "out/32_e0_long_LM_AA.txt")

# Death rates
write.table(
  summary_lmx_long %>% 
    select(country, sex, simulation_23=simulation, age, year, fore_23=lmx),
  row.names = F,
  "out/32_lmx_long_LM_AA.txt")






