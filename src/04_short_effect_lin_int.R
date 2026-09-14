
#'Evaluating the Impact of Mortality Crises on Life Expectancy: 
#Period and Cohort Perspectives'
# Short-term effect and linear Interpolation with WPP

# This code is to be used to replicate Figure S2 of the paper "Evaluating the 
#Impact of Mortality Crises on Life Expectancy: Period and Cohort Perspectives". 
#Specifically, it can be used to obtain cohort life expectancy at birth 
#estimates considering a linear interpolation of WPP rates and a short-term 
#disturbance using as example the COVID-19 pandemic.



#-------------------------------------------------------------------------------

# Cleaning the workspace
rm(list=ls(all=TRUE))

# Packages
require(tidyverse)
require(data.table)

# Functions
source("fun/functions.R")
# List of functions includes 'lifetable.e0' for calculating life expectancy at 
#birth and 'lin_inter' to obtain linear interpolation between two years.




#-------------------------------------------------------------------------------
# Obtaining data from the WPP 2024

# Exposures - "Population on 01 July, by single age. 1950-2023"
dat_exposure <- 
  fread('dat/WPP2024_PopulationBySingleAgeSex_Medium_1950-2023.csv.gz')

# Deaths - "Deaths, by single age. 1950-2023"
dat_death <- 
  fread('dat/WPP2024_DeathsBySingleAgeSex_Medium_1950-2023.csv.gz')

#### Forecasts - medium scenario ####
# Forecasts Exposures - "Population on 01 July, by single age. 2014-2100, medium"
dat_fore_exposure <- 
  fread('dat/WPP2024_PopulationBySingleAgeSex_Medium_2024-2100.csv.gz')

# Forecasts Deaths - "Deaths, by single age. 2024-2100"
dat_fore_death <- 
  fread('dat/WPP2024_DeathsBySingleAgeSex_Medium_2024-2100.csv.gz')


dat_exposure <- 
  dat_exposure %>% 
  full_join(dat_fore_exposure)

dat_death <- 
  dat_death %>% 
  full_join(dat_fore_death)





#-------------------------------------------------------------------------------
# Filtering countries for analyses

# List of selected countries
filter_country <- c('United States of America','Brazil','New Zealand',
                    'Italy','Mexico','Japan')


# Table to manually determine the last year for interpolation
# PS: Just in case of considering different returning years per country
table_interp <- 
  data.table(
    country = c('United States of America','Brazil','New Zealand',
                'Italy','Mexico','Japan'),
    last_year = c(2023, 2023, 2023, 2023, 2023, 2023))




#-------------------------------------------------------------------------------
# Data manipulation

# Deaths, Exposures and Death rates
mx1dt <-
  dat_death %>% 
  select(country=Location, year=Time, age=AgeGrpStart, 
         DeathMale, DeathFemale) %>% 
  
  left_join(dat_exposure %>% 
              select(country=Location , year=Time, age=AgeGrpStart , 
                     PopMale, PopFemale),
            by = c('country', 'year', 'age'), relationship ="many-to-many") %>% 
  mutate(mxM = DeathMale/PopMale,
         mxF = DeathFemale/PopFemale) %>% 
  select(c(country, year, age, Male=mxM, Female=mxF)) %>% 
  
  #  Remove in case of not selecting countries
  filter(country %in% filter_country) %>% 
  pivot_longer(-c(country, year, age), names_to = 'sex', values_to = 'mx') %>% 
  mutate(lmx = log(mx)) %>% 
  select(-mx)


# mx1dt: Log age-specific death rates for selected countries from 1950 to 2100.
#Variables: 'country', 'year', 'age', 'sex' and 'lmx'




#-------------------------------------------------------------------------------
# Obtaining cohort life expectancy for all countries using linear interpolation
# Short-term disturbance (2020-2023) interpolating WPP rates from 2020 to 2023

# List of countries to apply the loop
countries <- unique(mx1dt$country)

filter_sex <- c('Male', 'Female')

summary_e0 <- tibble()
summary_lmx <- tibble()


# Loop for all countries
for (c in 1:length(countries)) {
  
  aux_mx <- 
    mx1dt %>% 
    filter(country == countries[c])
  
  # Years to be interpolated by country
  interp_years <- 
    as.character(2020:table_interp[country==countries[c],]$last_year)
  
  # Considering cohorts finishing at age 100
  min_coh <- min(aux_mx$year)
  max_coh <- max(aux_mx$year)-100
  
  a <- as.integer(seq(0,100,1))
  
  # Loop by sex
  for (s in 1:length(filter_sex)) {
    
    # lmx without linear interpolation
    lmx_period <- 
      aux_mx %>% 
      filter(sex==filter_sex[s]) %>% 
      select(-c(country, sex)) %>% 
      pivot_wider(names_from = year, values_from = lmx) %>% 
      arrange(age) %>% 
      select(-age) %>% 
      mutate(`2101` = `2100`) %>% 
      as.matrix()
    
    
    # lmx with linear interpolation
    lmx_period_inter <- 
      aux_mx %>% 
      filter(sex==filter_sex[s]) %>% 
      select(-c(country, sex)) %>% 
      pivot_wider(names_from = year, values_from = lmx) %>% 
      arrange(age) %>%
      # Changing years to be interpolated
      mutate_at(interp_years, ~ na_if(., .)) %>%
      pivot_longer(-age, names_to = 'year', values_to = 'lmx') %>% 
      arrange(year, age) %>% 
      group_by(age) %>% 
      # Imputing the years of covid with the linear interpolation
      group_modify(~ lin_inter(data = .x)) %>% 
      ungroup() %>%
      pivot_wider(names_from = year,
                  values_from = lmx) %>% 
      arrange(age) %>% 
      select(-age) %>% 
      mutate(`2101` = `2100`)%>% 
      as.matrix()
    
    
    #### COHORT ####
    # From Period to Cohort using function from "functions.R"
    lmx_cohort <- 
      period_cohort(lmx=lmx_period,
                    a=a,first_cohort=min_coh,last_cohort=max_coh)
    
    lmx_cohort_inter <- 
      period_cohort(lmx=lmx_period_inter,
                    a=a,first_cohort=min_coh,last_cohort=max_coh)
    
    
    e0_cohort <- 
      apply(exp(lmx_cohort), 2, lifetable.e0, x=a, sex=filter_sex[s])
    
    e0_cohort_inter <- 
      apply(exp(lmx_cohort_inter), 2, lifetable.e0, x=a, sex=filter_sex[s])
    
    
    # Saving results cohort e0
    summary_e0_sex <- 
      data.frame(country = countries[c],
                 sex = filter_sex[s],
                 cohort = seq(min_coh, max_coh),
                 e0 = e0_cohort,
                 e0_inter = e0_cohort_inter) %>% 
      as_tibble()
    
    
    # Saving results period lmx (with interpolation)
    summary_lmx_sex <-
      lmx_period_inter %>% 
      as.data.frame() %>% 
      mutate(age = a,
             sex = filter_sex[s],
             country = countries[c]) %>% 
      pivot_longer(-c(age,sex,country), names_to='year', values_to='lmx_inter') %>% 
      mutate(year = as.integer(year))
    
    
    summary_e0 <- 
      summary_e0 %>%
      bind_rows(summary_e0_sex)
    
    summary_lmx <- 
      summary_lmx %>%
      bind_rows(summary_lmx_sex)
    
  }
  
  
  summary_e0 <- 
    summary_e0 %>%
    bind_rows(summary_e0_sex)
  
  summary_lmx <- 
    summary_lmx %>%
    bind_rows(summary_lmx_sex)
  
}

# Adding lmx without interpolation
summary_lmx <-
  summary_lmx %>% 
  left_join(mx1dt, by = c('age', 'sex', 'country', 'year'))

# summary_e0: Cohort life expectancy at birth by birth cohort and sex for  
#selected countries. Cohorts vary from 1950 to 2000. 'e0' uses data until 2023 
#and forecasts. 'e0_inter' uses linear interpolation.

# summary_lmx: Period log age-specific death rates for selected countries.
#Years vary from 1950 to 2100. 'lmx' uses data until 2023 and forecasts.
#'lmx_inter' uses linear interpolation.




#-------------------------------------------------------------------------------
#Saving results of short-term effect with linear interpolation of WPP rates

# Cohort life expectancy estimates
write.table(summary_e0 %>% 
              rename(e0_short=e0, e0_short_baseline=e0_inter),
            row.names = F,
            'out/04_e0_short_lin_int.txt')

# Death rates
write.table(summary_lmx,
            row.names = F,
            'out/04_lmx_short_lin_int.txt')








