
#'How period life expectancy can distort our interpretation of mortality crises'
# Data cleaning

# This code is to be used to replicate the results of the paper "How period life
#expectancy can distort our interpretation of mortality crises". Specifically, 
#it can be used to obtain cleaned data using death rates from World Population 
#prospects 2024 and Human Mortality Database.


#-------------------------------------------------------------------------------

# Cleaning the workspace
rm(list=ls(all=TRUE))

# Packages
require(tidyverse)
require(data.table)
require(countrycode)


#-------------------------------------------------------------------------------
# Obtaining data from the WPP 2024 and HMD
# WPP
# Exposures - "Population on 01 July, by single age. 1950-2023"
dat_exposure <- 
  fread("dat/WPP2024_PopulationBySingleAgeSex_Medium_1950-2023.csv.gz")

# Deaths - "Deaths, by single age. 1950-2023"
dat_death <- 
  fread("dat/WPP2024_DeathsBySingleAgeSex_Medium_1950-2023.csv.gz")

# HMD
# Death rates Mx_1x1 for USA, ITA, JPN and NZL_NP
hmd <- 
  read.table("dat/HMD_Mx_1x1_acessed01102025.txt", header = T)



#-------------------------------------------------------------------------------
# Filtering countries for analyses

# List of selected countries (names should equal those in WPP database)
filter_country <- c("United States of America","Brazil","New Zealand",
                    "Italy", "Mexico", "Japan")



#-------------------------------------------------------------------------------
# Data manipulation

## WPP
# Deaths, Exposures and Death rates
mx1dt_wpp <- 
  dat_death %>% 
  select(country=Location , year=Time, age=AgeGrpStart , 
         Male=DeathMale, Female=DeathFemale) %>% 
  pivot_longer(-c(country, year, age), 
               names_to = "sex", 
               values_to = "death") %>% 
  
  left_join(dat_exposure %>% 
              select(country=Location , year=Time, age=AgeGrpStart , 
                     Male=PopMale, Female=PopFemale) %>% 
              pivot_longer(-c(country, year, age), 
                           names_to="sex", 
                           values_to="exp"),
            
            by = c("country", "year", "age", "sex")) %>%
  
  #  Remove in case of not selecting countries
  filter(country %in% filter_country)


# Filling 0 deaths and exposure with mean last four years by age
# Otherwise, lmx = -Inf and can't use svd for Lee-Miller
mx1dt_wpp <- 
  mx1dt_wpp %>% 
  group_by(country, sex, age) %>% 
  mutate(
    #Deaths
    aux = round(
      # Numerator
      (lag(death,1L)+lag(death,2L)+lag(death,3L)+lag(death,4L))/
        # Denominator
        ((lag(death,1L)!=0) + (lag(death,2L)!=0) +
           (lag(death,3L)!=0) + (lag(death,4L)!=0)),
      3), 
    death = case_when(death == 0 ~ aux,
                      is.na(death) == TRUE ~ aux,
                      TRUE ~ death),
    #Exposure
    aux2 = round(
      # Numerator
      (lag(exp,1L)+lag(exp,2L)+lag(exp,3L)+lag(exp,4L))/
        # Denominator
        ((lag(exp,1L)!=0) + (lag(exp,2L)!=0) +
           (lag(exp,3L)!=0) + (lag(exp,4L)!=0)),
      3), 
    exp = case_when(exp == 0 ~ aux2,
                    is.na(exp) == TRUE ~ aux2,
                    TRUE ~ exp)
    
  ) %>% 
  ungroup() %>% 
  select(-c(aux, aux2)) %>% 
  mutate(mx = death/exp)

# mx1dt_wpp: Age-specific death rates for selected countries with WPP data.
#Years from 1950 to 2023.
# Variables: "country", "year", "age", "sex", "mx", "death" and "exp"


## WPP and HMD
mx1dt_hmd <- 
  hmd %>% 
  mutate(CNTRY = case_when(CNTRY == "GBR_NP" ~ "GBR",
                           CNTRY == "NZL_NP" ~ "NZL",
                           CNTRY == "DEUTNP" ~ "DEU",
                           CNTRY == "FRATNP" ~ "FRA",
                           TRUE ~ CNTRY),
         country = countrycode(CNTRY, origin = "iso3c",
                               destination = "un.name.en")) %>% 
  select(country, year=Year, age=Age, Male, Female) %>% 
  pivot_longer(-c(country, year, age), names_to = "sex", values_to = "mx") %>% 
  filter(age <= 100, year < 1950, year > 1919)

mx1dt_all <- 
  mx1dt_hmd %>% 
  full_join(mx1dt_wpp, 
            by = c("country","year","age","sex","mx"))

# mx1dt_hmd: Age-specific death rates for selected countries combining HMD and WPP.
#Years vary from beginning of HMD until 2100. For "death" and "exp", data starts
#in 1950, as this information comes from WPP. NA for years before that.
# Variables: "country", "year", "age", "sex", "mx", "death" and "exp"




#-------------------------------------------------------------------------------
# Saving the data

# WPP data
write.table(mx1dt_wpp, row.names = F, "dat/01_wpp_cleaned.txt")

# WPP and HMD
write.table(mx1dt_all, row.names = F, "dat/02_wpp_hmd_cleaned.txt")











































