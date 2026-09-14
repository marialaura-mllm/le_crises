
#'Evaluating the Impact of Mortality Crises on Life Expectancy: 
#Period and Cohort Perspectives'
# Long-term effect -  Product Ratio

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
require(demography)

# Functions
source("fun/functions.R")
# List of functions includes 'lifetable.e0' for calculating life expectancy at 
#birth and 'lin_inter' to obtain linear interpolation between two years.



#-------------------------------------------------------------------------------
# Reading data from WPP already cleaned and filtered by country

mx1dt_wpp <- 
  read.table("dat/01_wpp_cleaned.txt", header = T)

# mx1dt_wpp: Age-specific death rates for selected countries with WPP data.
#Years from 1950 to 2023.
# Variables: "country", "year", "age", "sex", "mx", "death" and "exp"



#-------------------------------------------------------------------------------
# Obtaining cohort life expectancy for all countries using the Product Ratio
# Long-term disturbance

summary_e0 <- tibble()
summary_lmx <- tibble()
period_e0 <- tibble()
filter_country <- mx1dt_wpp %>% distinct(country) %>% pull()

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
    dplyr::select(-c(country, death, mx)) %>% 
    pivot_wider(names_from = year, values_from = exp) %>% 
    arrange(sex, age)
  
  # Deaths
  aux_D <-
    aux_all %>% 
    dplyr::select(-c(country, exp, mx)) %>% 
    pivot_wider(names_from = year, values_from = death) %>% 
    arrange(sex, age)
  
  # Death rates
  aux_MX <-
    aux_all %>% 
    dplyr::select(-c(country, exp, death)) %>% 
    pivot_wider(names_from = year, values_from = mx) %>% 
    arrange(sex, age)
  
  # Initial cohort
  aux_ini_cohort <- aux_all %>%  dplyr::select(year) %>%  min
  
  ## time periods
  t1 <- as.integer(seq(1950,2023,1)) # Baseline period
  tF <- as.integer(seq(2024,2101,1)) # Forecast horizon
  nS <- 100 # Number of simulations
  
  ## data for females
  LMX1 <- log(as.matrix(aux_MX %>% 
                        filter(sex == "Female") %>% 
                          dplyr::select(-c(age, sex))))
  e01.obs <- apply(exp(LMX1),2,lifetable.e0,x=a,sex="Female")
  E1 <- as.matrix(aux_E %>% 
                  filter(sex == "Female") %>% 
                    dplyr::select(-c(age, sex)))
  
  data.fem <- demogdata(data=exp(LMX1), pop=E1, ages=a, years=t1, 
                        type="mortality", label=filter_country[c], 
                        name="female", lambda=0)
  
  
  ## data for males
  LMX2 <- log(as.matrix(aux_MX %>% 
                       filter(sex == "Male") %>% 
                         dplyr::select(-c(age, sex))))
  e02.obs <- apply(exp(LMX2),2,lifetable.e0,x=a,sex="Male")
  
  E2 <- as.matrix(aux_E %>% 
                    filter(sex == "Male") %>% 
                    dplyr::select(-c(age, sex)))
  
  data.mal <- demogdata(data=exp(LMX2), pop=E2, ages=a, years=t1, 
                        type="mortality", label=filter_country[c], 
                        name="male", lambda=0)
  
  ## combine
  data.comb <- list(
    type   = data.fem$type,
    label  = data.fem$label,
    lambda = data.fem$lambda,
    year   = data.fem$year,
    age    = data.fem$age,
    rate   = list(female = data.fem$rate[[1]], 
                  male   = data.mal$rate[[1]]),
    pop    = list(female = data.fem$pop[[1]],  
                  male   = data.mal$pop[[1]])
  )
  class(data.comb) <- "demogdata"
  
  
  ## fitting and forecasting with Product-Ratio Li-Lee
  fit.product.ratio <- coherentfdm(data.comb)
  fore.product.ratio <- forecast(fit.product.ratio,h = length(tF))
  plot(fore.product.ratio$male)

  LMX1_fore <- log(fore.product.ratio$female$rate$female)
  LMX2_fore <- log(fore.product.ratio$male$rate$male)
  
  e01.fore <- apply(exp(LMX1_fore),2,lifetable.e0,x=a,sex="Female")
  e02.fore <- apply(exp(LMX2_fore),2,lifetable.e0,x=a,sex="Male")
  
  # Death rates from 1950 to end of forecast 
  lmx1_all <- cbind(LMX1, LMX1_fore)
  lmx2_all <- cbind(LMX2, LMX2_fore)
  
  # Death rates from 1950 to end of forecast with observed 2020-23 
  lmx1_all_obs <- cbind(LMX1, LMX1_fore)
  lmx2_all_obs <- cbind(LMX2, LMX2_fore)
  
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
               # simulation = paste0("sim_", sim),
               cohort = seq(aux_ini_cohort, 2000),
               e0 = e01_cohort,
               e0_obs = e01_cohort_obs) %>% 
    as_tibble()
  
  
  # Saving results period lmx (with forecasts)
  summary_lmx_fem <- 
    cbind(country = filter_country[c],
          sex = "Female",
          # simulation = paste0("sim_", sim),
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
    as_tibble()
  
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
               # simulation = paste0("sim_", sim),
               cohort = seq(aux_ini_cohort, 2000),
               e0 = e02_cohort,
               e0_obs = e02_cohort_obs) %>% 
    as_tibble()
  
  
  # Saving results period lmx (with forecasts)
  summary_lmx_mal <- 
    cbind(country = filter_country[c],
          sex = "Male",
          # simulation = paste0("sim_", sim),
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
    as_tibble()
  
  
  ## combine
  summary_e0 <- 
    summary_e0 %>%
    bind_rows(summary_e0_fem) %>% 
    bind_rows(summary_e0_mal) 
  
  summary_lmx <- 
    summary_lmx %>%
    bind_rows(summary_lmx_fem) %>% 
    bind_rows(summary_lmx_mal)
  
  ## saving central period forecasts 
  period_e0_obs_fem <- 
    tibble(country = filter_country[c],
           sex = "Female",
           year = t1,
           e0 = e01.obs,ProdRatio=NA)
  period_e0_fore_fem <- 
    tibble(country = filter_country[c],
           sex = "Female",
           year = tF,e0=NA,
           ProdRatio = e01.fore)
  period_e0_fem <- period_e0_obs_fem %>% 
    bind_rows(period_e0_fore_fem)
  
  period_e0_obs_mal <- 
    tibble(country = filter_country[c],
           sex = "Male",
           year = t1,
           e0 = e02.obs,ProdRatio=NA)
  period_e0_fore_mal <- 
    tibble(country = filter_country[c],
           sex = "Male",
           year = tF,e0=NA,
           ProdRatio = e02.fore)
  period_e0_mal <- period_e0_obs_mal %>% 
    bind_rows(period_e0_fore_mal)
  
  period_e0_temp <- period_e0_fem %>% 
    bind_rows(period_e0_mal)
  
  period_e0 <- period_e0 %>% 
    bind_rows(period_e0_temp)
  
  
}

# summary_e0: Cohort life expectancy at birth by birth cohort and sex for  
#selected countries. Cohorts vary from 1950 until 2000. Results 
#were obtained forecasting with the Product Ratio using data until 2019.

# summary_lmx: Period log age-specific death rates for selected countries.
#Years vary from 1950 until 2100. Results were obtained forecasting 
#with  the Product Ratio using data until 2019.




#-------------------------------------------------------------------------------
# Saving results for Long-term disturbance

# Getting results from short to include information on baseline mortality
summary_e0_short <- 
  read.table("out/26_e0_short_product_ratio.txt", header = T)


# Cohort life expectancy estimates
write.table(
  # Baseline from Short-term
  summary_e0_short %>% 
    dplyr::select(country, sex, cohort, 
                  fore_19=e0_short_baseline) %>% 
    # Results from Long-term
    left_join(summary_e0 %>% 
                dplyr::select(country, sex, cohort, 
                              fore_23=e0),
              by = c("country", "sex", "cohort"), relationship= "many-to-many"),
  row.names = F,
  "out/36_e0_long_product_ratio.txt")

# Death rates
write.table(
  summary_lmx %>% 
    dplyr::select(country, sex, age, year, fore_23=lmx),
  row.names = F,
  "out/36_lmx_long_product_ratio.txt")



