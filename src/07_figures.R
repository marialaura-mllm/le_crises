
#'Evaluating the Impact of Mortality Crises on Life Expectancy: 
#Period and Cohort Perspectives'
# Figures

# This code is to be used to replicate the figures of the paper "Evaluating the 
#Impact of Mortality Crises on Life Expectancy: Period and Cohort Perspectives".



#-------------------------------------------------------------------------------

# Cleaning the workspace
rm(list=ls(all=TRUE))

# Packages
require(tidyverse)
require(ggpubr)
require(xtable)
require(countrycode)
require(LexisPlotR)

# Functions
source("fun/functions.R")
# List of functions includes 'lifetable.e0' for calculating life expectancy at 
#birth and 'lin_inter' to obtain linear interpolation between two years.





#-------------------------------------------------------------------------------
# Figure 1 - Trends in period life expectancy
# Needs codes src/01_ and 02_

# Age-specific death rates
e0_period <- 
  read.table("out/02_lmx_short_LL_mod.txt", header = T)


aux_trend <- 
  e0_period %>% 
  filter(year >=2017 & year <= 2023) %>% 
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  group_by(country, sex, year) %>% 
  summarise(e0c_med = median(e0_fore2019, na.rm=T),
            e0c_low = quantile(e0_fore2019,prob=(1-0.95)/2, na.rm=T),
            e0c_upp = quantile(e0_fore2019,prob=1-(1-0.95)/2, na.rm=T),
            e0_observed = unique(e0_observed)) %>%
  ungroup() %>% 
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12,
         change_e0_period_low = (e0_observed-e0c_low)*12,
         change_e0_period_upp = (e0_observed-e0c_upp)*12,
         country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil","Mexico")),
         label_med = case_when(year == 2021 & country %in% 
                                 c("United States of America","Brazil","Mexico") 
                               ~ change_e0_period_med,
                               year == 2022 & 
                                 country %in% c("New Zealand","Japan","Italy") 
                               ~ change_e0_period_med,
                               TRUE ~ NA))

# Plot
aux_trend %>% 
  ggplot(aes(x=year, y=e0_observed, group=sex, color=sex, 
             label=round(label_med,0)))+
  
  # Data until 2019
  geom_line(data=subset(aux_trend,year<=2019), lwd=1.3, linetype="solid")+
  
  # 95% PI
  geom_ribbon(data=subset(aux_trend,year>=2019), color="grey85", fill="black", 
              aes(ymin=e0c_low, ymax=e0c_upp), alpha=0.1)+
  
  # Expected trend
  geom_line(aes(y=e0c_med), data=subset(aux_trend,year>=2019), 
            linetype="twodash", lwd=0.8, alpha = 0.5, color="black")+
  # Observed after 2019
  geom_line(data=subset(aux_trend,year>=2019), lwd=1.3, alpha = 0.6)+
  geom_point(data=subset(aux_trend,year>=2020), size =2, alpha = 0.5)+
  
  #Labels
  geom_text(hjust = 0,vjust = 1.1, size =3.5, angle=0, aes(y=e0_observed),
            show.legend = F, color = "black")+
  
  facet_wrap(~country, scales = "free_x")+
  
  # Extra changes
  theme_classic()+
  labs(y = expression(""~e[0]^P),
       x = "Year",
       color = "Sex")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_color_manual("Sex", 
                     values = c("Male"="steelblue3", "Female"="firebrick2"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand=c(0, 0))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), expand=c(0, 0),
                     limits = c(63,90))


#ggsave(path = "out/",filename = "fig1.pdf", width = 12, height = 9)
#ggsave(path = "out/",filename = "fig1.png", width = 12, height = 9)




#-------------------------------------------------------------------------------
# Figure 2 - Short-term disturbance
# Needs codes src/01_, and 02_

# Cohort life expectancy from short-term disturbance scenario
e0_short <- 
  read.table("out/02_e0_short_LL_mod.txt", header = T) %>% 
  distinct()


aux <-
  e0_short %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12,
         country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil", "Mexico"))) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T),
            e0c_low = quantile(change_e0_temp,prob=(1-0.95)/2, na.rm=T),
            e0c_upp = quantile(change_e0_temp,prob=1-(1-0.95)/2, na.rm=T)) %>%
  ungroup()



aux %>%
  ggplot(aes(x = cohort,
             group = sex, 
             color = sex))+
  geom_line(aes(y=e0c_med), lwd=1.1)+
  # 95% PI
  geom_ribbon(color = "grey85",
              aes(ymin=e0c_low, ymax=e0c_upp, fill = sex), alpha=0.2)+
  
  geom_hline(yintercept = 0)+
  
  facet_wrap(~country, ncol = 3, scales = "free_x")+
  
  # Extra changes
  theme_classic()+
  
  labs(x = "Birth cohort",
       y = expression("Change in months "~e[0]^C),
       color = "Sex")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_color_manual("Sex", 
                     values = c("Male"="steelblue3", "Female"="firebrick2"))+
  scale_fill_manual("Sex", 
                    values = c("Male"="steelblue3", "Female"="firebrick2"))+
  
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), expand=c(0, 0))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))


#ggsave(path = "out/",filename = "fig2.pdf", width = 12, height = 9)
#ggsave(path = "out/",filename = "fig2.png", width = 12, height = 9)




#-------------------------------------------------------------------------------
# Figure 3 - Ratio Period/Cohort
# Needs codes src/01_, and 02_

## Period estimates
# Age-specific death rates
e0_period <- 
  read.table("out/02_lmx_short_LL_mod.txt", header = T)


aux_period <- 
  e0_period %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  summarise(e0c_med = median(e0_fore2019, na.rm=T),
            e0_observed = unique(e0_observed)) %>%
  ungroup() %>% 
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort <- 
  read.table("out/02_e0_short_LL_mod.txt", header = T) %>% 
  distinct()


aux_cohort <-
  e0_cohort %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))


aux <- 
  aux_period %>% 
  left_join(aux_cohort, by = c("sex", "country")) %>% 
  mutate(ratio = period/cohort)


# Plot
aux %>% 
  mutate(country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil","Mexico")),
         position = (period - cohort)/2) %>% 
  ggplot(aes(y=period, x=country, group=sex, color=sex,
             label=paste0(round(ratio, 0),"x")))+
  geom_hline(yintercept = 0)+
  
  # Line connecting period and cohort estimates
  geom_segment(aes(yend=cohort), lwd=0.8, alpha = 0.6, 
               position=position_nudge(x = ifelse(aux$sex=="Male", 0.2, -0.2)))+
  
  # Period and cohort estimates
  geom_point(aes(shape = "Period"), size =3,
             position = position_nudge(x = ifelse(aux$sex=="Male", 0.2, -0.2)))+
  geom_point(aes(y=cohort, shape = "Cohort"), size =3,
             position=position_nudge(x = ifelse(aux$sex == "Male", 0.2, -0.2)))+
  
  # Labels
  ## Ratio
  geom_text(vjust = 0,hjust = 1.2, size =6, angle=0, aes(y=position),
            show.legend = F,fontface = "bold",
            position=position_nudge(x = ifelse(aux$sex == "Male", 0.2, -0.2)))+
  ## Period
  geom_text(vjust = 0.9,hjust = 1.4, size =4.5, angle=0,
            aes(label=round(period, 0)), show.legend = F, color = "Black", 
            position=position_nudge(x = ifelse(aux$sex == "Male", 0.2, -0.2)))+
  
  ## Cohort
  geom_text(vjust = 0.9,hjust = 1.4, size =4.5, angle=0, 
            aes(y=cohort,label=round(cohort, 1)), show.legend=F, color="Black",
            position=position_nudge(x = ifelse(aux$sex == "Male", 0.2, -0.2)))+
  
  
  # Extra changes
  theme_classic()+
  
  labs(y = expression("Change in months " ~e[0]^P~" and "~e[0]^C),
       x = "",
       color = "Sex")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_color_manual("Sex", 
                     values = c("Male"="steelblue3", "Female"="firebrick2"))+
  scale_shape_manual("", 
                     values = c("Period"=17, "Cohort"=16))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), expand=c(0, 0.8),
                     limits = c(-72,0))

#ggsave(path = "out/",filename = "fig3.pdf", width = 12, height = 7)
#ggsave(path = "out/",filename = "fig3.png", width = 12, height = 7)




#-------------------------------------------------------------------------------
# Figure 4 - Age-specific death rates for impact measurement
# Needs codes src/01_, 02_, and 03_

lmx_short <- 
  read.table("out/02_lmx_short_LL_mod.txt", header = T) 

lmx_long <- 
  read.table("out/03_lmx_long_LL_mod.txt", header = T) 


aux_lmx60 <- 
  lmx_short %>% 
  filter(country == "Brazil",
         sex == "Male",
         age == 60,
         year >= 2017 & year <= 2027) %>%
  group_by(country, sex, year, age) %>% 
  mutate(lmx_med = median(lmx_short_baseline, na.rm=T),
         lmx_low = quantile(lmx_short_baseline,prob=(1-0.95)/2, na.rm=T),
         lmx_upp = quantile(lmx_short_baseline,prob=1-(1-0.95)/2, na.rm=T)) %>% 
  # Create manually 2024 point for better plotting
  mutate(lmx_short = case_when(year == 2024 ~ lmx_med,
                               TRUE ~ lmx_short))

aux_lmx60_long <- 
  lmx_long %>% 
  filter(country == "Brazil",
         sex == "Male",
         age == 60,
         year >= 2017 & year <= 2027) %>%
  group_by(country, sex, year, age) %>% 
  mutate(lmx_med_last = median(fore_23, na.rm=T),
         lmx_low_last = quantile(fore_23,prob=(1-0.95)/2, na.rm=T),
         lmx_upp_last = quantile(fore_23,prob=1-(1-0.95)/2, na.rm=T))

aux_lmx60_all <-
  aux_lmx60 %>% 
  left_join(aux_lmx60_long,
            by = c("country", "sex", "age", "year"))


# Short-term

# Plot
g1 <- 
  aux_lmx60  %>% 
  ggplot(aes(x=year))+
  
  
  geom_line(data=subset(aux_lmx60,year>=2019),
            aes(y=lmx_short_baseline, linetype="Simulations", group=simulation),
            lwd=1.1, color="grey50", alpha = 0.1)+
  
  # Impact
  geom_ribbon(data=subset(aux_lmx60,year<=2024),
              aes(ymin=lmx_med, ymax=lmx_short, fill = "Impact"), alpha=0.4)+
  # Observed
  geom_line(data=subset(aux_lmx60,year<=2024),
            aes(y=lmx_short, linetype="Observed"), lwd=1.2, color="#0F4C5C")+
  
  #Medium estimate
  geom_line(data=subset(aux_lmx60,year>=2019),
            aes(y=lmx_med, linetype="Medium estimate"), lwd=1.2, 
            color = "red")+
  
  theme_classic()+
  labs(title = "a.",
       y=expression("log("~m[60]^P~")"),
       x="Year")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey85",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_linetype_manual("Type of data",
                        values = c("WPP 2024"="solid", "Observed"="solid",
                                   "Simulations"="solid",
                                   "Medium estimate"="solid"))+
  scale_fill_manual("", values = c("Impact"="#E8871E"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand=c(0, 0),
                     limits = c(2017,2027))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8),
                     limits = c(-4.68, -4))




# Long-term

# Plot
g2 <-
  aux_lmx60_all  %>% 
  ggplot(aes(x=year))+
  
  
  # Simulations
  geom_line(data=subset(aux_lmx60_all,year>=2019),
            aes(y=lmx_short_baseline, linetype="Simulations", group=simulation),
            lwd=1.1, color="grey50", alpha = 0.1)+
  # Fore 1050-2023
  geom_line(data=subset(aux_lmx60_all,year>=2023),
            aes(y=fore_23, linetype="Simulations", group = simulation_23),
            lwd=1.1, color="grey50", alpha = 0.1)+
  
  # Impact
  geom_ribbon(aes(ymin=lmx_med_last, ymax=lmx_med, fill = "Impact"), alpha=0.4)+
  
  # Observed
  geom_line(data=subset(aux_lmx60_all,year<=2023),
            aes(y=lmx_short, linetype="Observed"), lwd=1.2, color="#0F4C5C")+
  
  # Medium estimates
  geom_line(data=subset(aux_lmx60_all,year>=2019),
            aes(y=lmx_med, linetype="Medium estimate"), 
            lwd=1.2, color="red")+
  
  geom_line(data=subset(aux_lmx60_all,year>=2023),
            aes(y=lmx_med_last, linetype="Medium estimate"), 
            lwd=1.2, color="red")+
  
  
  
  theme_classic()+
  labs(title = "b.",
       y=expression("log("~m[60]^P~")"),
       x="Year")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey85",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_linetype_manual("Type of data",
                        values = c("WPP 2024"="solid", "Observed"="solid",
                                   "Simulations"="solid",
                                   "Medium estimate"="solid"))+
  scale_fill_manual("", values = c("Impact"="#E8871E"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand=c(0, 0),
                     limits = c(2017,2027))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8),
                     limits = c(-4.68, -4))


ggarrange(g1, g2,ncol = 2, common.legend = TRUE, legend="bottom") 

#ggsave(path = "out/",filename = "fig4.pdf", width = 14, height = 6)
#ggsave(path = "out/",filename = "fig4.png", width = 14, height = 6)




#-------------------------------------------------------------------------------
# Table S1 - Summary of peak losses
# Needs codes src/01_, and 02_

## Period estimates
# Age-specific death rates
e0_period <- 
  read.table("out/02_lmx_short_LL_mod.txt", header = T)

aux_period <- 
  e0_period %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  summarise(e0c_med = median(e0_fore2019, na.rm=T),
            e0_observed = unique(e0_observed)) %>%
  ungroup() %>% 
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  slice_min(change_e0_period_med) %>%
  ungroup()


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort <- 
  read.table("out/02_e0_short_LL_mod.txt", header = T) %>% 
  distinct()


aux_cohort <-
  e0_cohort %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  slice_min(e0c_med) %>%
  ungroup()


aux <- 
  aux_period %>% 
  left_join(aux_cohort, by = c("sex", "country")) %>% 
  select(Country=country, Sex=sex, 
         'Year Peak'=year, 'Change Period'=change_e0_period_med,
         'Cohort Peak'=cohort, 'Change Cohort'=e0c_med.y) %>% 
  arrange(factor(Country, c("New Zealand","Japan","Italy",
                            "United States of America",
                            "Brazil", "Mexico")))


#print(xtable(aux, type="latex"),include.rownames=FALSE,file = "out/tabS1.tex")


#-------------------------------------------------------------------------------
# Figure S1 - Short-term disturbance with HMD data
# Needs codes src/01_, and 05_

# Cohort life expectancy from short-term disturbance scenario including HMD
e0_hmd <- 
  read.table("out/05_e0_short_hmd_LL_mod.txt", header = T) %>% 
  distinct()


aux <-
  e0_hmd %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12,
         country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil", "Mexico"))) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T),
            e0c_low = quantile(change_e0_temp,prob=(1-0.95)/2, na.rm=T),
            e0c_upp = quantile(change_e0_temp,prob=1-(1-0.95)/2, na.rm=T)) %>%
  ungroup()



aux %>%
  ggplot(aes(x = cohort,
             group = sex, 
             color = sex))+
  geom_line(aes(y=e0c_med), lwd=1.1)+
  # 95% PI
  geom_ribbon(color = "grey85",
              aes(ymin=e0c_low, ymax=e0c_upp, fill = sex), alpha=0.2)+
  theme_classic()+
  facet_wrap(~country, ncol = 3, scales = "free_x")+
  labs(x = "Birth cohort",
       y = expression("Change in months "~e[0]^C),
       color = "Sex")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_color_manual("Sex", 
                     values = c("Male"="steelblue3", "Female"="firebrick2"))+
  scale_fill_manual("Sex", 
                    values = c("Male"="steelblue3", "Female"="firebrick2"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), expand=c(0, 0), 
                     limits = c(1920, 2000))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  geom_hline(yintercept = 0)

#ggsave(path = "out/",filename = "figS1.pdf", width = 12, height = 9)
#ggsave(path = "out/",filename = "figS1.png", width = 12, height = 9)




#-------------------------------------------------------------------------------
# Figure S2 - Long-term disturbance
# Needs codes src/01_, 02_, and 03_

# Cohort life expectancy from long-term disturbance scenario
e0_long <-
  read.table("out/03_e0_long_LL_mod.txt", header = T) %>% 
  distinct()


aux <-
  e0_long %>% 
  mutate(change_e0_last = (fore_23-fore_19)*12,
         country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil", "Mexico"))) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_last, na.rm=T),
            e0c_low = quantile(change_e0_last,prob=(1-0.95)/2, na.rm=T),
            e0c_upp = quantile(change_e0_last,prob=1-(1-0.95)/2, na.rm=T)) %>%
  ungroup()

aux %>%
  ggplot(aes(x = cohort,
             group = sex, 
             color = sex))+
  geom_line(aes(y=e0c_med), lwd=1.1)+
  # 95% PI
  geom_ribbon(color = "grey85",
              aes(ymin=e0c_low, ymax=e0c_upp, fill = sex), alpha=0.2)+
  geom_hline(yintercept = 0)+
  
  facet_wrap(~country, ncol = 3, scales = "free_x")+
  
  # Extra changes
  theme_classic()+
  
  labs(x = "Birth cohort",
       y = expression("Change in months "~e[0]^C),
       color = "Sex")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_color_manual("Sex", 
                     values = c("Male"="steelblue3", "Female"="firebrick2"))+
  scale_fill_manual("Sex", 
                    values = c("Male"="steelblue3", "Female"="firebrick2"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), expand=c(0, 0), 
                     limits = c(1950, 2000))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), 
                     limits = c(-70, 63))

#ggsave(path = "out/",filename = "figS2.pdf", width = 12, height = 9)
#ggsave(path = "out/",filename = "figS2.png", width = 12, height = 9)




#-------------------------------------------------------------------------------
# Figure S3 - Short-term disturbance with Linear Interpolation WPP data
# Needs codes src/01_, and 04_

# Cohort life expectancy from short-term disturbance scenario using linear 
#interpolation from WPP 2024
e0_lin_int <-
  read.table("out/04_e0_short_lin_int.txt", header = T) %>% 
  distinct()


# Temporary Effect
e0_lin_int %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12,
         country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil", "Mexico"))) %>%
  
  ggplot(aes(x = cohort,
             group = sex, 
             color = sex))+
  geom_line(aes(y=change_e0_temp), lwd=1.1)+
  theme_classic()+
  facet_wrap(~country, ncol = 3, scales = "free_x")+
  labs(x = "Birth cohort",
       y = expression("Change in months "~e[0]^C),
       color = "Sex")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_color_manual("Sex", 
                     values = c("Male"="steelblue3", "Female"="firebrick2"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), expand=c(0, 0))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  geom_hline(yintercept = 0)

#ggsave(path = "out/",filename = "figS3.pdf", width = 12, height = 9)
#ggsave(path = "out/",filename = "figS3.png", width = 12, height = 9)




#-------------------------------------------------------------------------------
# Figure S4 - Lexis diagram on the different uses of age-specific death rates


# Creating scales
dta <- data.frame(
  year = rep(seq(1920,2100,1),length(seq(0,100,1))),
  age = rep(seq(0,100,1), length(seq(1920,2100,1)))
)

# For lexis
segment_data = data.frame(
  x = seq(1920,2020, by=10),
  xend = c(seq(2020,2100, by=10), 2100, 2100), 
  y = rep(0, length(seq(1920,2020, by=10))),
  yend = c(rep(100, length(seq(1920,2000, by=10))),90, 80)
)


dta %>% 
  ggplot(aes(x=year, y=age))+
  geom_point(alpha=0)+
  geom_vline(xintercept=c(1949.8, 2021), linetype="dashed", size=1)+
  # Lexis
  geom_hline(yintercept=seq(0, 100, by = 10), size=0.1, alpha = 0.2) +
  geom_vline(xintercept=seq(1920, 2100, by = 10), size=0.1, alpha = 0.2) +
  geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),
               , size=0.1, alpha = 0.2)+
  # Texts
  annotate(geom = "text", x = 2063, y = 94, size = 15/.pt,
           label = "Assumptions", vjust=-4.5)+
  annotate(geom = "text", x = 1983, y = 94, size = 15/.pt,
           label = "WPP", vjust=-4.5)+
  annotate(geom = "text", x = 1935, y = 94, size = 15/.pt,
           label = "HMD", vjust=-4.5)+
  annotate(geom = "text", x = 2020, y = 81, size = 13/.pt,
           label = expression(e[0]^P ~ "(2020)"), vjust=-4.5)+
  annotate(geom = "text", x = 2050, y = 81, size = 13/.pt,
           label = expression(e[0]^C ~ "(1950)"), vjust=-4.5)+
  
  # Shaded area
  annotate(geom = "polygon", x = c(1920,1950,1950), y = c(0,0,30), 
           fill = "black", alpha = 0.2)+ #Light gray
  annotate(geom = "polygon", x = c(1950,2020,2020,1950), y = c(0,0,100,30), 
           fill = "black", alpha = 0.6)+ #Dark gray
  annotate("segment", x = 1950, xend = 2020, y = 0, yend = 70, 
           colour = "black", size=2)+ #Cohort observed
  annotate("segment", x = 2021, xend = 2050, y = 71, yend = 100, 
           colour = "black", size=2, alpha=0.8, linetype="longdash")+ #Coh forec
  geom_vline(xintercept=2020,size=2, colour = "firebrick2")+ #Period
  
  # Extra changes
  labs(y = "Age",
       x = "Year")+
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey30"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15), expand=c(0, 0))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 13), expand=c(0, 0))+
  coord_cartesian(clip="off")


#ggsave(path = "out/",filename = "figS4.pdf", width = 11, height = 7)
#ggsave(path = "out/",filename = "figS4.png", width = 11, height = 7)



#-------------------------------------------------------------------------------
# Figure S5 - Comparison of model performances - Short-term scenario
# Needs codes src/01_, 02_, and src/11_ ...


### Lee-Miller(2001)
## Period estimates
# Age-specific death rates
e0_period_LM <- 
  read.table("out/21_lmx_short_LM.txt", header = T)


aux_period_LM <- 
  e0_period_LM %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  summarise(e0c_med = median(e0_fore2019, na.rm=T),
            e0_observed = unique(e0_observed)) %>%
  ungroup() %>% 
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_LM <- 
  read.table("out/21_e0_short_LM.txt", header = T) %>% 
  distinct()

aux_cohort_LM <-
  e0_cohort_LM %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_LM <- 
  aux_period_LM %>% 
  left_join(aux_cohort_LM, 
            by = c("sex", "country")) %>%
  mutate(model = "Lee-Miller")




### Lee-Miller(2001) Auto Arima
## Period estimates
# Age-specific death rates
e0_period_LM_AA <- 
  read.table("out/22_lmx_short_LM_AA.txt", header = T)


aux_period_LM_AA <- 
  e0_period_LM_AA %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  summarise(e0c_med = median(e0_fore2019, na.rm=T),
            e0_observed = unique(e0_observed)) %>%
  ungroup() %>% 
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_LM_AA <- 
  read.table("out/22_e0_short_LM_AA.txt", header = T) %>% 
  distinct()

aux_cohort_LM_AA <-
  e0_cohort_LM_AA %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_LM_AA <- 
  aux_period_LM_AA %>% 
  left_join(aux_cohort_LM_AA, 
            by = c("sex", "country")) %>%
  mutate(model = "Lee-Miller AA")




## CP-splines (Camarda 2019)
## Period estimates
# Age-specific death rates
e0_period_CPSplines <- 
  read.table("out/23_lmx_short_CPsplines.txt", header = T)


aux_period_CPSplines <- 
  e0_period_CPSplines %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  reframe(e0c_med = median(e0_fore2019, na.rm=T),
          e0_observed = unique(e0_observed)) %>%
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_CPSplines <- 
  read.table("out/23_e0_short_CPsplines.txt", header = T) %>% 
  distinct()

aux_cohort_CPSplines <-
  e0_cohort_CPSplines %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_CPSplines <- 
  aux_period_CPSplines %>% 
  left_join(aux_cohort_CPSplines, 
            by = c("sex", "country")) %>%
  mutate(model = "CP Splines")




## Li-Lee
## Period estimates
# Age-specific death rates
e0_period_LL <- 
  read.table("out/24_lmx_short_LL.txt", header = T)


aux_period_LL <- 
  e0_period_LL %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  reframe(e0c_med = median(e0_fore2019, na.rm=T),
          e0_observed = unique(e0_observed)) %>%
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_LL <- 
  read.table("out/24_e0_short_LL.txt", header = T) %>% 
  distinct()

aux_cohort_LL <-
  e0_cohort_LL %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_LL <- 
  aux_period_LL %>% 
  left_join(aux_cohort_LL, 
            by = c("sex", "country")) %>%
  mutate(model = "Li-Lee")


## Li-Lee Modified
## Period estimates
# Age-specific death rates
e0_period_LL_mod <- 
  read.table("out/02_lmx_short_LL_mod.txt", header = T)


aux_period_LL_mod <- 
  e0_period_LL_mod %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  reframe(e0c_med = median(e0_fore2019, na.rm=T),
          e0_observed = unique(e0_observed)) %>%
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_LL_mod <- 
  read.table("out/02_e0_short_LL_mod.txt", header = T) %>% 
  distinct()

aux_cohort_LL_mod <-
  e0_cohort_LL_mod %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_LL_mod <- 
  aux_period_LL_mod %>% 
  left_join(aux_cohort_LL_mod, 
            by = c("sex", "country")) %>%
  mutate(model = "Li-Lee Modified")


## Li-Lee Modified Auto Arima
## Period estimates
# Age-specific death rates
e0_period_LL_mod_AA <- 
  read.table("out/25_lmx_short_LL_mod_AA.txt", header = T)


aux_period_LL_mod_AA <- 
  e0_period_LL_mod_AA %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  reframe(e0c_med = median(e0_fore2019, na.rm=T),
          e0_observed = unique(e0_observed)) %>%
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_LL_mod_AA <- 
  read.table("out/25_e0_short_LL_mod_AA.txt", header = T) %>% 
  distinct()

aux_cohort_LL_mod_AA <-
  e0_cohort_LL_mod_AA %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_LL_mod_AA <- 
  aux_period_LL_mod_AA %>% 
  left_join(aux_cohort_LL_mod_AA, 
            by = c("sex", "country")) %>%
  mutate(model = "Li-Lee Modified AA")


## Product Ratio
## Period estimates
# Age-specific death rates
e0_period_PR <- 
  read.table("out/26_lmx_short_product_ratio.txt", header = T)


aux_period_PR <- 
  e0_period_PR %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country) %>% 
  summarise(e0c_med = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
            e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_PR <- 
  read.table("out/26_e0_short_product_ratio.txt", header = T) %>% 
  distinct()

aux_cohort_PR <-
  e0_cohort_PR %>% 
  mutate(change_e0_temp = (e0_short-e0_short_baseline)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(change_e0_temp))

aux_PR <- 
  aux_period_PR %>% 
  left_join(aux_cohort_PR, 
            by = c("sex", "country")) %>%
  mutate(model = "Product Ratio")

## ADD OTHER MODELS HERE



aux <- 
  # Lee-Miller results
  aux_LM %>% 
  
  # Lee-Miller Auto Arima
  full_join(aux_LM_AA) %>% 
  
  # CP Splines
  full_join(aux_CPSplines) %>% 
  
  # Li-Lee
  full_join(aux_LL) %>% 
  
  # Li-Lee Modified
  full_join(aux_LL_mod) %>% 
  
  # Li-Lee Modified Auto Arima
  full_join(aux_LL_mod_AA) %>% 
  
  # Product ratio
  full_join(aux_PR)

## ADD OTHER MODELS HERE





# Plot
aux %>% 
  mutate(country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil","Mexico")),
         model = factor(model, c("Li-Lee Modified","Li-Lee Modified AA",
                                 "Li-Lee","Product Ratio","Lee-Miller",
                                 "Lee-Miller AA","CP Splines")),
         position = (period - cohort)/2) %>% 
  
  ggplot(aes(y=period, x=sex, group=model, color=model))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 1.5, linetype = "dashed")+
  annotate("rect", fill="grey72", alpha=0.3,
           ymin=-Inf, ymax=0, xmin=1.5, xmax=2.6)+
  
  # Line connecting period and cohort estimates
  geom_segment(aes(yend=cohort), lwd=0.8, alpha = 0.6, 
               position=
                 position_nudge(x=ifelse(aux$model=="Li-Lee Modified", -0.39,
                                  ifelse(aux$model=="Li-Lee Modified AA", -0.26,
                                  ifelse(aux$model=="Li-Lee", -0.13,
                                  ifelse(aux$model=="Product Ratio", 0,
                                  ifelse(aux$model=="Lee-Miller", 0.13,
                                  ifelse(aux$model=="Lee-Miller AA", 0.26,
                                  ifelse(aux$model=="CP Splines", 0.39, 0)))))))))+
  
  # Period and cohort estimates
  geom_point(aes(shape = "Period"), size =3,
             position=
               position_nudge(x=ifelse(aux$model=="Li-Lee Modified", -0.39,
                                ifelse(aux$model=="Li-Lee Modified AA", -0.26,
                                ifelse(aux$model=="Li-Lee", -0.13,
                                ifelse(aux$model=="Product Ratio", 0,
                                ifelse(aux$model=="Lee-Miller", 0.13,
                                ifelse(aux$model=="Lee-Miller AA", 0.26,
                                ifelse(aux$model=="CP Splines", 0.39, 0)))))))))+
  
  geom_point(aes(y=cohort, shape = "Cohort"), size =3,
             position=
               position_nudge(x=ifelse(aux$model=="Li-Lee Modified", -0.39,
                                ifelse(aux$model=="Li-Lee Modified AA", -0.26,
                                ifelse(aux$model=="Li-Lee", -0.13,
                                ifelse(aux$model=="Product Ratio", 0,
                                ifelse(aux$model=="Lee-Miller", 0.13,
                                ifelse(aux$model=="Lee-Miller AA", 0.26,
                                ifelse(aux$model=="CP Splines", 0.39, 0)))))))))+
  
  geom_text(vjust = 1.7,hjust = 0.5, size =2.7, angle=0,
            aes(label=round(period, 0)), show.legend = F, color = "Black", 
            position=
              position_nudge(x=ifelse(aux$model=="Li-Lee Modified", -0.39,
                               ifelse(aux$model=="Li-Lee Modified AA", -0.26,
                               ifelse(aux$model=="Li-Lee", -0.13,
                               ifelse(aux$model=="Product Ratio", 0,
                               ifelse(aux$model=="Lee-Miller", 0.13,
                               ifelse(aux$model=="Lee-Miller AA", 0.26,
                               ifelse(aux$model=="CP Splines", 0.39, 0)))))))))+
  
  geom_text(vjust = -1.5,hjust = 0.5, size =2.7, angle=0, 
            aes(y=cohort,label=round(cohort, 1)), show.legend=F, color="Black",
            position=
              position_nudge(x=ifelse(aux$model=="Li-Lee Modified", -0.39,
                               ifelse(aux$model=="Li-Lee Modified AA", -0.26,
                               ifelse(aux$model=="Li-Lee", -0.13,
                               ifelse(aux$model=="Product Ratio", 0,
                               ifelse(aux$model=="Lee-Miller", 0.13,
                               ifelse(aux$model=="Lee-Miller AA", 0.26,
                               ifelse(aux$model=="CP Splines", 0.39, 0)))))))))+
  facet_wrap(~country)+
  
  
  # Extra changes
  theme_classic()+
  
  labs(y = expression("Change in months " ~e[0]^P~" and "~e[0]^C),
       x = "",
       color = "Sex")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_color_manual("Model", 
                     values = c("Lee-Miller"="#0F4C5C", 
                                "Lee-Miller AA"="#73956F",
                                "Li-Lee"="#E8871E",
                                "Li-Lee Modified" = "#5F0F40",
                                "Li-Lee Modified AA" = "#4E6766",
                                "Product Ratio" = "#81B29A",
                                "CP Splines"="#9A031E"))+
  scale_shape_manual("", values = c("Period"=17, "Cohort"=16))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), expand=c(0, 6),
                     limits = c(-110,5))

#ggsave(path = "out/",filename = "figS5.pdf", width = 14, height = 8)
#ggsave(path = "out/",filename = "figS5.png", width = 14, height = 8)




#-------------------------------------------------------------------------------
# Figure S6 - Comparison of model performances - Lasting effects 
# Needs codes src/01_, 02_, and src/11_


### Lee-Miller(2001)
## Period estimates
# Age-specific death rates
e0_period_LM <- 
  read.table("out/21_lmx_short_LM.txt", header = T)


aux_period_LM <- 
  e0_period_LM %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  summarise(e0c_med = median(e0_fore2019, na.rm=T),
            e0_observed = unique(e0_observed)) %>%
  ungroup() %>% 
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_LM <- 
  read.table("out/31_e0_long_LM.txt", header = T) %>% 
  distinct()

aux_cohort_LM <-
  e0_cohort_LM %>% 
  mutate(change_e0_temp = (fore_23-fore_19)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_LM <- 
  aux_period_LM %>% 
  left_join(aux_cohort_LM, 
            by = c("sex", "country")) %>%
  mutate(model = "Lee-Miller")




### Lee-Miller(2001) Auto Arima
## Period estimates
# Age-specific death rates
e0_period_LM_AA <- 
  read.table("out/22_lmx_short_LM_AA.txt", header = T)


aux_period_LM_AA <- 
  e0_period_LM_AA %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  summarise(e0c_med = median(e0_fore2019, na.rm=T),
            e0_observed = unique(e0_observed)) %>%
  ungroup() %>% 
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_LM_AA <- 
  read.table("out/32_e0_long_LM_AA.txt", header = T) %>% 
  distinct()

aux_cohort_LM_AA <-
  e0_cohort_LM_AA %>% 
  mutate(change_e0_temp = (fore_23-fore_19)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_LM_AA <- 
  aux_period_LM_AA %>% 
  left_join(aux_cohort_LM_AA, 
            by = c("sex", "country")) %>%
  mutate(model = "Lee-Miller AA")




## CP-splines (Camarda 2019)
## Period estimates
# Age-specific death rates
e0_period_CPSplines <- 
  read.table("out/23_lmx_short_CPsplines.txt", header = T)


aux_period_CPSplines <- 
  e0_period_CPSplines %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  reframe(e0c_med = median(e0_fore2019, na.rm=T),
          e0_observed = unique(e0_observed)) %>%
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_CPSplines <- 
  read.table("out/33_e0_long_CPsplines.txt", header = T) %>% 
  distinct()

aux_cohort_CPSplines <-
  e0_cohort_CPSplines %>% 
  mutate(change_e0_temp = (fore_23-fore_19)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_CPSplines <- 
  aux_period_CPSplines %>% 
  left_join(aux_cohort_CPSplines, 
            by = c("sex", "country")) %>%
  mutate(model = "CP Splines")




## Li-Lee
## Period estimates
# Age-specific death rates
e0_period_LL <- 
  read.table("out/24_lmx_short_LL.txt", header = T)


aux_period_LL <- 
  e0_period_LL %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  reframe(e0c_med = median(e0_fore2019, na.rm=T),
          e0_observed = unique(e0_observed)) %>%
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_LL <- 
  read.table("out/34_e0_long_LL.txt", header = T) %>% 
  distinct()

aux_cohort_LL <-
  e0_cohort_LL %>% 
  mutate(change_e0_temp = (fore_23-fore_19)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_LL <- 
  aux_period_LL %>% 
  left_join(aux_cohort_LL, 
            by = c("sex", "country")) %>%
  mutate(model = "Li-Lee")


## Li-Lee Modified
## Period estimates
# Age-specific death rates
e0_period_LL_mod <- 
  read.table("out/02_lmx_short_LL_mod.txt", header = T)


aux_period_LL_mod <- 
  e0_period_LL_mod %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  reframe(e0c_med = median(e0_fore2019, na.rm=T),
          e0_observed = unique(e0_observed)) %>%
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_LL_mod <- 
  read.table("out/03_e0_long_LL_mod.txt", header = T) %>% 
  distinct()

aux_cohort_LL_mod <-
  e0_cohort_LL_mod %>% 
  mutate(change_e0_temp = (fore_23-fore_19)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_LL_mod <- 
  aux_period_LL_mod %>% 
  left_join(aux_cohort_LL_mod, 
            by = c("sex", "country")) %>%
  mutate(model = "Li-Lee Modified")


## Li-Lee Modified Auto Arima
## Period estimates
# Age-specific death rates
e0_period_LL_mod_AA <- 
  read.table("out/25_lmx_short_LL_mod_AA.txt", header = T)


aux_period_LL_mod_AA <- 
  e0_period_LL_mod_AA %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country, simulation) %>% 
  mutate(e0_fore2019 = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
         e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  # Getting medium estimates
  group_by(country, sex, year) %>% 
  reframe(e0c_med = median(e0_fore2019, na.rm=T),
          e0_observed = unique(e0_observed)) %>%
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_LL_mod_AA <- 
  read.table("out/35_e0_long_LL_mod_AA.txt", header = T) %>% 
  distinct()

aux_cohort_LL_mod_AA <-
  e0_cohort_LL_mod_AA %>% 
  mutate(change_e0_temp = (fore_23-fore_19)*12) %>% 
  group_by(country, sex, cohort) %>% 
  summarise(e0c_med = median(change_e0_temp, na.rm=T)) %>%
  ungroup() %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(e0c_med))

aux_LL_mod_AA <- 
  aux_period_LL_mod_AA %>% 
  left_join(aux_cohort_LL_mod_AA, 
            by = c("sex", "country")) %>%
  mutate(model = "Li-Lee Modified AA")


## Product Ratio
## Period estimates
# Age-specific death rates
e0_period_PR <- 
  read.table("out/26_lmx_short_product_ratio.txt", header = T)


aux_period_PR <- 
  e0_period_PR %>% 
  filter(year >=2020 & year <= 2023) %>% 
  # Life Expectancy calculations
  group_by(year, sex, country) %>% 
  summarise(e0c_med = lifetable.e0(x= seq(0,100), mx= exp(lmx_short_baseline)),
            e0_observed = lifetable.e0(x= seq(0,100), mx= exp(lmx_short))) %>% 
  ungroup() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  arrange(change_e0_period_med) %>% 
  slice_head(n=1) %>% 
  dplyr::select(country, sex, year, period=change_e0_period_med)


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort_PR <- 
  read.table("out/36_e0_long_product_ratio.txt", header = T) %>% 
  distinct()

aux_cohort_PR <-
  e0_cohort_PR %>% 
  mutate(change_e0_temp = (fore_23-fore_19)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(cohort = min(change_e0_temp))

aux_PR <- 
  aux_period_PR %>% 
  left_join(aux_cohort_PR, 
            by = c("sex", "country")) %>%
  mutate(model = "Product Ratio")

## ADD OTHER MODELS HERE



aux <- 
  # Lee-Miller results
  aux_LM %>% 
  
  # Lee-Miller Auto Arima
  full_join(aux_LM_AA) %>% 
  
  # CP Splines
  full_join(aux_CPSplines) %>% 
  
  # Li-Lee
  full_join(aux_LL) %>% 
  
  # Li-Lee Modified
  full_join(aux_LL_mod) %>% 
  
  # Li-Lee Modified Auto Arima
  full_join(aux_LL_mod_AA) %>% 
  
  # Product ratio
  full_join(aux_PR)

## ADD OTHER MODELS HERE





# Plot
aux %>% 
  mutate(country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil","Mexico")),
         model = factor(model, c("Li-Lee Modified","Li-Lee Modified AA",
                                 "Li-Lee","Product Ratio","Lee-Miller",
                                 "Lee-Miller AA","CP Splines")),
         position = (period - cohort)/2) %>% 
  
  ggplot(aes(y=period, x=sex, group=model, color=model))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 1.5, linetype = "dashed")+
  annotate("rect", fill="grey72", alpha=0.3,
           ymin=-Inf, ymax=0, xmin=1.5, xmax=2.6)+
  
  # Line connecting period and cohort estimates
  geom_segment(aes(yend=cohort), lwd=0.8, alpha = 0.6, 
               position=
                 position_nudge(x=ifelse(aux$model=="Li-Lee Modified", -0.39,
                                  ifelse(aux$model=="Li-Lee Modified AA", -0.26,
                                  ifelse(aux$model=="Li-Lee", -0.13,
                                  ifelse(aux$model=="Product Ratio", 0,
                                  ifelse(aux$model=="Lee-Miller", 0.13,
                                  ifelse(aux$model=="Lee-Miller AA", 0.26,
                                  ifelse(aux$model=="CP Splines", 0.39, 0)))))))))+
  
  # Period and cohort estimates
  geom_point(aes(shape = "Period"), size =3,
             position=
               position_nudge(x=ifelse(aux$model=="Li-Lee Modified", -0.39,
                                ifelse(aux$model=="Li-Lee Modified AA", -0.26,
                                ifelse(aux$model=="Li-Lee", -0.13,
                                ifelse(aux$model=="Product Ratio", 0,
                                ifelse(aux$model=="Lee-Miller", 0.13,
                                ifelse(aux$model=="Lee-Miller AA", 0.26,
                                ifelse(aux$model=="CP Splines", 0.39, 0)))))))))+
  
  geom_point(aes(y=cohort, shape = "Cohort"), size =3,
             position=
               position_nudge(x=ifelse(aux$model=="Li-Lee Modified", -0.39,
                                ifelse(aux$model=="Li-Lee Modified AA", -0.26,
                                ifelse(aux$model=="Li-Lee", -0.13,
                                ifelse(aux$model=="Product Ratio", 0,
                                ifelse(aux$model=="Lee-Miller", 0.13,
                                ifelse(aux$model=="Lee-Miller AA", 0.26,
                                ifelse(aux$model=="CP Splines", 0.39, 0)))))))))+
  
  geom_text(vjust = 1.7,hjust = 0.5, size =2.7, angle=0,
            aes(label=round(period, 0)), show.legend = F, color = "Black", 
            position=
              position_nudge(x=ifelse(aux$model=="Li-Lee Modified", -0.39,
                               ifelse(aux$model=="Li-Lee Modified AA", -0.26,
                               ifelse(aux$model=="Li-Lee", -0.13,
                               ifelse(aux$model=="Product Ratio", 0,
                               ifelse(aux$model=="Lee-Miller", 0.13,
                               ifelse(aux$model=="Lee-Miller AA", 0.26,
                               ifelse(aux$model=="CP Splines", 0.39, 0)))))))))+
  
  geom_text(vjust = -1.5,hjust = 0.5, size =2.7, angle=0, 
            aes(y=cohort,label=round(cohort, 0)), show.legend=F, color="Black",
            position=
              position_nudge(x=ifelse(aux$model=="Li-Lee Modified", -0.39,
                               ifelse(aux$model=="Li-Lee Modified AA", -0.26,
                               ifelse(aux$model=="Li-Lee", -0.13,
                               ifelse(aux$model=="Product Ratio", 0,
                               ifelse(aux$model=="Lee-Miller", 0.13,
                               ifelse(aux$model=="Lee-Miller AA", 0.26,
                               ifelse(aux$model=="CP Splines", 0.39, 0)))))))))+
  facet_wrap(~country)+
  
  
  # Extra changes
  theme_classic()+
  
  labs(y = expression("Change in months " ~e[0]^P~" and "~e[0]^C),
       x = "",
       color = "Sex")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_color_manual("Model", 
                     values = c("Lee-Miller"="#0F4C5C", 
                                "Lee-Miller AA"="#73956F",
                                "Li-Lee"="#E8871E",
                                "Li-Lee Modified" = "#5F0F40",
                                "Li-Lee Modified AA" = "#4E6766",
                                "Product Ratio" = "#81B29A",
                                "CP Splines"="#9A031E"))+
  scale_shape_manual("", values = c("Period"=17, "Cohort"=16))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), expand=c(0, 6),
                     limits = c(-110,13))

#ggsave(path = "out/",filename = "figS6.pdf", width = 14, height = 8)
#ggsave(path = "out/",filename = "figS6.png", width = 14, height = 8)



#-------------------------------------------------------------------------------
# Figure S7 - Comparison of period le across models
# Needs codes src/01_, 02_, and src/11_

## loading data - Li-Lee and Lee-Miller
load("out/24_period_e0_LL.Rdata")

## renaming models
plot_data <- period_e0 %>%
  rename(`Lee-Miller`=LM, `Li-Lee`=LL)


## loading data - Li-Lee-Modified
load("out/02_period_e0_LL_mod.Rdata")
period_e0 <- period_e0 %>% 
  rename(`Li-Lee Modified`=LL) %>% 
  select(-LM)

## merging
plot_data <- plot_data %>% 
  left_join(period_e0)


## loading data - Li-Miller AA
load("out/22_period_e0_LM_AA.Rdata")
period_e0 <- period_e0 %>% 
  rename(`Lee-Miller AA`=LM_AA) %>% 
  select(-LM) %>% 
  mutate(country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil","Mexico")))

## merging
plot_data <- plot_data %>% 
  left_join(period_e0)



plot_data <- plot_data %>%
  pivot_longer(
    cols = c(`Lee-Miller`, `Lee-Miller AA`, `Li-Lee`, `Li-Lee Modified`),
    names_to = "model",
    values_to = "forecast"
  ) %>% 
  mutate(country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil","Mexico")),
         model = factor(model, c("Lee-Miller","Lee-Miller AA",
                                 "Li-Lee","Li-Lee Modified")))


## plotting
ggplot() +
  
  # Observed life expectancy
  geom_point(
    data = period_e0,
    aes(x = year, y = e0, shape = sex),
    size = 1
  ) +
  
  # Forecasts
  geom_line(
    data = plot_data,
    aes(
      x = year,
      y = forecast,
      colour = model,
      group = interaction(sex, model)
    ),
    linewidth = 0.8
  ) +
  
  facet_wrap(~country, scales = "free_y") +
  
  scale_shape_manual(
    values = c("Female" = 21, "Male" = 22)
  ) +
  
  scale_colour_manual(
    values = c(
      `Lee-Miller` = "#0F4C5C",
      `Lee-Miller AA` = "#73956F",
      `Li-Lee` = "#E8871E",
      `Li-Lee Modified` = "#5F0F40"
    )
  ) +
  
  geom_vline(
    xintercept = 2019.5,
    linetype = "dashed",
    colour = "grey50"
  ) +
  
  # Extra changes
  theme_classic()+
  
  labs(x = "Year",
       y = "Life expectancy at birth",
       shape = "Sex",
       colour = "Model")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))

#ggsave(path = "out/",filename = "figS7.pdf", width = 12, height = 8)
#ggsave(path = "out/",filename = "figS7.png", width = 12, height = 8)




#-------------------------------------------------------------------------------
# Figure S8 - Simulations
# Needs codes src/01_, 02_, and src/11_

e0_period <- 
  read.table('out/06_e0_period_sim.txt', header = T)

e0_cohort <- 
  read.table('out/06_e0_cohort_sim.txt', header = T)


g1 <- 
  e0_period %>% 
  mutate(n_sim = factor(n_sim),
         country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil","Mexico"))) %>% 
  ggplot(aes(x=sex, group=n_sim, color=n_sim))+
  # Period
  geom_point(aes(y = change_p_med, shape = "Period"), size =2,
             show.legend = F,
             position=position_nudge(x = ifelse(e0_period$n_sim=="100", -0.3,
                                                ifelse(e0_period$n_sim=="1000", 0,
                                                       ifelse(e0_period$n_sim=="10000", 0.3, 0)))))+
  # Period CI
  geom_errorbar(aes(y = change_p_med, ymin = change_p_low, ymax = change_p_upp), 
                width = 0.2, lwd=0.7,
                position=position_nudge(x = ifelse(e0_period$n_sim=="100", -0.3,
                                                   ifelse(e0_period$n_sim=="1000", 0,
                                                          ifelse(e0_period$n_sim=="10000", 0.3, 0)))))+
  
  facet_wrap(~country, labeller = labeller(country = label_wrap_gen(width = 17)))+
  # Extra changes
  theme_classic()+
  
  labs(y = expression("Change in months " ~e[0]^P),
       x = "", title = "a. Period")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_color_manual("Number of simulations", 
                     values = c("100"="#0F4C5C", 
                                "1000"="#9A031E",
                                "10000"="#E8871E"))+
  scale_shape_manual("", values = c("Period"=17, "Cohort"=16))



g2 <- e0_cohort %>% 
  mutate(n_sim = factor(n_sim),
         country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil","Mexico"))) %>% 
  ggplot(aes(x=sex, group=n_sim, color=n_sim))+
  # Cohort
  geom_point(aes(y = change_c_med, shape = "Cohort"), size =2,
             show.legend = F,
             position=position_nudge(x = ifelse(e0_cohort$n_sim=="100", -0.3,
                                                ifelse(e0_cohort$n_sim=="1000", 0,
                                                       ifelse(e0_cohort$n_sim=="10000", 0.3, 0)))))+
  # Cohort CI
  geom_errorbar(aes(y = change_c_med, ymin = change_c_low, ymax = change_c_upp), 
                width = 0.2, lwd=0.7,
                position=position_nudge(x = ifelse(e0_cohort$n_sim=="100", -0.3,
                                                   ifelse(e0_cohort$n_sim=="1000", 0,
                                                          ifelse(e0_cohort$n_sim=="10000", 0.3, 0)))))+
  
  facet_wrap(~country, labeller = labeller(country = label_wrap_gen(width = 17)))+
  # Extra changes
  theme_classic()+
  
  labs(y = expression("Change in months " ~e[0]^C),
       x = "", title = "b. Cohort")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_color_manual("Number of simulations", 
                     values = c("100"="#0F4C5C", 
                                "1000"="#9A031E",
                                "10000"="#E8871E"))+
  scale_shape_manual("", values = c("Period"=17, "Cohort"=16))



ggarrange(g1, g2, ncol =2,common.legend = TRUE, legend="bottom")

#ggsave(path = "out/",filename = "figS8.pdf", width = 14, height = 8)
#ggsave(path = "out/",filename = "figS8.png", width = 14, height = 8)




#-------------------------------------------------------------------------------
# Figure S9 - Relative difference between WPP and HMD

# Data WPP
mx1dt_wpp <- 
  read.table("dat/01_wpp_cleaned.txt", header = T) %>% 
  filter(year == 1950) %>% 
  select(country, age, sex, mx_wpp=mx)


# Data HMD
# HMD
# Death rates Mx_1x1 for USA, ITA, JPN and NZL_NP
mx1dt_hmd <- 
  read.table("dat/HMD_Mx_1x1_acessed01102025.txt", header = T) %>% 
  mutate(CNTRY = case_when(CNTRY == "GBR_NP" ~ "GBR",
                           CNTRY == "NZL_NP" ~ "NZL",
                           CNTRY == "DEUTNP" ~ "DEU",
                           CNTRY == "FRATNP" ~ "FRA",
                           TRUE ~ CNTRY),
         country = countrycode(CNTRY, origin = "iso3c",
                               destination = "un.name.en")) %>% 
  select(country, year=Year, age=Age, Male, Female) %>% 
  pivot_longer(-c(country, year, age), names_to = "sex", values_to = "mx") %>% 
  filter(age <= 100, year == 1950) %>% 
  select(country, age ,sex, mx_hmd=mx)

mx1dt_all <- 
  mx1dt_hmd %>% 
  left_join(mx1dt_wpp, 
            by = c("country","age","sex")) %>% 
  mutate(ratio = (mx_hmd-mx_wpp)/mx_wpp)

mx1dt_all %>% 
  mutate(country = factor(country, c("New Zealand","Japan","Italy",
                                     "United States of America",
                                     "Brazil","Mexico"))) %>% 
  ggplot(aes(x=age, y=ratio, group=sex, color=sex))+
  
  geom_hline(yintercept = 0, alpha = 0.7)+
  
  geom_line(lwd=1.1)+
  facet_wrap(~country, ncol=2)+
  
  
  # Extra changes
  theme_classic()+
  
  labs(x = "Age",
       y = "Relative Difference",
       color = "Sex")+
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        axis.title.x = element_text(vjust = -1),
        axis.line = element_line(colour = "grey70"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(color = "grey90",linewidth = 0.1), 
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.7,"cm"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 
                           "inches"))+
  scale_color_manual("Sex", 
                     values = c("Male"="steelblue3", "Female"="firebrick2"))+
  
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand=c(0, 0))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  
  
#ggsave(path = "out/",filename = "figS9.pdf", width = 12, height = 9)
#ggsave(path = "out/",filename = "figS9.png", width = 12, height = 9)

















