
#'How period life expectancy can distort our interpretation of mortality crises'
# Figures

# This code is to be used to replicate the figures of the paper "How period life
#expectancy can distort our interpretation of mortality crises".



#-------------------------------------------------------------------------------

# Cleaning the workspace
rm(list=ls(all=TRUE))

# Packages
require(tidyverse)
require(ggpubr)

# Functions
source("fun/functions.R")
# List of functions includes 'lifetable.e0' for calculating life expectancy at 
#birth and 'lin_inter' to obtain linear interpolation between two years.




#-------------------------------------------------------------------------------
# Figure 1 - Trends in period life expectancy
# Needs codes src/01_ and 02_

# Age-specific death rates
e0_period <- 
  read.table("out/lmx_short.txt", header = T)


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
            e0_observed = e0_observed) %>%
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




#-------------------------------------------------------------------------------
# Figure 2 - Short-term disturbance
# Needs codes src/01_, and 02_

# Cohort life expectancy from short-term disturbance scenario
e0_short <- 
  read.table("out/e0_short.txt", header = T) %>% 
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




#-------------------------------------------------------------------------------
# Figure 3 - Ratio Period/Cohort
# Needs codes src/01_, and 02_

## Period estimates
# Age-specific death rates
e0_period <- 
  read.table("out/lmx_short.txt", header = T)


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
            e0_observed = e0_observed) %>%
  ungroup() %>% 
  distinct() %>% 
  mutate(change_e0_period_med = (e0_observed-e0c_med)*12) %>% 
  group_by(sex,country) %>% 
  # Most substantial impact
  summarise(period = min(change_e0_period_med))


## Cohort estimates
# Cohort life expectancy from short-term disturbance scenario
e0_cohort <- 
  read.table("out/e0_short.txt", header = T) %>% 
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
  geom_text(vjust = 0,hjust = 1.2, size =4, angle=0, aes(y=position),
            show.legend = F,fontface = "bold",
            position=position_nudge(x = ifelse(aux$sex == "Male", 0.2, -0.2)))+
  
  geom_text(vjust = 0.9,hjust = 1.4, size =2.7, angle=0,
            aes(label=round(period, 0)), show.legend = F, color = "Black", 
            position=position_nudge(x = ifelse(aux$sex == "Male", 0.2, -0.2)))+
  
  geom_text(vjust = 0.9,hjust = 1.4, size =2.7, angle=0, 
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




#-------------------------------------------------------------------------------
# Figure 4 - Age-specific death rates for impact measurement
# Needs codes src/01_, 02_, and 03_

lmx_short <- 
  read.table("out/lmx_short.txt", header = T) 

lmx_long <- 
  read.table("out/lmx_long.txt", header = T) 


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
  
  geom_line(data=subset(aux_lmx60,year<=2024),
            aes(y=lmx_short, linetype="Observed"), lwd=1.1, color="steelblue3")+
  geom_line(data=subset(aux_lmx60,year>=2019),
            aes(y=lmx_short_baseline, linetype="Simulations", group=simulation),
            lwd=1.1, color="grey50", alpha = 0.2)+
  geom_line(data=subset(aux_lmx60,year>=2019),
            aes(y=lmx_med, linetype="Medium estimate"), lwd=1.1, 
            color = "red")+
  # Impact
  geom_ribbon(data=subset(aux_lmx60,year<=2024),
              aes(ymin=lmx_med, ymax=lmx_short, fill = "Impact"), alpha=0.3)+
  
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
  scale_fill_manual("", values = c("Impact"="orangered1"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand=c(0, 0),
                     limits = c(2017,2027))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8),
                     limits = c(-4.68, -4))




# Long-term

# Plot
g2 <-
  aux_lmx60_all  %>% 
  ggplot(aes(x=year))+
  
  # Observed
  geom_line(data=subset(aux_lmx60_all,year<=2023),
            aes(y=lmx_short, linetype="Observed"), lwd=1.1, color="steelblue3")+
  # Simulations
  geom_line(data=subset(aux_lmx60_all,year>=2019),
            aes(y=lmx_short_baseline, linetype="Simulations", group=simulation),
            lwd=1.1, color="grey50", alpha = 0.2)+
  # Fore 1050-2023
  geom_line(data=subset(aux_lmx60_all,year>=2023),
            aes(y=fore_23, linetype="Simulations", group = simulation_23),
            lwd=1.1, color="grey50", alpha = 0.2)+
  geom_line(data=subset(aux_lmx60_all,year>=2019),
            aes(y=lmx_med, linetype="Medium estimate"), 
            lwd=1.1, color="red")+
  
  geom_line(data=subset(aux_lmx60_all,year>=2023),
            aes(y=lmx_med_last, linetype="Medium estimate"), 
            lwd=1.1, color="red")+
  # Impact
  geom_ribbon(aes(ymin=lmx_med_last, ymax=lmx_med, fill = "Impact"), alpha=0.3)+
  
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
  scale_fill_manual("", values = c("Impact"="orangered1"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand=c(0, 0),
                     limits = c(2017,2027))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8),
                     limits = c(-4.68, -4))


ggarrange(g1, g2,ncol = 2, common.legend = TRUE, legend="bottom") 

#ggsave(path = "out/",filename = "fig4.pdf", width = 14, height = 6)




#-------------------------------------------------------------------------------
# Figure S1 - Long-term disturbance
# Needs codes src/01_, 02_, and 03_

# Cohort life expectancy from long-term disturbance scenario
e0_long <-
  read.table("out/e0_long.txt", header = T) %>% 
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
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))

#ggsave(path = "out/",filename = "figS1.pdf", width = 12, height = 9)




#-------------------------------------------------------------------------------
# Figure S2 - Short-term disturbance with Linear Interpolation WPP data
# Needs codes src/01_, and 04_

# Cohort life expectancy from short-term disturbance scenario using linear 
#interpolation from WPP 2024
e0_lin_int <-
  read.table("out/e0_short_lin_int.txt", header = T) %>% 
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

#ggsave(path = "out/",filename = "figS2.pdf", width = 12, height = 9)




#-------------------------------------------------------------------------------
# Figure S3 - Short-term disturbance with HMD data
# Needs codes src/01_, and 05_

# Cohort life expectancy from short-term disturbance scenario including HMD
e0_hmd <- 
  read.table("out/e0_short_hmd.txt", header = T) %>% 
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

#ggsave(path = "out/",filename = "figS3.pdf", width = 12, height = 9)










