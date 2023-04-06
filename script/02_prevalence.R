###################################
###   Prevalence of infections  ###
###      Created: 07/11/21      ###
###      Updated: 07/11/21      ###
###################################

library(dplyr)
# library(readr)
# library(haven)
# library(stringr)
library(lubridate)
# library(gmodels)
# library(reshape2)
library(ggplot2)
library(tidyr)
library(patchwork)

# load("data/ftd_h120622.rda")
# load("data/ftd_c120622.rda")
load("data/ftd_h122222.rda")
load("data/ftd_c020722.rda")

source("script/ftd_functions.R")

# Create datasets of viral prevalence for plotting
site_ftd_c <- prev_func(ftd_c,"country_id")
site_ftd_h <- prev_func(ftd_h,"country_id")
site_sex_ftd_c <- prev_func(ftd_c,"country_id", "p_sex")
site_sex_ftd_h <- prev_func(ftd_h,"country_id", "p_sex")
site_age_ftd_c <- prev_func(ftd_c,"country_id", "hosp_infage_cat")%>%
  #add to data_proc
  mutate(age_cat = case_when(
    hosp_infage_cat == 0 ~ "<14 weeks",
    hosp_infage_cat == 1 ~ "14-26 weeks",
    hosp_infage_cat == 2 ~ "27-39 weeks",
    hosp_infage_cat == 3 ~ "40-52 weeks",
  ))
site_age_ftd_h <- prev_func(ftd_h,"country_id", "hosp_infage_cat")%>%
  mutate(age_cat = case_when(
    hosp_infage_cat == 0 ~ "<14 weeks",
    hosp_infage_cat == 1 ~ "14-26 weeks",
    hosp_infage_cat == 2 ~ "27-39 weeks",
    hosp_infage_cat == 3 ~ "40-52 weeks",
  ))

# Create tables of viral prevalence
top_five_h <- site_ftd_h%>%
  group_by(country_id)%>%
  arrange(country_id, -prevalence)%>%
  slice_head(n=5)%>%
  mutate(cumm_prev = cumsum(prevalence))

top_five_c <- site_ftd_c%>%
  group_by(country_id)%>%
  arrange(country_id, -prevalence)%>%
  slice_head(n=5)%>%
  mutate(cumm_prev = cumsum(prevalence))

# Plot prevalence by country and faceted by sex and age
site_c <- prev_bar(site_ftd_c)+
  ggtitle("B")+
  theme(plot.title = element_text(size = 20))
site_h <- prev_bar(site_ftd_h)+
  ggtitle("A")+
  theme(plot.title = element_text(size = 20))
prev_bar(site_sex_ftd_c, "p_sex")
prev_bar(site_sex_ftd_h, "p_sex")
prev_bar(site_age_ftd_c, "age_cat")
prev_bar(site_age_ftd_h, "age_cat")

prev_panel <- (site_h|site_c)+
  plot_layout(guides = "collect")
  
png(file="figures/prev_panelV2.png", width=6000,height=3000,res = 300)
prev_panel
dev.off()

png(file="figures/prev_age_hV2.png", width=5000,height=3000,res = 300)
prev_bar(site_age_ftd_h, "age_cat")
dev.off()

png(file="figures/prev_age_cV2.png", width=5000,height=3000,res = 300)
prev_bar(site_age_ftd_c, "age_cat")
dev.off()

ggplot(site_age_ftd_h, aes(fill=country_id, y=prevalence, x=virus2)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 40, 10))+
  coord_flip()+
  labs(y = "Prevalence", x = " ", size=1)+
  scale_fill_manual(values = cols4, breaks = country_id, labels = country_id, name=NULL)+
  theme_classic()+
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        legend.text = element_text(size=20))+
  # facet_grid(cols = vars(age_cat))
  facet_wrap(as.formula(sprintf('~%s',"age_cat")),nrow = 1)

# get mean age and iqr for overall cohort
ftd_h_sub <- select(ftd_h, study_id, p_dob, p_date, rsv_result, rv_ev,
                    hadv, infect_d, p_sex, )%>%
  ungroup()%>%
  mutate(age_m = interval(p_dob, p_date)/months(1))

ftd_c_sub <- select(ftd_c, study_id, p_dob, p_date, rsv_result, rv_ev,
                    hadv, infect_d)%>%
  ungroup()%>%
  mutate(age_m = interval(p_dob, p_date)/months(1))

ftd_tot <- bind_rows(ftd_h_sub, ftd_c_sub)%>%
  mutate(age_wk = interval(p_dob, p_date)/weeks(1))

summary(ftd_tot$age_m)
7.032-1.207

ftd_tot%>%
  filter(rsv_result==1)%>%
  ggplot(aes(age_wk))+
geom_histogram()

table(ftd_tot$infect_d)
#### Extra code
# top_ten_h <- site_ftd_h%>%
#   group_by(country_id)%>%
#   arrange(country_id, -prevalence)%>%
#   slice_head(n=10)
# 
# top_ten_c <- site_ftd_c%>%
#   group_by(country_id)%>%
#   arrange(country_id, -prevalence)%>%
#   slice_head(n=10)
