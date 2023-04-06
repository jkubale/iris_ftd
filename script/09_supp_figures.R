# Supplemental figures
## Create and save as .svg for pasting into quarto doc as resizing too finicky

library(reshape2)
library(scico)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggupset)
library(stringr)
library(kableExtra)

load("data/ftd_h120622.rda") # ../ goes up a level since qmd. root is different
load("data/ftd_c120622.rda")

source("script/ftd_functions.R")


# Supplemental Figure 1 ----
df2 <- suppressMessages(ftd_h %>%
                          as_tibble() %>%
                          select(study_id, all_of(vir_ord))%>%
                          rename_with(~  vir_labs, .cols = which(colnames(.)%in%vir_ord))%>%
                          pivot_longer(-study_id,
                                       names_to = "infection", 
                                       values_to = "value")%>%
                          filter(value==1) %>%
                          select(-value)%>%
                          group_by(study_id)%>%
                          summarise(infection = list(infection))%>%
                          filter(lengths(infection)>1))



# create upset plot
svg(file = "figures/supp_fig1.svg", width = 12, height = 10)
suppressMessages(df2%>%
                   ggplot(aes(x=infection)) + # specify list variable (combinations to count)
                   geom_bar() + # makes bar plot of combination freq, specifies color (here using hex code), can remove for default 
                   scale_y_continuous(limits=c(0, 175), breaks=seq(0, 175, 25))+
                   ylab("Count")+ # specify y axis label
                   xlab(" ")+ # specify x axis label
                   scale_x_upset(n_intersections = 10)+ # how many of most frequent combinations to display (here displays top 15)
                   theme_classic()+ 
                   theme(axis.title.y = element_text(size = 26), # specify size of y axis title
                         axis.text.y = element_text(size = 26))+ # specify size of x axis title
                   theme_combmatrix(combmatrix.label.text = element_text(size=26)))
dev.off()

# Supplemental Figure 2 ----
ftd_tot_sub <- subset_merge()


ftd_tot_sub <- ftd_tot_sub%>%
  mutate(study_wk = as.numeric(floor(interval("2015-06-26", e_specimen_resp_date)/weeks(1)+1)),
         study_m = as.numeric(floor(interval("2015-06-26", e_specimen_resp_date)/months(1)+1)))

m_dates <- data.frame("date" = seq(as.Date("2015-06-21", format="%Y-%m-%d"), 
                                   as.Date("2017-04-20", format="%Y-%m-%d"), 
                                   by="months"))

m_dates <- m_dates%>%
  mutate(study_m = row_number())

dates <- m_dates[,1]
study_m <- m_dates[,2]
country_id = c("Albania","Jordan","Nicaragua","Philippines")
m_dates2 <- expand.grid(study_m, country_id, vir_ord)
colnames(m_dates2) <- c("study_m", "country_id", "variable")
m_dates3 <- left_join(m_dates2, m_dates, by = "study_m")

ftd_tot_long <- ftd_tot_sub%>%
  group_by(country_id, study_m)%>%
  summarise(across(.cols = all_of(vir_ord),sum),
            tot_samp = n())%>%
  select(-tot_samp)%>%
  melt(id.vars = c("country_id", "study_m"))%>%
  ungroup()%>%
  group_by(country_id, variable)%>%
  mutate(tot_num = sum(value),
         max_wk_num = max(value),
         prop_tot = round(value/tot_num, 3),
         prop_max = round(value/max_wk_num, 3))%>% # proportion of max monthly positives of virus by site
  ungroup()

ftd_tot_long2 <- left_join(m_dates3, ftd_tot_long, by=c("country_id", "study_m", "variable"))%>%
  mutate(
    value_cat = factor(case_when(
      is.na(value) ~ "None Collected",
      value == 0 ~ "0",
      value < 5 ~ "1-4",
      value < 10 & value > 4 ~ "5-9",
      value < 15 & value > 9 ~ "10-14",
      value < 20 & value > 14 ~ "15-19",
      value < 25 & value > 19 ~ "20-24",
      value >= 25 ~ "25+"
    ),
    levels = c("None Collected","0","1-4","5-9","10-14","15-19","20-24","25+")),
    variable2 = factor(variable, 
                       levels = vir_ord,
                       labels = vir_labs))

svg(file = "figures/supp_fig2.svg", width = 22, height = 15)
ggplot(ftd_tot_long2,aes(date,variable2,fill=prop_max))+
  geom_tile(color= "white",linewidth=0.1, size=2) + 
  scale_fill_gradient(low = "#023A7A",
                      high = "#9F711B",
                      na.value = "white",
                      name = NULL)+
  scale_x_date(date_labels="%b-%y", date_breaks="2 months", expand=c(0,0)) + 
  theme_bw()+
  theme(axis.text.y = element_text(size = 26),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 26, angle=90, vjust = 0.5),
        # strip.text.x = element_text(size = 30),
        legend.text=element_text(size=26)
  )+
  guides(guide_colorsteps(frame.colour = "black"))+
  facet_grid(~country_id)+
  theme(strip.text = element_text(size=30))
dev.off()

## Supplemental Figure 3 ----
site_eti <- read.csv("data/site_estimates.csv")%>%
  mutate(eti_prop = round(eti_mean, 3)*100,
         lcl_prop = round(ci_025, 3)*100,
         ucl_prop = round(ci_975, 3)*100)

svg(file = "figures/supp_fig3V2.svg", width = 22, height = 15)
site_eti%>%
  mutate(cause = factor(case_when(
    cause=="RV_EV" ~ "RV/EV",
    cause=="other" ~ "Other",
    T ~ cause
  ), ordered = F, levels=c(all_of(vir_labs3))))%>%
  ggplot(aes(y=eti_prop, ymin = lcl_prop, ymax = ucl_prop, x=cause, color=site)) + 
  geom_pointrange(position=position_dodge(width = 0.30) , size=2, linewidth=2)+
  scale_y_continuous(limits=c(0, 100), breaks=seq(0, 80, 20))+
  coord_flip()+
  labs(y = "Population Etiologic Fraction", x = " ", size=1)+
  scale_color_manual(values = cols4, breaks = c("Albania", "Jordan", "Nicaragua", "Philippines"), labels = c("Albania", "Jordan", "Nicaragua", "Philippines"), name=NULL)+
  theme_bw()+
  theme(axis.text = element_text(size = 26),
        axis.title.x = element_text(size = 28),
        legend.position = "none")+
  facet_grid(~site)+
  theme(strip.text = element_text(size=30))
dev.off()

## Sensitivity analysis plot of ETI ----
site_eti <- read.csv("data/site_estimates_sub.csv")%>%
  mutate(eti_prop = round(eti_mean, 3)*100,
         lcl_prop = round(ci_025, 3)*100,
         ucl_prop = round(ci_975, 3)*100)

svg(file = "figures/sens_site_etiV2.svg", width = 22, height = 15)
site_eti%>%
  mutate(cause = factor(case_when(
    cause=="RV_EV" ~ "RV/EV",
    cause=="other" ~ "Other",
    T ~ cause
  ), ordered = F, levels=c(all_of(vir_labs3))))%>%
  ggplot(aes(y=eti_prop, ymin = lcl_prop, ymax = ucl_prop, x=cause, color=site)) + 
  geom_pointrange(position=position_dodge(width = 0.30), size=2, linewidth=2 )+
  scale_y_continuous(limits=c(0, 100), breaks=seq(0, 80, 20))+
  coord_flip()+
  labs(y = "Population Etiologic Fraction", x = " ", size=1)+
  scale_color_manual(values = cols4, breaks = c("Albania", "Jordan", "Nicaragua", "Philippines"), labels = c("Albania", "Jordan", "Nicaragua", "Philippines"), name=NULL)+
  theme_bw()+
  theme(axis.text = element_text(size = 26),
        axis.title.x = element_text(size = 28),
        legend.position = "none")+
  facet_grid(~site)+
  theme(strip.text = element_text(size=30))
dev.off()
