###############################
###   ALRI ETI figures      ###
###   Created: 9/19/22      ###
###   Updated: 9/19/22      ###
###############################

library(dplyr)
library(readr)
library(here)
library(ggplot2)
library(scico)
library(patchwork)


eti_site <- read_csv("data/site_estimates.csv")
eti_age <- read_csv("data/age_estimates.csv")


eti_age <- eti_age%>%
  mutate(eti_round = round(eti_mean*100,2),
         lcl_round = round(ci_025*100,2),
         ucl_round = round(ci_975*100,2),
         cause2 = case_when(
           cause=="Metapneumovirus" ~ "HMPV",
           cause=="Adenovirus" ~ "HAdV",
           cause=="other"~ "Other",
           T~cause
         ))

top_five_age <- eti_age%>%
  group_by(age)%>%
  arrange(age, -eti_round)%>%
  slice_head(n=5)%>%
  mutate(cum_eti = cumsum(eti_round))
# Dot whisker of marginal estimates



# Dot whisker of marginal estimates by site
eti_site <- eti_site%>%
  mutate(eti_round = round(eti_mean*100,2),
         lcl_round = round(ci_025*100,2),
         ucl_round = round(ci_975*100,2))

top_five_site <- eti_site%>%
  group_by(site)%>%
  arrange(site, -eti_round)%>%
  slice_head(n=5)%>%
  mutate(cum_eti = cumsum(eti_round),
         cause2 = case_when(
           cause=="Metapneumovirus" ~ "HMPV",
           cause=="Adenovirus" ~ "HAdV",
           cause=="other"~ "Other",
          T~cause
         ))
# change Metapneumovirus, other, Adenovirus
  # scale_color_scico_d(palette = "vik", name = "", end=.8)
# scico(8, end=.8, palette = "vik")
# class(eti_site$cause)
# 
# top_ten_eti_panel <- top_ten_site

eti_site_panel <- function(dat,loc, color){
  dat%>%
    filter(site==loc)%>%
    ggplot(aes(x = site, y = eti_round, ymin = lcl_round, ymax = ucl_round, color=as.factor(cause2)))+
    geom_pointrange(position=position_dodge(width = .8), size=2 , linewidth = 2)+
    # scale_color_scico_d(palette = "vik", name = "", end=.8)+
    scale_color_manual(values = color, name=NULL)+
    labs(x = "", y="Population Etiological Fraction")+
    scale_y_continuous(limits = c(0,80),breaks = seq(from=0, to=80, by=20))+
    theme_bw()+
    theme(axis.text.y = element_text(size = 28),
          axis.title.y = element_text(size = 30),
          legend.title = element_text(size = 28),
          legend.text = element_text(size = 28),
          plot.title = element_text(size = 30),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())+
    ggtitle(loc)
}

cols_alb <- c("Influenza A" = "#77ACC6", "Influenza B" = "#D3E3EA", 
              "HMPV" = "#1F719E", "RSV" = "#C2A56A", "RV_EV" = "#9F711B")
cols_jord <- c("Influenza B" = "#D3E3EA",  "HMPV" = "#1F719E", "Other" = "#E4D7BD", 
               "RSV" = "#C2A56A", "RV_EV" = "#9F711B")
cols_nic <- c("HAdV" = "#001260", "HPIV3" = "#023F7E",  "HMPV" = "#1F719E",
              "RSV" = "#C2A56A", "RV_EV" = "#9F711B")
cols_phil <- c("HPIV3" = "#023F7E", "HMPV" = "#1F719E", "Other" = "#E4D7BD", 
               "RSV" = "#C2A56A", "RV_EV" = "#9F711B")

## need to fix the colors so that they aren't reused for different pathogens
## still need to fix hmpv and hpiv3 order in Nicaragua panel?
alb_top5 <- eti_site_panel(top_five_site, "Albania", cols_alb)
jord_top5 <- eti_site_panel(top_five_site, "Jordan", cols_jord)
nic_top5 <- eti_site_panel(top_five_site, "Nicaragua", cols_nic)
phil_top5 <- eti_site_panel(top_five_site, "Philippines", cols_phil)

# change Metapneumovirus, other, Adenovirus
top5_panel <- (alb_top5|jord_top5)/(nic_top5|phil_top5)

png(file="figures/top5_panel022423.png", width=6000,height=6000,res = 300)
top5_panel
dev.off()




scico(3, end=0.8, palette = "vik")
# cols_age <- c("Under 3 months" = "#001260", "3-5 months" = "#A5C9D9", "6-11 months"= "#9F711B")

## RSV by age
eti_age_panel <- function(dat, age_cat, color){
dat%>%
  filter(age==age_cat)%>%
ggplot(aes(x = age, y = eti_round, ymin = lcl_round, ymax = ucl_round, color=as.factor(cause2)))+
geom_pointrange(position=position_dodge(width = .8), size=2 , linewidth=2)+
  # scale_color_scico_d(palette = "vik", name = "", end=.8)+
  scale_color_manual(values = color, name=NULL)+
  labs(x = "", y="Population Etiological Fraction")+
  scale_y_continuous(limits = c(0,80),breaks = seq(from=0, to=80, by=20))+
  theme_bw()+
  theme(axis.text.y = element_text(size = 28),
        axis.title.y = element_text(size = 30),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 30),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
    ggtitle(age_cat)
}

cols_un3m <- c("HPIV3" = "#023F7E", "HMPV" = "#1F719E", "Other" = "#E4D7BD", 
               "RSV" = "#C2A56A", "RV_EV" = "#9F711B")
cols_35m <- c("Influenza A" = "#77ACC6", "HPIV3" = "#023F7E", 
              "HMPV" = "#1F719E", "RSV" = "#C2A56A", "RV_EV" = "#9F711B")
cols_611m <- c("Influenza A" = "#77ACC6", "Influenza B" = "#D3E3EA", 
               "HMPV" = "#1F719E", "RSV" = "#C2A56A", "RV_EV" = "#9F711B")

eti_un3 <- eti_age_panel(top_five_age, "Under 3 months", cols_un3m)
eti_35 <- eti_age_panel(top_five_age, "3-5 months", cols_35m)
eti_611 <- eti_age_panel(top_five_age, "6-11 months", cols_611m)

top5_age_panel <- (eti_un3|eti_35|eti_611)


png(file="figures/top5_age_panel022423.png", width=8000,height=4000,res = 300)
top5_age_panel
dev.off()
