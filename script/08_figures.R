##############################
###   IRIS data figures    ###
###   Created: 12/18/22    ###
###   Updated: 12/18/22    ###
##############################

library(dplyr)
library(ggplot2)
library(ggupset)
library(lubridate)
library(tidyr)
library(RColorBrewer)

# Supplemental figure 1:----
# table with top 5 viruses and cumulative prop as plot by site

# Supplemental figure 2:----
# Upset plots of co-infection combinations --------------------------------
vir3 <- c("flua_univ_result","flub_univ_result","rsv_result","cor63","cor229",
          "cor43","corhku1","hpiv2","hpiv3","hpiv4","hpiv1","hmpvab","hbov", "mpneu",
          "hadv","rv_ev","hpev")

ftd_h_upset <- ftd_h%>%select(study_id, all_of(vir3), infect)
ftd_c_upset <- ftd_c%>%select(study_id, all_of(vir3), infect)
ftd_tot <- bind_rows(ftd_h_upset, ftd_c_upset)

# ftd_tot <- ftd_tot%>%
#   ungroup()%>%
#   mutate(rv2 = case_when(
#     rv==0 ~0,
#     rv==1 & ev==0~1,
#     rv==1 & ev==1~0,
#     T~99
#   ))
table(ftd_tot$rv,useNA = "always")

## need to update gather to pivot
ftd_tot2 <- ftd_tot %>%
  as_tibble() %>%
  filter(infect >1)%>%
  select(study_id, all_of(vir3))%>%
  rename(Flu_A = flua_univ_result,
         Flu_B = flub_univ_result)%>%
  tidyr::gather(virus, Member, -study_id) %>%
  filter(Member==1) %>%
  select(- Member)%>%
  group_by(study_id)%>%
  summarise(virus = list(virus))

cols2 <- c("0" = "#78C679","1" ="#006837")

## tweak labels
tot <- ftd_tot2%>%
  # distinct(title, year, length, .keep_all=TRUE) %>%
  ggplot(aes(x=virus)) +
  geom_bar(fill="#006837") +
  # scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2))+
  ylab("Count")+
  xlab(" ")+
  scale_x_upset(n_intersections = 10)+
  theme_classic()+
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  theme_combmatrix(combmatrix.panel.point.color.fill = "#006837",
                   combmatrix.panel.line.size = 0,
                   combmatrix.label.text = element_text(size=16))

tot

## Hospitalized

## Non-ill

# Supplemental figure 3:----
# seasonality
ftd_tot_sub <- subset_merge()

# get number of samples per week to standardize across pathogens
min(ftd_tot_sub$e_specimen_resp_date) # 6/26/2015 -- will use 6/21/2015 as week start since it was a Sunday

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
  # mutate(across(.cols = all_of(vir_ord), ~.x/tot_samp))%>%
  select(-tot_samp)%>%
  melt(id.vars = c("country_id", "study_m"))


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

summary(ftd_tot_long2$value)
chk <- filter(ftd_tot_long2, value>0)
hist(chk$value)


scales::viridis_pal(option = "C")(6)
cols_heat <- c( "0" = "light gray","1-4" = "#0D0887FF","5-9" = "#6A00A8FF","10-14" = "#B12A90FF",
                "15-19" = "#E16462FF","20-24" = "#FCA636FF","25+" = "#F0F921FF")


seas_heat <- ggplot(ftd_tot_long2,aes(date,variable2,fill=value_cat))+
  geom_tile(color= "white",size=0.1) + 
  # scale_fill_gradient()+
  # scale_fill_manual(na.value = "gray",breaks = c(0.03, 0.1, 0.2, 0.3))+
  # scale_fill_viridis(na.value = "gray",discrete=T, name = NULL, na.translate=F, option = "C")+
  scale_fill_manual(na.value = "white",values = cols_heat, name = "Num positives \n per month")+
  scale_x_date(date_labels="%b-%y", date_breaks="2 months", expand=c(0,0)) + # sets format for month-year ("%b-%y") and displays only every 2 months
  theme_bw()+
  theme(axis.text.y = element_text(size = 16),
        axis.title = element_blank(),
        # axis.text.x = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle=90, vjust = 0.5),
        strip.text.x = element_text(size = 16),
        legend.text=element_text(size=16),
        legend.title = element_text(size=18)
  )+
  guides(guide_colorsteps(frame.colour = "black"))+
  facet_grid(~country_id)

png(file="figures/seas_panel82322.png", width=8200,height=4500,res = 300)
seas_heat
dev.off()