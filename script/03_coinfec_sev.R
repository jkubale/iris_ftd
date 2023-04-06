############################################
###   Co-infection and severity models   ###
###           Created: 04/19/21          ###
###           Updated: 07/11/21          ###
############################################


library(dplyr)
library(readr)
library(haven)
library(stringr)
library(lubridate)
library(gmodels)
library(reshape2)
library(ggplot2)
library(mgcv)
library(gtsummary)
library(performance)


load("data/ftd_h120622.rda")
load("data/ftd_c120622.rda")
source("script/ftd_functions.R")

# Move code creating these vars to 01_data_proc script
ftd_h$site <- as.factor(ftd_h$country_id)
ftd_c$site <- as.factor(ftd_c$country_id)
ftd_h$sex <- as.factor(ftd_h$p_sex)
ftd_c$sex <- as.factor(ftd_c$p_sex)
ftd_h$hosp <- 1
ftd_c$hosp <- 0

ftd_h$age_cat <- factor(ftd_h$hosp_infage_cat, labels = c("<14 weeks", "14-26 weeks",
                                                          "27-39 weeks", "40+ weeks"))
ftd_h$infect_cat <- factor(ftd_h$infect, labels = c("None detected",
                                                    "Single infection",
                                                    "Co-infection"))
ftd_c$age_cat <- factor(ftd_c$hosp_infage_cat, labels = c("<14 weeks", "14-26 weeks",
                                                          "27-39 weeks", "40+ weeks"))
ftd_c$infect_cat <- factor(ftd_c$infect, labels = c("None detected",
                                                    "Single infection",
                                                    "Co-infection"))

ftd_h$any_chronic <- factor(ftd_h$any_chronic, labels = c("No chronic condition", 
                                                          "Chronic condition"))

ftd_c$any_chronic <- ifelse(ftd_c$e_chronic=="No",0,
                            ifelse(ftd_c$e_chronic=="Yes",1,NA))

ftd_h$any_chronic <- factor(ftd_h$any_chronic, labels = c("No chronic condition", 
                                                          "Chronic condition"))

ftd_c$any_chronic <- factor(ftd_c$any_chronic, labels = c("No chronic condition", 
                                                          "Chronic condition"))

table(ftd_c$any_chronic)
table(ftd_h$any_chronic)
# Days hospitalized
# dh_mod_p <- gam(days_hosp ~ sex + age_cat + infect_cat + any_chronic + 
#                 s(site, bs="re"), data=ftd_h, family = poisson())
# summary(dh_mod_p)
# 
# # model checks
# check_overdispersion(dh_mod_p)
# check_zeroinflation(dh_mod_p)
ftd_h$infect_cat2 <- factor(ftd_h$infect_cat, levels = c("Single infection",
                                          "None detected",
                                          "Co-infection"))

dh_mod_nb <- gam(days_hosp ~ sex + age_cat + infect_cat2 + any_chronic + 
                   s(site, bs="re"), data=ftd_h, family = nb())
summary(dh_mod_nb)

ftd_h_norsv <- ftd_h%>%
  filter(rsv_result !=1, rsv2!=1)

dh_mod_nb2 <- gam(days_hosp ~ sex + age_cat + infect_cat2 + any_chronic + 
                    s(site, bs="re"), data=ftd_h_norsv, family = nb())
summary(dh_mod_nb2)
exp(0.07294)
# model checks
check_zeroinflation(dh_mod_nb) # how to handle overfitting of zeros?
#look at residuals and see if
boxplot(resid(dh_mod_nb, type = "pearson"))

nb_pred <- predict.gam(dh_mod_nb, exclude = "site")
plot(nb_pred, resid(dh_mod_nb, type = "pearson"))

# Create regression table
vars <- c("sex","age_cat","infect_cat2", "any_chronic")

dh_table <- tbl_regression(dh_mod_nb, exponentiate = T,
                           label = list(
                             "sex" ~ "Sex",
                             "age_cat" ~ "Age",
                             "infect_cat2" ~ "Infection",
                             "any_chronic" ~ "Chronic Condition"),
                           include = all_of(vars)
)%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value

dh_table2 <- tbl_regression(dh_mod_nb2, exponentiate = T,
                            label = list(
                              "sex" ~ "Sex",
                              "age_cat" ~ "Age",
                              "infect_cat2" ~ "Infection",
                              "any_chronic" ~ "Chronic Condition"),
                            include = all_of(vars))%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value

# save table
gt::gtsave(as_gt(dh_table), expand = 10, file = "Box/IRIS_/figures/dhsco_tab.png")
gt::gtsave(as_gt(dh_table2), expand = 10, file = "Box/IRIS_/figures/dhsco_tab_norsv.png")

dh_table

exp(-0.14273)

# Resp vs. non-resp
dhrnr_mod_p <- gam(days_hosp ~ sex + age_cat + as.factor(m_clx_primary) + any_chronic + 
                     s(site, bs="re"), data=ftd_h, family = poisson())
summary(dhrnr_mod_p)
check_overdispersion(dhrnr_mod_p)
check_zeroinflation(dhrnr_mod_p)

dhrnr_mod_nb <- gam(days_hosp ~ sex + age_cat + as.factor(m_clx_primary) + any_chronic + 
                      s(site, bs="re"), data=ftd_h, family = nb())
summary(dhrnr_mod_nb) ## participants with non-resp illness ~13% shorter average days hosp
## ^need to check residuals

# RSV and RV models maybe HADV and pooled coronaviruses----
spec_path <- data_merge()
spec_path <- spec_path%>%
  dplyr::mutate(hosp = factor(hosp, levels = c(0,1), labels = c("No", "Yes")))

rsv <- spec_path%>%
  dplyr::mutate(sex = as.factor(p_sex))%>%
  filter(rsv_result==1)
rv <- spec_path%>%
  dplyr::mutate(sex = as.factor(p_sex))%>%
  filter(rv_ev==1)
# hadv <- spec_path%>%
#   filter(hadv==1)
# table(hadv$site)

# cor_tot <- spec_path%>%
#   mutate(cor_any = case_when(
#     cor43==0 & cor63==0 & cor229==0 & corhku1==0 ~ 0,
#     cor43==1 | cor63==1 | cor229==1 | corhku1==1 ~ 1
#   ))%>%
#   filter(cor_any==1)
# table(cor_tot$site)

rsv_co1 <- gam(hosp ~ age_cat + sex + infect_d + s(site, bs = "re") + any_chronic,
               data=rsv, family = binomial("logit"))
summary(rsv_co1)

co_vars1 <- c("sex","age_cat","infect_d", "any_chronic")
rsv_table <- tbl_regression(rsv_co1, exponentiate = T,
                            label = list(
                              "sex" ~ "Sex",
                              "age_cat" ~ "Age",
                              "infect_d" ~ "Co-infection",
                              "any_chronic" ~ "Chronic condition"),
                            include = all_of(co_vars1))%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value
gt::gtsave(as_gt(rsv_table), expand = 10, file = "Box/IRIS_/figures/rsv_co_tab.png")                             

rsv <- rsv%>%
  ungroup()%>%
  dplyr::mutate(rv_co = as.factor(case_when(
    rv_ev==1~1,
    T~0
  )))
table(rsv$rv_co, useNA = "always")

rsv_co2  <- gam(hosp ~ age_cat + sex + rv_co + s(site, bs = "re") + any_chronic,
                data=rsv, family = binomial("logit"))
summary(rsv_co2)

co_vars2 <- c("sex","age_cat","rv_co", "any_chronic")
rsv_table2 <- tbl_regression(rsv_co2, exponentiate = T,
                            label = list(
                              "sex" ~ "Sex",
                              "age_cat" ~ "Age",
                              "rv_co" ~ "RV Co-infection",
                            "any_chronic" ~ "Chronic condition"),
                            include = all_of(co_vars2))%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value
gt::gtsave(as_gt(rsv_table2), expand = 10, file = "Box/IRIS_/figures/rsv_co_tab2.png")


rv_co1 <- gam(hosp ~ age_cat + sex + infect_d + s(site, bs = "re")+ any_chronic,
              data=rv, family = binomial("logit"))
summary(rv_co1)


rv_table <- tbl_regression(rv_co1, exponentiate = T,
                           label = list(
                             "sex" ~ "Sex",
                             "age_cat" ~ "Age",
                             "infect_d" ~ "Co-infection",
                             "any_chronic" ~ "Chronic condition"),
                           include = all_of(co_vars1))%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value
gt::gtsave(as_gt(rv_table), expand = 10, file = "Box/IRIS_/figures/rv_co_tab.png")

rv <- rv%>%
  ungroup()%>%
  dplyr::mutate(rsv_co = as.factor(case_when(
    rsv_result==1~1,
    T~0
  )))
table(rsv$rv_co, useNA = "always")

rv_co2  <- gam(hosp ~ age_cat + sex + rsv_co + s(site, bs = "re"),
                data=rv, family = binomial("logit"))
summary(rv_co2)

co_vars3 <- c("sex","age_cat","rsv_co")
rv_table2 <- tbl_regression(rv_co2, exponentiate = T,
                           label = list(
                             "sex" ~ "Sex",
                             "age_cat" ~ "Age",
                             "rsv_co" ~ "RSV Co-Infection"),
                           include = all_of(co_vars3))%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value
gt::gtsave(as_gt(rv_table2), expand = 10, file = "Box/IRIS_/figures/rv_co_tab2.png")

### Drop RV_EV/RSV co-infections
rv_sub <- rv%>%
  filter(rsv_co==0)

rv_co3  <- gam(hosp ~ age_cat + sex + infect_d + s(site, bs = "re")+ any_chronic,
               data=rv_sub, family = binomial("logit"))
summary(rv_co3)


rv_table3 <- tbl_regression(rv_co3, exponentiate = T,
                           label = list(
                             "sex" ~ "Sex",
                             "age_cat" ~ "Age",
                             "infect_d" ~ "Co-infection",
                             "any_chronic" ~ "Chronic condition"),
                           include = all_of(co_vars1))%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value
gt::gtsave(as_gt(rv_table3), expand = 10, file = "Box/IRIS_/figures/rv_co_tab2.png")

chk <- filter(ftd_h, rv_ev==1, infect_d==1)
table(chk$rv_ev, chk$rsv_result, useNA = "always") ##205
table(chk$rv_ev, chk$hadv, useNA = "always") ##103
(103+205)/482

# limited power
hadv_co1 <- gam(infect_d ~ age_cat + sex + hosp ,
                data=hadv, family = binomial("logit"))
summary(hadv_co1)

cor_co1 <- gam(infect_d ~ age_cat + sex + hosp ,
               data=cor_tot, family = binomial("logit"))
summary(cor_co1)
exp(1.03355)
### Unexpected results -- among these models
# ICU admission
table(ftd_h$location_crf, useNA = "always", ftd_h$country_id)

icu_mod1 <- gam(location_crf ~ sex + age_cat + infect_cat + any_chronic + 
                  s(site, bs="re"), data=ftd_h, family = binomial())
summary(icu_mod1)

# Death 
ftd_h$death <- factor(ftd_h$death_indicator, ordered = F)
class(ftd_h$infect_d)
dth_mod1 <- gam(death ~ as.factor(infect_cat) + sex + as.factor(hosp_infage_dic) + s(site, bs ="re"),
                data=ftd_h, family = binomial())
summary(dth_mod1)
