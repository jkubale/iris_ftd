# Co-infection severity (cleaned script)


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

ftd_h$infect_cat2 <- factor(ftd_h$infect_cat, levels = c("Single infection",
                                                         "None detected",
                                                         "Co-infection"))

ftd_c$infect_cat2 <- factor(ftd_c$infect_cat, levels = c("Single infection",
                                                         "None detected",
                                                         "Co-infection"))

# Odds of hospitalization
tot_dat <- data_merge()
  filter(rsv_result==0)%>%
  mutate(infect_cat3 = factor(infect, labels = c("None detected",
                                                  "Single infection",
                                                  "Co-infection")))
tot_dat <- tot_dat%>%
  dplyr::mutate(hosp2 = factor(hosp, levels = c(0,1), labels = c("No", "Yes")),
                sex = as.factor(p_sex))

any_hosp  <- gam(hosp ~ age_cat + sex + infect_cat2 + s(site, bs = "re")+ any_chronic,
                  data=tot_dat, family = binomial("logit"))
summary(any_hosp)
table(tot_dat$hosp)
vars <- c("sex","age_cat","infect_cat2", "any_chronic")
any_hosp_table <- tbl_regression(any_hosp, exponentiate = T,
                                  label = list(
                                    "sex" ~ "Sex",
                                    "age_cat" ~ "Age",
                                    "infect_cat2" ~ "Infection type",
                                    "any_chronic" ~ "Chronic condition"),
                                  include = all_of(vars))%>%
  modify_table_body(dplyr::select,-p.value)

any_hosp_table

# Days hospitalized
## any co-infection vs. any single infection



dh_mod_nb <- gam(days_hosp ~ sex + age_cat + infect_cat2 + any_chronic + 
                   s(site, bs="re"), data=ftd_h, family = nb())
summary(dh_mod_nb)

## any co-infection vs. any single infection (both excluding RSV)
ftd_h_norsv <- ftd_h%>%
  filter(rsv_result !=1, rsv2!=1)

dh_mod_nb2 <- gam(days_hosp ~ sex + age_cat + infect_cat2 + any_chronic + 
                    s(site, bs="re"), data=ftd_h_norsv, family = nb())
summary(dh_mod_nb2)

# Create regression tables


## comparison with RSV (regardless of pathogen combination)
dh_table <- tbl_regression(dh_mod_nb, exponentiate = T,
                           label = list(
                             "sex" ~ "Sex",
                             "age_cat" ~ "Age",
                             "infect_cat2" ~ "Infection type",
                             "any_chronic" ~ "Chronic Condition"),
                           include = all_of(vars)
)%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value

## comparison without RSV (regardless of pathogen combination)
dh_table2 <- tbl_regression(dh_mod_nb2, exponentiate = T,
                            label = list(
                              "sex" ~ "Sex",
                              "age_cat" ~ "Age",
                              "infect_cat2" ~ "Infection",
                              "any_chronic" ~ "Chronic Condition"),
                            include = all_of(vars))%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value

# save table
gt::gtsave(as_gt(dh_table), expand = 10, file = "figures/dhsco_tab.png")
gt::gtsave(as_gt(dh_table2), expand = 10, file = "figures/dhsco_tab_norsv.png")

# RSV and RV_EV specific models - take 1----
spec_path <- data_merge()
spec_path <- spec_path%>%
  dplyr::mutate(hosp = factor(hosp, levels = c(0,1), labels = c("No", "Yes")))

rsv <- spec_path%>%
  dplyr::mutate(sex = as.factor(p_sex))%>%
  filter(rsv_result==1)
rv <- spec_path%>%
  dplyr::mutate(sex = as.factor(p_sex))%>%
  filter(rv_ev==1)

## RSV single vs co-infection severity
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
gt::gtsave(as_gt(rsv_table), expand = 10, file = "figures/rsv_co_tab.png")                             

## RSV/RV_EV co-infections vs. RSV single infections
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

## RV_EV co-infections vs. RV_EV single infections
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
gt::gtsave(as_gt(rv_table), expand = 10, file = "figures/rv_co_tab.png")

rv <- rv%>%
  ungroup()%>%
  dplyr::mutate(rsv_co = as.factor(case_when(
    rsv_result==1~1,
    T~0
  )),
  hadv_co = as.factor(case_when(
    hadv==1~1,
    T~0
  )))
table(rsv$rv_co, useNA = "always")

rv_co2  <- gam(hosp ~ age_cat + sex + rsv_co + s(site, bs = "re") + any_chronic,
               data=rv, family = binomial("logit"))

co_vars3 <- c("sex","age_cat","rsv_co")
rv_table2 <- tbl_regression(rv_co2, exponentiate = T,
                            label = list(
                              "sex" ~ "Sex",
                              "age_cat" ~ "Age",
                              "rsv_co" ~ "RSV Co-Infection"),
                            include = all_of(co_vars3))%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value
gt::gtsave(as_gt(rv_table2), expand = 10, file = "figures/rv_co_tab2.png")

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
gt::gtsave(as_gt(rv_table3), expand = 10, file = "figures/rv_co_tab3.png")

# RSV and RV_EV specific models - take 2----
## need to drop any non RSV/RV_EV co-infections 

vir_ord2 <- c("hpev",  "mpneu", "hpiv1", "hpiv2",  "hpiv3", "hpiv4", 
             "cor229", "cor63", "corhku1", "cor43", "flub_univ_result",
             "flua_univ_result", "hbov", "hmpvab", "hadv")

ftd_h_sub <- ftd_h%>%
  filter(rsv_result==1|rv_ev==1)%>%
  filter(if_all(.cols=all_of(vir_ord2), ~ . ==0))%>%
  mutate(rsv_rv_cat = factor(case_when(
    rsv_result==1 & rv_ev==1 ~ "Co-infection",
    rsv_result==1 & rv_ev==0 ~ "RSV infection",
    rsv_result==0 & rv_ev==1 ~ "RV/EV infection"
  ), ordered = F, levels = c("Co-infection", "RSV infection", "RV/EV infection")))



dh_mod_nb_rsvrvev <- gam(days_hosp ~ sex + age_cat + rsv_rv_cat + any_chronic + 
                   s(site, bs="re"), data=ftd_h_sub, family = nb())
summary(dh_mod_nb_rsvrvev)


# check_overdispersion(dh_mod_nb_rsvrvev)
# check_zeroinflation(dh_mod_nb_rsvrvev)


# Create regression tables
vars2 <- c("sex","age_cat","rsv_rv_cat", "any_chronic")

## comparison with RSV (regardless of pathogen combination)
dh_table_rsvrvev <- tbl_regression(dh_mod_nb_rsvrvev, exponentiate = T,
                           label = list(
                             "sex" ~ "Sex",
                             "age_cat" ~ "Age",
                             "rsv_rv_cat" ~ "Infection",
                             "any_chronic" ~ "Chronic Condition"),
                           include = all_of(vars2)
)%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value


## odds of hospitalization
spec_path <- data_merge()
spec_path <- spec_path%>%
  dplyr::mutate(hosp = factor(hosp, levels = c(0,1), labels = c("No", "Yes")))

rsv_rvev <- spec_path%>%
  dplyr::mutate(sex = as.factor(p_sex))%>%
  filter(rsv_result==1|rv_ev==1)%>%
  filter(if_all(.cols=all_of(vir_ord2), ~ . ==0))%>%
  mutate(rsv_rv_cat = factor(case_when(
    rsv_result==1 & rv_ev==1 ~ "Co-infection",
    rsv_result==1 & rv_ev==0 ~ "RSV infection",
    rsv_result==0 & rv_ev==1 ~ "RV/EV infection"
  ), ordered = F, levels = c("Co-infection", "RSV infection", "RV/EV infection")))

rsv_rvev2  <- gam(hosp ~ age_cat + sex + rsv_rv_cat + s(site, bs = "re")+ any_chronic,
               data=rsv_rvev, family = binomial("logit"))
summary(rsv_rvev2)

rsv_rv_co_vars2 <- c("sex","age_cat","rsv_rv_cat", "any_chronic")
rsv_rvev_table2 <- tbl_regression(rsv_rvev2, exponentiate = T,
                            label = list(
                              "sex" ~ "Sex",
                              "age_cat" ~ "Age",
                              "rsv_rv_cat" ~ "Infection type",
                              "any_chronic" ~ "Chronic condition"),
                            include = all_of(rsv_rv_co_vars2))%>%
  modify_table_body(dplyr::select,-p.value)

# RV_EV and HAdV specific models ----
## need to drop any non HAdV/RV_EV co-infections 

vir_ord3 <- c("hpev",  "mpneu", "hpiv1", "hpiv2",  "hpiv3", "hpiv4", 
              "cor229", "cor63", "corhku1", "cor43", "flub_univ_result",
              "flua_univ_result", "hbov", "hmpvab", "rsv_result")

ftd_h_sub2 <- ftd_h%>%
  filter(hadv==1|rv_ev==1)%>%
  filter(if_all(.cols=all_of(vir_ord3), ~ . ==0))%>%
  mutate(hadv_rv_cat = factor(case_when(
    hadv==1 & rv_ev==1 ~ "Co-infection",
    hadv==1 & rv_ev==0 ~ "HAdV infection",
    hadv==0 & rv_ev==1 ~ "RV/EV infection"
  ), ordered = F, levels = c("Co-infection", "HAdV infection", "RV/EV infection")))



dh_mod_nb_hadvrvev <- gam(days_hosp ~ sex + age_cat + hadv_rv_cat + any_chronic + 
                           s(site, bs="re"), data=ftd_h_sub2, family = nb())
summary(dh_mod_nb_hadvrvev)


# check_overdispersion(dh_mod_nb_rsvrvev)
# check_zeroinflation(dh_mod_nb_rsvrvev)


# Create regression tables
vars3 <- c("sex","age_cat","hadv_rv_cat", "any_chronic")

## comparison with RSV (regardless of pathogen combination)
dh_table_hadrvev <- tbl_regression(dh_mod_nb_hadvrvev, exponentiate = T,
                                   label = list(
                                     "sex" ~ "Sex",
                                     "age_cat" ~ "Age",
                                     "hadv_rv_cat" ~ "Infection",
                                     "any_chronic" ~ "Chronic Condition"),
                                   include = all_of(vars3)
)%>%
  modify_table_body(dplyr::select,-p.value)# remove p-value


## odds of hospitalization
spec_path <- data_merge()
spec_path <- spec_path%>%
  dplyr::mutate(hosp = factor(hosp, levels = c(0,1), labels = c("No", "Yes")))

rsv_rvev <- spec_path%>%
  dplyr::mutate(sex = as.factor(p_sex))%>%
  filter(rsv_result==1|rv_ev==1)%>%
  filter(if_all(.cols=all_of(vir_ord2), ~ . ==0))%>%
  mutate(rsv_rv_cat = factor(case_when(
    rsv_result==1 & rv_ev==1 ~ "Co-infection",
    rsv_result==1 & rv_ev==0 ~ "RSV infection",
    rsv_result==0 & rv_ev==1 ~ "RV/EV infection"
  ), ordered = F, levels = c("Co-infection", "RSV infection", "RV/EV infection")))

rsv_rvev2  <- gam(hosp ~ age_cat + sex + rsv_rv_cat + s(site, bs = "re")+ any_chronic,
                  data=rsv_rvev, family = binomial("logit"))
summary(rsv_rvev2)

rsv_rv_co_vars2 <- c("sex","age_cat","rsv_rv_cat", "any_chronic")
rsv_rvev_table2 <- tbl_regression(rsv_rvev2, exponentiate = T,
                                  label = list(
                                    "sex" ~ "Sex",
                                    "age_cat" ~ "Age",
                                    "rsv_rv_cat" ~ "Infection type",
                                    "any_chronic" ~ "Chronic condition"),
                                  include = all_of(rsv_rv_co_vars2))%>%
  modify_table_body(dplyr::select,-p.value)
