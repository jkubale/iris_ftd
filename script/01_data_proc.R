###############################
###   Processing IRIS data  ###
###   Created: 11/18/20     ###
###   Updated: 11/25/20     ###
###############################

library(dplyr)
library(readr)
library(haven)
library(stringr)
library(lubridate)
library(gmodels)
library(reshape2)
library(ggplot2)
library(tidyr)

source("script/ftd_functions.R")

# Quick ref ---------------------------------------------------------------
# country_id:
# 1=Nicaragua
# 2=Philippines
# 3=Jordan
# 4=Albania

# Pre-term cutoffs:
# Late preterm 34-36 weeks
# Moderately preterm 32-34 weeks
# Very preterm <32 weeks

# Low birthweight cutoffs:
# Low birthweight <2500 g
# Very low birthweight <1500 g
# Extremely low birthweight <1000 g

# Notes -------------------------------------------------------------------
#1)
## will use univ_result for influenza a and b as there is currently a discrepancy
## use rsv_result vs. rsv2
## may check raw data for inconsistency
## confirm before final

#2)
## how to (or if to) use obs with no IC or with ct >33?

# Create vectors of variable names ----------------------------------------

# create vector of virus result variables
vir <- c("flua_univ_result","flub_univ_result","rsv_result","flua",
         "flub","flua_h1n1swl","rv","cor63","cor229","cor43","corhku1",
         "hpiv2","hpiv3","hpiv4","ic","ic_ct","hpiv1","hmpvab","hbov",
         "mpneu","rsv2","hadv","ev","EVD68","hpev")

hosp_vars <- c("study_id", "country_id", "study_year", "p_dob", "p_calc_age", "hosp_infage_cat", "hosp_infage_dic",
               "p_sex", "s_onset", "p_date", "p_admitdate", "death_indicator", "resp_crf", "temp_crf", "m_oxygen",
               "m_oxygen_mode___1", "m_oxygen_mode___2", "m_oxygen_mode___3", "m_oxygen_mode___4", "m_oxygen_mode___5",
               "m_oxygen_mode___6", "m_oxygen_mode___7", "m_oxygen_mode___8", "m_oxygen_mode___9",
               "o2sat_crf", "location_crf", "o2rec_crf", "oxygenmode1_crf", "oxygenmode2_crf", "oxygenmode3_crf",
               "oxygenmode4_crf", "oxygenmode5_crf", "oxygenmode6_crf", "oxygenmode7_crf", "oxygenmode8_crf",
               "oxygenmode9_crf", "m_dx_r_pneu___1","m_dx_r_pneu___3", "e_chronic", "e_priormeds", "e_priormeds_abx", "m_clx_primary", "d_clx_primary",
               all_of(vir), "e_specimen_resp_type", "e_specimen_resp_id", "e_specimen_resp_date", "e_mat_delivery", 
               "e_mat_preglength_wks", "e_mat_prem", "e_mat_icu", "e_mat_vent", "e_mat_birthweight", "e_bf_ever", "e_bf_exclusive",
               "e_bf_current", "e_bf_agestop", "e_hh_care", "e_smoke", "e_rank", "e_rank_edu", "e_rank_job", "e_rank_money", "d_date",
               "e_mother_edu", "e_father_edu", "f_birthweight", "e_mat_birthweight", 
               "e_mat_birthweight_verify")

non_ill_vars <- c("study_id", "country_id", "study_year", "p_dob", "p_calc_age", "hosp_infage_cat", "hosp_infage_dic",
                  "p_sex", "s_onset", "p_date", "p_admitdate", "death_indicator", all_of(vir), "e_chronic",
                  "e_specimen_resp_type", "e_specimen_resp_id", "e_specimen_resp_date", "e_smoke",
                  "e_mat_delivery", "e_mat_preglength_wks", "e_mat_prem", "e_mat_icu",
                  "e_mat_vent", "e_mat_birthweight", "e_rank", "e_rank_edu", "e_rank_job", "e_rank_money",
                  "e_mother_edu", "e_father_edu", "f_birthweight", "e_mat_birthweight", 
                  "e_mat_birthweight_verify")

labs <- c("flua_univ_result","flub_univ_result","rsv_result","flua",
          "flub","rv","cor63","cor229","cor43","corhku1",
          "hpiv2","hpiv3","hpiv4","ic","hpiv1","hmpvab","hbov",
          "mpneu","rsv2","hadv","ev","hpev")

# same as labs but doesn't include ic
vir2<- c("flua_univ_result","flub_univ_result","rsv_result","rv_ev","cor63","cor229",
         "cor43","corhku1","hpiv2","hpiv3","hpiv4","hpiv1","hmpvab","hbov", "mpneu",
         "hadv","hpev")## add EVD68? -- 3 pos that are ev neg, continue using rsv_result?


# "e_hh_electricity", "e_hh_toilet", "e_hh_floor"
# Create vectors for various variable combinations (lab results, diagnoses, severity)

# Read in data ------------------------------------------------------------
# Way too many variables in overall dataset
## start with basics and build out from there
## rename vars so don't have to reference data dictionary
## harmonize format of vars between hosp and controls

iris_basic_h <- read_dta("data/IRIS_hosp_clean_24Jun2020.dta")%>%
  zap_label()%>%
  dplyr::select(all_of(hosp_vars))

iris_basic_c <- read_dta("data/IRIS_FTD_nonill_24Jun2020.dta")%>%
  zap_label()%>%
  dplyr::select(all_of(non_ill_vars), starts_with("cf"))

# harmonize variable formats ----------------------------------------------
## Need harmonizing: country_id, study_year, e_specimen_resp_type, p_sex, all labs
## use hosp format: study_year, all labs, e_specimen_resp_type, e_mat_delivery, e_mat_prem, e_mat_icu, e_mat_vent, e_mat_birthweight
## use control format: country_id, p_sex 

con_class <- data.frame(lapply(iris_basic_c, class))
hosp_class <- data.frame(lapply(iris_basic_h, class))

iris_basic_h <- harm_ftd_h(iris_basic_h)

iris_basic_c2 <- iris_basic_c
iris_basic_c2[,which(colnames(iris_basic_c) %in% labs)] <- apply(iris_basic_c2[,which(colnames(iris_basic_c) %in% labs)],2, function(x) {ifelse(x == "Negative", 0,
                                                                                                                                                ifelse(x == "Positive",1,NA))})

ftd_c <- harm_ftd_c(iris_basic_c2)

# Data transformation and basic analysis ----------------------------------
# Create variable for co-infection vs. single infection
ftd_c$del <- as.numeric(apply(ftd_c[,which(colnames(ftd_c) %in% labs)], 1, function(x) {ifelse(any(is.na(x)), 1, 0)}))
# table(ftd_c$ic, useNA = "always")
# chk <- filter(ftd_c, del==1) ## one sample has missing value for hpiv4 only
# will recode as 0 for now to distinguish between single and co-infections
# chk2 <- filter(ftd_c, ic_ct >33|ic==0|is.na(ic))
# table(ftd_c$hpiv4, useNA = "always")
ftd_c$hpiv4 <- ifelse(ftd_c$study_id=="JC-0008", 0, ftd_c$hpiv4)

ftd_c <- ftd_c%>%
  ungroup()%>%
  mutate(rv_ev = case_when(
    rv==0 & ev==0 ~0,
    rv==1 | ev==1~1))

# get number of pathogens detected per participant
# community controls
ftd_c$num_in <- as.numeric(apply(ftd_c[,which(colnames(ftd_c) %in% vir2)], 1, sum))
# table(ftd_c$num_in, useNA = "always")

ftd_c <- ftd_c%>%
  dplyr::ungroup()%>%
  dplyr::mutate(infect = case_when(
    num_in == 0 ~0,
    num_in == 1 ~1,
    num_in >1 ~2,
    T~99
  ),
  infect_d = case_when(
    infect>1 ~1,
    T~0
  ))
  

# table(ftd_c$num_in) 
#  0   1   2    3 
# 553 396 111   8 


# hospitalized
ftd_h <- filter(iris_basic_h, is.na(rv)==F)
ftd_h$del <- as.numeric(apply(ftd_h[,which(colnames(ftd_h) %in% labs)], 1, function(x) {ifelse(any(is.na(x)), 1, 0)}))
# table(ftd_h$del, useNA = "always")
# chk2 <- filter(ftd_h, del==1) ## IC is missing from all 358
# table(chk2$country_id, useNA = "always")

# Will use all (even without IC) for now
# get number of pathogens per participant
ftd_h <- ftd_h%>%
  dplyr::ungroup()%>%
  dplyr::mutate(rv2 = case_when(
    rv==0 ~0,
    rv==1 & ev==0~1,
    rv==1 & ev==1~0
  ),
  rv_ev = case_when(
    rv==0 & ev==0 ~0,
    rv==1 | ev==1~1))

# table(ftd_h$rv2, useNA = "always")
ftd_h$num_in <- as.numeric(apply(ftd_h[,which(colnames(ftd_h) %in% vir2)], 1, sum))
# table(ftd_h$num_in, useNA = "always")
#     0    1    2    3    4 
#   1080 1830  599  112   11    

# table(ftd_h$num_in,ftd_h$rsv_result, useNA = "always")

ftd_h <- ftd_h%>%
  dplyr::ungroup()%>%
  dplyr::mutate(infect = case_when(
    num_in == 0 ~0,
    num_in == 1 ~1,
    num_in >1 ~2,
    T~99
  ))

ftd_h$infect_d <- ifelse(ftd_h$infect >1, 1,0)

# colnames(ftd_h)
# table(ftd_h$infect_d, useNA = "always")
# table(ftd_c$infect_d, useNA = "always")

# co_inf <- matrix(c(722, 119, 2190, 949), ncol=2,byrow = T)
# chisq.test(co_inf) # p < 0.001


ftd_h$day_dif <- as.numeric(interval(start = ftd_h$s_onset, end = ftd_h$e_specimen_resp_date)%/% days(1))
ftd_h$days_hosp <- as.numeric(interval(start = ftd_h$p_admitdate, end = ftd_h$d_date)%/% days(1))

# table(ftd_h$infect, useNA = "always")
# table(ftd_h$days_hosp, useNA = "always")
ftd_h <- filter(ftd_h, days_hosp >=0) #1 obs from philippines with earlier discharge date, dropping for now
# table(ftd_h$hadv, useNA = "always")

# table(ftd_h$num_in, useNA = "always") #29.7% had 0 infections, 70.3% had at least 1
# table(ftd_c$num_in, useNA = "always") #51.8% had 0 infections, 48.2% had at least 1
# table(ftd_h$rsv_result,ftd_h$num_in, useNA = "always")
# table(ftd_c$rsv_result,ftd_c$num_in, useNA = "always")

ftd_h$any_chronic <- ifelse(ftd_h$e_chronic==0,0,
                            ifelse(ftd_h$e_chronic==1,1,NA))
ftd_c$any_chronic <- ifelse(ftd_c$e_chronic=="No",0,
                            ifelse(ftd_c$e_chronic=="Yes",1,NA))

ftd_h <- ftd_h%>%
  mutate(age_cat = case_when(
  hosp_infage_cat == 0 ~ "<14 weeks",
  hosp_infage_cat == 1 ~ "14-26 weeks",
  hosp_infage_cat == 2 ~ "27-39 weeks",
  hosp_infage_cat == 3 ~ "40-52 weeks",
))

ftd_c <- ftd_c%>%
  mutate(age_cat = case_when(
    hosp_infage_cat == 0 ~ "<14 weeks",
    hosp_infage_cat == 1 ~ "14-26 weeks",
    hosp_infage_cat == 2 ~ "27-39 weeks",
    hosp_infage_cat == 3 ~ "40-52 weeks",
  ))

# table(ftd_c$any_chronic, useNA = "always")

save(ftd_h, file = "data/ftd_h122222.rda")
save(ftd_c, file = "data/ftd_c020722.rda")

load("data/ftd_h122222.rda")
load("data/ftd_c020722.rda")

############################  old  ################################################ 
# table(ftd_c$infect, useNA = "always")
# ftd_c$infect_d <- ifelse(ftd_c$infect >1, 1,0)
# table(ftd_c$infect_d, useNA = "always")


#table(iris_basic_c2$rv2, useNA = "always")## all were missing rv results


colnames(iris_basic_c2)
ftd_c <- filter(iris_basic_c2, is.na(rv)==F)
