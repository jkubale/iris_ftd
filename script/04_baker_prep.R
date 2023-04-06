############################################
###      Prep data for NPLCM models      ###
###           Created: 07/11/21          ###
###           Updated: 09/01/21          ###
############################################

library(dplyr)
library(haven)
# library(R2WinBUGS)
# library(baker)
library(reshape2)
library(lubridate)
# library(rjags)
library(ggplot2)

source("script/ftd_functions.R")

# Create vector of variables to import
vir <- c("flua_univ_result","flub_univ_result","rsv_result","flua",
         "flub","flua_h1n1swl","rv","cor63","cor229","cor43","corhku1",
         "hpiv2","hpiv3","hpiv4","hpiv1","hmpvab","hbov",
         "mpneu","rsv2","hadv", "ev","EVD68","hpev")

vir2 <- c("flua_univ_result","flub_univ_result","rsv_result",
          "cor63","cor229","cor43","corhku1",
          "hpiv2","hpiv3","hpiv4","hpiv1","hmpvab","hbov",
          "mpneu","hadv","rv_ev", "hpev")

vir3 <- c("flua_univ_result","flub_univ_result","rsv_result",
          "cor63","cor229","cor43","corhku1",
          "hpiv2","hpiv3","hpiv4","hpiv1","hmpvab","hbov",
          "mpneu","hadv","rv_ev","hpev")


labs <- c("flua_univ_result","flub_univ_result","rsv_result","flua",
          "flub","rv","cor63","cor229","cor43","corhku1",
          "hpiv2","hpiv3","hpiv4", "hpiv1","hmpvab","hbov",
          "mpneu","rsv2","hadv","ev", "hpev")

labs2 <- c("flua_univ_result","flub_univ_result","rsv_result",
          "rv","cor63","cor229","cor43","corhku1",
          "hpiv2","hpiv3","hpiv4", "hpiv1","hmpvab","hbov",
          "mpneu","hadv","ev", "hpev")

hosp_vars <- c("study_id", "country_id", "study_year", "p_dob", "hosp_infage_cat", 
               "hosp_infage_dic", "p_sex", "s_onset", "p_date", "p_admitdate", 
               "death_indicator", "m_dx_r_pneu___1","m_dx_r_pneu___3", "m_dx_r_ards___1","m_dx_r_ards___3",
               "m_dx_r_bronchio___1","m_dx_r_bronchio___3", "m_dx_r_respfail___1","m_dx_r_respfail___3", 
               "e_chronic", "e_priormeds", "e_priormeds_abx", "m_clx_primary", "d_clx_primary",
               all_of(vir), "e_specimen_resp_id", "e_specimen_resp_date", "e_mat_delivery", 
               "e_mat_preglength_wks", "e_mat_prem", "e_mat_icu", "e_mat_vent", 
               "e_mat_birthweight", "e_bf_ever", "e_bf_exclusive", "e_bf_current", "e_bf_agestop",
               "e_hh_care", "e_smoke", "d_date", "oxygenmode4_crf", "oxygenmode5_crf", 
               "oxygenmode8_crf", "oxygenmode9_crf", "m_oxygen_mode___4", "m_oxygen_mode___5",
               "m_oxygen_mode___8", "m_oxygen_mode___9")

non_ill_vars <- c("study_id", "country_id", "study_year", "p_dob", "hosp_infage_cat", "hosp_infage_dic",
                  "p_sex", "s_onset", "p_date", "p_admitdate", "death_indicator", all_of(vir), "e_chronic",
                  "e_specimen_resp_type", "e_specimen_resp_id", "e_specimen_resp_date", "e_smoke",
                  "e_mat_delivery", "e_mat_preglength_wks", "e_mat_prem", "e_mat_icu",
                  "e_mat_vent", "e_mat_birthweight", "e_mother_edu", "e_father_edu")

cause_list_iris <- c("Influenza A", "Influenza B", "RSV", "Cor63",
                     "Cor229", "Cor43", "CorHKU1", "HPIV2", "HPIV3", "HPIV4", "HPIV1",
                     "Metapneumovirus", "Bocavirus", "M. pneumoniae", "Adenovirus", "RV_EV", "Parechovirus")

iris_hlca <- read_dta("data/IRIS_hosp_clean_24Jun2020.dta")%>%
  zap_label()%>%
  dplyr::select(all_of(hosp_vars))

iris_clca <- read_dta("data/IRIS_FTD_nonill_24Jun2020.dta")%>%
  zap_label()%>%
  dplyr::select(all_of(non_ill_vars))

con_class <- data.frame(lapply(iris_clca, class))
hosp_class <- data.frame(lapply(iris_hlca, class))

iris_hlca <- harm_ftd_h(iris_hlca)

iris_clca2 <- iris_clca
iris_clca2[,which(colnames(iris_clca) %in% labs)] <- apply(iris_clca2[,which(colnames(iris_clca2) %in% labs)],2, function(x) {ifelse(x == "Negative", 0,
                                                                                                                                     ifelse(x == "Positive",1,NA))})
iris_clca2 <- harm_ftd_c(iris_clca2)

# table(iris_clca2$rv2, useNA = "always")## all were missing rv results
iris_clca3 <- filter(iris_clca2, is.na(rv)==F)
iris_clca3$del <- as.numeric(apply(iris_clca3[,which(colnames(iris_clca3) %in% labs)], 1, function(x) {ifelse(any(is.na(x)), 1, 0)}))
chk <- filter(iris_clca3, del==1) ## one sample has missing value for hpiv4 only
# will recode as 0 for now to distinguish between single and co-infections
# chk2 <- filter(ftd_c, ic_ct >33|ic==0|is.na(ic))
table(iris_clca3$hpiv4, useNA = "always")
iris_clca3$hpiv4 <- ifelse(iris_clca3$study_id=="JC-0008", 0,  iris_clca3$hpiv4)
iris_clca3$case_pneu <- 0
iris_clca3$case_alri <- 0
iris_clca3$case_sev_alri <- 0
iris_clca3$control <- 1

iris_clca3$age_m <- floor(interval(iris_clca3$p_dob,iris_clca3$p_date)/months(1))
table(iris_clca3$age_m, useNA = "always")
iris_clca3 <- iris_clca3%>%
  dplyr::ungroup()%>%
  dplyr::mutate(rv_ev = case_when(
    rv==0 & ev==0 ~0,
    rv==1 | ev==1~1
  ))

## hospitalized
iris_hlca2 <- filter(iris_hlca, is.na(rv)==F)
iris_hlca2$del <- as.numeric(apply(iris_hlca2[,which(colnames(iris_hlca2) %in% labs2)], 1, function(x) {ifelse(any(is.na(x)), 1, 0)}))
table(iris_hlca2$del, useNA = "always")
chk <- filter(iris_hlca2, del==1)# have singleplex results for flu, don't need to drop

# Will use all (even without IC) for now
# get number of pathogens per participant
iris_hlca3 <- iris_hlca2%>%
  dplyr::ungroup()%>%
  dplyr::mutate(rv_ev = case_when(
    rv==0 & ev==0 ~0,
    rv==1 | ev==1~1
  ))
iris_hlca3$days_hosp <- as.numeric(interval(start = iris_hlca3$p_admitdate, end = iris_hlca3$d_date)%/% days(1))
iris_hlca3 <- filter(iris_hlca3, days_hosp >=0) #1 obs from philippines with earlier discharge date, dropping for now
table(iris_hlca3$m_dx_r_pneu___1, useNA = "always")
table(iris_hlca3$m_dx_r_pneu___3, useNA = "always")
iris_hlca3 <- filter(iris_hlca3, is.na(m_dx_r_pneu___1)==F) ## 1 missing response for m_dx_r_pneu___1 and m_dx_r_pneu___3


iris_hlca3 <- iris_hlca3%>%
  ungroup()%>%
  mutate(case_pneu = case_when(
    m_dx_r_pneu___1==1 | m_dx_r_pneu___3==1 ~1,
    m_dx_r_pneu___1==0 ~0
  ),
  case_alri = case_when(
    case_pneu==1| m_dx_r_ards___1==1 | m_dx_r_ards___3==1|
      m_dx_r_bronchio___1==1 | m_dx_r_bronchio___3==1|
      m_dx_r_respfail___1==1 | m_dx_r_respfail___3==1 ~1,
    T~0
  ),
  case_sev_alri = case_when(
    ((case_pneu==1| m_dx_r_ards___1==1 | m_dx_r_ards___3==1|
      m_dx_r_bronchio___1==1 | m_dx_r_bronchio___3==1|
      m_dx_r_respfail___1==1 | m_dx_r_respfail___3==1) & 
      (oxygenmode4_crf==1|oxygenmode5_crf==1|
         oxygenmode8_crf==1|oxygenmode9_crf==1|
         m_oxygen_mode___4==1|m_oxygen_mode___5==1|
         m_oxygen_mode___8==1|m_oxygen_mode___9==1
         ))~1,
    T~0
  ),
  age_m = floor(interval(p_dob, p_date)/months(1)))%>%
  filter(del==0)
iris_hlca3$control <- 0
table(iris_hlca3$case_pneu, useNA = "always") ## 1148/3630 = 31.6% hospitalized for pneumonia -- compare to PERCH
table(iris_hlca3$case_alri, useNA = "always") ## 1743/3630 = 48.0% hospitalized for ALRI -- compare to PERCH
table(iris_hlca3$case_sev_alri, useNA = "always") ## 77 
table(iris_hlca3$age_m, useNA = "always")
table(iris_hlca3$del, useNA = "always")

colnames(iris_hlca3)
# merge basic data (demo, labs, and case status)
hlca_sub <- iris_hlca3%>%
  select(study_id, country_id, control, age_m, p_sex, p_date, case_alri, case_pneu, all_of(vir2))

clca_sub <- iris_clca3%>%
  select(study_id, country_id, control, age_m, p_sex, p_date, case_alri, case_pneu, all_of(vir2))
tot_lca_sub <- bind_rows(hlca_sub, clca_sub)


tot_lca_sub$oth_cause <-  as.numeric(apply(tot_lca_sub[,which(colnames(tot_lca_sub) %in% vir2)], 1, function(x) {ifelse(any(x==1), 0, 1)}))

tot_lca_sub <- tot_lca_sub%>%
  rename_with(~  cause_list_iris, .cols = which(colnames(.)%in%vir3))

tot_lca_sub_pneu <- tot_lca_sub%>%
  dplyr::mutate(del = case_when(
    case_pneu==0 & control==0 ~1,
    T~0
  ))%>%
  filter(del==0)%>%
  select(-del)

tot_lca_sub_alri <- tot_lca_sub%>%
  dplyr::mutate(del = case_when(
    case_alri==0 & control==0 ~1,
    T~0
  ))%>%
  filter(del==0)%>%
  select(-del)
table(tot_lca_sub_alri$case_alri)

# save basic dataset
save(tot_lca_sub_pneu, file = "data/basic_lca_dat_pneu110321.rda")
save(tot_lca_sub_alri, file = "data/basic_lca_dat_alri110321.rda")

load("data/basic_lca_dat_alri081721.rda")


# no covariates dataset -- all sites
# arrange so cases on top
tot_lca_sub <- tot_lca_sub%>%
  ungroup()%>%
  arrange(-case_pneu)

MBS <- data.frame(select(tot_lca_sub, all_of(cause_list_iris)))
MBS1 <- list("MBS1"=MBS)
Mobs <- list("MBS"=MBS1,"MSS"=NULL,"MGS"=NULL)
Y = tot_lca_sub$case_pneu
data_iris_noreg <- list("Mobs"=Mobs, "Y"=Y, "X"=NULL)

save(data_iris_noreg, file = "data/data_iris_noreg081721.rda")


# site only dataset for discrete only model (1 variable)
# arrange so that cases on top, but also grouped by site
tot_lca_sub <- tot_lca_sub%>%
  ungroup()%>%
  arrange(-case_pneu, country_id)

Y <- tot_lca_sub$case_pneu
SITE <- tot_lca_sub$country_id
X <- list("SITE" = SITE)
data_iris_discrete <- list("Mobs"=Mobs, "Y"=Y, "X"=X)

save(data_iris_discrete, file = "data/data_iris_disc081721.rda")

# ALRI model adjusted for site, age, and sex

load("data/basic_lca_dat_alri110321.rda")

tot_lca_sub_alri_mult3 <- tot_lca_sub_alri%>%
  ungroup()%>%
  mutate(age_cat = as.factor(case_when(
    age_m <3 ~1,
    age_m >=3 & age_m <6 ~2,
    age_m >=6~3
  )),
sex = as.factor(p_sex))%>%
  arrange(-case_alri, country_id)%>%
  select(-`Other cause`)%>%
  rename(RV_EV = `RV/EV`)

save(tot_lca_sub_alri_mult3, file = "data/tot_lca_sub_alri_mult3fix.rda")


MBS <- data.frame(select(tot_lca_sub_alri_mult3, all_of(cause_list_iris)))
MBS1 <- list("MBS1"=MBS)
Mobs <- list("MBS"=MBS1,"MSS"=NULL,"MGS"=NULL)
Y <- tot_lca_sub_alri_mult3$case_alri
SITE <- tot_lca_sub_alri_mult3$country_id
X <- list("SITE" = SITE, "AGE" = tot_lca_sub_alri_mult3$age_cat, "SEX" = tot_lca_sub_alri_mult3$sex)

tot_lca_sub_alri_mult3fix <- list("Mobs"=Mobs, "Y"=Y, "X"=X)

save(tot_lca_sub_alri_mult3fix, file = "data/tot_lca_sub_alri_mult3_121621.rda")


# ALRI model adjusted for site, age, and sex -- without periods where cases/controls mis-matched
load("data/tot_lca_sub_alri_mult3fix.rda")

as.numeric(floor(interval("2015-06-29", "2017-04-20")/weeks(1)))
study_wks <- as.numeric(c(0:94))
country_id <- c("Albania", "Jordan", "Nicaragua","Philippines")
case_alri <- as.numeric(c(0,1))

temp_dat2 <- expand.grid("country_id" = country_id, "case_alri" = case_alri, "study_wk" = study_wks)

wk_enroll <- tot_lca_sub_alri%>%
  mutate(study_wk = as.numeric(floor(interval("2015-06-29", p_date)/weeks(1))))%>%
  ungroup()%>%
  group_by(country_id, case_alri, study_wk)%>%
  summarise(num = n())

tot_wk_enroll <- left_join(temp_dat2, wk_enroll, by=c("country_id", "case_alri", "study_wk"))
tot_wk_enroll$num <- ifelse(is.na(tot_wk_enroll$num), 0, tot_wk_enroll$num)

tot_wk_enroll <- tot_wk_enroll%>%
  ungroup()%>%
  arrange(country_id, study_wk)%>%
  mutate(prop_con = case_when(
    case_alri == 0 ~ num/(num + lead(num, n=1)),
    T~ NA_real_))

# start with week 15 in Nicaragua and pre 34 and post 53 in Philippines

tot_lca_sub_alri_mult3a <- tot_lca_sub_alri_mult3%>%
  ungroup()%>%
  mutate(study_wk = as.numeric(floor(interval("2015-06-29", p_date)/weeks(1))),
  #        age_cat = as.factor(case_when(
  #   age_m <3 ~1,
  #   age_m >=3 & age_m <6 ~2,
  #   age_m >=6~3
  # )),
  # RV_EV = case_when(
  #   RV==1 | EV==1 ~1,
  #   T~0
  # ), 
  # sex = as.factor(p_sex),
  del = case_when(
    country_id == "Nicaragua" & study_wk <15 ~1,
    country_id == "Philippines" & study_wk >33 & study_wk <54 ~ 1,
    T~0
  ))%>%
  # select(-c(RV, EV))%>%
  filter(del==0)%>%
  arrange(-case_alri, country_id)


MBS <- data.frame(select(tot_lca_sub_alri_mult3a, all_of(cause_list_iris)))
MBS1 <- list("MBS1"=MBS)
Mobs <- list("MBS"=MBS1,"MSS"=NULL,"MGS"=NULL)
Y <- tot_lca_sub_alri_mult3a$case_alri
SITE <- tot_lca_sub_alri_mult3a$country_id
X <- list("SITE" = SITE, "AGE" = tot_lca_sub_alri_mult3a$age_cat, "SEX" = tot_lca_sub_alri_mult3a$sex)

tot_lca_sub_alri_mult3afix <- list("Mobs"=Mobs, "Y"=Y, "X"=X)

save(tot_lca_sub_alri_mult3afix, file = "data/tot_lca_sub_alri_mult3a_121621.rda")

# site only dataset for discrete only model (1 variable) -- stratified by age
# Under 3 months
# pneumonia
tot_lca_sub_u3mp <- tot_lca_sub_pneu%>%
  filter(age_m <3)%>%
  ungroup()%>%
  arrange(-case_pneu, country_id)


Y <- tot_lca_sub_u3mp$case_pneu
SITE <- tot_lca_sub_u3mp$country_id
X <- list("SITE" = SITE)
data_iris_disc_u3p <- list("Mobs"=Mobs, "Y"=Y, "X"=X)

save(data_iris_disc_u3p, file = "data/data_iris_disc_u3p081721.rda")

#alri
tot_lca_sub_u3ma <- tot_lca_sub_alri%>%
  filter(age_m <3)%>%
  ungroup()%>%
  arrange(-case_alri, country_id)


Y <- tot_lca_sub_u3ma$case_alri
SITE <- tot_lca_sub_u3ma$country_id
X <- list("SITE" = SITE)
data_iris_disc_u3a <- list("Mobs"=Mobs, "Y"=Y, "X"=X)

save(data_iris_disc_u3a, file = "data/data_iris_disc_u3a081721.rda")

# 3-5 months
# pneumonia
tot_lca_sub_35mp <- tot_lca_sub_pneu%>%
  filter(age_m >=3, age_m <6)%>%
  ungroup()%>%
  arrange(-case_pneu, country_id)


Y <- tot_lca_sub_35mp$case_pneu
SITE <- tot_lca_sub_35mp$country_id
X <- list("SITE" = SITE)
data_iris_disc_35p <- list("Mobs"=Mobs, "Y"=Y, "X"=X)

save(data_iris_disc_35p, file = "data/data_iris_disc_35p081721.rda")

# alri
tot_lca_sub_35ma <- tot_lca_sub_alri%>%
  filter(age_m >=3, age_m<6)%>%
  ungroup()%>%
  arrange(-case_alri, country_id)


Y <- tot_lca_sub_35ma$case_alri
SITE <- tot_lca_sub_35ma$country_id
X <- list("SITE" = SITE)
data_iris_disc_35a <- list("Mobs"=Mobs, "Y"=Y, "X"=X)

save(data_iris_disc_35a, file = "data/data_iris_disc_35a081721.rda")

#6-11 months
# pneumonia
tot_lca_sub_o6mp <- tot_lca_sub_pneu%>%
  filter(age_m >=6)%>%
  ungroup()%>%
  arrange(-case_pneu, country_id)


Y <- tot_lca_sub_o6mp$case_pneu
SITE <- tot_lca_sub_o6mp$country_id
X <- list("SITE" = SITE)
data_iris_disc_o6p <- list("Mobs"=Mobs, "Y"=Y, "X"=X)

save(data_iris_disc_o6p, file = "data/data_iris_disc_o6p081721.rda")


# alri
tot_lca_sub_o6ma <- tot_lca_sub_alri%>%
  filter(age_m >=6)%>%
  ungroup()%>%
  arrange(-case_alri, country_id)


Y <- tot_lca_sub_o6ma$case_alri
SITE <- tot_lca_sub_o6ma$country_id
X <- list("SITE" = SITE)
data_iris_disc_o6a <- list("Mobs"=Mobs, "Y"=Y, "X"=X)

save(data_iris_disc_o6a, file = "data/data_iris_disc_o6a081721.rda")

# discrete only model with site, age, and sex as covariates
tot_lca_sub_multd <-  tot_lca_sub_alri%>%
  ungroup()%>%
  mutate(age_d = factor(case_when(
    age_m < 6 ~ 1,
    T ~ 2
  ), ordered = F),
  sex = as.factor(p_sex))%>%
  arrange(-case_pneu, country_id)

Y <- tot_lca_sub_multd$case_pneu
SITE <- tot_lca_sub_multd$country_id
X <- list("SITE" = SITE, "AGE" = tot_lca_sub_multd$age_d, "SEX" = tot_lca_sub_multd$sex)
Z <- list()
J <- list()
data_iris_disc_multd <- list("Mobs"=Mobs, "Y"=Y, "X"=X)

save(data_iris_disc_multd, file = "data/data_iris_disc_multd072121.rda")
