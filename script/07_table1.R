# Table 1

library(gtsummary)
library(flextable)
load("data/ftd_h120622.rda")
load("data/ftd_c120622.rda")

h_table1 <- table_vars(ftd_h)
c_table1 <- table_vars(ftd_c)

table(ftd_c$e_father_edu, useNA = "always")

## To do:
# relevel/reorder household smoking variable levels
# fix house_smoke and house_smk both showing up in table
## try improving autofit?


table1_tot <- bind_rows(h_table1, c_table1)%>%
  mutate(house_smk = factor(house_smoke,ordered = F,
                            levels = c("Never",
                                       "< Monthly",
                                       "Monthly",
                                       "Weekly",
                                       "Daily")))

table(table1_tot$male, useNA = "always")
table(table1_tot$house_smk, useNA = "always")

class(table1_tot$house_smoke)

table1 <- table1_tot%>%
tbl_strata(
  strata = country_id,
  .tbl_fun =
    ~.x %>%
    tbl_summary(
      by = cohort,
  statistic = list(c("p_calc_age") ~ "{mean} ({sd})",
                   all_categorical() ~ "{n} ({p}%)"),
  digits = all_continuous() ~ 2,
  label = list(p_calc_age ~ "Age (weeks)",
               male ~ "Male",
               parent_sec_edu ~ "Parent education \n (secondary or higher)",
               low_bw ~ "Low birthweight",
               house_smk ~ "Smoking in household"),
  missing = "no"
))%>%
  as_flex_table()
  # add_n(),
  # .header = "**{strata}**")

flextable::save_as_docx(table1, path = "figures/table1.docx")


?chisq.test

## see if proportion with at least 1 pathogen detected varied by cohort
x <- matrix(c(2547,1085, 513,555), ncol = 2)
x
chisq.test(x)

## Extra code

# age in months - p_calc_age
# male - p_sex -> create male ind 
# parent secondary edu+ -> create parent edu e_mother_edu, e_father_edu
# low birth weight -> f_birthweight, e_mat_birthweight, e_mat_birthweight_verify
# household smoking -> e_smoke
# by cohort (hosp or non-ill) and site

# ftd_h <- ftd_h%>%
#   ungroup()%>%
#   mutate(birthweight = case_when(
#     is.na(e_mat_birthweight)==F & is.na(f_birthweight)==T ~ e_mat_birthweight,
#     is.na(e_mat_birthweight)==F & is.na(f_birthweight)==F ~ f_birthweight,
#     is.na(e_mat_birthweight)==T & is.na(f_birthweight)==F ~ f_birthweight,
#     T ~ NA_real_
#   ),
#   low_bw = factor(case_when(
#     is.na(birthweight)==F & birthweight < 2.5 ~ "1",
#     is.na(birthweight)==F & birthweight >= 2.5 ~ "0",
#     is.na(birthweight)==T ~ NA_character_
#   ), ordered = F),
#   male = factor(case_when(
#     p_sex=="Male"~1,
#     T~0
#   ), ordered = F),
#   parent_sec_edu = factor(case_when(
#     e_mother_edu %in% c(5:9) | e_father_edu %in% c(5:9)~"1",
#     e_mother_edu %in% c(1:4) & e_father_edu %in% c(1:4)~"0",
#     e_mother_edu %in% c(77,88,99)|e_father_edu %in% c(77,88,99)~NA_character_
#   ), ordered = F))
# 
# table(ftd_h$male, useNA = "always")
# table(ftd_h$parent_sec_edu, useNA = "always")
# 
# class(ftd_h$e_mother_edu)
# 
# ftd_c <- ftd_c%>%
#   ungroup()%>%
#   mutate(birthweight = case_when(
#     is.na(e_mat_birthweight)==F & is.na(f_birthweight)==T ~ e_mat_birthweight,
#     is.na(e_mat_birthweight)==F & is.na(f_birthweight)==F ~ f_birthweight,
#     is.na(e_mat_birthweight)==T & is.na(f_birthweight)==F ~ f_birthweight,
#     T ~ NA_real_
#   ),
#   low_bw = factor(case_when(
#     is.na(birthweight)==F & birthweight < 2.5 ~ "1",
#     is.na(birthweight)==F & birthweight >= 2.5 ~ "0",
#     is.na(birthweight)==T ~ NA_character_
#   ), ordered = F))
# 
# 
# table(ftd_h$low_bw, useNA = "always")
# 
# table(ftd_h$e_mat_birthweight_verify, useNA = "always")
# table(ftd_h$e_mat_birthweight, useNA = "always")
# 
# chk <- filter(ftd_h, is.na(e_mat_birthweight))
# table(chk$f_birthweight, useNA = "always")
# 
# chk2 <- filter(ftd_c, is.na(e_mat_birthweight))
# table(chk2$f_birthweight, useNA = "always")
