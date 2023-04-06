##############################
###   IRIS data figures    ###
###   Created: 12/03/20    ###
###   Updated: 12/03/20    ###
##############################

library(here)
library(dplyr)
library(ggplot2)
library(ggupset)
library(lubridate)
library(tidyr)
library(RColorBrewer)


# Load data ---------------------------------------------------------------


# Sampling of cases and controls over time --------------------------------
h_sub <- ftd_h%>%
  dplyr::select(study_id, country_id, study_year, p_dob, hosp_infage_cat,
                p_sex, e_specimen_resp_date)%>%
  dplyr::mutate(study = "hosp")

c_sub <- ftd_c%>%
  dplyr::select(study_id, country_id, study_year, p_dob, hosp_infage_cat,
                p_sex, e_specimen_resp_date)%>%
  dplyr::mutate(study = "non-ill")

tot_sub <- bind_rows(c_sub, h_sub)
min(tot_sub$e_specimen_resp_date) # will start study week on Sunday 2015-06-21
max(tot_sub$e_specimen_resp_date) # end dates (see below) on Saturday 2017-04-22 
tot_sub$study_wk <- as.numeric(interval(start = "2015-06-21", end = tot_sub$e_specimen_resp_date)%/% weeks(1)+1)

tot_sum <- tot_sub%>%
  ungroup()%>%
  group_by(study_wk, country_id, study)%>%
  summarise(num = n())


wk_dates <- data.frame("date" = seq(as.Date("2015-06-21", format="%Y-%m-%d"), as.Date("2017-04-22", format="%Y-%m-%d"), by="weeks"))
wk_dates <- wk_dates%>%
  mutate(study_wk = row_number())

dates <- seq(as.Date("2015-06-21", format="%Y-%m-%d"), as.Date("2017-04-22", format="%Y-%m-%d"), by="weeks")

country_id <- c("Albania", "Jordan", "Nicaragua", "Philippines")
study <- c("hosp", "non-ill")
dates2 <- expand.grid(dates, country_id, study)
names(dates2) <- c("date", "country_id", "study")

samp_wk1 <- left_join(dates2, wk_dates, by="date")
samp_wk2 <- left_join(samp_wk1, tot_sum, by=c("study_wk", "country_id", "study"))
samp_wk2$num <- ifelse(is.na(samp_wk2$num),0,samp_wk2$num)
table(samp_wk2$num, useNA = "always")
samp_wk2$studyf <- factor(samp_wk2$study)
## need to tweak colors
scico(2, palette = "vik", end=0.8)

cols2 <- c("hosp"="red", "non-ill"= "blue")
cols2a <- c("hosp"="#B2182B", "non-ill"= "#2166AC")
cols2b <- c("hosp"="#001260", "non-ill"= "#9F711B")
table(samp_wk1$study_wk, useNA = "always")

samp1<-ggplot(samp_wk2, aes(x=date, y=num, group=study)) +
  geom_bar(aes(fill=as.factor(study)), stat = "identity")+ # specifies plot type (geom_line = line) and varies color by case/control status
  labs(x = " ", y = "Number of samples", size = 1)+ # Specify axis labels
  scale_x_date(date_labels="%b-%y", date_breaks="2 months", expand=c(0,0))+  # sets format for months-year ("%b-%y") and displays only every 2 months
  scale_fill_manual(values = cols2b, labels =c("Hospitalized", "Non-ill"), name = NULL)+ # assigns pre specified colors to case/control status
  theme_bw()+ # removes some background elements for cleaner look
  theme(axis.text.y = element_text(size = 20))+ # sets rotation and location of y axis labels
  theme(axis.text.x = element_text(size = 20,angle = 45, vjust = 1, hjust=1))+ # sets rotation and location of x axis labels
  # facet_grid(rows = vars(country_id))+# switch moves country label from top to bottom, scales and space make month/yr readable
  # facet_wrap(~ country_id, ncol=1)+
  facet_grid(rows=vars(country_id))+
  theme(strip.text = element_text(size = 20))+
  theme(legend.text=element_text(size=20))
samp1

png(file="figures/sample_timeV3b.png",width=3600,height=4000,res = 300)
samp1+theme(axis.title.y = element_text(size=24))
dev.off()

# shows potential case/control mismatch in Nica at beginning of study and 
# in middle of study in the Philippines

# Explore data ------------------------------------------------------------

# Look at histograms of days between onset and sampling across infection types
chk5 <- filter(ftd_h, infect==1)
chk6 <- filter(ftd_h, infect==2)
chk7 <- filter(ftd_h, infect==0)

hist(chk7$day_dif)
hist(chk5$day_dif)
hist(chk6$day_dif)


test <- glm(infect_d ~ day_dif, family = "binomial", data = chk4)
summary(test)
ftd_c$case <- 0
ftd_h$case <- 1

ftd_tot <- bind_rows(ftd_h, ftd_c)
ftd_tot$infect_d <- ifelse(ftd_tot$infect >1, 1,0)
table(ftd_tot$infect, useNA = "always")
table(ftd_tot$infect_d, useNA = "always")

coinf1 <- glm(case ~ as.factor(infect_d), family = "binomial", data=ftd_tot)
summary(coinf1)
vir2


# Upset plots of co-infection combinations --------------------------------
ftd_tot <- ftd_tot%>%
  ungroup()%>%
  mutate(rv2 = case_when(
    rv==0 ~0,
    rv==1 & ev==0~1,
    rv==1 & ev==1~0,
    T~99
  ))
table(ftd_tot$rv,useNA = "always")
vir3 <- c("flua_univ_result","flub_univ_result","rsv_result","rv2","cor63","cor229",
          "cor43","corhku1","hpiv2","hpiv3","hpiv4","hpiv1","hmpvab","hbov", "mpneu",
          "hadv","ev","hpev")

ftd_tot2 <- ftd_tot %>%
  as_tibble() %>%
  filter(infect >1)%>%
  select(study_id, all_of(vir3))%>%
  rename(Flu_A = flua_univ_result,
         Flu_B = flub_univ_result)%>%
  gather(virus, Member, -study_id) %>%
  filter(Member==1) %>%
  select(- Member)%>%
  group_by(study_id)%>%
  summarise(virus = list(virus))

cols2 <- c("0" = "#78C679","1" ="#006837")

test <- ftd_tot2%>%
  # distinct(title, year, length, .keep_all=TRUE) %>%
  ggplot(aes(x=virus)) +
  geom_bar(fill="#006837") +
  # scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2))+
  ylab("Count")+
  xlab(" ")+
  scale_x_upset(n_intersections = 15)+
  theme_classic()+
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  theme_combmatrix(combmatrix.panel.point.color.fill = "#006837",
                   combmatrix.panel.line.size = 0,
                   combmatrix.label.text = element_text(size=16))

test

# Seasonal coronaviruses----
cor_dat_h <-  ftd_h%>%
  dplyr::select(cor63, cor229, cor43, corhku1, infect_d)%>%
  dplyr::group_by(infect_d)%>%
  summarise(across(.cols = all_of(cor),sum))%>%
  melt(id.vars = "infect_d")%>%
  dplyr::mutate(infect = factor(infect_d, levels = c("0","1"), 
                                labels = c("Single infection", "Co-infection")))%>%
  ungroup()%>%
  dplyr::group_by(variable)%>%
  dplyr::mutate(tot = sum(value),
                prop = round(value/tot,2)*100)

cor_dat_hc <-  ftd_h%>%
  dplyr::select(cor63, cor229, cor43, corhku1, infect_d, country_id)%>%
  dplyr::group_by(country_id, infect_d)%>%
  summarise(across(.cols = all_of(cor),sum))%>%
  melt(id.vars = c("country_id","infect_d"))%>%
  dplyr::mutate(infect = factor(infect_d, levels = c("0","1"), 
                                labels = c("Single infection", "Co-infection")))%>%
  ungroup()%>%
  dplyr::group_by(country_id, variable)%>%
  dplyr::mutate(tot = sum(value),
                prop = round(value/tot,2)*100)

cor_dat_cc <-  ftd_c%>%
  dplyr::select(cor63, cor229, cor43, corhku1, infect_d, country_id)%>%
  dplyr::group_by(country_id, infect_d)%>%
  summarise(across(.cols = all_of(cor),sum))%>%
  melt(id.vars = c("country_id","infect_d"))%>%
  dplyr::mutate(infect = factor(infect_d, levels = c("0","1"), 
                                labels = c("Single infection", "Co-infection")))%>%
  ungroup()%>%
  dplyr::group_by(country_id, variable)%>%
  dplyr::mutate(tot = sum(value),
                prop = round(value/tot,2)*100)

# proportion of coronavirus pos that were single vs. co-inf marginalized across site
png(file="Box/IRIS_/figures/cor_proph.png", width=5000,height=3000,res = 300)
cor_bplot(cor_dat_h, prop,"Proportion")+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

# number of coronavirus pos that were single vs. co-inf marginalized across site
png(file="Box/IRIS_/figures/cor_numh.png", width=5000,height=3000,res = 300)
cor_bplot(cor_dat_h, value,"Number")+
  scale_y_continuous(limits = c(0,80),breaks = c(0,20,40,60,80), labels = c(0,20,40,60,80))+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

# Barplot of proportion of coronaviruses that were single vs. co-infections among hosp
png(file="Box/IRIS_/figures/cor_hpanprop.png", width=6000,height=3000,res = 300)
cor_bplot(cor_dat_hc, prop,"Proportion")+
  facet_grid(~country_id)+
  theme(axis.text.x = element_text(size = 20),
        strip.text.x = element_text(size = 16))
dev.off()

# Barplot of proportion of coronaviruses that were single vs. co-infections among non-ill
png(file="Box/IRIS_/figures/cor_cpanprop.png", width=6000,height=3000,res = 300)
cor_bplot(cor_dat_cc, prop,"Proportion")+
  facet_grid(~country_id)+
  theme(axis.text.x = element_text(size = 20),
        strip.text.x = element_text(size = 16))
dev.off()

# Barplot of number of coronaviruses that were single vs. co-infections among hosp
png(file="Box/IRIS_/figures/cor_hpannum.png", width=6000,height=3000,res = 300)
cor_bplot(cor_dat_hc, value,"Number")+
  facet_grid(~country_id)+
  scale_y_continuous(limits = c(0,30),breaks = c(0,10,20,30), labels = c(0,10,20,30))+
  theme(axis.text.x = element_text(size = 20),
        strip.text.x = element_text(size = 16))
dev.off()

# Barplot of number of coronaviruses that were single vs. co-infections among non-ill
png(file="Box/IRIS_/figures/cor_cpannum.png", width=6000,height=3000,res = 300)
cor_bplot(cor_dat_cc, value,"Number")+
  facet_grid(~country_id)+
  scale_y_continuous(breaks = c(0,3,6,9))+
  theme(axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 16))
dev.off()

# What to coronaviruses tend to co-infect with?
colnames(ftd_tot_sub)

cor_pos <- ftd_tot_sub%>%
  filter(cor43==1|cor63==1|cor229==1|corhku1==1)%>%
  mutate(num_in = as.numeric(apply(.[,which(colnames(.) %in% vir_ord)], 1, sum)))%>%
  filter(num_in>1)

# upset plots
cor_upset(cor_upmunge(cor_pos))## still not working need to figure out renaming


cor_pos <- ftd_tot_sub%>%
  filter(cor43==1|cor63==1|cor229==1|corhku1==1)%>%
  mutate(num_in = as.numeric(apply(.[,which(colnames(.) %in% vir_ord)], 1, sum)))%>%
  filter(num_in>1)

cor_posh <- ftd_h%>%
  filter(cor43==1|cor63==1|cor229==1|corhku1==1)%>%
  # mutate(num_in = as.numeric(apply(.[,which(colnames(.) %in% vir_ord)], 1, sum)))%>%
  filter(num_in>1)

cor_posc <- ftd_c%>%
  filter(cor43==1|cor63==1|cor229==1|corhku1==1)%>%
  # mutate(num_in = as.numeric(apply(.[,which(colnames(.) %in% vir_ord)], 1, sum)))%>%
  filter(num_in>1)

cor_pos <- cor_pos %>%
  as_tibble() %>%
  mutate(fid = row_number())%>%
  select(fid, all_of(vir_ord))%>%
  rename_with(~ newnames, .cols = colnames(.))

colnames(cor_pos) <- c("FID",vir_labs)

cor_pos <- cor_pos%>%
  gather(virus, Member, -FID) %>%
  filter(Member==1) %>%
  select(- Member)%>%
  group_by(FID)%>%
  summarise(virus = list(virus))


tot_corup <- cor_pos%>%
  # distinct(title, year, length, .keep_all=TRUE) %>%
  ggplot(aes(x=virus)) +
  geom_bar(fill="#B2182B") +
  #scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2))+
  ylab("Number")+
  xlab(" ")+
  scale_x_upset(n_intersections = 11)+
  theme_classic()+
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  theme_combmatrix(combmatrix.panel.point.color.fill = "#B2182B",
                   combmatrix.panel.line.size = 0,
                   combmatrix.label.text = element_text(size=16))

tot_corup


cor_posh <- cor_posh %>%
  as_tibble() %>%
  mutate(fid = row_number())%>%
  select(fid, all_of(vir_ord))

colnames(cor_posh) <- c("FID",vir_labs)

cor_posh <- cor_posh%>%
  gather(virus, Member, -FID) %>%
  filter(Member==1) %>%
  select(- Member)%>%
  group_by(FID)%>%
  summarise(virus = list(virus))

#cols2 <- c("0" = "#78C679","1" ="#006837")
library(ggupset)
h_corup <- cor_posh%>%
  # distinct(title, year, length, .keep_all=TRUE) %>%
  ggplot(aes(x=virus)) +
  geom_bar(fill="#B2182B") +
  #scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2))+
  ylab("Number")+
  xlab(" ")+
  scale_x_upset(n_intersections = 11)+
  theme_classic()+
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  theme_combmatrix(combmatrix.panel.point.color.fill = "#B2182B",
                   combmatrix.panel.line.size = 0,
                   combmatrix.label.text = element_text(size=16))

h_corup


cor_posc <- cor_posc %>%
  as_tibble() %>%
  mutate(fid = row_number())%>%
  select(fid, all_of(vir_ord))

colnames(cor_posc) <- c("FID",vir_labs)

cor_posc <- cor_posc%>%
  gather(virus, Member, -FID) %>%
  filter(Member==1) %>%
  select(- Member)%>%
  group_by(FID)%>%
  summarise(virus = list(virus))

c_corup <- cor_posc%>%
  # distinct(title, year, length, .keep_all=TRUE) %>%
  ggplot(aes(x=virus)) +
  geom_bar(fill="#B2182B") +
  #scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2))+
  ylab("Number")+
  xlab(" ")+
  scale_x_upset(n_intersections = 11)+
  theme_classic()+
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  theme_combmatrix(combmatrix.panel.point.color.fill = "#B2182B",
                   combmatrix.panel.line.size = 0,
                   combmatrix.label.text = element_text(size=16))

c_corup
png(file="Box/IRIS_/figures/cortot_upset.png",width=3600,height=3000,res = 300)
tot_corup
dev.off()

png(file="Box/IRIS_/figures/corh_upset.png",width=3600,height=3000,res = 300)
h_corup
dev.off()

png(file="Box/IRIS_/figures/corc_upset.png",width=3600,height=3000,res = 300)
c_corup
dev.off()

### Single vs. co-infections----
ftd_tot <- subset_merge()

ftd_tot <- ftd_tot%>%
  dplyr::mutate(num_in = as.numeric(apply(ftd_tot[,which(colnames(ftd_tot) %in% vir_ord)], 1, sum)))%>%
  rename_with(~  vir_labs, .cols = which(colnames(.)%in%vir_ord))%>%
  filter(num_in>0)%>%
  dplyr::mutate(infect_d = case_when(
    num_in == 1 ~ 0,
    num_in >1 ~1
  ))%>%
  dplyr::group_by(infect_d)%>%
  summarise(across(.cols = all_of(vir_labs),sum))%>%
  melt(id.vars = "infect_d")%>%
  dplyr::mutate(infect = factor(infect_d, levels = c("0","1"), 
                                labels = c("Single infection", "Co-infection")))%>%
  ungroup()%>%
  dplyr::group_by(variable)%>%
  dplyr::mutate(tot = sum(value),
                prop = round(value/tot,2)*100)

ftd_h1 <- ftd_h%>%
  select(country_id, e_specimen_resp_date, all_of(vir_ord), num_in)%>%
  rename_with(~  vir_labs, .cols = which(colnames(.)%in%vir_ord))%>%
  filter(num_in>0)%>%
  dplyr::mutate(infect_d = case_when(
    num_in == 1 ~ 0,
    num_in >1 ~1
  ))%>%
  dplyr::group_by(infect_d)%>%
  summarise(across(.cols = all_of(vir_labs),sum))%>%
  melt(id.vars = "infect_d")%>%
  dplyr::mutate(infect = factor(infect_d, levels = c("0","1"), 
                                labels = c("Single infection", "Co-infection")))%>%
  ungroup()%>%
  dplyr::group_by(variable)%>%
  dplyr::mutate(tot = sum(value),
                prop = round(value/tot,2)*100)

ftd_c1 <- ftd_c%>%
  select(country_id, e_specimen_resp_date, all_of(vir_ord), num_in)%>%
  rename_with(~  vir_labs, .cols = which(colnames(.)%in%vir_ord))%>%
  filter(num_in>0)%>%
  dplyr::mutate(infect_d = case_when(
    num_in == 1 ~ 0,
    num_in >1 ~1
  ))%>%
  dplyr::group_by(infect_d)%>%
  summarise(across(.cols = all_of(vir_labs),sum))%>%
  melt(id.vars = "infect_d")%>%
  dplyr::mutate(infect = factor(infect_d, levels = c("0","1"), 
                                labels = c("Single infection", "Co-infection")))%>%
  ungroup()%>%
  dplyr::group_by(variable)%>%
  dplyr::mutate(tot = sum(value),
                prop = round(value/tot,2)*100)

# Hosp
png(file="Box/IRIS_/figures/hsing_v_co.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_h1, prop, "Proportion")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

png(file="Box/IRIS_/figures/hsing_v_conum.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_h1, value, "Number")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

ftd_h_sub <- ftd_h1%>%
  filter(!variable %in% c("RSV","Rhinovirus"))

png(file="Box/IRIS_/figures/hsing_v_co_sub.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_h_sub, prop, "Proportion")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

png(file="Box/IRIS_/figures/hsing_v_conum_sub.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_h_sub, value, "Number")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

# Non-ill
png(file="Box/IRIS_/figures/csing_v_co.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_c1, prop, "Proportion")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

png(file="Box/IRIS_/figures/csing_v_conum.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_c1, value, "Number")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

ftd_c_sub <- ftd_c1%>%
  filter(!variable %in% c("RSV","Rhinovirus"))

png(file="Box/IRIS_/figures/csing_v_co_sub.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_c_sub, prop, "Proportion")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

png(file="Box/IRIS_/figures/csing_v_conum_sub.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_c_sub, value, "Number")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

# Total
png(file="Box/IRIS_/figures/sing_v_co.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_tot, prop, "Proportion")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

png(file="Box/IRIS_/figures/sing_v_conum.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_tot, value, "Number")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

#w/o RSV and RV
ftd_tot_sub <- ftd_tot%>%
  filter(!variable %in% c("RSV","Rhinovirus"))

png(file="Box/IRIS_/figures/sing_v_co_sub.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_tot_sub, prop, "Proportion")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

png(file="Box/IRIS_/figures/sing_v_conum_sub.png",width=4500,height=2000,res = 300)
cor_bplot(ftd_tot_sub, value, "Number")+
  coord_flip()+
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24))
dev.off()

# Heatmaps ----------------------------------------------------------------

# tot <- filter(ftd_tot, country_id=="Nicaragua")
heat_mprep <- function(input, outname){
  output <- filter(input, infect_d==1)
  output <- output[,which(colnames(output) %in% vir2)]
  output <- output%>%
    rename(flu_a = flua_univ_result,
           flu_b = flub_univ_result,
           rsv = rsv_result)
  assign(x = outname, value = output, envir = globalenv())
}

heat_munge <- function (input, outname){
  names = colnames(input)
  output =  matrix(, nrow = length(names), ncol = length(names))
  rownames(output) = toupper(names)
  colnames(output) = toupper(names)
  
  for (i in 1:ncol(output))
  {for (j in 1:nrow(output)) 
  {output[i,j] <- ifelse(sum(input[,i]&input[,j])>0,sum(input[,i]&input[,j]),0) }}
  
  diag(output) <- NA
  # outname <<- output
  assign(x = outname, value = output, envir = globalenv())
}

png(file="Box/IRIS_/nica_hosp_heat.png",width=3600,height=2100,res = 300)
ggplot(aes(x=Var1, y=Var2), data=(melt(cor(m2))))+
  geom_tile(aes(fill=value))+
  xlab(NULL)+
  ylab(NULL)+
  theme(axis.text.x = element_text(size=20, angle=45, hjust = 1),
        axis.text.y = element_text(size=20),
        legend.title = element_text( size=24),
        legend.text = element_text(size = 16))+
  labs(fill="Number")+
  scale_fill_gradient('Proportion', limit = c(0, 0.3), low = "blue", high = "red") 
dev.off()



### Viral prevalence----
# Total
h_tot <- ftd_h%>%
  group_by(country_id)%>%
  summarise(tot_n = n(),
            across(.cols = all_of(vir_ord),sum))%>%
  rename_with(~  vir_labs, .cols = which(colnames(.)%in%vir_ord))%>%
  mutate(across(.cols = all_of(vir_labs), ~ round(.x/tot_n*100,2)))%>%
  pivot_longer(!c(country_id, tot_n), names_to = "virus", values_to = "prevalence")%>%
  mutate(virus2 = factor(virus, levels = vir_labs))

h_prev_tot <- ggplot(h_tot, aes(fill=country_id, y=prevalence, x=virus2)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 40, 10))+
  coord_flip()+
  labs(y = "Prevalence", x = " ", size=1)+
  scale_fill_manual(values = cols4, breaks = country_id, labels = country_id, name=NULL)+
  theme_classic()+
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        legend.text = element_text(size=20))

c_tot <- ftd_c%>%
  group_by(country_id)%>%
  summarise(tot_n = n(),
            across(.cols = all_of(vir_ord),sum))%>%
  rename_with(~  vir_labs, .cols = which(colnames(.)%in%vir_ord))%>%
  mutate(across(.cols = all_of(vir_labs), ~ round(.x/tot_n*100,2)))%>%
  pivot_longer(!c(country_id, tot_n), names_to = "virus", values_to = "prevalence")%>%
  mutate(virus2 = factor(virus, levels = vir_labs))

c_prev_tot <- ggplot(c_tot, aes(fill=country_id, y=prevalence, x=virus2)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 40, 10))+
  coord_flip()+
  labs(y = "Prevalence", x = " ", size=1)+
  scale_fill_manual(values = cols4, breaks = country_id, labels = country_id, name=NULL)+
  theme_classic()+
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        legend.text = element_text(size=20))

# By age
c_age <- ftd_c%>%
  group_by(country_id, hosp_infage_dic)%>%
  summarise(tot_n = n(),
            across(.cols = all_of(vir_ord),sum))%>%
  rename_with(~  vir_labs, .cols = which(colnames(.)%in%vir_ord))%>%
  mutate(across(.cols = all_of(vir_labs), ~ round(.x/tot_n*100,2)))%>%
  pivot_longer(!c(country_id, hosp_infage_dic, tot_n), names_to = "virus", values_to = "prevalence")%>%
  mutate(virus2 = factor(virus, levels = vir_labs),
         age = factor(hosp_infage_dic, labels = c("<6 months","6-12 months")))

c_prev_age <- ggplot(c_age, aes(fill=country_id, y=prevalence, x=virus2)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 40, 10))+
  coord_flip()+
  labs(y = "Prevalence", x = " ", size=1)+
  scale_fill_manual(values = cols4, breaks = country_id, labels = country_id, name=NULL)+
  theme_classic()+
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        legend.text = element_text(size=20))+
  facet_wrap(~age)

h_age <- ftd_h%>%
  group_by(country_id, hosp_infage_dic)%>%
  summarise(tot_n = n(),
            across(.cols = all_of(vir_ord),sum))%>%
  rename_with(~  vir_labs, .cols = which(colnames(.)%in%vir_ord))%>%
  mutate(across(.cols = all_of(vir_labs), ~ round(.x/tot_n*100,2)))%>%
  pivot_longer(!c(country_id, hosp_infage_dic, tot_n), names_to = "virus", values_to = "prevalence")%>%
  mutate(virus2 = factor(virus, levels = vir_labs),
         age = factor(hosp_infage_dic, labels = c("<6 months","6-12 months")))

h_prev_age <- ggplot(h_age, aes(fill=country_id, y=prevalence, x=virus2)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 40, 10))+
  coord_flip()+
  labs(y = "Prevalence", x = " ", size=1)+
  scale_fill_manual(values = cols4, breaks = country_id, labels = country_id, name=NULL)+
  theme_classic()+
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        legend.text = element_text(size=20))+
  facet_wrap(~age)

# By sex
h_sex <- ftd_h%>%
  group_by(country_id,p_sex)%>%
  summarise(tot_n = n(),
            across(.cols = all_of(vir_ord),sum))%>%
  rename_with(~  vir_labs, .cols = which(colnames(.)%in%vir_ord))%>%
  mutate(across(.cols = all_of(vir_labs), ~ round(.x/tot_n*100,2)))%>%
  pivot_longer(!c(country_id,p_sex, tot_n), names_to = "virus", values_to = "prevalence")%>%
  mutate(virus2 = factor(virus, levels = vir_labs))

c_sex <- ftd_c%>%
  group_by(country_id,p_sex)%>%
  summarise(tot_n = n(),
            across(.cols = all_of(vir_ord),sum))%>%
  rename_with(~  vir_labs, .cols = which(colnames(.)%in%vir_ord))%>%
  mutate(across(.cols = all_of(vir_labs), ~ round(.x/tot_n*100,2)))%>%
  pivot_longer(!c(country_id,p_sex, tot_n), names_to = "virus", values_to = "prevalence")%>%
  mutate(virus2 = factor(virus, levels = vir_labs))

hosp_prev_s <- ggplot(h_sex, aes(fill=country_id, y=prevalence, x=virus2)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 40, 10))+
  coord_flip()+
  labs(y = "Prevalence", x = " ", size=1)+
  scale_fill_manual(values = cols4, breaks = country_id, labels = country_id, name=NULL)+
  theme_classic()+
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        legend.text = element_text(size=20))+
  facet_wrap(~p_sex)

c_prev_s <- ggplot(c_sex, aes(fill=country_id, y=prevalence, x=virus2)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 40, 10))+
  coord_flip()+
  labs(y = "Prevalence", x = " ", size=1)+
  scale_fill_manual(values = cols4, breaks = country_id, labels = country_id, name=NULL)+
  theme_classic()+
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        legend.text = element_text(size=20))+
  facet_wrap(~p_sex)

png(file="Box/IRIS_/figures/hosp_prev_age.png",width=3600,height=2000,res = 300)
h_prev_age
dev.off()

png(file="Box/IRIS_/figures/hosp_prev_sex.png",width=3600,height=2000,res = 300)
hosp_prev_s
dev.off()

png(file="Box/IRIS_/figures/hosp_prev_tot.png",width=3600,height=2000,res = 300)
h_prev_tot
dev.off()

png(file="Box/IRIS_/figures/nonill_prev_sex.png",width=3600,height=2000,res = 300)
c_prev_s
dev.off()

png(file="Box/IRIS_/figures/nonill_prev_age.png",width=3600,height=2000,res = 300)
c_prev_age
dev.off()

png(file="Box/IRIS_/figures/nonill_prev_tot.png",width=3600,height=2000,res = 300)
c_prev_tot
dev.off()
