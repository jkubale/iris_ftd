# Functions and vectors for FTD analysis

### Vectors----
vir_ord <- c("hpev",  "mpneu", "hpiv1", "hpiv2",  "hpiv3", "hpiv4", 
             "cor229", "cor63", "corhku1", "cor43", "flub_univ_result",
             "flua_univ_result", "hbov", "hmpvab", "hadv", "rv_ev", "rsv_result")

vir_labs <- c("Parechovirus", "M. pneumoniae", "HPIV1", "HPIV2", "HPIV3", "HPIV4", 
              "Cor229", "Cor63", "CorHKU1", "Cor43", "Influenza B",
               "Influenza A", "Bocavirus",
              "HMPV", "Adenovirus", "RV/EV", "RSV")

vir_ord2 <- c("hpev",  "mpneu", "EVD68", "hpiv1", "hpiv2",  "hpiv3", "hpiv4", 
             "cor229", "cor63", "corhku1", "cor43", "flub_univ_result",
             "flua_univ_result", "hbov", "hmpvab", "hadv", "rv_ev2", "rsv_result")

vir_labs2 <- c("Parechovirus", "M. pneumoniae", "EV-D68", "HPIV1", "HPIV2", "HPIV3", "HPIV4", 
              "Cor229", "Cor63", "CorHKU1", "Cor43", "Influenza B",
              "Influenza A", "Bocavirus",
              "HMPV", "Adenovirus", "RV/EV", "RSV")

vir_labs3 <- c("Other", "Parechovirus", "M. pneumoniae", "HPIV1", "HPIV2", "HPIV3", "HPIV4", 
              "Cor229", "Cor63", "CorHKU1", "Cor43", "Influenza B",
              "Influenza A", "Bocavirus",
              "Metapneumovirus", "Adenovirus", "RV/EV", "RSV")

newnames <- c("FID",all_of(vir_labs))

country_id <- c("Albania","Jordan","Nicaragua","Philippines")

plotnames <- c("country","RSV","Rhinovirus","Adenovirus","HMPV","Bocavirus","Enterovirus",
               "Influenza A","HPIV3","Cor43","Influenza B","CorHKU1","Cor63","HPIV4",
               "HPIV2","Cor229","M. pneumoniae","HPIV1","Parechovirus")

# cols4 <- c("Albania"="#B2182B", "Jordan"="#F4A582", "Nicaragua"="#92C5DE", "Philippines"= "#2166AC")
cols4 <- c("Albania"="#001260", "Jordan"="#70A7C3", "Nicaragua"="#C7AD78", "Philippines"= "#601200")


### Harmonize data before merge----
harm_ftd_h  <- function(dat){
  
  dat%>%
    dplyr::mutate(country_id1 = case_when(
      country_id==1 ~ "Nicaragua",
      country_id==2 ~ "Philippines",
      country_id==3 ~ "Jordan",
      country_id==4 ~ "Albania",
      T~"unk"),
      p_sex1 = case_when(
        p_sex == 0 ~ "Male",
        p_sex == 1 ~ "Female",
        T ~ "unk"
      )
    )%>%
    dplyr::select(-c(country_id,p_sex))%>%
    dplyr::rename(country_id = country_id1,
                  p_sex = p_sex1)
}

harm_ftd_c <- function(dat){
  dat%>%
    dplyr::mutate(study_year1 = case_when(
      study_year=="Year 1" ~1,
      study_year=="Year 2" ~2,
      T~99
    ),
    spec_type = case_when(
      e_specimen_resp_type=="Both nasal and OP"~1,
      e_specimen_resp_type=="Nasal only"~2,
      e_specimen_resp_type=="OP only"~3,
      e_specimen_resp_type=="Aspirate only"~4,
      T~99
    ),
    flua_h1n1swl1 = case_when(
      flua_h1n1swl=="Negative"|flua_h1n1swl=="Indeterminate"~0,
      flua_h1n1swl=="Positive"~1,
      flua_h1n1swl==""~NA_real_,
      T~99
    ),
    e_mat_delivery2 = case_when(
      e_mat_delivery == "Unknown" ~88,
      e_mat_delivery == "By pushing through the birth canal (also called vaginal delivery)"~1,
      e_mat_delivery == "By a caesarean section (also called a C-section)"~2
    ),
    e_mat_prem2 = case_when(
      e_mat_prem == "Yes"~1,
      e_mat_prem == "No" ~0
    ),
    e_mat_icu2 = case_when(
      e_mat_icu == "No"~0,
      e_mat_icu == "Yes" ~1,
      e_mat_icu == "Unknown"~88
    ),
    e_mat_vent2 = case_when(
      e_mat_vent == "No"~0,
      e_mat_vent == "Yes" ~1,
      e_mat_vent == "Unknown"~88,
      e_mat_vent == "Refused"~99
    ),
    e_smoke2 = case_when(
      e_smoke == "Never" ~ 0,
      e_smoke == "Daily" ~ 1,
      e_smoke == "Weekly" ~ 2,
      e_smoke == "Monthly" ~ 3,
      e_smoke == "Less than monthly" ~ 4,
      e_smoke == "Unknown" ~ 88,
      e_smoke == "Refused" ~ 99
    ),
    e_mother_edu2 = case_when(
      e_mother_edu== "No formal schooling" ~ 1,
      e_mother_edu== "Some primary school" ~ 2,
      e_mother_edu== "Completed primary school" ~ 3,
      e_mother_edu== "Some secondary/high school" ~ 4,
      e_mother_edu== "Completed secondary/high school" ~ 5,
      e_mother_edu== "Technical or associate degree" ~ 6,
      e_mother_edu== "Some college or university" ~ 7,
      e_mother_edu== "Completed college or university" ~ 8,
      e_mother_edu== "Master's degree or other higher degree" ~ 9,
      e_mother_edu== "Not applicable" ~ 77,
      e_mother_edu== "Unknown" ~ 88,
      e_mother_edu== "Refused" ~ 99,
    ),
    e_father_edu2 = case_when(
      e_father_edu== "No formal schooling" ~ 1,
      e_father_edu== "Some primary school" ~ 2,
      e_father_edu== "Completed primary school" ~ 3,
      e_father_edu== "Some secondary/high school" ~ 4,
      e_father_edu== "Completed secondary/high school" ~ 5,
      e_father_edu== "Technical or associate degree" ~ 6,
      e_father_edu== "Some college or university" ~ 7,
      e_father_edu== "Completed college or university" ~ 8,
      e_father_edu== "Master's degree or other higher degree" ~ 9,
      e_father_edu== "Not applicable" ~ 77,
      e_father_edu== "Unknown" ~ 88,
      e_father_edu== "Refused" ~ 99,
    ),
    rv2 = case_when(
      rv==0 ~0,
      rv==1 & ev==0~1,
      rv==1 & ev==1~0))%>%
    dplyr::select(-c(study_year,e_specimen_resp_type,flua_h1n1swl, e_mat_delivery, 
                     e_smoke, e_mat_prem, e_mat_icu, e_mat_vent, e_mother_edu, 
                     e_father_edu))%>%
    dplyr::rename(study_year = study_year1,
                  e_specimen_resp_type = spec_type,
                  flua_h1n1swl = flua_h1n1swl1,
                  e_mat_delivery = e_mat_delivery2,
                  e_mat_prem = e_mat_prem2,
                  e_mat_icu = e_mat_icu2,
                  e_mat_vent = e_mat_vent2,
                  e_smoke = e_smoke2,
                  e_mother_edu = e_mother_edu2,
                  e_father_edu = e_father_edu2)%>%
    dplyr::filter(is.na(rv)==F)
}



### Data munge functions----
data_merge <- function(){
  ftd_h_sub <- ftd_h%>%
    select(site, e_specimen_resp_date, p_sex, age_cat, infect_d, all_of(vir_ord), num_in, hosp, any_chronic)
  
  ftd_c_sub <- ftd_c%>%
    select(site, e_specimen_resp_date, p_sex, age_cat, infect_d, all_of(vir_ord), num_in, hosp, any_chronic)
  
  bind_rows(ftd_h_sub, ftd_c_sub)
}

subset_merge <- function(){
  ftd_h_sub <- ftd_h%>%
    select(country_id, e_specimen_resp_date, all_of(vir_ord))
  
  ftd_c_sub <- ftd_c%>%
    select(country_id, e_specimen_resp_date, all_of(vir_ord))
  
  bind_rows(ftd_h_sub, ftd_c_sub)
}

melt_clean <- function(dat, ...){
  
  vars <- c(...)
  
  dat%>%
    # mutate(study_m = as.numeric(floor(interval(as.Date("2015-06-21"), 
    #                                            e_specimen_resp_date)/months(1)))+1)%>%
    group_by(across(all_of(vars)))%>%
    summarise(across(.cols = all_of(vir_ord),sum))%>%
    melt(id.vars = across(all_of(vars)))%>%
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
}

add_dates <- function(){
  m_dates <- data.frame("date" = seq(as.Date(min(ftd_tot_sub$e_specimen_resp_date),
                                             format="%Y-%m-%d"), 
                                     as.Date(max(ftd_tot_sub$e_specimen_resp_date),
                                             format="%Y-%m-%d"), 
                                     by="months"))
  
  m_dates <- m_dates%>%
    mutate(study_m = row_number())
  
  left_join(ftd_tot_long, m_dates, by="study_m")
}


### Table 1 functions----
table_vars <- function(dat){

  dat_nm <- deparse(substitute(dat))
  
  dat%>%
    ungroup()%>%
    mutate(birthweight = case_when(
      is.na(e_mat_birthweight)==F & is.na(f_birthweight)==T ~ e_mat_birthweight,
      is.na(e_mat_birthweight)==F & is.na(f_birthweight)==F ~ f_birthweight,
      is.na(e_mat_birthweight)==T & is.na(f_birthweight)==F ~ f_birthweight,
      T ~ NA_real_
    ),
    low_bw = case_when(
      is.na(birthweight)==F & birthweight < 2.5 ~ 1,
      is.na(birthweight)==F & birthweight >= 2.5 ~ 0,
      is.na(birthweight)==T ~ NA_real_
    ),
    male = case_when(
      p_sex=="Male"~1,
      T~0
    ),
    parent_sec_edu = case_when(
      e_mother_edu %in% c(5:9) | e_father_edu %in% c(5:9)~1,
      e_mother_edu %in% c(1:4) & e_father_edu %in% c(1:4)~0,
      e_mother_edu %in% c(77,88,99)|e_father_edu %in% c(77,88,99)~NA_real_
    ),
    house_smoke = factor(case_when(
      e_smoke==0 ~"Never",
      e_smoke==1 ~"Daily",
      e_smoke==2 ~"Weekly",
      e_smoke==3 ~"Monthly",
      e_smoke==4 ~"< Monthly",
      T~NA_character_
    ), ordered=F),
    cohort = case_when(
      dat_nm=="ftd_h" ~"Hospitalized",
      T~"Non-ill"
    ))%>%
    select(country_id, cohort, p_calc_age, male, parent_sec_edu, low_bw, house_smoke)
}

table_vars_s <- function(dat, gvar){
  dat_nm <- deparse(substitute(dat))
  gvar <- deparse(substitute(gvar))
  
  suppressMessages(dat%>%
    ungroup()%>%
    mutate( EVD68 = case_when(
      is.na(EVD68)==T ~0,
      T ~ EVD68
    ),
    rv_ev2 = case_when(
      rv_ev==1 & EVD68 != 1 ~ 1,
      rv_ev==1 & EVD68==1 ~0,
      rv_ev==0 ~0
    )
   )%>%
    select(all_of(gvar), all_of(vir_ord2))%>%
    rename_with(~  vir_labs2, .cols = which(colnames(.)%in%vir_ord2))%>%
    # rename(virus = virus2)%>%
    mutate(cohort = case_when(
      dat_nm=="ftd_h" ~"Hospitalized",
      T~"Non-ill"
    )))
}

prev_tables <- function(dat, gvar){
  gvar <- deparse(substitute(gvar))
  
  dat%>%
  gtsummary::tbl_strata(
    strata = gvar,
    .tbl_fun =
      ~.x %>%
      gtsummary::tbl_summary(
        by = cohort,
        statistic = list(gtsummary::all_continuous() ~ "{n} ({p}%)"),
        digits = gtsummary::all_continuous() ~ 2
      ))%>%
    gtsummary::as_flex_table()
}

### table of mean ETI and 95% Credible interval by site
table_site_eti <- function(site_nm){
  
  if(site_nm=="Albania"){
    site_eti%>%
      filter(site==site_nm)%>%
      mutate(tot_eti = paste0(eti_prop, " ","(",lcl_prop, ",", ucl_prop,")"))%>%
      select(cause, tot_eti)
  }
  else{
    site_eti%>%
      filter(site==site_nm)%>%
      mutate(tot_eti = paste0(eti_prop, " ","(",lcl_prop, ",", ucl_prop,")"))%>%
      select(tot_eti)
  }
}

### Single vs. co-infection for all viruses----

### Prevalence of viral infections----
prev_func <- function(dat, ...){
  df <- dat
  
  vars <- c(...)
  
  df%>%
    group_by(across(all_of(vars)))%>%
    summarise(tot_n = n(),
              across(.cols = all_of(vir_ord),sum))%>%
    ungroup()%>%
    rename_with(~  vir_labs, .cols = which(colnames(.)%in%vir_ord))%>%
    mutate(across(.cols = all_of(vir_labs), ~ round(.x/tot_n*100,2)))%>%
    select(-tot_n)%>%
    pivot_longer(!all_of(vars), names_to = "virus", values_to = "prevalence")%>%
    mutate(virus2 = factor(virus, levels = vir_labs))
    
}

### Plotting prevalence with facets for strata
prev_bar <- function(dat, facets){
  df <- dat
  if(missing(facets)) {
    ggplot(df, aes(fill=country_id, y=prevalence, x=virus2)) + 
      geom_bar(position="dodge", stat="identity")+
      scale_y_continuous(limits=c(0, 45), breaks=seq(0, 40, 10))+
      coord_flip()+
      labs(y = "Prevalence", x = " ", size=1)+
      scale_fill_manual(values = cols4, breaks = country_id, labels = country_id, name=NULL)+
      theme_classic()+
      theme(axis.text = element_text(size = 20),
            axis.title.x = element_text(size = 24),
            legend.text = element_text(size=20))
  } else {
    ggplot(df, aes(fill=country_id, y=prevalence, x=virus2)) + 
      geom_bar(position="dodge", stat="identity")+
      scale_y_continuous(limits=c(0, 45), breaks=seq(0, 40, 10))+
      coord_flip()+
      labs(y = "Prevalence", x = " ", size=1)+
      scale_fill_manual(values = cols4, breaks = country_id, labels = country_id, name=NULL)+
      theme_classic()+
      theme(axis.text = element_text(size = 20),
            axis.title.x = element_text(size = 24),
            legend.text = element_text(size=20))+
      facet_wrap(as.formula(sprintf('~%s',facets)))
  }
  
}



##----------------------------
### Function to plot stratified by a discrete covariate using ggplot 
## Requres the ggplot, dplyr, and reshape2 packages to be loaded 
## The output is a list of plots, one for each stratum 
## The last plot in the list is of the marginal posterior etiology probabiltiies computed across all strata 
## The plot results are similar to plot_etiology_strat but users have more flexibility to customize with ggplot functions (e.g. plot themes)
##----------------------------

gg_plot_etiology_strat <- function(DIR_NPLCM,strata_weights=NULL,truth=NULL, category){
  
  # Check if result folder filepath is correctly specified 
  if (!is_jags_folder(DIR_NPLCM)){
    stop("==[baker] oops, not a folder baker recognizes. 
         Try a folder generated by baker, e.g., a temporary folder?==")
  }
  
  # Read in the nplcm data and model options 
  data_nplcm <- dget(file.path(DIR_NPLCM,"data_nplcm.txt"))  
  model_options <- dget(file.path(DIR_NPLCM,"model_options.txt"))  
  
  if (is.null(data_nplcm)|is.null(model_options)){
    stop("==[baker] oops, nplcm data and model options settings not found in the result folder!=")
  }
  
  ## get the list of causes 
  cause_list = model_options$likelihood$cause_list
  
  ## read in the posterior samples 
  new_env <- new.env()
  source(file.path(DIR_NPLCM,"jagsdata.txt"),local=new_env)
  bugs.dat <- as.list(new_env)
  rm(new_env)
  
  res_nplcm1 <- coda::read.coda(file.path(DIR_NPLCM,"CODAchain1.txt"),
                                file.path(DIR_NPLCM,"CODAindex.txt"),
                                quiet=TRUE)
  res_nplcm2 <-  coda::read.coda(file.path(DIR_NPLCM,"CODAchain2.txt"),
                                 file.path(DIR_NPLCM,"CODAindex.txt"),
                                 quiet=TRUE)
  
  res_nplcm <- mcmc.list(res_nplcm1, res_nplcm2)
  res_nplcm <- as.matrix(res_nplcm)
  
  print_res <- function(x) plot(res_nplcm[,grep(x,colnames(res_nplcm))])
  get_res   <- function(x) res_nplcm[,grep(x,colnames(res_nplcm))]
  
  # structure the posterior samples:
  JBrS_1        <- ncol(bugs.dat$MBS_1) # number of pathogens in MBS1
  n_samp_kept   <- nrow(res_nplcm)      # number of posterior sample after burn-in
  Jcause        <- bugs.dat$Jcause      # number of causes
  Nd            <- bugs.dat$Nd          # case size
  Nu            <- bugs.dat$Nu          # control size
  
  likelihood <- model_options$likelihood
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  
  # generate design matrix for etiology regression:
  Eti_formula <- likelihood$Eti_formula
  Z_Eti <- stats::model.matrix(Eti_formula, data.frame(X, Y)[Y == 
                                                               1, , drop = FALSE])
  ncol_dm_Eti <- ncol(Z_Eti)
  Eti_colname_design_mat <- attributes(Z_Eti)$dimnames[[2]]
  attributes(Z_Eti)[names(attributes(Z_Eti)) != "dim"] <- NULL
  unique_Eti_level <- unique(Z_Eti)
  n_unique_Eti_level <- nrow(unique_Eti_level)
  
  
  ## Check if model has been fit with a discrete covariate
  if (is.null(n_unique_Eti_level)){
    stop("==[baker] Number of stratums is null! ? ==")
  }
  
  # etiology:
  plotid_Eti <- which(data_nplcm$Y==1) #want the etiology for cases only 
  #  format eti probabilities into an array by stratum and cause
  Eti_prob_scale <- array(get_res("pEti"),c(n_samp_kept,n_unique_Eti_level,Jcause))
  
  # weight to marginalize posterior etiology distributions
  user_weight <- rep(1/n_unique_Eti_level,n_unique_Eti_level)
  if(!is.null(strata_weights) && length(strata_weights==n_unique_Eti_level)){
    user_weight <- strata_weights
  }
  
  # marginalized posterior etiology over all sites using user-defined weights
  Eti_overall_usr_weight <- apply(Eti_prob_scale,1,function(S) t(S)%*%matrix(user_weight,ncol=1))
  
  
  plot_list <- list()
  for(site in 1:n_unique_Eti_level){
    
    # shape probabilities array into dataset for ggplot 
    etiData = data.frame(Eti_prob_scale[,site,])
    names(etiData) = cause_list
    plotData = melt(etiData,id.vars = NULL)
    names(plotData) = c("cause", "prob")
    
    #plotData <- plotData[plotData$cause=="RSV",]
    # compute posterior means 
    summaryData <- plotData %>%
      dplyr::group_by(cause) %>%
      dplyr::summarise(eti_mean = mean(prob), 
                       ci_025 = quantile(prob, c(0.025)), 
                       ci_975 = quantile(prob, c(0.975)))
    
    ## plot histograms of the posterior probabilities for each stratum 
    plot_list[[site]] <- ggplot(plotData, aes(x=prob)) +
      geom_histogram(fill="#ffc04d") + 
      geom_vline(data=summaryData, 
                 mapping = aes(xintercept=eti_mean), colour="#005b96") +
      geom_vline(data=summaryData, 
                 mapping = aes(xintercept=ci_025), colour="green", linetype="dashed") +
      geom_vline(data=summaryData, 
                 mapping = aes(xintercept=ci_975), colour="green", linetype="dashed") +
      facet_wrap(~cause) +
      labs(title=paste0("Plot of posterior etiology probabilities for stratum: ", site),
           subtitle = "Mean posterior probability displayed as solid line \n 95% CrIs displayed as dashed lines",
           x = "Etiology", y ="Samples")
  }
  
  # shape the overall marginal etiology probabilities as a data frame 
  overallEti = data.frame(Eti_overall_usr_weight)
  overallEti$cause = cause_list
  overallEti = melt(overallEti, id.vars = "cause", value.name = "prob")
  
  # compute the posterior means 
  summaryEti <- overallEti %>%
    dplyr::group_by(cause) %>%
    dplyr::summarise(overall_mean = mean(prob),
                     ci_025 = quantile(prob, c(0.025)), 
                     ci_975 = quantile(prob, c(0.975)))
  
  ## add the plot of overall marginal eti probabiltiies (across sites)
  plot_list[[n_unique_Eti_level+1]] <-  ggplot(overallEti, aes(x=prob)) +
    geom_histogram(fill="#ffc04d") +
    geom_vline(data=summaryEti, 
               mapping = aes(xintercept=overall_mean), colour="#005b96") +
    geom_vline(data=summaryEti, 
               mapping = aes(xintercept=ci_025), colour="green", linetype="dashed") +
    geom_vline(data=summaryEti, 
               mapping = aes(xintercept=ci_975), colour="green", linetype="dashed") +
    facet_wrap(~cause) +
    facet_wrap(~cause)+
    labs(title=paste("Marginal posterior etiology probabilities (across all sites)",category),
         subtitle= "Overall marginal posterior mean displayed as solid line",
         x = "Etiology", y ="Samples")
  
  return(plot_list)
  
}

#### Functions to extract marginal and site-specific posterior means by pathogen
extract_post_means <- function(DIR_NPLCM,strata_weights=NULL,truth=NULL, category){
  # Check if result folder filepath is correctly specified 
  if (!is_jags_folder(DIR_NPLCM)){
    stop("==[baker] oops, not a folder baker recognizes. 
         Try a folder generated by baker, e.g., a temporary folder?==")
  }
  
  # Read in the nplcm data and model options 
  data_nplcm <- dget(file.path(DIR_NPLCM,"data_nplcm.txt"))  
  model_options <- dget(file.path(DIR_NPLCM,"model_options.txt"))  
  
  if (is.null(data_nplcm)|is.null(model_options)){
    stop("==[baker] oops, nplcm data and model options settings not found in the result folder!=")
  }
  
  ## get the list of causes 
  cause_list = model_options$likelihood$cause_list
  
  ## read in the posterior samples 
  new_env <- new.env()
  source(file.path(DIR_NPLCM,"jagsdata.txt"),local=new_env)
  bugs.dat <- as.list(new_env)
  rm(new_env)
  
  res_nplcm1 <- coda::read.coda(file.path(DIR_NPLCM,"CODAchain1.txt"),
                                file.path(DIR_NPLCM,"CODAindex.txt"),
                                quiet=TRUE)
  res_nplcm2 <-  coda::read.coda(file.path(DIR_NPLCM,"CODAchain2.txt"),
                                 file.path(DIR_NPLCM,"CODAindex.txt"),
                                 quiet=TRUE)
  
  res_nplcm <- mcmc.list(res_nplcm1, res_nplcm2)
  res_nplcm <- as.matrix(res_nplcm)
  
  print_res <- function(x) plot(res_nplcm[,grep(x,colnames(res_nplcm))])
  get_res   <- function(x) res_nplcm[,grep(x,colnames(res_nplcm))]
  
  # structure the posterior samples:
  JBrS_1        <- ncol(bugs.dat$MBS_1) # number of pathogens in MBS1
  n_samp_kept   <- nrow(res_nplcm)      # number of posterior sample after burn-in
  Jcause        <- bugs.dat$Jcause      # number of causes
  Nd            <- bugs.dat$Nd          # case size
  Nu            <- bugs.dat$Nu          # control size
  
  likelihood <- model_options$likelihood
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  
  # generate design matrix for etiology regression:
  Eti_formula <- likelihood$Eti_formula
  Z_Eti <- stats::model.matrix(Eti_formula, data.frame(X, Y)[Y == 
                                                               1, , drop = FALSE])
  ncol_dm_Eti <- ncol(Z_Eti)
  Eti_colname_design_mat <- attributes(Z_Eti)$dimnames[[2]]
  attributes(Z_Eti)[names(attributes(Z_Eti)) != "dim"] <- NULL
  unique_Eti_level <- unique(Z_Eti)
  n_unique_Eti_level <- nrow(unique_Eti_level)
  
  
  ## Check if model has been fit with a discrete covariate
  if (is.null(n_unique_Eti_level)){
    stop("==[baker] Number of stratums is null! ? ==")
  }
  
  # etiology:
  plotid_Eti <- which(data_nplcm$Y==1) #want the etiology for cases only 
  #  format eti probabilities into an array by stratum and cause
  Eti_prob_scale <- array(get_res("pEti"),c(n_samp_kept,n_unique_Eti_level,Jcause))
  
  # weight to marginalize posterior etiology distributions
  user_weight <- rep(1/n_unique_Eti_level,n_unique_Eti_level)
  if(!is.null(strata_weights) && length(strata_weights==n_unique_Eti_level)){
    user_weight <- strata_weights
  }
  
  # marginalized posterior etiology over all sites using user-defined weights
  Eti_overall_usr_weight <- apply(Eti_prob_scale,1,function(S) t(S)%*%matrix(user_weight,ncol=1))
  
  
  site_list <<- list()
  for(site in 1:n_unique_Eti_level){
    
    # shape probabilities array into dataset for ggplot 
    etiData = data.frame(Eti_prob_scale[,site,])
    names(etiData) = cause_list
    plotData = melt(etiData,id.vars = NULL)
    names(plotData) = c("cause", "prob")
    
    #plotData <- plotData[plotData$cause=="RSV",]
    # compute posterior means 
    site_list[[site]] <<- plotData %>%
      dplyr::group_by(cause) %>%
      dplyr::summarise(eti_mean = mean(prob), 
                       ci_025 = quantile(prob, c(0.025)), 
                       ci_975 = quantile(prob, c(0.975)))
   
  }}

#### Functions to extract marginal and site-specific posterior means by pathogen
extract_post_margmeans <- function(DIR_NPLCM,strata_weights=NULL,truth=NULL, category){
  # Check if result folder filepath is correctly specified 
  if (!is_jags_folder(DIR_NPLCM)){
    stop("==[baker] oops, not a folder baker recognizes. 
         Try a folder generated by baker, e.g., a temporary folder?==")
  }
  
  # Read in the nplcm data and model options 
  data_nplcm <- dget(file.path(DIR_NPLCM,"data_nplcm.txt"))  
  model_options <- dget(file.path(DIR_NPLCM,"model_options.txt"))  
  
  if (is.null(data_nplcm)|is.null(model_options)){
    stop("==[baker] oops, nplcm data and model options settings not found in the result folder!=")
  }
  
  ## get the list of causes 
  cause_list = model_options$likelihood$cause_list
  
  ## read in the posterior samples 
  new_env <- new.env()
  source(file.path(DIR_NPLCM,"jagsdata.txt"),local=new_env)
  bugs.dat <- as.list(new_env)
  rm(new_env)
  
  res_nplcm1 <- coda::read.coda(file.path(DIR_NPLCM,"CODAchain1.txt"),
                                file.path(DIR_NPLCM,"CODAindex.txt"),
                                quiet=TRUE)
  res_nplcm2 <-  coda::read.coda(file.path(DIR_NPLCM,"CODAchain2.txt"),
                                 file.path(DIR_NPLCM,"CODAindex.txt"),
                                 quiet=TRUE)
  
  res_nplcm <- mcmc.list(res_nplcm1, res_nplcm2)
  res_nplcm <- as.matrix(res_nplcm)
  
  print_res <- function(x) plot(res_nplcm[,grep(x,colnames(res_nplcm))])
  get_res   <- function(x) res_nplcm[,grep(x,colnames(res_nplcm))]
  
  # structure the posterior samples:
  JBrS_1        <- ncol(bugs.dat$MBS_1) # number of pathogens in MBS1
  n_samp_kept   <- nrow(res_nplcm)      # number of posterior sample after burn-in
  Jcause        <- bugs.dat$Jcause      # number of causes
  Nd            <- bugs.dat$Nd          # case size
  Nu            <- bugs.dat$Nu          # control size
  
  likelihood <- model_options$likelihood
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  
  # generate design matrix for etiology regression:
  Eti_formula <- likelihood$Eti_formula
  Z_Eti <- stats::model.matrix(Eti_formula, data.frame(X, Y)[Y == 
                                                               1, , drop = FALSE])
  ncol_dm_Eti <- ncol(Z_Eti)
  Eti_colname_design_mat <- attributes(Z_Eti)$dimnames[[2]]
  attributes(Z_Eti)[names(attributes(Z_Eti)) != "dim"] <- NULL
  unique_Eti_level <- unique(Z_Eti)
  n_unique_Eti_level <- nrow(unique_Eti_level)
  
  
  ## Check if model has been fit with a discrete covariate
  if (is.null(n_unique_Eti_level)){
    stop("==[baker] Number of stratums is null! ? ==")
  }
  
  # etiology:
  plotid_Eti <- which(data_nplcm$Y==1) #want the etiology for cases only 
  #  format eti probabilities into an array by stratum and cause
  Eti_prob_scale <- array(get_res("pEti"),c(n_samp_kept,n_unique_Eti_level,Jcause))
  
  # weight to marginalize posterior etiology distributions
  user_weight <- rep(1/n_unique_Eti_level,n_unique_Eti_level)
  if(!is.null(strata_weights) && length(strata_weights==n_unique_Eti_level)){
    user_weight <- strata_weights
  }
  
  # marginalized posterior etiology over all sites using user-defined weights
  Eti_overall_usr_weight <- apply(Eti_prob_scale,1,function(S) t(S)%*%matrix(user_weight,ncol=1))
  
  
  
  
    
    # shape the overall marginal etiology probabilities as a data frame 
    overallEti = data.frame(Eti_overall_usr_weight)
    overallEti$cause = cause_list
    overallEti = melt(overallEti, id.vars = "cause", value.name = "prob")
    
    # compute the posterior means 
    summaryEti <- overallEti %>%
      dplyr::group_by(cause) %>%
      dplyr::summarise(overall_mean = mean(prob),
                       ci_025 = quantile(prob, c(0.025)), 
                       ci_975 = quantile(prob, c(0.975)))
    
  }

############### old functions ----
### Coronavirus plots----

# cor_bplot <- function(dat, y_val, ax_lab){  
#   dat%>% 
#     ggplot(aes(x=variable, y={{y_val}}, fill=infect)) +
#     geom_bar(aes(fill=as.factor(infect)), stat = "identity", position = "dodge")+ 
#     labs(x = " ", y = ax_lab, size = 1)+ 
#     scale_fill_manual(values = cols2cor, name = NULL)+ 
#     theme_bw()+ 
#     theme(axis.text.y = element_text(size = 20), 
#           axis.title.y = element_text(size = 24), 
#           axis.text.x = element_text(size = 20), 
#           legend.text=element_text(size=20))
# }
# Upset plots of coronavirus co-infections -- need to fix the renaming...
# cor_upmunge <- function(dat){
#   
#   df <- dat %>%
#     as_tibble() %>%
#     mutate(fid = row_number())%>%
#     select(fid, all_of(vir_ord))%>%
#     rename_with(~ .data[[newnames]], .cols = colnames(.))
#   
#   
#   df%>%
#     gather(virus, Member, -FID) %>%
#     filter(Member==1) %>%
#     select(- Member)%>%
#     group_by(FID)%>%
#     summarise(virus = list(virus))
#   
# }

# cor_upset <- function(dat){
#   dat%>%
#     ggplot(aes(x=virus)) +
#     geom_bar(fill="#B2182B") +
#     ylab("Number")+
#     xlab(" ")+
#     scale_x_upset(n_intersections = 11)+
#     theme_classic()+
#     theme(axis.title.y = element_text(size = 20),
#           axis.text.y = element_text(size = 20))+
#     theme_combmatrix(combmatrix.panel.point.color.fill = "#B2182B",
#                      combmatrix.panel.line.size = 0,
#                      combmatrix.label.text = element_text(size=16))
# }


# prev_function <- function(p_dat, p_var){
#   p_dat%>%
#     group_by(country_id,p_sex)%>%
#     summarise(cases = sum({{ p_var }}),
#               tot_n = n())%>%
#     mutate("prev_{{p_var}}" := round((cases/tot_n)*100, digits = 2))%>%
#     dplyr::select(last_col())
# }
## not currently working
# prev_function <- function(p_dat,g_vars, p_var){
#   p_dat%>%
#     dplyr::group_by({{g_vars}})%>%
#     dplyr::summarise("cases_{{p_var}}" = sum({{ p_var }}),
#               "tot_{{p_var}}" = n())%>%
#     dplyr::mutate("prev_{{p_var}}" := round(("cases_{{p_var}}"/"tot_{{p_var}}")*100, digits = 2))
#     # dplyr::select(last_col())
# }
# prev_function(ftd_h, c("country_id", "p_sex"), vir_ord)
# g_vars <- c("country_id", "p_sex")

# cor <- c("cor63", "cor229", "cor43", "corhku1")

# cols2cor <- c("Co-infection"="#B2182B", "Single infection"= "#2166AC")

# iris_full <- read_csv("Box/IRIS_0620/finaldataset.csv")
# import_dat <- function(data_out, path, var_imp){
#   data_out <-read_dta(path)%>%
#     zap_label()%>%
#     dplyr::select(all_of({{var_imp}}))
# }
# import_dat

# ,
# e_rank_edu = case_when(
#   e_rank_edu==44 ~4,
#   e_rank_edu==55 ~5,
#   T~e_rank_edu
# )