# Run BAKER models

library(dplyr)
library(R2WinBUGS)
library(baker)
library(reshape2)
library(lubridate)
library(rjags)
library(ggplot2)


# load("data/data_iris_disc_u6071221.rda")
# load("data/data_iris_disc_o6071221.rda")
# load("data/data_iris_disc_u3081721.rda")
# load("data/data_iris_disc_o3081721.rda")
# load("data/data_iris_disc_u3a081721.rda")
# load("data/data_iris_disc_35a081721.rda")
# load("data/data_iris_disc_o6a081721.rda")
# load("data/data_iris_disc_multd072121.rda")
# load("data/tot_lca_sub_alri_mult3_110321.rda")
# load("data/tot_lca_sub_alri_mult3a_120221.rda")
load("data/tot_lca_sub_alri_mult3_121621.rda")

source("script/ftd_functions.R")

cause_list_iris <- c("Influenza A", "Influenza B", "RSV", "Cor63",
                     "Cor229", "Cor43", "CorHKU1", "HPIV2", "HPIV3", "HPIV4", "HPIV1",
                     "Metapneumovirus", "Bocavirus", "M. pneumoniae", "Adenovirus", "RV_EV", 
                     "Parechovirus", "other")
K <- 1 
N.CAUSE <- 18

tot_lca_sub_alri_mult3fix$X <- as.data.frame(tot_lca_sub_alri_mult3fix$X)
# tot_lca_sub_alri_mult3a$X <- as.data.frame(tot_lca_sub_alri_mult3a$X)
# tot_lca_sub_alri_mult3$X <- as.data.frame(tot_lca_sub_alri_mult3$X)
# data_iris_disc_u3$X <- as.data.frame(data_iris_disc_u3$X) 
# data_iris_disc_o3$X <- as.data.frame(data_iris_disc_o3$X)
# data_iris_disc_u3a$X <- as.data.frame(data_iris_disc_u3a$X)
# data_iris_disc_35a$X <- as.data.frame(data_iris_disc_35a$X)
# data_iris_disc_o6a$X <- as.data.frame(data_iris_disc_o6a$X)
# data_iris_disc_o3a$X <- as.data.frame(data_iris_disc_o3a$X)
# data_iris_disc_u6$X <- as.data.frame(data_iris_disc_u6$X) 
# data_iris_disc_o6$X <- as.data.frame(data_iris_disc_o6$X) 
# data_iris_disc_multd$X <- as.data.frame(data_iris_disc_multd$X) 
# N.SITE <- length(unique(tot_lca_sub_alri_mult3$X$SITE)) 
# N.SITE <- length(unique(tot_lca_sub_alri_mult3a$X$SITE)) 
N.SITE <- length(unique(tot_lca_sub_alri_mult3fix$X$SITE)) 

# Under 3m--pneumonia----
table(data_iris_disc_u3$Y)
Nd <- 517
Nu <- 429


BrS_object_1 <- make_meas_object(
  patho = cause_list_iris, 
  specimen = "MBS", test = "1",
  quality = "BrS", 
  cause_list = cause_list_iris)

# also save the data clean options:

# parent directory for testing code (LOCAL):
data_dir    <- "Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd/data"
working_dir <- "/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd"

clean_options <- list(
  BrS_objects = make_list(BrS_object_1),         # all bronze-standard measurements.
  # SS_objects  = make_list(SS_object_1),  ## SS measurements 
  X_extra = c("SITE"), ##ZW add this to record the covariates considered for analysis; here only "SITE" for now.
  patho_taxo_dir = file.path(data_dir, "pathogen_category_IRIS.csv"), ##ZW would change the name to "pathogen_category_IRIS.csv" moving forward.
  allow_missing = FALSE)

m_opt1 <- list(likelihood   = list(cause_list = cause_list_iris,                # <---- fitted causes.
                                   k_subclass = c(1),                      # <---- no. of subclasses.
                                   Eti_formula = ~ -1+as.factor(SITE),                    # <---- etiology regression formula; only for cases.
                                   FPR_formula = list(
                                     MBS1 =  ~ -1+as.factor(SITE))),
               use_measurements = c("BrS"),
               prior = list(Eti_prior   = t(sapply(1:N.SITE, function(i) overall_uniform(1, cause_list_iris))), 
                            TPR_prior   = list(
                              BrS  = list(info  = "informative",
                                          input = "match_range",
                                          val   = list(
                                            MBS1 = list(up = list(rep(0.99,length(BrS_object_1$patho))),
                                                        low = list(c(rep(0.5,length(BrS_object_1$patho)-2),0.001,0.5))
                                            ))))))     

model_options_discrete_reg <- m_opt1

assign_model(model_options_discrete_reg, data_iris_disc_u3)

working_dir = getwd()

Date     <- gsub("-", "_", Sys.Date())

set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) 

dated_strat_name_discrete_u3   <- file.path(working_dir,paste0("results/site_u3r"))
#dir.create(dated_strat_name_discrete_u6)

result_folder_discrete_u3 <- dated_strat_name_discrete_u3

mcmc_options_discrete_u3 <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(10000), #10000
  n.burnin = as.integer(5000), #5000
  n.thin = 1,
  individual.pred = FALSE, 
  ppd = !TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_u3,
  bugsmodel.dir = result_folder_discrete_u3,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_u3$result.folder, "data_clean_options.txt"))
dput(data_iris_disc_u3,file.path(mcmc_options_discrete_u3$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
start.time <- Sys.time()
nplcm_reg_discrete_only_u3 <- nplcm(data_iris_disc_u3,model_options_discrete_reg,mcmc_options_discrete_u3)
end.time <- Sys.time()
print(end.time-start.time)

gg_plot_etiology_strat(result_folder_discrete_u3, category= "<3 months")

# Over 3m--pneumonia----
table(data_iris_disc_o3$Y)
Nd <- 631
Nu <- 639


BrS_object_1 <- make_meas_object(
  patho = cause_list_iris, 
  specimen = "MBS", test = "1",
  quality = "BrS", 
  cause_list = cause_list_iris)

# also save the data clean options:

# parent directory for testing code (LOCAL):
data_dir    <- "Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd/data"
working_dir <- "/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd"

clean_options <- list(
  BrS_objects = make_list(BrS_object_1),         # all bronze-standard measurements.
  # SS_objects  = make_list(SS_object_1),  ## SS measurements 
  X_extra = c("SITE"), ##ZW add this to record the covariates considered for analysis; here only "SITE" for now.
  patho_taxo_dir = file.path(data_dir, "pathogen_category_IRIS.csv"), ##ZW would change the name to "pathogen_category_IRIS.csv" moving forward.
  allow_missing = FALSE)

m_opt1 <- list(likelihood   = list(cause_list = cause_list_iris,                # <---- fitted causes.
                                   k_subclass = c(1),                      # <---- no. of subclasses.
                                   Eti_formula = ~ -1+as.factor(SITE),                    # <---- etiology regression formula; only for cases.
                                   FPR_formula = list(
                                     MBS1 =  ~ -1+as.factor(SITE))),
               use_measurements = c("BrS"),
               prior = list(Eti_prior   = t(sapply(1:N.SITE, function(i) overall_uniform(1, cause_list_iris))), 
                            TPR_prior   = list(
                              BrS  = list(info  = "informative",
                                          input = "match_range",
                                          val   = list(
                                            MBS1 = list(up = list(rep(0.99,length(BrS_object_1$patho))),
                                                        low = list(c(rep(0.5,length(BrS_object_1$patho)-2),0.001,0.5))
                                            ))))))     

model_options_discrete_reg <- m_opt1

assign_model(model_options_discrete_reg, data_iris_disc_o3)

working_dir = getwd()

Date     <- gsub("-", "_", Sys.Date())

set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) 

dated_strat_name_discrete_o3   <- file.path(working_dir,paste0("results/site_o3r"))
#dir.create(dated_strat_name_discrete_u6)

result_folder_discrete_o3 <- dated_strat_name_discrete_o3

mcmc_options_discrete_o3 <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(10000), #10000
  n.burnin = as.integer(5000), #5000
  n.thin = 1,
  individual.pred = FALSE, 
  ppd = !TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_o3,
  bugsmodel.dir = result_folder_discrete_o3,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_o3$result.folder, "data_clean_options.txt"))
dput(data_iris_disc_o3,file.path(mcmc_options_discrete_o3$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
start.time <- Sys.time()
nplcm_reg_discrete_only_o3 <- nplcm(data_iris_disc_o3,model_options_discrete_reg,mcmc_options_discrete_o3)
end.time <- Sys.time()
print(end.time-start.time)

gg_plot_etiology_strat(result_folder_discrete_o3, category= "3+ months")






# Under 3m alri----
table(data_iris_disc_u3a$Y)
Nd <- 808
Nu <- 429


BrS_object_1 <- make_meas_object(
  patho = cause_list_iris, 
  specimen = "MBS", test = "1",
  quality = "BrS", 
  cause_list = cause_list_iris)

# also save the data clean options:

# parent directory for testing code (LOCAL):
data_dir    <- "Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd/data"
working_dir <- "/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd"

clean_options <- list(
  BrS_objects = make_list(BrS_object_1),         # all bronze-standard measurements.
  # SS_objects  = make_list(SS_object_1),  ## SS measurements 
  X_extra = c("SITE"), ##ZW add this to record the covariates considered for analysis; here only "SITE" for now.
  patho_taxo_dir = file.path(data_dir, "pathogen_category_IRIS.csv"), ##ZW would change the name to "pathogen_category_IRIS.csv" moving forward.
  allow_missing = FALSE)

m_opt1 <- list(likelihood   = list(cause_list = cause_list_iris,                # <---- fitted causes.
                                   k_subclass = c(1),                      # <---- no. of subclasses.
                                   Eti_formula = ~ -1+as.factor(SITE),                    # <---- etiology regression formula; only for cases.
                                   FPR_formula = list(
                                     MBS1 =  ~ -1+as.factor(SITE))),
               use_measurements = c("BrS"),
               prior = list(Eti_prior   = t(sapply(1:N.SITE, function(i) overall_uniform(1, cause_list_iris))), 
                            TPR_prior   = list(
                              BrS  = list(info  = "informative",
                                          input = "match_range",
                                          val   = list(
                                            MBS1 = list(up = list(rep(0.99,length(BrS_object_1$patho))),
                                                        low = list(c(rep(0.5,length(BrS_object_1$patho)-2),0.001,0.5))
                                            ))))))     

model_options_discrete_reg <- m_opt1

assign_model(model_options_discrete_reg, data_iris_disc_u3a)

working_dir = getwd()

Date     <- gsub("-", "_", Sys.Date())

set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) 

dated_strat_name_discrete_u3a   <- file.path(working_dir,paste0("results/site_u3ra"))
#dir.create(dated_strat_name_discrete_u6)

result_folder_discrete_u3a <- dated_strat_name_discrete_u3a

mcmc_options_discrete_u3a <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(10000), #10000
  n.burnin = as.integer(5000), #5000
  n.thin = 1,
  individual.pred = FALSE, 
  ppd = !TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_u3a,
  bugsmodel.dir = result_folder_discrete_u3a,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_u3a$result.folder, "data_clean_options.txt"))
dput(data_iris_disc_u3a,file.path(mcmc_options_discrete_u3a$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
start.time <- Sys.time()
nplcm_reg_discrete_only_u3a <- nplcm(data_iris_disc_u3a,model_options_discrete_reg,mcmc_options_discrete_u3a)
end.time <- Sys.time()
print(end.time-start.time)

gg_plot_etiology_strat(result_folder_discrete_u3a, category= "<3 months")
u3a <- extract_post_margmeans(result_folder_discrete_u3a)
# 3-5m alri----
table(data_iris_disc_35a$Y)
Nd <- 397
Nu <- 247


BrS_object_1 <- make_meas_object(
  patho = cause_list_iris, 
  specimen = "MBS", test = "1",
  quality = "BrS", 
  cause_list = cause_list_iris)

# also save the data clean options:

# parent directory for testing code (LOCAL):
data_dir    <- "Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd/data"
working_dir <- "/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd"

clean_options <- list(
  BrS_objects = make_list(BrS_object_1),         # all bronze-standard measurements.
  # SS_objects  = make_list(SS_object_1),  ## SS measurements 
  X_extra = c("SITE"), ##ZW add this to record the covariates considered for analysis; here only "SITE" for now.
  patho_taxo_dir = file.path(data_dir, "pathogen_category_IRIS.csv"), ##ZW would change the name to "pathogen_category_IRIS.csv" moving forward.
  allow_missing = FALSE)

m_opt1 <- list(likelihood   = list(cause_list = cause_list_iris,                # <---- fitted causes.
                                   k_subclass = c(1),                      # <---- no. of subclasses.
                                   Eti_formula = ~ -1+as.factor(SITE),                    # <---- etiology regression formula; only for cases.
                                   FPR_formula = list(
                                     MBS1 =  ~ -1+as.factor(SITE))),
               use_measurements = c("BrS"),
               prior = list(Eti_prior   = t(sapply(1:N.SITE, function(i) overall_uniform(1, cause_list_iris))), 
                            TPR_prior   = list(
                              BrS  = list(info  = "informative",
                                          input = "match_range",
                                          val   = list(
                                            MBS1 = list(up = list(rep(0.99,length(BrS_object_1$patho))),
                                                        low = list(c(rep(0.5,length(BrS_object_1$patho)-2),0.001,0.5))
                                            ))))))     

model_options_discrete_reg <- m_opt1

assign_model(model_options_discrete_reg, data_iris_disc_35a)

working_dir = getwd()

Date     <- gsub("-", "_", Sys.Date())

set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) 

dated_strat_name_discrete_35a   <- file.path(working_dir,paste0("results/site_35ra"))
#dir.create(dated_strat_name_discrete_u6)

result_folder_discrete_35a <- dated_strat_name_discrete_35a

mcmc_options_discrete_35a <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(10000), #10000
  n.burnin = as.integer(5000), #5000
  n.thin = 1,
  individual.pred = FALSE, 
  ppd = !TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_35a,
  bugsmodel.dir = result_folder_discrete_35a,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_35a$result.folder, "data_clean_options.txt"))
dput(data_iris_disc_35a,file.path(mcmc_options_discrete_35a$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
start.time <- Sys.time()
nplcm_reg_discrete_only_35a <- nplcm(data_iris_disc_35a,model_options_discrete_reg,mcmc_options_discrete_35a)
end.time <- Sys.time()
print(end.time-start.time)

# gg_plot_etiology_strat(result_folder_discrete_u3a, category= "<3 months")
gg_plot_etiology_strat(result_folder_discrete_35a, category= "3-5 months")
marg35a <- extract_post_margmeans(result_folder_discrete_35a)
# 6-11m alri----
table(data_iris_disc_o6a$Y)
Nd <- 538
Nu <- 392


BrS_object_1 <- make_meas_object(
  patho = cause_list_iris, 
  specimen = "MBS", test = "1",
  quality = "BrS", 
  cause_list = cause_list_iris)

# also save the data clean options:

# parent directory for testing code (LOCAL):
data_dir    <- "Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd/data"
working_dir <- "/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd"

clean_options <- list(
  BrS_objects = make_list(BrS_object_1),         # all bronze-standard measurements.
  # SS_objects  = make_list(SS_object_1),  ## SS measurements 
  X_extra = c("SITE"), ##ZW add this to record the covariates considered for analysis; here only "SITE" for now.
  patho_taxo_dir = file.path(data_dir, "pathogen_category_IRIS.csv"), ##ZW would change the name to "pathogen_category_IRIS.csv" moving forward.
  allow_missing = FALSE)

m_opt1 <- list(likelihood   = list(cause_list = cause_list_iris,                # <---- fitted causes.
                                   k_subclass = c(1),                      # <---- no. of subclasses.
                                   Eti_formula = ~ -1+as.factor(SITE),                    # <---- etiology regression formula; only for cases.
                                   FPR_formula = list(
                                     MBS1 =  ~ -1+as.factor(SITE))),
               use_measurements = c("BrS"),
               prior = list(Eti_prior   = t(sapply(1:N.SITE, function(i) overall_uniform(1, cause_list_iris))), 
                            TPR_prior   = list(
                              BrS  = list(info  = "informative",
                                          input = "match_range",
                                          val   = list(
                                            MBS1 = list(up = list(rep(0.99,length(BrS_object_1$patho))),
                                                        low = list(c(rep(0.5,length(BrS_object_1$patho)-2),0.001,0.5))
                                            ))))))     

model_options_discrete_reg <- m_opt1

assign_model(model_options_discrete_reg, data_iris_disc_o6a)

working_dir = getwd()

Date     <- gsub("-", "_", Sys.Date())

set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) 

dated_strat_name_discrete_o6a   <- file.path(working_dir,paste0("results/site_o6ra"))
#dir.create(dated_strat_name_discrete_u6)

result_folder_discrete_o6a <- dated_strat_name_discrete_o6a

mcmc_options_discrete_o6a <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(10000), #10000
  n.burnin = as.integer(5000), #5000
  n.thin = 1,
  individual.pred = FALSE, 
  ppd = !TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_o6a,
  bugsmodel.dir = result_folder_discrete_o6a,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_o6a$result.folder, "data_clean_options.txt"))
dput(data_iris_disc_o6a,file.path(mcmc_options_discrete_o6a$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
start.time <- Sys.time()
nplcm_reg_discrete_only_o6a <- nplcm(data_iris_disc_o6a,model_options_discrete_reg,mcmc_options_discrete_o6a)
end.time <- Sys.time()
print(end.time-start.time)

gg_plot_etiology_strat(result_folder_discrete_o6a, category= "6-11 months")
margo6a <- extract_post_margmeans(result_folder_discrete_o6a)
# Under 6m----
table(data_iris_disc_u6$Y)
Nd <- 767
Nu <- 2400


BrS_object_1 <- make_meas_object(
  patho = cause_list_iris, 
  specimen = "MBS", test = "1",
  quality = "BrS", 
  cause_list = cause_list_iris)

# also save the data clean options:

# parent directory for testing code (LOCAL):
data_dir    <- "Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd/data"
working_dir <- "/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd"

clean_options <- list(
  BrS_objects = make_list(BrS_object_1),         # all bronze-standard measurements.
  # SS_objects  = make_list(SS_object_1),  ## SS measurements 
  X_extra = c("SITE"), ##ZW add this to record the covariates considered for analysis; here only "SITE" for now.
  patho_taxo_dir = file.path(data_dir, "pathogen_category_simulation2.csv"), ##ZW would change the name to "pathogen_category_IRIS.csv" moving forward.
  allow_missing = FALSE)

m_opt1 <- list(likelihood   = list(cause_list = cause_list_iris,                # <---- fitted causes.
                                   k_subclass = c(1),                      # <---- no. of subclasses.
                                   Eti_formula = ~ -1+as.factor(SITE),                    # <---- etiology regression formula; only for cases.
                                   FPR_formula = list(
                                     MBS1 =  ~ -1+as.factor(SITE))),
               use_measurements = c("BrS"),
               prior = list(Eti_prior   = t(sapply(1:N.SITE, function(i) overall_uniform(1, cause_list_iris))), 
                    TPR_prior   = list(
                      BrS  = list(info  = "informative",
                                  input = "match_range",
                                  val   = list(
                                    MBS1 = list(up = list(rep(0.99,length(BrS_object_1$patho))),
                                                low = list(c(rep(0.5,length(BrS_object_1$patho)-2),0.001,0.5))
                                    ))))))     

model_options_discrete_reg <- m_opt1

assign_model(model_options_discrete_reg, data_iris_disc_u6)

working_dir = getwd()

Date     <- gsub("-", "_", Sys.Date())

set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) 

dated_strat_name_discrete_u6   <- file.path(working_dir,paste0("results/site_u6r"))
#dir.create(dated_strat_name_discrete_u6)

result_folder_discrete_u6 <- dated_strat_name_discrete_u6

mcmc_options_discrete_u6 <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(10000), #10000
  n.burnin = as.integer(5000), #5000
  n.thin = 1,
  individual.pred = FALSE, 
  ppd = !TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_u6,
  bugsmodel.dir = result_folder_discrete_u6,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_u6$result.folder, "data_clean_options.txt"))
dput(data_iris_disc_u6,file.path(mcmc_options_discrete_u6$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
start.time <- Sys.time()
nplcm_reg_discrete_only_u6 <- nplcm(data_iris_disc_u6,model_options_discrete_reg,mcmc_options_discrete_u6)
end.time <- Sys.time()
print(end.time-start.time)

gg_plot_etiology_strat(result_folder_discrete_u6, category= "<6 months")


# 6+ months----
table(data_iris_disc_o6$Y)
Nd <- 381
Nu <- 1150


assign_model(model_options_discrete_reg, data_iris_disc_o6)

#working_dir = tempdir()
#Date     <- gsub("-", "_", Sys.Date())

set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) ##ZW:  create results folder.

dated_strat_name_discrete_o6   <- file.path(working_dir,paste0("results/site_o6r"))
# dir.create(dated_strat_name_discrete_o6)

result_folder_discrete_o6 <- dated_strat_name_discrete_o6

mcmc_options_discrete_o6 <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(2000), #10000
  n.burnin = as.integer(1000), #5000
  n.thin = 1,
  individual.pred = FALSE, 
  ppd = !TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_o6,
  bugsmodel.dir = result_folder_discrete_o6,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_o6$result.folder, "data_clean_options.txt"))
dput(data_iris_disc_u6,file.path(mcmc_options_discrete_o6$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
# start.time <- Sys.time()
# nplcm_reg_discrete_only_o6 <- nplcm(data_iris_disc_o6,model_options_discrete_reg,mcmc_options_discrete_o6)
# end.time <- Sys.time()
# print(end.time-start.time)

gg_plot_etiology_strat(result_folder_discrete_o6, category = "6-11 months")


# model diagnostics----
## Under 6 months

pool_chains <- function(folder, num){
  coda::read.coda(file.path(folder,paste0("CODAchain",num,".txt")),
                  file.path(folder,"CODAindex.txt"),
                  quiet=TRUE)
}
u6_chainlist <- mcmc.list()
for (k in 1:mcmc_options_discrete_u6$n.chains){
 u6_chainlist[[k]] <- pool_chains(result_folder_discrete_u6, k)
  
}

# coda::read.coda(file.path(folder,paste0("CODAchain",k,".txt")),
#                 file.path(folder,"CODAindex.txt"),
#                 quiet=TRUE)
# 
# res_nplcm1 <- coda::read.coda(file.path(result_folder_discrete_u6,"CODAchain1.txt"),
#                               file.path(result_folder_discrete_u6,"CODAindex.txt"),
#                               quiet=TRUE)
# res_nplcm2 <-  coda::read.coda(file.path(result_folder_discrete_u6,"CODAchain2.txt"),
#                                file.path(result_folder_discrete_u6,"CODAindex.txt"),
#                                quiet=TRUE)
# 
# res_nplcm <- mcmc.list(res_nplcm1, res_nplcm2)
library(MCMCvis)
mcmc_vis_res <- MCMCchains(
  res_nplcm,
  params = "all",
  excl = NULL,
  ISB = TRUE,
  exact = TRUE,
  mcmc.list = TRUE,
  chain_num = NULL
)


MCMCtrace(
  u6_chainlist,
  params = c('thetaBS_1', 'psiBS_1', 'pEti'),
  pdf = T,
  ISB = F,
  exact = F,
  open_pdf = F,
  filename = "under6_traceplots",
  wd = "C:/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd/results/site_u6r")


# MCMCtrace(
#   res_nplcm,
#   params = c('psiBS_1'),
#   pdf = T,
#   ISB = F,
#   exact = F,
#   open_pdf = F,
#   filename = "psi_trace",
#   wd = "C:/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd")
# 
# MCMCtrace(
#   res_nplcm,
#   params = c('pEti'),
#   pdf = T,
#   ISB = F,
#   exact = F,
#   open_pdf = F,
#   filename = "pEti_trace",
#   wd = "C:/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd")


# res_nplcm_discrete_u6 <- coda::read.coda(file.path(result_folder_discrete_u6,"CODAchain1.txt"),
#                                       file.path(result_folder_discrete_u6,"CODAindex.txt"),
#                                       quiet=TRUE)
# 
# print_res_reg <- function(x) for (i in grep(x,colnames(res_nplcm_discrete_u6))) plot(res_nplcm_discrete_u6[,i],main=paste0('Trace of ',colnames(res_nplcm_discrete_u6)[i]))
# 
# print_res_reg("thetaBS") # check capitalization if not running

## 6+ months



# site, age_cat, and sex discrete model----
table(tot_lca_sub_alri_mult3fix$X$SITE)
table(tot_lca_sub_alri_mult3fix$X$SITE, tot_lca_sub_alri_mult3fix$X$AGE)
table(tot_lca_sub_alri_mult3fix$X$SITE, tot_lca_sub_alri_mult3fix$X$SEX)

table(tot_lca_sub_alri_mult3fix$Y)
Nd <- 1743
Nu <- 1068

BrS_object_1 <- make_meas_object(
  patho = cause_list_iris, 
  specimen = "MBS", test = "1",
  quality = "BrS", 
  cause_list = cause_list_iris)

# also save the data clean options:

# parent directory for testing code (LOCAL):
data_dir    <- "Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd/data"
working_dir <- "/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd"

clean_options <- list(
  BrS_objects = make_list(BrS_object_1),      
  X_extra = c("SITE", "AGE","SEX"), 
  patho_taxo_dir = file.path(data_dir, "pathogen_category_IRIS.csv"), 
  allow_missing = FALSE)

m_opt1 <- list(likelihood   = list(cause_list = cause_list_iris,               
                                   k_subclass = c(1),                      
                                   Eti_formula = ~ -1+as.factor(SITE)+ as.factor(AGE) + as.factor(SEX),                    # <---- etiology regression formula; only for cases.
                                   FPR_formula = list(
                                     MBS1 =  ~ -1+as.factor(SITE) + as.factor(AGE) + as.factor(SEX))),
               use_measurements = c("BrS"),
               prior = list(Eti_prior   = t(sapply(1:(N.SITE*3*2), function(i) overall_uniform(1, cause_list_iris))), 
                            TPR_prior   = list(
                              BrS  = list(info  = "informative",
                                          input = "match_range",
                                          val   = list(
                                            MBS1 = list(up = list(rep(0.99,length(BrS_object_1$patho))),
                                                        low = list(c(rep(0.5,2), 0.75, rep(0.5, length(BrS_object_1$patho)-5),0.001,0.5))
                                            ))))))     

model_options_discrete_reg <- m_opt1

assign_model(model_options_discrete_reg, tot_lca_sub_alri_mult3fix)

working_dir = getwd()

Date     <- gsub("-", "_", Sys.Date())

set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) 

dated_strat_name_discrete_mult3fix   <- file.path(working_dir,paste0("results/site_mult3_121621"))
dir.create(dated_strat_name_discrete_mult3fix)

result_folder_discrete_mult3 <- dated_strat_name_discrete_mult3fix

mcmc_options_discrete_mult3 <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(10000), #10000
  n.burnin = as.integer(5000), #5000
  n.thin = 1,
  individual.pred = TRUE, 
  ppd = !TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_mult3,
  bugsmodel.dir = result_folder_discrete_mult3,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_mult3$result.folder, "data_clean_options.txt"))
dput(tot_lca_sub_alri_mult3fix,file.path(mcmc_options_discrete_mult3$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
start.time <- Sys.time()
nplcm_reg_discrete_only_mult3fix <- nplcm(tot_lca_sub_alri_mult3fix,model_options_discrete_reg,mcmc_options_discrete_mult3)
end.time <- Sys.time()
print(end.time-start.time)


pool_chains <- function(folder, num){
  coda::read.coda(file.path(folder,paste0("CODAchain",num,".txt")),
                  file.path(folder,"CODAindex.txt"),
                  quiet=TRUE)
}

mult3_chainlist <- mcmc.list()
for (k in 1:mcmc_options_discrete_mult3$n.chains){
  mult3_chainlist[[k]] <- pool_chains(result_folder_discrete_mult3, k)
  
}

# mult3_chainlist(result_folder_discrete_mult3)
library(MCMCvis)
MCMCtrace(
  mult3_chainlist,
  params = c('thetaBS_1', 'psiBS_1', 'pEti'),
  pdf = T,
  ISB = F,
  exact = F,
  open_pdf = F,
  filename = "mult3_traceplots")

extract_post_means(result_folder_discrete_mult3)
gg_plot_etiology_strat(result_folder_discrete_mult3)


# res_nplcm_discrete_mult3 <- coda::read.coda(file.path(result_folder_discrete_mult3,"CODAchain1.txt"),
#                                                                                file.path(result_folder_discrete_mult3,"CODAindex.txt"),
#                                                                                quiet=TRUE)
# 
# print_res_reg <- function(x) for (i in grep(x,colnames(res_nplcm_discrete_mult3))) plot(res_nplcm_discrete_mult3[,i],main=paste0('Trace of ',colnames(res_nplcm_discrete_mult3)[i]))
# # 
# print_res_reg("thetaBS")

# site, age_cat, and sex discrete model -- without periods with mismatched cases/controls----
table(tot_lca_sub_alri_mult3a$X$SITE)
table(tot_lca_sub_alri_mult3a$X$SITE, tot_lca_sub_alri_mult3a$X$AGE)
table(tot_lca_sub_alri_mult3a$X$SITE, tot_lca_sub_alri_mult3a$X$SEX)

table(tot_lca_sub_alri_mult3a$Y)
Nd <- 1611
Nu <- 1063

BrS_object_1 <- make_meas_object(
  patho = cause_list_iris, 
  specimen = "MBS", test = "1",
  quality = "BrS", 
  cause_list = cause_list_iris)

# also save the data clean options:

# parent directory for testing code (LOCAL):
data_dir    <- "Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd/data"
working_dir <- "/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd"

clean_options <- list(
  BrS_objects = make_list(BrS_object_1),      
  X_extra = c("SITE", "AGE","SEX"), 
  patho_taxo_dir = file.path(data_dir, "pathogen_category_IRIS.csv"), 
  allow_missing = FALSE)

m_opt1 <- list(likelihood   = list(cause_list = cause_list_iris,               
                                   k_subclass = c(1),                      
                                   Eti_formula = ~ -1+as.factor(SITE)+ as.factor(AGE) + as.factor(SEX),                    # <---- etiology regression formula; only for cases.
                                   FPR_formula = list(
                                     MBS1 =  ~ -1+as.factor(SITE) + as.factor(AGE) + as.factor(SEX))),
               use_measurements = c("BrS"),
               prior = list(Eti_prior   = t(sapply(1:(N.SITE*3*2), function(i) overall_uniform(1, cause_list_iris))), 
                            TPR_prior   = list(
                              BrS  = list(info  = "informative",
                                          input = "match_range",
                                          val   = list(
                                            MBS1 = list(up = list(rep(0.99,length(BrS_object_1$patho))),
                                                        low = list(c(rep(0.5,length(BrS_object_1$patho)-2),0.001,0.5))
                                            ))))))     

model_options_discrete_reg <- m_opt1

assign_model(model_options_discrete_reg, tot_lca_sub_alri_mult3a)

working_dir = getwd()

Date     <- gsub("-", "_", Sys.Date())

set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) 

dated_strat_name_discrete_mult3a   <- file.path(working_dir,paste0("results/site_mult3a_120221"))
# dir.create(dated_strat_name_discrete_mult3a)

result_folder_discrete_mult3a <- dated_strat_name_discrete_mult3a

mcmc_options_discrete_mult3a <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(10000), #10000
  n.burnin = as.integer(5000), #5000
  n.thin = 1,
  individual.pred = TRUE, 
  ppd = !TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_mult3a,
  bugsmodel.dir = result_folder_discrete_mult3a,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_mult3a$result.folder, "data_clean_options.txt"))
dput(tot_lca_sub_alri_mult3a,file.path(mcmc_options_discrete_mult3a$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
start.time <- Sys.time()
nplcm_reg_discrete_only_mult3a <- nplcm(tot_lca_sub_alri_mult3a,model_options_discrete_reg,mcmc_options_discrete_mult3a)
end.time <- Sys.time()
print(end.time-start.time)


pool_chains <- function(folder, num){
  coda::read.coda(file.path(folder,paste0("CODAchain",num,".txt")),
                  file.path(folder,"CODAindex.txt"),
                  quiet=TRUE)
}

mult3_chainlist <- mcmc.list()
for (k in 1:mcmc_options_discrete_mult3$n.chains){
  mult3_chainlist[[k]] <- pool_chains(result_folder_discrete_mult3, k)
  
}

# mult3_chainlist(result_folder_discrete_mult3)
library(MCMCvis)
MCMCtrace(
  mult3_chainlist,
  params = c('thetaBS_1', 'psiBS_1', 'pEti'),
  pdf = T,
  ISB = F,
  exact = F,
  open_pdf = F,
  filename = "mult3_traceplots")

extract_post_means(result_folder_discrete_mult3)
gg_plot_etiology_strat(result_folder_discrete_mult3)




# site, age, and sex discrete model with k = 2----
table(tot_lca_sub_alri_mult3$X$SITE)
table(tot_lca_sub_alri_mult3$X$SITE, tot_lca_sub_alri_mult3$X$AGE)
table(tot_lca_sub_alri_mult3$X$SITE, tot_lca_sub_alri_mult3$X$SEX)

table(tot_lca_sub_alri_mult3$Y)
Nd <- 1743
Nu <- 1068

BrS_object_1 <- make_meas_object(
  patho = cause_list_iris, 
  specimen = "MBS", test = "1",
  quality = "BrS", 
  cause_list = cause_list_iris)

# also save the data clean options:

# parent directory for testing code (LOCAL):
data_dir    <- "Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd/data"
working_dir <- "/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd"

clean_options <- list(
  BrS_objects = make_list(BrS_object_1),      
  X_extra = c("SITE", "AGE","SEX"), 
  patho_taxo_dir = file.path(data_dir, "pathogen_category_IRIS.csv"), 
  allow_missing = FALSE)

m_opt1 <- list(likelihood   = list(cause_list = cause_list_iris,               
                                   k_subclass = c(2),                      
                                   Eti_formula = ~ -1+as.factor(SITE)+ as.factor(AGE) + as.factor(SEX),
                                   # <---- etiology regression formula; only for cases.
                                   FPR_formula = list(
                                     MBS1 =  ~ -1+as.factor(SITE) + as.factor(AGE) + as.factor(SEX))),
               use_measurements = c("BrS"),
               prior = list(Eti_prior   = c(2,2),# <--- etiology prior; sd of a Gaussian for coefficient of (ps basis, non-ps basis).# <-- ZW edits.
                                    half_nu_s2    = c(1/2,100/2), # <--- scale2/2 parameter in the t-distn for the intercept. #<-- ZW edits. 
                                    FPR_coef_prior = c(2,2), # <--- sd of a Gaussian for coefficient of (ps basis, non-ps basis). # <-- ZW edits.
                            TPR_prior   = list(
                              BrS  = list(info  = "informative",
                                          input = "match_range",
                                          val   = list(
                                            MBS1 = list(up = list(rep(0.99,length(BrS_object_1$patho))),
                                                        low = list(c(rep(0.5,length(BrS_object_1$patho)-2),0.001,0.5))
                                            ))))))     

model_options_discrete_reg <- m_opt1

assign_model(model_options_discrete_reg, tot_lca_sub_alri_mult3)

working_dir = getwd()

Date     <- gsub("-", "_", Sys.Date())

set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) 

dated_strat_name_discrete_mult3k2   <- file.path(working_dir,paste0("results/site_mult3k2_110921"))
#dir.create(dated_strat_name_discrete_mult3k2)

result_folder_discrete_mult3k2 <- dated_strat_name_discrete_mult3k2

mcmc_options_discrete_mult3k2 <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(10000), #10000
  n.burnin = as.integer(5000), #5000
  n.thin = 1,
  individual.pred = TRUE, 
  ppd = !TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_mult3k2,
  bugsmodel.dir = result_folder_discrete_mult3k2,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_mult3k2$result.folder, "data_clean_options.txt"))
dput(tot_lca_sub_alri_mult3,file.path(mcmc_options_discrete_mult3k2$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
start.time <- Sys.time()
nplcm_reg_discrete_only_mult3k2 <- nplcm(tot_lca_sub_alri_mult3,model_options_discrete_reg,mcmc_options_discrete_mult3k2)
end.time <- Sys.time()
print(end.time-start.time)



# site, age, and sex discrete model__old----
table(data_iris_disc_multd$Y)
Nd <- 1148
Nu <- 3550

BrS_object_1 <- make_meas_object(
  patho = cause_list_iris, 
  specimen = "MBS", test = "1",
  quality = "BrS", 
  cause_list = cause_list_iris)

# also save the data clean options:

# parent directory for testing code (LOCAL):
data_dir    <- "Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd/data"
working_dir <- "/Users/jkubale/Dropbox (University of Michigan)/Gordon Lab Data/Project Data/Personal Project Workspace/jkubale_gord/iris_ftd"

clean_options <- list(
  BrS_objects = make_list(BrS_object_1),      
    X_extra = c("SITE", "AGE","SEX"), 
  patho_taxo_dir = file.path(data_dir, "pathogen_category_IRIS.csv"), 
  allow_missing = FALSE)

m_opt1 <- list(likelihood   = list(cause_list = cause_list_iris,               
                                   k_subclass = c(1),                      
                                   Eti_formula = ~ -1+as.factor(SITE)+ as.factor(AGE) + as.factor(SEX),                    # <---- etiology regression formula; only for cases.
                                   FPR_formula = list(
                                     MBS1 =  ~ -1+as.factor(SITE) + as.factor(AGE) + as.factor(SEX))),
               use_measurements = c("BrS"),
               prior = list(Eti_prior   = t(sapply(1:(N.SITE*2*2), function(i) overall_uniform(1, cause_list_iris))), 
                            TPR_prior   = list(
                              BrS  = list(info  = "informative",
                                          input = "match_range",
                                          val   = list(
                                            MBS1 = list(up = list(rep(0.99,length(BrS_object_1$patho))),
                                                        low = list(c(rep(0.5,length(BrS_object_1$patho)-2),0.001,0.5))
                                            ))))))     

model_options_discrete_reg <- m_opt1

assign_model(model_options_discrete_reg, data_iris_disc_multd)

working_dir = getwd()

Date     <- gsub("-", "_", Sys.Date())

set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) 

dated_strat_name_discrete_multd   <- file.path(working_dir,paste0("results/site_mult"))
#dir.create(dated_strat_name_discrete_multd)

result_folder_discrete_multd <- dated_strat_name_discrete_multd

mcmc_options_discrete_multd <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(10000), #10000
  n.burnin = as.integer(5000), #5000
  n.thin = 1,
  individual.pred = FALSE, 
  ppd = !TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_multd,
  bugsmodel.dir = result_folder_discrete_multd,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_multd$result.folder, "data_clean_options.txt"))
dput(data_iris_disc_multd,file.path(mcmc_options_discrete_multd$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
start.time <- Sys.time()
nplcm_reg_discrete_only_multd <- nplcm(data_iris_disc_multd,model_options_discrete_reg,mcmc_options_discrete_multd)
end.time <- Sys.time()
print(end.time-start.time)

# gg_plot_etiology_strat(result_folder_discrete_mltd)

# With ppd----
set.seed(894)
# create folders to store the MCMC output and model results 

# dir.create(file.path(working_dir,"results")) 

# dated_strat_name_discrete_multdppd   <- file.path(working_dir,paste0("results/site_mult_ppd"))
#dir.create(dated_strat_name_discrete_multd)

result_folder_discrete_multdppd <- dated_strat_name_discrete_multdppd

mcmc_options_discrete_multdppd <-  list(
  debugstatus = TRUE,
  n.chains   = 4,
  n.itermcmc = as.integer(10000), #10000
  n.burnin = as.integer(5000), #5000
  n.thin = 1,
  individual.pred = FALSE, 
  ppd = TRUE, #ZW: set to false now, this is for posterior predictive checking; can set to true later.
  get.pEti = FALSE,
  result.folder = result_folder_discrete_multdppd,
  bugsmodel.dir = result_folder_discrete_multdppd,
  jags.dir = "",
  use_jags = TRUE
)

dput(clean_options, file.path(mcmc_options_discrete_multdppd$result.folder, "data_clean_options.txt"))
dput(data_iris_disc_multd,file.path(mcmc_options_discrete_multdppd$result.folder,"data_nplcm.txt")) 

# Commented out so I don't accidentally re-run
start.time <- Sys.time()
nplcm_reg_discrete_only_multdppd <- nplcm(data_iris_disc_multd,model_options_discrete_reg,mcmc_options_discrete_multdppd)
end.time <- Sys.time()
print(end.time-start.time)