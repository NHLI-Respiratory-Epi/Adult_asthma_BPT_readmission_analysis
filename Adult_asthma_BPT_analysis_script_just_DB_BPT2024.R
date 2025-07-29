# Adult asthma BPT analysis script

library(tidyverse)
library(dplyr)
library(tidyr)
library(lme4)
library(stringr)
source("C:/Users/aadamson/Documents/R/My R functions/tidyoutput.R")
source("C:/Users/aadamson/Documents/R/My R functions/lintestOR.R")

dat <- readRDS("C:/Alex Harley/Audit_2023_onwards/2022-2023/AA/Data/tidyData/linked_audit_HES_PEDW_ONS_data_AA_2022-23.RDS")

dat %>% filter(CCIweighted == 0) %>% select(DIAG_01:DIAG_20) %>% print(n=1000)

dat %>% select(RSR_24hour, RSR_BPT) %>% table(useNA = "ifany")

dat %>% filter(is.na(RSR_BPT)) %>% select(life_status, discharge_bundle_yes_no) %>% summary()

dat %>% colnames()

# get rid of people 

# how many index admissions?

dat %>% filter(index_admission_flag == 1) %>% nrow()

# remove admissions ineligible for the BPT (those who died in hospital or who were transferred)
# and those who were in Wales because we didn't request access to Welsh data.

# number to exclude:
dat %>% filter(index_admission_flag == 1) %>% 
  filter(discharge_bundle == "Patient transferred to another hospital" | life_status == "Died as inpatient" |
           country == "Wales") %>% nrow()

dat <- dat %>% filter(discharge_bundle != "Patient transferred to another hospital") %>%
  filter(life_status == "Alive") %>% filter(country == "England")


# And just keep those who were male or female and relevel gender to remove the other options

dat <- dat %>% filter(gender %in% c("Male", "Female")) %>% 
  mutate(gender = factor(gender, levels  = c("Male", "Female")))

summary(dat$gender)
# let's just see how many inappropriate admissions we have according to the asthma diagnosis in HES

nrow(dat)
dat %>% select(hosp_code) %>% unique() %>% nrow()

dat %>% filter(index_admission_flag == 1) %>% select(adult_asthma_coded_admission) %>% table()

# wow...2227 inappropriately included admissions??

dat %>% filter(index_admission_flag == 1) %>% filter(adult_asthma_coded_admission == 0) %>% 
  select(DIAG_01) %>% table()

dat %>% filter(index_admission_flag == 1) %>% filter(adult_asthma_coded_admission == 0) %>% 
  mutate(DIAG_01_3d = as.character(DIAG_01_3d)) %>% select(DIAG_01_3d) %>% table() %>% sort()

dat %>% filter(index_admission_flag == 1) %>% filter(adult_asthma_coded_admission == 0) %>% 
  select(DIAG_04) %>% table() %>% sort()

dat %>% select(adult_asthma_coded_admission, discharge_bundle) %>% table()

colnames(dat)

# Independent variables:


# calculate asthma readmissions

table(dat$read30_row_flag, dat$index_admission_flag, useNA = "ifany")
table(dat$read30_row_flag, dat$index_admission_flag, useNA = "ifany")
table(dat$adult_asthma_coded_admission)

# keep this variable just for the index admissions
dat$adult_asthma_coded_admission_including_index <- dat$adult_asthma_coded_admission

# now recode them so they are only for the readmissions

dat$adult_asthma_coded_admission[dat$index_admission_flag == 1] <- 0

table(dat$adult_asthma_coded_admission)

dat$asthma_read30_row_flag <- dat$read30_row_flag
dat$asthma_read30_row_flag[dat$adult_asthma_coded_admission == 0] <- 0

dat$asthma_read90_row_flag <- dat$read90_row_flag
dat$asthma_read90_row_flag[dat$adult_asthma_coded_admission == 0] <- 0



dat <- dat %>% group_by(patient_ID) %>% mutate(read30 = max(read30_row_flag),
                                               read90 = max(read90_row_flag),
                                               asthma_read30 = max(asthma_read30_row_flag),
                                               asthma_read90 = max(asthma_read90_row_flag)) %>%
  ungroup()
                                                 

# keep just the index admissions

dat <- dat %>% filter(index_admission_flag == 1)

dat$RSR_BPT_bin <- NA
dat$RSR_BPT_bin[dat$RSR_BPT == "<24 hours"] <- 1
dat$RSR_BPT_bin[dat$RSR_BPT == ">=24 hours or no RSR"] <- 0

dat$RSR_bin <- NA
dat$RSR_bin[dat$RSR == "Yes"] <- 1
dat$RSR_bin[dat$RSR == "No"] <- 0

table(dat$RSR_BPT, dat$RSR_BPT_bin, useNA = "ifany")
table(dat$RSR, dat$RSR_bin, useNA = "ifany")

# RSR_BPT
# All DB variables
# Discharge bundle
# Combination: RSR_BPT + PAAP + inhaler + smoking cessation + DB_spec_review_4_weeks

dat <- dat %>% mutate(BPT2024 = ifelse(DB_inhaler == 1 & DB_maintenance == 1 & 
                                         DB_PAAP == 1 & (DB_smoke == 1 | is.na(DB_smoke)) &
                                         DB_spec_review_4_weeks == 1 & 
                                         discharge_bundle_yes_no == "Yes" & RSR_BPT_bin == 1, 1, 0))

dat <- dat %>% mutate(DB_BPT2024 = ifelse(DB_inhaler == 1 & DB_maintenance == 1 & 
                                         DB_PAAP == 1 & (DB_smoke == 1 | is.na(DB_smoke)) &
                                         DB_spec_review_4_weeks == 1 & 
                                         discharge_bundle_yes_no == "Yes", 1, 0))

dat <- dat %>% mutate(BPT2023 = ifelse(DB_inhaler == 1 & DB_maintenance == 1 & 
                                         DB_PAAP == 1 & (DB_smoke == 1 | is.na(DB_smoke)) & 
                                         discharge_bundle_yes_no == "Yes" & RSR_BPT_bin == 1, 1, 0))


dat %>% select(BPT2023, BPT2024) %>% table()

# and sort out some of the variables we need

summary(dat$CCIweightedcat)
levels(dat$CCIweightedcat)

dat$CCIweightedcat <- as.character(dat$CCIweightedcat) 
dat$CCIweightedcat[dat$CCIweightedcat %in% c("6", "7+")] <- "6+"

dat$CCIweightedcat <- factor(dat$CCIweightedcat, levels = c("0-1", "2", "3", "4", "5", "6+"))


# and we centre age

mean(dat$age)
table(dat$age)

# mean is ~50 so we subtract 50.

dat$age_orig <- dat$age

dat$age <- dat$age - 50

# splitting things into independent variables, confounders, and outcomes

# ID
ID_vars <- c("patient_ID", "hosp_code", "trust_code", "country",
             "hosp_name", "trust_name", "ICS", "region")

cluster_var <- c("hosp_code")



# covariates
covars <- c("IMD_quintile", "age", "gender", "smoke_status", "inhaled_steroids_dis",                        
            "oral_steroids_dis", "oral_steroids_rescue_history", "IMD_quintile",
            "asthma_sev", "CCIweightedcat")  

for (i in covars) {
  print(i)
    print(summary(dat[[i]]))
    print(levels(dat[[i]]))
  
}

covars_alt_smoking <- c("IMD_quintile", "age", "gender", "smoke_status_DB", "inhaled_steroids_dis",                        
            "oral_steroids_dis", "oral_steroids_rescue_history", "IMD_quintile",
            "asthma_sev", "CCIweightedcat")  


# make the covariates into a suitable format

lintestOR(xx = dat, depvar = "read30", indvar = "agecat", logodds = TRUE)

dat$agecat2 <- factor(dat$age)
dat$CCIweighted2 <- factor(dat$CCIweighted)

lintestOR(xx = dat, depvar = "read30", indvar = "agecat2", logodds = TRUE)
lintestOR(xx = dat, depvar = "read30", indvar = "CCIweightedcat", logodds = TRUE)
lintestOR(xx = dat, depvar = "read30", indvar = "CCIweighted2", logodds = TRUE)



# independent variables
independent_vars <- c("discharge_bundle_yes_no", "DB_inhaler", "DB_maintenance",
                      "DB_adherence", "DB_PAAP", "DB_triggers", "DB_smoke",                                    
                      "DB_comm_FU_2_days", "DB_FU_any", "DB_spec_review_4_weeks",
                      "DB_none", "RSR_bin", "RSR_BPT_bin", "BPT2024", "DB_BPT2024")

# other variables
other_vars <- c("arrival_day_of_week", "died30", "died90", "BPT2023")

# outcome variables
outcome_vars <- c("read30", "asthma_read30", "read90", "asthma_read90") 

dat <- dat %>% select(all_of(c(ID_vars, covars, independent_vars, outcome_vars)))

# And let's make a data frame for all the analyses we have to do



nrow(dat)
dat %>% select(hosp_code) %>% group_by(hosp_code) %>% unique() %>% nrow()

# unadjusted analyses

analyses_df_unadj <- data.frame(outcome_vars = rep(outcome_vars, each = length(independent_vars)), 
                                independent_vars = rep(independent_vars, 4), 
                                varcheck = NA, est = NA, lo = NA, hi = NA, pval = NA)

colnames(analyses_df_unadj)

model_output_list_unadj <- vector(nrow(analyses_df_unadj), mode = "list")


# for (i in 1:3) {
for (i in 1:nrow(analyses_df_unadj)) {

  model_output_list_unadj[i] <-  glmer(formula(paste(analyses_df_unadj$outcome_vars[i]," ~ ",
                                               analyses_df_unadj$independent_vars[i]," + ", 
                                               "(1 | ", cluster_var, ")", collapse = "")),
      family = binomial(link = "logit"), data = dat)

  
  print(tidyoutput(x = model_output_list_unadj[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2))
  analyses_df_unadj[i, 3:7] <- tidyoutput(x = model_output_list_unadj[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2)[2, ]
    
}

# adjusted analyses - needs to have smoking_DB in a seperate part

independent_vars_no_DB_smoke <- independent_vars[-7]
covars_no_smoke_status <- covars[-4]

analyses_df_adj <- data.frame(outcome_vars = rep(outcome_vars, each = length(independent_vars_no_DB_smoke)), 
                                independent_vars = rep(independent_vars_no_DB_smoke, 4), 
                                varcheck = NA, est = NA, lo = NA, hi = NA, pval = NA)

colnames(analyses_df_adj)

model_output_list_adj <- vector(nrow(analyses_df_adj), mode = "list")


# for (i in 1:3) {
for (i in 1:nrow(analyses_df_adj)) {
  
  model_output_list_adj[i] <-  glmer(formula(paste(analyses_df_adj$outcome_vars[i]," ~ ",
                                                     analyses_df_adj$independent_vars[i]," + ",
                                                     paste(covars, collapse = "+"),
                                                     "+ (1 | ", cluster_var, ")", collapse = "")),
                                       family = binomial(link = "logit"), data = dat)
  
  
  print(tidyoutput(x = model_output_list_adj[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2))
  analyses_df_adj[i, 3:7] <- tidyoutput(x = model_output_list_adj[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2)[2, ]
  
}

# smoking specific model
table(dat$smoke_status, dat$DB_smoke)



analyses_df_adj_smoking <- data.frame(outcome_vars = outcome_vars, 
                              independent_vars = rep("DB_smoke", 4), 
                              varcheck = NA, est = NA, lo = NA, hi = NA, pval = NA)

colnames(analyses_df_adj_smoking)

model_output_list_adj_smoking <- vector(nrow(analyses_df_adj_smoking), mode = "list")

dat %>% select(hosp_code, smoke_status) %>% table()

library(blme)



# for (i in 4) {
for (i in 1:nrow(analyses_df_adj_smoking)) {
  
  model_output_list_adj_smoking[i] <-  glmer(formula(paste(analyses_df_adj_smoking$outcome_vars[i]," ~ ",
                                                           analyses_df_adj_smoking$independent_vars[i]," + ",
                                                   paste(covars_no_smoke_status, collapse = "+"),
                                                   "+ (1 | ", cluster_var, ")", collapse = "")),
                                     family = binomial(link = "logit"), 
                                     control = glmerControl(optCtrl = list(maxfun = 100000)),
                                     data = (dat %>% filter(smoke_status %in% 
                                                              c("Current smoker", "Current smoker and current vaper")) %>%
                                               mutate(hosp_code = factor(hosp_code))))
  
  
  print(tidyoutput(x = model_output_list_adj_smoking[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2))
  analyses_df_adj_smoking[i, 3:7] <- tidyoutput(x = model_output_list_adj_smoking[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2)[2, ]
  
}


summary(model_output_list_adj_smoking[[1]])

saveRDS(analyses_df_unadj, "C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/analyses_df_unadj_DB_BPT2024.RDS")
saveRDS(analyses_df_adj, "C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/analyses_df_adj_DB_BPT2024.RDS")
saveRDS(analyses_df_adj_smoking, "C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/analyses_df_adj_smoking_DB_BPT2024.RDS")
saveRDS(model_output_list_unadj, "C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/model_output_list_unadj_DB_BPT2024.RDS")
saveRDS(model_output_list_adj, "C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/model_output_list_adj_DB_BPT2024.RDS")
saveRDS(model_output_list_adj_smoking, "C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/model_output_list_adj_smoking_DB_BPT2024.RDS")



# final model

# need to combine the smoking status variable. But, we also need to think about what we are trying to show.
# We have all the individual elements now. Now we move onto: 
# - individual elements adjusted for other elements but not double counted
# the overall BPT not adjusted for elements that make it up


# individual elements not double counted:

independent_vars

independent_vars_indiv_elements <- c("DB_inhaler",
                                     "DB_maintenance",
                                     "DB_adherence",
                                     "DB_PAAP",
                                     "DB_triggers",
                                     "DB_comm_FU_2_days",
                                     "DB_spec_review_4_weeks",
                                     "RSR_BPT_bin")
                                     

# Combination: RSR_BPT + PAAP + inhaler + smoking cessation + DB_spec_review_4_weeks

independent_vars_BPT_plus <- c("BPT2024",
                               "DB_adherence",
                               "DB_triggers",
                               "DB_comm_FU_2_days")

model_all_indiv <- vector(4, mode = "list")

for (i in 1:length(outcome_vars)) {
  
model_all_indiv[i] <-  glmer(formula(paste(outcome_vars[i]," ~ ",
                                  paste(independent_vars_indiv_elements, collapse = "+"),
                                  " + ", paste(covars, collapse = "+"),
                                  " + (1 | ", cluster_var, ")", collapse = "")),
                                                  family = binomial(link = "logit"), data = dat)

print(tidyoutput(x = model_all_indiv[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2))

}
library(dplyr)
dat %>% select(one_of(independent_vars_indiv_elements, covars))  %>% mutate_all(~as.numeric()) %>%  str() # psych::corr.test()
                      

saveRDS(model_all_indiv, "C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/model_all_indiv_DB_BPT2024.RDS")

model_all_BPT <- vector(4, mode = "list")

for (i in 1:length(outcome_vars)) {
  
  model_all_BPT[i] <-  glmer(formula(paste(outcome_vars[i]," ~ ",
                                             paste(independent_vars_BPT_plus, collapse = "+"),
                                             " + ", paste(covars, collapse = "+"),
                                             " + (1 | ", cluster_var, ")", collapse = "")),
                               family = binomial(link = "logit"), data = dat)
  
  print(tidyoutput(x = model_all_BPT[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2))
  
}

saveRDS(model_all_BPT, "C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/model_all_BPT_DB_BPT2024.RDS")


summary(model_output_list_adj[[1]])

analyses_df_unadj
analyses_df_adj

warnings(model_all_BPT[[1]])

output_folder <- "C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/"

analyses_df_unadj$pval <- round(analyses_df_unadj$pval, 3)
write.csv(analyses_df_unadj, file.path(output_folder, "analyses_df_unadj.csv"),
          row.names = FALSE)

analyses_df_adj$pval <- round(analyses_df_adj$pval, 3)
write.csv(analyses_df_adj, file.path(output_folder, "analyses_df_adj.csv"),
          row.names = FALSE)

analyses_df_adj_smoking$pval <- round(analyses_df_adj_smoking$pval, 3)
write.csv(analyses_df_adj_smoking, file.path(output_folder, "analyses_df_adj_smoking.csv"),
          row.names = FALSE)

library(stringr)
library(tibble)

analyses_df_adj_smoking$outcome_vars[1]

for (i in 1:4) {
  

tidyoutput(x = model_all_BPT[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2) %>%
    mutate(pval = round(pval, 3)) %>% 
    write.csv(file.path(output_folder, paste0("model_all_BPT_", analyses_df_adj_smoking$outcome_vars[i], ".csv")), row.names = FALSE)

  tidyoutput(x = model_all_indiv[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2) %>%
    mutate(pval = round(pval, 3)) %>% 
    write.csv(file.path(output_folder, paste0("model_all_indiv_", analyses_df_adj_smoking$outcome_vars[i], ".csv")), row.names = FALSE)
  
}




# New from James' comments

DB_BPT2024

"DB_adherence",
"DB_triggers",
"DB_comm_FU_2_days")

model_all_indiv <- readRDS("C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/model_all_indiv.RDS")

tidyoutput(x = model_all_indiv[[4]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2) %>%
  filter(Variable != "smoke_statusNever smoked and current vaper") %>% select(-pval)

simple_model_all_indiv <- vector(4, mode = "list")

for (i in 1:length(outcome_vars)) {
  
  simple_model_all_indiv[[i]] <-  glm(formula(paste(outcome_vars[i]," ~ ",
                                             paste(independent_vars_indiv_elements, collapse = "+"),
                                             " + ", paste(covars, collapse = "+"),
                                             # " + (1 | ", cluster_var, ")", 
                                             collapse = "")),
                               family = binomial(link = "logit"), data = dat)
  
print(tidyoutput(x = simple_model_all_indiv[[i]], meth = "Wald", MEM = FALSE, pval = TRUE, roundno = 2))
  
}



confint(coef(simple_model_all_indiv[[1]]))

library(dplyr)
dat %>% select(one_of(independent_vars_indiv_elements, covars))  %>% mutate_all(~as.numeric()) %>%  str() # psych::corr.test()


saveRDS(model_all_indiv, "C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/model_all_indiv_DB_BPT2024.RDS")

model_all_BPT <- vector(4, mode = "list")

for (i in 1:length(outcome_vars)) {
  
  model_all_BPT[i] <-  glmer(formula(paste(outcome_vars[i]," ~ ",
                                           paste(independent_vars_BPT_plus, collapse = "+"),
                                           " + ", paste(covars, collapse = "+"),
                                           " + (1 | ", cluster_var, ")", collapse = "")),
                             family = binomial(link = "logit"), data = dat)
  
  print(tidyoutput(x = model_all_BPT[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2))
  
}

saveRDS(model_all_BPT, "C:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/model_all_BPT_DB_BPT2024.RDS")




