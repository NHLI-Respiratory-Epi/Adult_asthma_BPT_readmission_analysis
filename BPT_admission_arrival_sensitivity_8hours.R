# Adult asthma BPT analysis script

# For this one, we will add on 8 hours to let people make the BPT.

library(tidyverse)
library(car)
library(lme4)
library(finalfit)
library(parallel)
library(MuMIn)
library(DHARMa)
library(domir)
source("G:/Alex Harley/Audit_2023_onwards/My R functions/tidyoutput.R")
source("G:/Alex Harley/Audit_2023_onwards/My R functions/lintestOR.R")

clean_up <- function(df) {
  df %>% mutate(coefficients = paste0(sprintf("%.2f", .$est), " (", 
                                      sprintf("%.2f", .$lo), " to ", 
                                      sprintf("%.2f", .$hi), ")"),
                pval = sprintf("%.3f", round(pval, 3)),
                est = NULL,
                lo = NULL,
                hi = NULL)
}


dat <- readRDS("G:/Alex Harley/Audit_2023_onwards/2022-2023/AA/Data/tidyData/linked_audit_HES_PEDW_ONS_data_AA_2022-23.RDS")




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


dat$read30_categorical <- "Not readmitted"
dat$read30_categorical[dat$read30 == 1] <- "Readmitted"
dat$read30_categorical <- factor(dat$read30_categorical, levels = c("Readmitted", "Not readmitted"))

dat$read90_categorical <- "Not readmitted"
dat$read90_categorical[dat$read90 == 1] <- "Readmitted"
dat$read90_categorical <- factor(dat$read90_categorical, levels = c("Readmitted", "Not readmitted"))

dat$asthma_read30_categorical <- "Not readmitted"
dat$asthma_read30_categorical[dat$asthma_read30 == 1] <- "Readmitted"
dat$asthma_read30_categorical <- factor(dat$asthma_read30_categorical, levels = c("Readmitted", "Not readmitted"))

dat$asthma_read90_categorical <- "Not readmitted"
dat$asthma_read90_categorical[dat$asthma_read90 == 1] <- "Readmitted"
dat$asthma_read90_categorical <- factor(dat$asthma_read90_categorical, levels = c("Readmitted", "Not readmitted"))


table(dat$asthma_read30, dat$asthma_read30_categorical)
table(dat$asthma_read90, dat$asthma_read90_categorical)

# As we are only using English data, we should use the English IMDs

summary(dat$LSOA)

imd <- read.csv("G:/Alex Harley/Audit_2023_onwards/General UK data/IMD/clean_IMD2019_England.csv")

imd <- imd %>% select(LSOA = LSOA_code_2011, IMD_decile)

dat <- dat %>% left_join(., imd, by = "LSOA")


table(dat$IMD_decile, dat$IMD_quintile)

dat$IMD_quintile <- NA

dat$IMD_quintile[dat$IMD_decile == 1] <- 1
dat$IMD_quintile[dat$IMD_decile == 2] <- 1
dat$IMD_quintile[dat$IMD_decile == 3] <- 2
dat$IMD_quintile[dat$IMD_decile == 4] <- 2
dat$IMD_quintile[dat$IMD_decile == 5] <- 3
dat$IMD_quintile[dat$IMD_decile == 6] <- 3
dat$IMD_quintile[dat$IMD_decile == 7] <- 4
dat$IMD_quintile[dat$IMD_decile == 8] <- 4
dat$IMD_quintile[dat$IMD_decile == 9] <- 5
dat$IMD_quintile[dat$IMD_decile == 10] <- 5
dat$IMD_quintile[is.na(dat$IMD_decile)] <- "Missing IMD quintile"

dat$IMD_quintile <- factor(dat$IMD_quintile, levels = c("1", "2", "3", "4", "5", "Missing IMD quintile"))

table(dat$IMD_decile, dat$IMD_quintile, useNA = "ifany")





# keep just the index admissions

# get rid of people 


# remove admissions ineligible for the BPT (those who died in hospital or who were transferred)
# and those who were in Wales because we didn't request access to Welsh data.

table(dat$index_admission_flag)
# start off with index admissions at English hospitals

dat <- dat %>% filter(index_admission_flag == 1)

# how many people at this point? 
nrow(dat)

table(dat$country)

dat <- dat %>% filter(country == "England")

nrow(dat)
# this many people died as inpatient: 
dat %>% filter(life_status == "Died as inpatient") %>% nrow()


# this many people died as inpatient: 
dat %>% filter(life_status == "Died as inpatient") %>% nrow()

# This many people were transferred to another hospital:
dat %>% filter(discharge_bundle == "Patient transferred to another hospital") %>% nrow()


dat <- dat %>% filter(discharge_bundle != "Patient transferred to another hospital") %>%
  filter(life_status == "Alive") 

nrow(dat)
# And just keep those who were male or female and relevel gender to remove the other options

dat %>% filter(!(gender %in% c("Male", "Female"))) %>% pull(gender) %>% table()
  
dat <- dat %>% filter(gender %in% c("Male", "Female")) %>% 
  mutate(gender = factor(gender, levels  = c("Male", "Female")))

summary(dat$gender)

# this leaves us with this many patients:
nrow(dat)

table(dat$died30_categorical, dat$read30_categorical)
table(dat$died90_categorical, dat$read90_categorical)

# only 26 and 52 people respectively died within 30- and 90-days and were not captured in the 
# readmission data. This number is so low that we do not think it will have a large effect 
# for censoring purposes. Less than 0.5% of patients died who were not captured in the readmissions data.

#=====================================================================================================#
# # let's just see how many inappropriate admissions we have according to the asthma diagnosis in HES
# # Come back to this section for sensitivity analysis
# 
# dat %>% select(hosp_code) %>% unique() %>% nrow()
# 
# dat %>% filter(index_admission_flag == 1) %>% select(adult_asthma_coded_admission) %>% table()
# 
# # wow...2227 inappropriately included admissions??
# 
# dat %>% filter(index_admission_flag == 1) %>% filter(adult_asthma_coded_admission == 0) %>% 
#   select(DIAG_01) %>% table()
# 
# dat %>% filter(index_admission_flag == 1) %>% filter(adult_asthma_coded_admission == 0) %>% 
#   mutate(DIAG_01_3d = as.character(DIAG_01_3d)) %>% select(DIAG_01_3d) %>% table() %>% sort()
# 
# dat %>% filter(index_admission_flag == 1) %>% filter(adult_asthma_coded_admission == 0) %>% 
#   select(DIAG_04) %>% table() %>% sort()
# 
# dat %>% select(adult_asthma_coded_admission, discharge_bundle) %>% table()
# 
# colnames(dat)
#=====================================================================================================#

# THESE ARE CHANGING NOW FOR THE +8 HOUR BPT SENSITIVITY ANALYSIS!!!!!

# Create independent variables:
colnames(dat)

summary(dat$arrival_to_RSR_hours)


dat$RSR_BPT_bin <- NA
dat$RSR_BPT_bin[dat$RSR_BPT == "<24 hours"] <- 1
dat$RSR_BPT_bin[dat$RSR_BPT == ">=24 hours or no RSR"] <- 0

table(dat$RSR_BPT_bin, useNA = "ifany")

dat$RSR_BPT_bin[dat$arrival_to_RSR_hours < 28] <- 1

table(dat$RSR_BPT_bin, useNA = "ifany")

# 6279 under arrival BPT, 7596 under + 8 hours leeway BPT

dat$RSR_bin <- NA
dat$RSR_bin[dat$RSR == "Yes"] <- 1
dat$RSR_bin[dat$RSR == "No"] <- 0

table(dat$RSR_BPT, dat$RSR_BPT_bin, useNA = "ifany")
table(dat$RSR, dat$RSR_bin, useNA = "ifany")


# Let's create each the medication element of the 'Asthma 4' discharge bundle
# All other elements are already created. 

dat <- dat %>% mutate(DB_A4_meds = ifelse(DB_adherence == 1 & DB_maintenance == 1 &
                                         DB_inhaler == 1, 1, 0))



# RSR_BPT
# All DB variables
# Discharge bundle
# Combination: RSR_BPT + PAAP + inhaler + smoking cessation + DB_spec_review_4_weeks



dat <- dat %>% mutate(DB_BPT = ifelse(DB_inhaler == 1 & DB_maintenance == 1 &
                                        DB_adherence == 1 & DB_PAAP == 1 &
                                        DB_FU_any == 1 &
                                        discharge_bundle == "Yes" & (DB_smoke == 1 | is.na(DB_smoke)),
                                      1, 0))


dat <- dat %>% mutate(BPT2024 = ifelse(DB_inhaler == 1 & DB_maintenance == 1 & 
                                         DB_adherence == 1 & DB_PAAP == 1 & (DB_smoke == 1 | is.na(DB_smoke)) &
                                         DB_FU_any == 1 & 
                                         discharge_bundle_yes_no == "Yes" & RSR_BPT_bin == 1, 1, 0))

dat$BPT2024_cat <- "BPT met"
dat$BPT2024_cat[dat$BPT2024 == 0] <- "BPT not met"

dat$BPT2024_cat <- factor(dat$BPT2024_cat, levels = c("BPT met", "BPT not met"))

table(dat$BPT2024, dat$BPT2024_cat)

# now, 4351 meet the BPT. Before, it was 3627.

# and sort out some of the variables we need

# # look at linearity
# 
# lintestOR(xx = dat, depvar = "read30", indvar = "agecat", logodds = TRUE)
# 
# dat$agecat2 <- factor(dat$age)
# dat$CCIweighted2 <- factor(dat$CCIweighted)
# 
# lintestOR(xx = dat, depvar = "read30", indvar = "agecat2", logodds = TRUE)
# lintestOR(xx = dat, depvar = "read30", indvar = "CCIweightedcat", logodds = TRUE)
# lintestOR(xx = dat, depvar = "read30", indvar = "CCIweighted2", logodds = TRUE)



summary(dat$CCIweightedcat)
levels(dat$CCIweightedcat)

dat$CCIweightedcat <- as.character(dat$CCIweightedcat) 
dat$CCIweightedcat[dat$CCIweightedcat %in% c("6", "7+")] <- "6+"

dat$CCIweightedcat <- factor(dat$CCIweightedcat, levels = c("0-1", "2", "3", "4", "5", "6+"))
summary(dat$CCIweightedcat)


# and we centre age

mean(dat$age)
table(dat$age)

# mean is ~51 so we subtract 50.

dat$age_orig <- dat$age

dat$age <- dat$age - 51
mean(dat$age)


# I'm recoding smoking to remove vaping, as it is not clear
# if being a smoker means they are not a vaper or if the question wasn't asked

dat$smoke_status[dat$smoke_status == "Never smoked and current vaper"] <- "Never smoked"
dat$smoke_status[dat$smoke_status == "Ex-smoker and current vaper"] <- "Ex-smoker"
dat$smoke_status[dat$smoke_status == "Current smoker and current vaper"] <- "Current smoker"
dat$smoke_status[dat$smoke_status == "Current vaper"] <- "Not recorded"

dat$smoke_status <- factor(dat$smoke_status, levels = c("Never smoked", "Ex-smoker", "Current smoker", "Not recorded"))

summary(dat$smoke_status)

dat$smoke_status_with_cessation_DB <- as.character(dat$smoke_status)

dat$smoke_status_with_cessation_DB[dat$smoke_status == "Current smoker" & dat$DB_smoke == 1] <- "Current smoker given cessation advice"
dat$smoke_status_with_cessation_DB[dat$smoke_status == "Current smoker" & dat$DB_smoke == 0] <- "Current smoker not given cessation advice"
dat$smoke_status_with_cessation_DB <- factor(dat$smoke_status_with_cessation_DB, 
                                             levels = c("Never smoked", "Ex-smoker", 
                                                        "Current smoker given cessation advice",
                                                        "Current smoker not given cessation advice",
                                                        "Not recorded"))

summary(dat$smoke_status_with_cessation_DB)


# splitting things into independent variables, confounders, and outcomes

# ID
ID_vars <- c("patient_ID", "hosp_code", "trust_code", "country",
             "hosp_name", "trust_name", "ICS", "region")

cluster_var <- c("hosp_code")




# covariates
covars <- c("IMD_quintile", "age", "gender", "smoke_status", "oral_steroids_rescue_history",
            "asthma_sev", "CCIweightedcat")  


for (i in covars) {
  print(i)
    print(summary(dat[[i]]))
    print(levels(dat[[i]]))
  
}

covars_for_smoke <- c("IMD_quintile", "age", "gender",  "oral_steroids_rescue_history",
            "asthma_sev", "CCIweightedcat")  

covars_for_all_elements_included <- c("IMD_quintile", "age", "gender", "smoke_status_with_cessation_DB",  "oral_steroids_rescue_history",
                                      "asthma_sev", "CCIweightedcat")  


# independent variables
independent_vars <- c("BPT2024", # overall BPT 
                      "RSR_BPT_bin", "RSR_bin", # RSR section 
                      "DB_BPT", "discharge_bundle_yes_no", # DB section 
                      "DB_A4_meds", "DB_inhaler", "DB_maintenance", "DB_adherence", # DB meds  
                      "DB_PAAP", # DB PAAP
                      "DB_FU_any", "DB_comm_FU_2_days", "DB_spec_review_4_weeks", # DB FU  
                      "DB_none")

independent_vars_smoke <- c("DB_smoke")


# We don't bother including RSR within 24 hours as part of the independent
# vars part as obviously not useful - instead, we just look at the discharge bundle elements

independent_vars_for_comb <- c("DB_inhaler", "DB_maintenance", "DB_adherence", # DB meds  
                               "DB_PAAP", "DB_comm_FU_2_days", "DB_spec_review_4_weeks")



# outcome variables
outcome_vars <- c("read30", "asthma_read30", "read90", "asthma_read90") 

outcome_vars_table1 <- c("read30_categorical", "asthma_read30_categorical", 
                         "read90_categorical", "asthma_read90_categorical") 


table1_explan <- c("IMD_quintile", "age_orig", "gender", "smoke_status", "inhaled_steroids_dis",                        
                 "oral_steroids_dis", "oral_steroids_rescue_history",
                 "asthma_sev", "CCIweightedcat", # (this is covars but with original age) 
                 independent_vars[-1], independent_vars_smoke, outcome_vars_table1)

table1_dep <- "BPT2024_cat"                        
                   

# look at DB_none

# recode DB_none to account for triggers not being included

table(dat$DB_none)

dat$DB_none[dat$DB_inhaler == 0 & dat$DB_maintenance == 0 & dat$DB_adherence == 0 &  
              dat$DB_PAAP == 0 & dat$DB_comm_FU_2_days == 0 & 
              dat$DB_spec_review_4_weeks == 0 &
              (dat$DB_smoke == 0 | is.na(dat$DB_smoke))] <- 1

table(dat$DB_none)
# 
# dat %>% filter(DB_inhaler == 0 & DB_maintenance == 0 & DB_adherence == 0 &  
#                DB_PAAP == 0 & DB_comm_FU_2_days == 0 & DB_spec_review_4_weeks == 0 &
#                  (DB_smoke == 0 | is.na(DB_smoke))) %>%
#   select(DB_none, DB_triggers) %>% table()
# 
# dat %>% filter(DB_inhaler == 0 & DB_maintenance == 0 & DB_adherence == 0 &  
#                  DB_PAAP == 0 & DB_comm_FU_2_days == 0 & DB_spec_review_4_weeks == 0 &
#                  (DB_smoke == 0 | is.na(DB_smoke))) %>%
#   filter(DB_triggers == 0 & DB_none == 0) %>% select(starts_with("DB")) 


dat <- dat %>% select(all_of(unique(c(ID_vars, cluster_var, covars, covars_for_smoke, 
                                      covars_for_all_elements_included, independent_vars, 
                                      independent_vars_smoke, outcome_vars, table1_explan, table1_dep))))


# label all the variables

dat <- dat %>% mutate(# patient_ID = ff_label(, ""),                      
                      # hosp_code = ff_label(, ""),
                      # trust_code = ff_label(, ""),
                      # country = ff_label(, ""),
                      # hosp_name = ff_label(, ""),
                      # trust_name = ff_label(, ""),
                      # ICS = ff_label(, ""),
                      # region = ff_label(, ""),
                      IMD_quintile = ff_label(IMD_quintile, "Index of Multiple Deprivation (IMD) quintile (1 = most deprived, 5 = least deprived)"),
                      age = ff_label(age, "Age (years)"),
                      gender = ff_label(gender, "Gender"),
                      smoke_status = ff_label(smoke_status, "Smoking status"),
                      inhaled_steroids_dis = ff_label(inhaled_steroids_dis, "Patient prescribed inhaled steroids at discharge"),
                      oral_steroids_dis = ff_label(oral_steroids_dis, "Patient prescribed oral steroids for at least 5 days during admission or at discharge"),
                      oral_steroids_rescue_history = ff_label(oral_steroids_rescue_history, "Patient prescribed three or more courses of rescue/emergency oral steroids in the 12 months prior to admission"),
                      asthma_sev = ff_label(asthma_sev, "Asthma attack severity"),
                      CCIweightedcat = ff_label(CCIweightedcat, "Charlson Comorbidity Index (CCI)"),
                      smoke_status_with_cessation_DB = ff_label(smoke_status_with_cessation_DB, "Smoking status"),
                      BPT2024 = ff_label(BPT2024, "Patient meets Best Practice Tariff (BPT) quality of care"),
                      RSR_BPT_bin = ff_label(RSR_BPT_bin, "Patient reviewed by a respiratory specialist within 24 hours of admission"),
                      RSR_bin = ff_label(RSR_bin, "Patient reviewed by a respiratory specialist during admission"),
                      DB_BPT = ff_label(DB_BPT, "Patient discharged with a discharge bundle including all elements required to meet the BPT"),
                      discharge_bundle_yes_no = ff_label(discharge_bundle_yes_no, "Patient marked as having received a discharge bundle"),
                      DB_A4_meds = ff_label(DB_A4_meds, "Patient received a review of their asthma medicine, adherence, and inhaler technique (i.e. the medicine section of the 'Asthma 4' discharge bundle element)"),
                      DB_inhaler = ff_label(DB_inhaler, "Patient inhaler technique checked (discharge bundle element)"),
                      DB_maintenance = ff_label(DB_maintenance, "Patient inhaler medication reviewed (discharge bundle element)"),
                      DB_adherence = ff_label(DB_adherence, "Patient medication adherence checked (discharge bundle element)"),
                      DB_PAAP = ff_label(DB_PAAP, "Patient received a Personalised Asthma Action Plan (PAAP) or patient's PAAP was reviewed (discharge bundle element)"),
                      DB_FU_any = ff_label(DB_FU_any, "Community follow-up requested within 2 working days or specialist review requested within 4 weeks (discharge bundle element)"),
                      DB_comm_FU_2_days = ff_label(DB_comm_FU_2_days, "Community follow-up requested within 2 working days (discharge bundle element)"),
                      DB_spec_review_4_weeks = ff_label(DB_spec_review_4_weeks, "Specialist review requested within 4 weeks (discharge bundle element)"),
                      DB_none = ff_label(DB_none, "Patient did not receive any discharge bundle elements"),
                      DB_smoke = ff_label(DB_smoke, "Patient received smoking cessation advice (discharge bundle element)"),
                      read30 = ff_label(read30, "Patient readmitted for any reason within 30 days of discharge"),
                      asthma_read30 = ff_label(asthma_read30, "Patient readmitted with asthma within 30 days of discharge"),
                      read90 = ff_label(read90, "Patient readmitted for any reason within 90 days of discharge"),
                      asthma_read90 = ff_label(asthma_read90, "Patient readmitted with asthma within 90 days of discharge"),
                      age_orig = ff_label(age_orig, "Age (years)"),
                      read30_categorical = ff_label(read30_categorical, "Patient readmitted for any reason within 30 days of discharge"),
                      asthma_read30_categorical = ff_label(asthma_read30_categorical, "Patient readmitted with asthma within 30 days of discharge"),
                      read90_categorical = ff_label(read90_categorical, "Patient readmitted for any reason within 90 days of discharge"),
                      asthma_read90_categorical = ff_label(asthma_read90_categorical, "Patient readmitted with asthma within 90 days of discharge"),
                      BPT2024_cat = ff_label(BPT2024_cat, "Patient meets Best Practice Tariff (BPT) quality of care"))

# keep these names to save them for use in other tables
keep_names <- data.frame(clean_name = extract_variable_label(dat)) %>% 
  rownames_to_column()  %>% rename(varname = rowname)

#=============================================#

# Table 1 !!!!!!



pop_table <-  dat %>% mutate(across(all_of(c(independent_vars[c(-1, -5)], independent_vars_smoke)),
                             ~fct_recode(factor(., levels = c("1", "0")), `No` = "0", `Yes` = "1"))) %>%
  mutate(discharge_bundle_yes_no = factor(discharge_bundle_yes_no, levels = c("Yes", "No"))) %>%
  ff_relabel_df(dat) %>% # used to relabel after tidyverse functions remove the labels
  summary_factorlist(., dependent = table1_dep,
                       explanatory = table1_explan, cont_nonpara = TRUE, total_col = TRUE,
                       add_col_totals = TRUE, p = TRUE)  

pop_table$label


# 


# write.csv(pop_table, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/Tidy results/Table1.csv",
#           row.names = FALSE)

# find all unique combinations of discharge bundle elements

table(dat$discharge_bundle_yes_no)

dat %>% filter(discharge_bundle_yes_no == "Yes") %>% select(starts_with("DB")) %>% summary()
dat %>% filter(discharge_bundle_yes_no == "No") %>% select(starts_with("DB")) %>% summary()
dat %>% select(DB_inhaler, DB_maintenance, DB_adherence, DB_PAAP, DB_comm_FU_2_days, DB_spec_review_4_weeks, DB_none) %>%
  unique() %>% nrow()

# 65 unique combinations...

dat %>% group_by(DB_inhaler, DB_maintenance, DB_adherence, DB_PAAP, DB_comm_FU_2_days, DB_spec_review_4_weeks, DB_none) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>% print(n=1000)

rough_heatmap <- dat %>% group_by(DB_inhaler, DB_maintenance, DB_adherence, DB_PAAP, DB_comm_FU_2_days, DB_spec_review_4_weeks, DB_none) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>% head(20)

heatmap_other_colnames <- dat %>% group_by(DB_inhaler, DB_maintenance, DB_adherence, DB_PAAP, DB_comm_FU_2_days, DB_spec_review_4_weeks, DB_none) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>% filter(n < 100) %>% colnames()  

heatmap_other_sum <- dat %>% group_by(DB_inhaler, DB_maintenance, DB_adherence, DB_PAAP, DB_comm_FU_2_days, DB_spec_review_4_weeks, DB_none) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>% filter(n < 100) %>% pull(n) %>% sum()  

heatmap2 <- data.frame(x1 = NA, x2 = NA, x3 = NA, x4 = NA, x5 = NA, x6 = NA, x7 = NA, n = heatmap_other_sum) 
colnames(heatmap2) <- heatmap_other_colnames

heatmap3 <- data.frame(x1 = sum(dat$DB_inhaler), x2 = sum(dat$DB_maintenance), x3 = sum(dat$DB_adherence), x4 = sum(dat$DB_PAAP),
                       x5 = sum(dat$DB_comm_FU_2_days), x6 = sum(dat$DB_spec_review_4_weeks), x7 = sum(dat$DB_none), n = nrow(dat))
colnames(heatmap3) <- heatmap_other_colnames

rough_heatmap <- bind_rows(rough_heatmap, heatmap2, heatmap3)

rough_heatmap$Combination <- c(rep(NA, 20), "Other combinations", "Total")

rough_heatmap <- rough_heatmap %>% relocate(Combination, .before = DB_inhaler)

clean_names <- keep_names %>% filter(varname %in% c("DB_inhaler", "DB_maintenance", "DB_adherence", "DB_PAAP", "DB_comm_FU_2_days", 
                                     "DB_spec_review_4_weeks", "DB_none")) %>% pull(clean_name)

clean_names # check name order matches up with heat map

colnames(rough_heatmap)[2:8] <- clean_names

rough_heatmap %>% print(n=100)

# write.csv(rough_heatmap, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/rough_heatmap_data.csv",
#           row.names = FALSE)

#========================================================================#


# And let's make a data frame for all the analyses we have to do


nrow(dat)
dat %>% select(hosp_code) %>% group_by(hosp_code) %>% unique() %>% nrow()

# unadjusted analyses

analyses_df_unadj <- data.frame(outcome_vars = rep(outcome_vars, each = length(independent_vars) + 1),
                                independent_vars = rep(c(independent_vars, independent_vars_smoke), 4),
                                varcheck = NA, est = NA, lo = NA, hi = NA, pval = NA)


model_output_list_unadj <- vector(nrow(analyses_df_unadj), mode = "list")
model_output_list_unadj_glm <- vector(nrow(analyses_df_unadj), mode = "list")



# print(tidyoutput(x = model_output_list_unadj[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2))
# analyses_df_unadj[i, 3:7] <- tidyoutput(x = model_output_list_unadj[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2)[2, ]

analyses_df_unadj_profile <- analyses_df_unadj
analyses_df_unadj_Wald <- analyses_df_unadj
analyses_df_unadj_glm <- analyses_df_unadj

# unhash this or the readRDS files

# FOR THE SENSITIVITY ANALYSIS, ONLY DOING RSR 24 HOURS AND BPT

# for (i in 1:3) {
# for (i in 1:nrow(analyses_df_unadj)) {
# for (i in c(14, 29, 44, 59)) { # for DB_none
for (i in c(1, 2, 16, 17, 31, 32, 46, 47)) { # for BPT/RSR 28 HOUR SENSITIVTY


  model_output_list_unadj[[i]] <-  glmer(formula(paste(analyses_df_unadj$outcome_vars[i]," ~ ",
                                               analyses_df_unadj$independent_vars[i]," + ",
                                               "(1 | ", cluster_var, ")", collapse = "")),
      family = binomial(link = "logit"), data = dat)

  model_output_list_unadj_glm[[i]] <-  glm(formula(paste(analyses_df_unadj$outcome_vars[i]," ~ ",
                                                     analyses_df_unadj$independent_vars[i], collapse = "")),
                                       family = binomial(link = "logit"), data = dat)

  try({  print(tidyoutput(x = model_output_list_unadj[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2)) })
  try({  analyses_df_unadj_Wald[i, 3:7] <- tidyoutput(x = model_output_list_unadj[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2)[2, ] })
  try({  analyses_df_unadj_profile[i, 3:7] <- tidyoutput(x = model_output_list_unadj[[i]], meth = "profile", MEM = TRUE, pval = TRUE, roundno = 2)[2, ]  })
  try({  analyses_df_unadj_glm[i, 3:7] <- tidyoutput(x = model_output_list_unadj_glm[[i]], meth = "profile", MEM = FALSE, pval = TRUE, roundno = 2)[2, ] })


}

# none_unadj <- analyses_df_unadj_Wald
BPT_RSR_28hour_unadj <- analyses_df_unadj_Wald

saveRDS(BPT_RSR_28hour_unadj, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/BPT_RSR_28hour_analyses_df_unadj_Wald.RDS")


#
#
# saveRDS(analyses_df_unadj_Wald, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/analyses_df_unadj_Wald.RDS")
# saveRDS(analyses_df_unadj_profile, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/analyses_df_profile.RDS")
# saveRDS(analyses_df_unadj_glm, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/analyses_df_glm.RDS")
# saveRDS(model_output_list_unadj, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/model_output_list_unadj.RDS")
# saveRDS(model_output_list_unadj_glm, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/model_output_list_unadj_glm.RDS")

# analyses_df_unadj_Wald <- readRDS("G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/analyses_df_unadj_Wald.RDS")
# analyses_df_unadj_profile <- readRDS("G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/analyses_df_profile.RDS")
# analyses_df_unadj_glm <- readRDS("G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/analyses_df_glm.RDS")
# model_output_list_unadj <- readRDS("G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/model_output_list_unadj.RDS")
# model_output_list_unadj_glm <- readRDS("G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/model_output_list_unadj_glm.RDS")





analyses_df_unadj_Wald <- analyses_df_unadj_Wald %>% clean_up() %>% select(outcome_vars, independent_vars, coefficients, pval)
analyses_df_unadj_profile <- analyses_df_unadj_profile %>% clean_up() %>% select(outcome_vars, independent_vars, coefficients, pval)
analyses_df_unadj_glm <- analyses_df_unadj_glm %>% clean_up() %>% select(outcome_vars, independent_vars, coefficients, pval)

analyses_df_unadj_compared <- analyses_df_unadj_Wald[ ,-4] %>% rename(Wald_coefficients = coefficients) %>%
  left_join(., rename(analyses_df_unadj_profile[ ,-4], profile_coefficients = coefficients), 
            by = c("outcome_vars", "independent_vars")) %>%
  left_join(., rename(analyses_df_unadj_glm[ ,-4], glm_coefficients = coefficients), 
            by = c("outcome_vars", "independent_vars"))


  
  # No major issues here. Use Wald confidence intervals! 
# This code uses the cleaned names used earlier when labelling the variables  
    
clean_unadj_analyses <- left_join(analyses_df_unadj_Wald, keep_names, by = c("independent_vars" = "varname")) %>%
  rename(clean_independent_variable_name = clean_name) %>% 
  left_join(., keep_names, by = c("outcome_vars" = "varname")) %>%
  rename(clean_output_variable_name = clean_name) %>% 
  select(outcome_vars, clean_output_variable_name, independent_vars, clean_independent_variable_name, 
         coefficients, pval)
#  select(clean_output_variable_name, clean_independent_variable_name, 
#         coefficients, pval)

# write.csv(clean_unadj_analyses, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/Tidy results/Supplemental_table_1.csv",
#           row.names = FALSE)



# Adjusted analyses

# multi-level
analyses_df_adj <- data.frame(outcome_vars = rep(outcome_vars, each = length(independent_vars)),
                                independent_vars = rep(independent_vars, 4),
                                varcheck = NA, est = NA, lo = NA, hi = NA, pval = NA)

# normal
analyses_df_adj_glm <- analyses_df_adj

# create the lists for model results
model_output_list_adj <- vector(nrow(analyses_df_adj), mode = "list")
model_output_list_adj_glm <- vector(nrow(analyses_df_adj), mode = "list")


 # for (i in 1:nrow(analyses_df_adj)) {
for (i in c(1, 2, 15, 16, 29, 30, 43, 44)) { # for BPT/RSR 28 HOUR SENSITIVTY
  
  model_output_list_adj[[i]] <-  glmer(formula(paste(analyses_df_adj$outcome_vars[i]," ~ ",
                                                     analyses_df_adj$independent_vars[i]," + ",
                                                     paste(covars, collapse = "+"),
                                                     "+ (1 | ", cluster_var, ")", collapse = "")),
                                       family = binomial(link = "logit"), data = dat)

  model_output_list_adj_glm[[i]] <-  glm(formula(paste(analyses_df_adj$outcome_vars[i]," ~ ",
                                                   analyses_df_adj$independent_vars[i]," + ",
                                                   paste(covars, collapse = "+"),
                                                   collapse = "")),
                                     family = binomial(link = "logit"), data = dat)


  print(summary(model_output_list_adj[[i]]))
  print(tidyoutput(x = model_output_list_adj[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2))
    analyses_df_adj[i, 3:7] <- tidyoutput(x = model_output_list_adj[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2)[2, ]
  analyses_df_adj_glm[i, 3:7] <- tidyoutput(x = model_output_list_adj_glm[[i]], meth = "Wald", MEM = FALSE, pval = TRUE, roundno = 2)[2, ]

}

BPT_RSR_28hour_adj <- analyses_df_adj

saveRDS(BPT_RSR_28hour_adj, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/BPT_RSR_28hour_analyses_df_adj.RDS")
summary(model_output_list_adj[[15]])
summary(model_output_list_adj_glm[[15]])

saveRDS(model_output_list_adj, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/model_output_list_adj_BPT_RSR_28hour.RDS")
saveRDS(model_output_list_adj_glm, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/model_output_list_adj_glm_BPT_RSR_28hour.RDS")

clean_adj_analyses <- BPT_RSR_28hour_adj %>% filter(!is.na(est))

clean_adj_analyses <- clean_adj_analyses %>% clean_up() %>% select(outcome_vars, independent_vars, coefficients, pval)
   

clean_adj_analyses <- left_join(clean_adj_analyses, keep_names, by = c("independent_vars" = "varname")) %>%
  rename(clean_independent_variable_name = clean_name) %>% 
  left_join(., keep_names, by = c("outcome_vars" = "varname")) %>%
  rename(clean_output_variable_name = clean_name) %>% 
  select(outcome_vars, clean_output_variable_name, independent_vars, clean_independent_variable_name, 
         coefficients, pval)
#  select(clean_output_variable_name, clean_independent_variable_name, 
#         coefficients, pval)

write.csv(clean_adj_analyses, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/Tidy results/Clean_adjusted_results_table_BPT_RSR_28hours.csv",
          row.names = FALSE)

# 
# # we probably want to keep the results with the confounders for the BPT2024 analysis.
# # We can do that here:
# 
# skeleton <- summary_factorlist(dat,
#                            dependent = "read30_categorical",
#                            explanatory = c(independent_vars[1], covars), fit_id = TRUE)
# 
# # the BPT2024 variable hasn't come across properly so need to do this:
# skeleton[2, 5] <- "BPT2024"
# skeleton[1:2, 2] <- c("No", "Yes")
# skeleton <- skeleton %>% select(label, levels, fit_id, index)
# 
# # and then we join them all up
# 
# #30 day read
# BPT1 <- fit2df(model_output_list_adj[[1]], digits = c(2, 2, 0), confint_sep = " to ") 
# 
# # don't want the p-values
# BPT1$OR <- str_remove(BPT1$OR, ", p<1")
# BPT1$OR <- str_remove(BPT1$OR, ", p=1")
# 
# BPT1 <- BPT1 %>% rename(`OR (95% CI) for 30-day all-cause readmission` = "OR")
# 
# # 30 day asthma read
# BPT2 <- fit2df(model_output_list_adj[[15]], digits = c(2, 2, 0), confint_sep = " to ") 
# 
# # don't want the p-values
# BPT2$OR <- str_remove(BPT2$OR, ", p<1")
# BPT2$OR <- str_remove(BPT2$OR, ", p=1")
# 
# BPT2 <- BPT2 %>% rename(`OR (95% CI) for 30-day asthma readmission` = "OR")
# 
# # 90 day read
# BPT3 <- fit2df(model_output_list_adj[[29]], digits = c(2, 2, 0), confint_sep = " to ") 
# 
# # don't want the p-values
# BPT3$OR <- str_remove(BPT3$OR, ", p<1")
# BPT3$OR <- str_remove(BPT3$OR, ", p=1")
# 
# BPT3 <- BPT3 %>% rename(`OR (95% CI) for 90-day all-cause readmission` = "OR")
# 
# # 90 day asthma read
# BPT4 <- fit2df(model_output_list_adj[[43]], digits = c(2, 2, 0), confint_sep = " to ") 
# 
# # don't want the p-values
# BPT4$OR <- str_remove(BPT4$OR, ", p<1")
# BPT4$OR <- str_remove(BPT4$OR, ", p=1")
# 
# BPT4 <- BPT4 %>% rename(`OR (95% CI) for 90-day asthma readmission` = "OR")
# 
# analyses_df_adj # this is used for the dot and whisker plot
# # test <- left_join(skeleton, BPT1, by = c("fit_id" = "explanatory"))
# 
# 
# clean_BPT <- ff_merge(skeleton, BPT2) %>% ff_merge(., BPT1) %>%
#   ff_merge(., BPT4) %>% ff_merge(., BPT3) %>% select(-index, -fit_id)
# 
# # write.csv(clean_BPT, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/Tidy results/BPT2024_full_results_read.csv",
# #           row.names = FALSE)
# 
# custom_labels = c(
# "BPT2024" = "Best Practice Tariff criteria met",
# "RSR_BPT_bin" =  "Respiratory specialist review within 24 hours of admission",
# "RSR_bin" =  "Respiratory specialist review during admission",
# "DB_BPT" =  "Discharge bundle elements meet BPT criteria",
# "discharge_bundle_yes_no" =  "Patient received a discharge bundle",
# "DB_A4_meds" =  "Review of asthma medicine, adherence, and inhaler technique\n('Asthma 4' bundle, section 1)",
# "DB_inhaler" =  "Inhaler technique checked",
# "DB_maintenance" =  "Inhaler medication reviewed",
# "DB_adherence" =  "Medication adherence checked",
# "DB_PAAP" =  "Receipt or review of Personalised Asthma Action Plan (PAAP)",
# "DB_FU_any" =  "Community follow-up requested within 2 working days or \nspecialist review requested within 4 weeks",
# "DB_comm_FU_2_days" =  "Community follow-up requested within 2 working days",
# "DB_spec_review_4_weeks" = "Specialist review requested within 4 weeks",
# "DB_none" = "No discharge bundle elements received",
# "DB_smoke" = "Smoking cessation advice received")
# # "read30" = "Patient readmitted for any reason within 30 days of discharge"
# # asthma_read30 = ff_label(asthma_read30, "Patient readmitted with asthma within 30 days of discharge",
# # read90 = ff_label(read90, "Patient readmitted for any reason within 90 days of discharge",
# # asthma_read90 = ff_label(asthma_read90, "Patient readmitted with asthma within 90 days of discharge")
# 
# str(analyses_df_adj)
# # analyses_df_adj_save <- analyses_df_adj
# 
# analyses_df_adj <- analyses_df_adj_save
# 
# analyses_df_adj <- analyses_df_adj %>% 
#   mutate(outcome_vars = fct_relevel(outcome_vars, "read30", "asthma_read30",
#                                                "read90", "asthma_read90")) %>%
#   mutate(outcome_vars = fct_recode(outcome_vars,
#                                    `30-day all-cause\nreadmission` = "read30", 
#                                    `30-day asthma\nreadmission` = "asthma_read30", 
#                                    `90-day all-cause\nreadmission` = "read90", 
#                                    `90-day asthma\nreadmission` = "asthma_read90")) 
# 
# analyses_df_adj <- analyses_df_adj %>% 
#   mutate(independent_vars = fct_relevel(independent_vars, 
#                                         rev(unique(analyses_df_adj$independent_vars))))
# 
# # Plot these results... 
# p <- ggplot(analyses_df_adj, mapping = aes(x = est, y = independent_vars, group = outcome_vars)) +
#   geom_point(size = 2) +
#   geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0) +
#   facet_wrap(~outcome_vars, nrow = 1) + # +
# #  facet_grid(rows = vars(independent_vars), cols = vars(outcome_vars), scales = "free_y", space = "free_y") +
# #  facet_grid(rows = vars(outcome_vars), cols = vars(outcome_vars), scales = "free_y", space = "free_y") +
#   coord_cartesian(xlim = c(0.5, 1.5)) +
#   geom_vline(xintercept = 1, linetype = "dashed", size = 0.5) +
#   theme_minimal() +
#   scale_y_discrete(labels = custom_labels) +
#   labs(x = "Adjusted odds ratio", y = NULL) +
#   theme(strip.placement = "outside",
#         strip.text.y = element_text(angle = 0, hjust = 1),
#         axis.ticks.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.y.right = element_blank(),
#         panel.spacing = unit(1, "lines")) 
# #        aspect.ratio = 1/2) 
#  
# p
# 
# ggsave("G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/Tidy results/adj_odds_ratios.pdf",
#        p, width = 12, height = 7, dpi = 300)
# 
# 
# 
# # Now, we identify what within the discharge bundle are the most important elements:
# 
# # need to combine the smoking status variable. But, we also need to think about what we are trying to show.
# # We have all the individual elements now. Now we move onto: 
# # - individual elements adjusted for other elements but not double counted
# # the overall BPT not adjusted for elements that make it up
# 
# 
# # individual elements not double counted:
# 
# independent_vars_for_comb 
# covars_for_all_elements_included
# 
# 
# 
# model_all_DB_indiv <- vector(4, mode = "list")
# model_all_DB_indiv_glm <- vector(4, mode = "list")
# 
# for (i in 1:length(outcome_vars)) {
#   
#   model_all_DB_indiv[[i]] <-  glmer(formula(paste(outcome_vars[i]," ~ ",
#                                              paste(independent_vars_for_comb, collapse = "+"),
#                                              " + ", paste(covars_for_all_elements_included, collapse = "+"),
#                                              " + (1 | ", cluster_var, ")", collapse = "")),
#                                family = binomial(link = "logit"), data = dat)
#   
#   print(tidyoutput(x = model_all_DB_indiv[[i]], meth = "Wald", MEM = TRUE, pval = TRUE, roundno = 2))
# 
#   model_all_DB_indiv_glm[[i]] <-  glm(formula(paste(outcome_vars[i]," ~ ",
#                                                   paste(independent_vars_for_comb, collapse = "+"),
#                                                   " + ", paste(covars_for_all_elements_included, collapse = "+"),
#                                                   collapse = "")),
#                                     family = binomial(link = "logit"), data = dat)
#   
# }
# 
# 
# # saveRDS(model_all_DB_indiv, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/model_all_DB_indiv.RDS")
# # saveRDS(model_all_DB_indiv_glm, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/model_all_DB_indiv_glm.RDS")
# 
# summary(model_all_DB_indiv[[1]])
# summary(model_all_DB_indiv[[2]])
# summary(model_all_DB_indiv[[3]])
# summary(model_all_DB_indiv[[4]])
# 
# # simulateResiduals(model_all_DB_indiv[[1]], plot = TRUE) # check residuals
# 
# 
# 
# # we are going to use the theoretical (pi/3) variance of the marginal model.
# # This means we are going to take the fixed variance of binomial models that exists
# # independently of the data collected, and we are going to use the marginal
# # R2 (predicting using the fixed effects) instead of the conditional R2 
# # (adding in hospital to the prediction)
# 
# # We are using the domir package for dominance analysis
# # it's a pretty complicated package but it means it's very flexible
# 
# # first thing we need to do is create a function that will extract
# # the R squared and return it in a list along with its name
# # Uses the MuMin package to extract it
# 
# extract.R2 <- function(x) {r.squaredGLMM(x)[1,1] %>% data.frame("R.square" = .) }
# 
# # then we need to create a wrapper function to use within the domin function
# # this will take the first argument of the domir function as its formula
# # everything else is supplied as normal
# 
# wrapper <- function(formula) {
#   glmer_model <- glmer(formula = formula, 
#                        family = binomial(link = "logit"), 
#                        data = dat)
# }
# 
# # and here is the function we use
# # each one of these is going to need to run the model 64 times.
# # Oh well! We have the weekend to run it...
# 
# # # this is the way to do it - don't use 'all', use 'consmodel'
# # dom_test <- domin(read30 ~ DB_inhaler + DB_maintenance, # random effects must be supplied in consmodel section, not here
# #       wrapper, # glmer wrapper function given above
# #       fitstat = list("extract.R2", "R.square"), # extraction function plus name of what is to be extracted
# #       consmodel = c("CCIweightedcat", "(1 | hosp_code)")) # specify random effects in this constant model section
# 
# 
# dom_read30 <- domin(read30 ~ DB_inhaler + DB_maintenance + DB_adherence + # formula. when using consmodel, take out covars.
#                       DB_PAAP + DB_comm_FU_2_days + DB_spec_review_4_weeks, # random effects MUST be supplied in consmodel section
#                     wrapper, # glmer wrapper function given above
#                     fitstat = list("extract.R2", "R.square"), # extraction function plus name of what is to be extracted
#                     consmodel = c("IMD_quintile", "age", "gender", "smoke_status_with_cessation_DB",
#                             "inhaled_steroids_dis", "oral_steroids_dis", 
#                             "oral_steroids_rescue_history", "asthma_sev", "CCIweightedcat", "(1 | hosp_code)")) 
# # specify the covavariates you don't want to do the dominance analysis for, and the random effects,
# # in the constant model (consmodel) section
# 
# saveRDS(dom_read30, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/dom_read30.RDS")
# 
# dom_read90 <- domin(read90 ~ DB_inhaler + DB_maintenance + DB_adherence + # formula. when using consmodel, take out covars.
#                       DB_PAAP + DB_comm_FU_2_days + DB_spec_review_4_weeks, # random effects MUST be supplied in consmodel section
#                     wrapper, # glmer wrapper function given above
#                     fitstat = list("extract.R2", "R.square"), # extraction function plus name of what is to be extracted
#                     consmodel = c("IMD_quintile", "age", "gender", "smoke_status_with_cessation_DB",
#                                   "inhaled_steroids_dis", "oral_steroids_dis", 
#                                   "oral_steroids_rescue_history", "asthma_sev", "CCIweightedcat", "(1 | hosp_code)")) 
# # specify the covavariates you don't want to do the dominance analysis for, and the random effects,
# # in the constant model (consmodel) section
# 
# saveRDS(dom_read90, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/dom_read90.RDS")
# 
# dom_asthma_read30 <- domin(asthma_read30 ~ DB_inhaler + DB_maintenance + DB_adherence + # formula. when using consmodel, take out covars.
#                              DB_PAAP + DB_comm_FU_2_days + DB_spec_review_4_weeks, # random effects MUST be supplied in consmodel section
#                            wrapper, # glmer wrapper function given above
#                     fitstat = list("extract.R2", "R.square"), # extraction function plus name of what is to be extracted
#                     consmodel = c("IMD_quintile", "age", "gender", "smoke_status_with_cessation_DB",
#                                   "inhaled_steroids_dis", "oral_steroids_dis", 
#                                   "oral_steroids_rescue_history", "asthma_sev", "CCIweightedcat", "(1 | hosp_code)")) 
# # specify the covavariates you don't want to do the dominance analysis for, and the random effects,
# # in the constant model (consmodel) section
# 
# saveRDS(dom_asthma_read30, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/dom_asthma_read30.RDS")
# 
# dom_asthma_read90 <- domin(asthma_read90 ~ DB_inhaler + DB_maintenance + DB_adherence + # formula. when using consmodel, take out covars.
#                              DB_PAAP + DB_comm_FU_2_days + DB_spec_review_4_weeks, # random effects MUST be supplied in consmodel section
#                            wrapper, # glmer wrapper function given above
#                     fitstat = list("extract.R2", "R.square"), # extraction function plus name of what is to be extracted
#                     consmodel = c("IMD_quintile", "age", "gender", "smoke_status_with_cessation_DB",
#                                   "inhaled_steroids_dis", "oral_steroids_dis", 
#                                   "oral_steroids_rescue_history", "asthma_sev", "CCIweightedcat", "(1 | hosp_code)")) 
# # specify the covavariates you don't want to do the dominance analysis for, and the random effects,
# # in the constant model (consmodel) section
# 
# saveRDS(dom_asthma_read90, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/dom_asthma_read90.RDS")
# 
# 
# 
# sort(dom_asthma_read30$Conditional_Dominance[ ,1], decreasing = TRUE)
# sort(dom_asthma_read30$Conditional_Dominance[ ,2], decreasing = TRUE)
# sort(dom_asthma_read30$Conditional_Dominance[ ,3], decreasing = TRUE)
# sort(dom_asthma_read30$Conditional_Dominance[ ,4], decreasing = TRUE)
# sort(dom_asthma_read30$Conditional_Dominance[ ,5], decreasing = TRUE)
# sort(dom_asthma_read30$Conditional_Dominance[ ,6], decreasing = TRUE)
# 
# sort(dom_asthma_read90$Conditional_Dominance[ ,1], decreasing = TRUE)
# sort(dom_asthma_read90$Conditional_Dominance[ ,2], decreasing = TRUE)
# sort(dom_asthma_read90$Conditional_Dominance[ ,3], decreasing = TRUE)
# sort(dom_asthma_read90$Conditional_Dominance[ ,4], decreasing = TRUE)
# sort(dom_asthma_read90$Conditional_Dominance[ ,5], decreasing = TRUE)
# sort(dom_asthma_read90$Conditional_Dominance[ ,6], decreasing = TRUE)
