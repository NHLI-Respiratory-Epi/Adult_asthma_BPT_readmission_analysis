# Adult asthma BPT analysis script

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


# start off with index admissions at English hospitals

dat <- dat %>% filter(index_admission_flag == 1) %>% filter(country == "England")

# how many people at this point? 
nrow(dat)

# this many people died as inpatient: 
dat %>% filter(life_status == "Died as inpatient") %>% nrow()


# this many people died as inpatient: 
dat %>% filter(life_status == "Died as inpatient") %>% nrow()

# This many people were transferred to another hospital:
dat %>% filter(discharge_bundle == "Patient transferred to another hospital") %>% nrow()


dat <- dat %>% filter(discharge_bundle != "Patient transferred to another hospital") %>%
  filter(life_status == "Alive") 

# And just keep those who were male or female and relevel gender to remove the other options

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

table(dat$adult_asthma_coded_admission, dat$adult_asthma_coded_admission_including_index)



# Create independent variables:
colnames(dat)

dat$RSR_BPT_bin <- NA
dat$RSR_BPT_bin[dat$RSR_BPT == "<24 hours"] <- 1
dat$RSR_BPT_bin[dat$RSR_BPT == ">=24 hours or no RSR"] <- 0

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
covars <- c("IMD_quintile", "age", "gender", "smoke_status", "inhaled_steroids_dis",                        
            "oral_steroids_dis", "oral_steroids_rescue_history",
            "asthma_sev", "CCIweightedcat")  


for (i in covars) {
  print(i)
    print(summary(dat[[i]]))
    print(levels(dat[[i]]))
  
}

covars_for_smoke <- c("IMD_quintile", "age", "gender", "inhaled_steroids_dis",                        
            "oral_steroids_dis", "oral_steroids_rescue_history",
            "asthma_sev", "CCIweightedcat")  

covars_for_all_elements_included <- c("IMD_quintile", "age", "gender", "smoke_status_with_cessation_DB", "inhaled_steroids_dis",                        
                                      "oral_steroids_dis", "oral_steroids_rescue_history",
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

summary(dat$discharge_bundle)






sex_table1_explan <- c("IMD_quintile", "age_orig",  "smoke_status", "inhaled_steroids_dis",                        
                 "oral_steroids_dis", "oral_steroids_rescue_history",
                 "asthma_sev", "CCIweightedcat", # (this is covars but with original age) 
                 independent_vars, independent_vars_smoke, outcome_vars_table1)




table1_dep <- "BPT2024_cat"                        


# look at DB_none

# recode DB_none to account for triggers not being included

table(dat$DB_none)

dat$DB_none[dat$DB_inhaler == 0 & dat$DB_maintenance == 0 & dat$DB_adherence == 0 &  
              dat$DB_PAAP == 0 & dat$DB_comm_FU_2_days == 0 & 
              dat$DB_spec_review_4_weeks == 0 &
              (dat$DB_smoke == 0 | is.na(dat$DB_smoke))] <- 1

table(dat$DB_none)





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

# Asthma primary first diagnosis vs not table

# # first we do everything we are interested in without p values
# sex_table <-  dat %>% mutate(across(all_of(c(independent_vars[c(-1, -5)], independent_vars_smoke)),
#                              ~fct_recode(factor(., levels = c("1", "0")), `No` = "0", `Yes` = "1"))) %>%
#   mutate(discharge_bundle_yes_no = factor(discharge_bundle_yes_no, levels = c("Yes", "No"))) %>%
#   ff_relabel_df(dat) %>% # used to relabel after tidyverse functions remove the labels
#   summary_factorlist(., dependent = sex_table_dep,
#                        explanatory = sex_table1_explan, cont_nonpara = TRUE, total_col = FALSE,
#                        add_col_totals = TRUE, p = FALSE, fit_id = TRUE)  

summary(dat$gender)

dat$gender <- relevel(dat$gender, ref = "Female")

# then we remove the sex primary codes variable and fit the p values (as it will crash otherwise)
sex_table <-  dat %>% mutate(across(all_of(c(independent_vars[c(-1, -5)], independent_vars_smoke)),
                                      ~fct_recode(factor(., levels = c("1", "0")), `No` = "0", `Yes` = "1"))) %>%
  mutate(discharge_bundle_yes_no = factor(discharge_bundle_yes_no, levels = c("Yes", "No"))) %>%
  ff_relabel_df(dat) %>% # used to relabel after tidyverse functions remove the labels
  summary_factorlist(., dependent = "gender",
                     explanatory = sex_table1_explan, cont_nonpara = TRUE, total_col = FALSE,
                     add_col_totals = TRUE, p = TRUE, fit_id = FALSE)  
# 
# sex_table_ps <- sex_table_ps %>% select(p, fit_id, index)
# 
# sex_table <- left_join(sex_table, sex_table_ps, by = c("fit_id", "index"))
# glimpse(sex_table)
# glimpse(sex_table_ps)
#
# sex_table <- sex_table %>% select(-fit_id, -index)



write.csv(sex_table, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/Tidy results/Supplemental Table 17 women vs men table.csv",
          row.names = FALSE)


