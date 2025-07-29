# try and calculate the absolute risk reduction etc.

library(tidyverse)
library(lme4)

model_output_list_adj <- readRDS("G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/model_output_list_adj.RDS")

summary(model_output_list_adj[[19]]) # this is the discharge bundle one for asthma read 30
summary(model_output_list_adj[[47]]) # this is the discharge bundle one for asthma read 90

model <- model_output_list_adj[[19]] # run for whichever one you want
# model <- model_output_list_adj[[47]] # run for whichever one you want

fulld <- model@frame
nrow(fulld)

d_exposed <- fulld
d_exposed$discharge_bundle_yes_no <- "Yes" 
d_unexposed <- fulld
d_exposed$discharge_bundle_yes_no <- "No" 




 # 

dat <- data.frame(
  DB = fulld$discharge_bundle_yes_no,
  orig = predict(model, fulld, type = "response"),
  exp = predict(model, d_exposed, type = "response"),
  unexp = predict(model, d_unexposed, type = "response"))

summary(dat)

mean_exp <- mean(dat$exp)
mean_unexp <- mean(dat$unexp)

risk_ratio <- mean_unexp/mean_exp
risk_reduction <- mean_exp - mean_unexp

NNT <- 1/risk_reduction
NNT

colnames(fulld)

model <- model_output_list_adj[[27]] # run for whichever one you want
