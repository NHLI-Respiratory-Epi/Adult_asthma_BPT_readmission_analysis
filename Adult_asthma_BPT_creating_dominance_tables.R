library(tidyverse)
library(domir)


dom_read30 <- readRDS("G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/dom_read30.RDS")

dom_read90 <- readRDS("G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/dom_read90.RDS")

dom_asthma_read30 <- readRDS("G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/dom_asthma_read30.RDS")

dom_asthma_read90 <- readRDS("G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/dom_asthma_read90.RDS")


dom_read30
dom_asthma_read30

str(dom_read30)

dom_read30$Standardized
dom_read30$Ranks
dom_read30$Conditional_Dominance

dom30 <- data.frame(dom_read30$Conditional_Dominance)
dom30$`General dominance` <- dom_read30$General_Dominance
dom30$`Standardised dominance` <- dom_read30$Standardized
dom30$`General dominance rank` <- dom_read30$Ranks
dom30 <- round(dom30, 5)
dom30$model <- "All-cause 30-day readmission"
dom30 <- rownames_to_column(dom30, var = "Discharge bundle element")
dom30 <- dom30 %>% select(model, `Discharge bundle element`, `General dominance rank`, `General dominance`,
                          `Standardised dominance`, IVs_1:IVs_6)

dom90 <- data.frame(dom_read90$Conditional_Dominance)
dom90$`General dominance` <- dom_read90$General_Dominance
dom90$`Standardised dominance` <- dom_read90$Standardized
dom90$`General dominance rank` <- dom_read90$Ranks
dom90 <- round(dom90, 5)
dom90$model <- "All-cause 90-day readmission"
dom90 <- rownames_to_column(dom90, var = "Discharge bundle element")
dom90 <- dom90 %>% select(model, `Discharge bundle element`, `General dominance rank`, `General dominance`,
                          `Standardised dominance`, IVs_1:IVs_6)

dom30_asthma <- data.frame(dom_asthma_read30$Conditional_Dominance)
dom30_asthma$`General dominance` <- dom_asthma_read30$General_Dominance
dom30_asthma$`Standardised dominance` <- dom_asthma_read30$Standardized
dom30_asthma$`General dominance rank` <- dom_asthma_read30$Ranks
dom30_asthma <- round(dom30_asthma, 5)
dom30_asthma$model <- "Asthma 30-day readmission"
dom30_asthma <- rownames_to_column(dom30_asthma, var = "Discharge bundle element")
dom30_asthma <- dom30_asthma %>% select(model, `Discharge bundle element`, `General dominance rank`, `General dominance`,
                          `Standardised dominance`, IVs_1:IVs_6)

dom90_asthma <- data.frame(dom_asthma_read90$Conditional_Dominance)
dom90_asthma$`General dominance` <- dom_asthma_read90$General_Dominance
dom90_asthma$`Standardised dominance` <- dom_asthma_read90$Standardized
dom90_asthma$`General dominance rank` <- dom_asthma_read90$Ranks
dom90_asthma <- round(dom90_asthma, 5)
dom90_asthma$model <- "Asthma 90-day readmission"
dom90_asthma <- rownames_to_column(dom90_asthma, var = "Discharge bundle element")
dom90_asthma <- dom90_asthma %>% select(model, `Discharge bundle element`, `General dominance rank`, `General dominance`,
                          `Standardised dominance`, IVs_1:IVs_6)

dom_all <- bind_rows(dom30_asthma, dom30, dom90_asthma, dom90)

write.csv(dom_all, "G:/Alex Harley/Audit_2023_onwards/2022-2023/Research/Adult_asthma_BPT/Analysis/Output/Analysis for paper/Tidy results/dominance_analysis.csv",
          row.names = FALSE)

row.names(dom)

str(summary(dom_read30))

sort(dom_asthma_read30$Conditional_Dominance[ ,1], decreasing = TRUE)
sort(dom_asthma_read30$Conditional_Dominance[ ,2], decreasing = TRUE)
sort(dom_asthma_read30$Conditional_Dominance[ ,3], decreasing = TRUE)
sort(dom_asthma_read30$Conditional_Dominance[ ,4], decreasing = TRUE)
sort(dom_asthma_read30$Conditional_Dominance[ ,5], decreasing = TRUE)
sort(dom_asthma_read30$Conditional_Dominance[ ,6], decreasing = TRUE)

sort(dom_asthma_read90$Conditional_Dominance[ ,1], decreasing = TRUE)
sort(dom_asthma_read90$Conditional_Dominance[ ,2], decreasing = TRUE)
sort(dom_asthma_read90$Conditional_Dominance[ ,3], decreasing = TRUE)
sort(dom_asthma_read90$Conditional_Dominance[ ,4], decreasing = TRUE)
sort(dom_asthma_read90$Conditional_Dominance[ ,5], decreasing = TRUE)
sort(dom_asthma_read90$Conditional_Dominance[ ,6], decreasing = TRUE)