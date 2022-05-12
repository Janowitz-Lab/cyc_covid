# ------------------------------------------------------------
# 
# Script name: cyc_analysis_finngen_cordioli.R
#
# Purpose: run validation analyses for association between PGS for Cystatin C and Covid19-related outcomes in FinnGen 
#
# Author: Mattia Cordioli - Insitute for Molecular Medicine Finland, University of Helsinki, Helsinki, FI
#
# Date created: 2021-11-05
#
# Email: mattia.cordioli@helsinki.fi
#
# ------------------------------------------------------------

rm(list=ls())
library(data.table)
library(dplyr)

# # # Read in infectious disease register and covariates files
covid <- fread('/finngen/library-red/finngen_R7/infectious_disease_register_corona_4.0/data/finngen_R7_infectious_disease_register_corona.txt')
covs <- fread('/finngen/library-red/finngen_R7/phenotype_4.0/data/finngen_R7_cov_1.0.txt.gz')

# # # Read in and standardize CyC scores, join with covariates
scores <- fread('/home/ivm/cyc_prs/FinnGen_R7.UKB380_PGS_LDPRED2_inner.hm3.sscore') %>% 
  mutate(score_std = as.numeric(scale(SCORE1_AVG))) %>% 
  inner_join(covs, by = c("IID" = "FINNGENID")) %>% 
  mutate(FINNGENID = IID)
  
covid <- covid %>% 
  inner_join(scores, by = "FINNGENID")

scores_bmi <- fread('/finngen/library-red/finngen_R7/prs_1.0/data/finngen_R7_SNP_gwas_mc_merge_nogc.tbl.uniq.sscore') %>% 
  mutate(score_bmi_std = as.numeric(scale(SCORE1_AVG))) %>% 
  select(FINNGENID = IID, score_bmi_std)

covid <- covid %>% 
  inner_join(scores_bmi, by = "FINNGENID")


# # # Define cases

# ICU OR any respiratory support:
# Respiratory support	0	No information about the given respiratory support
# Respiratory support	1	Extracorporeal membrane oxygenation
# Respiratory support	2	No given respiratory support
# Respiratory support	3	other oxygen therapy
# Respiratory support	4	ventilation (including non-invasive ventilation)

icu <- covid %>% 
  filter( Intensive_care == 1 | Respiratory_support %in% c(1,3,4)) %>% 
  mutate(case = 1)

# Hospitalized
hosp <- covid %>% 
  filter(Hospital_treatment == 1) %>% 
  mutate(case = 1)

# # # Define controls

# Covid+ controls
controls <- covid %>% 
  filter(!FINNGENID %in% c(icu$FINNGENID, hosp$FINNGENID)) %>% 
  mutate(case = 0)

icu1 <- bind_rows(icu, controls)
hosp1 <- bind_rows(hosp, controls)

# Population controls
pop_controls <- scores %>% 
  filter(!FINNGENID %in% covid$FINNGENID) %>% 
  mutate(case = 0)

icu2 <- bind_rows(icu,pop_controls)
hosp2 <- bind_rows(hosp,pop_controls)


# # # Run regression models

# ICU vs COVID19+
g_icu <- glm(formula = case ~ score_std + PC1 + PC2 + PC3 + PC4 + AGE_AT_DEATH_OR_END_OF_FOLLOWUP + SEX + batch,
        family = binomial, data = icu1)

# Hospitalized vs COVID19+
g_hosp <- glm(formula = case ~ score_std + PC1 + PC2 + PC3 + PC4 + AGE_AT_DEATH_OR_END_OF_FOLLOWUP + SEX + batch,
             family = binomial, data = hosp1)

# ICU vs COVID19+, BMI adjusted
g_icu_bmi <- glm(formula = case ~ score_std + BMI + PC1 + PC2 + PC3 + PC4 + AGE_AT_DEATH_OR_END_OF_FOLLOWUP + SEX + batch,
             family = binomial, data = icu1)

# Hospitalized vs COVID19+, BMI adjusted
g_hosp_bmi <- glm(formula = case ~ score_std + BMI + PC1 + PC2 + PC3 + PC4 + AGE_AT_DEATH_OR_END_OF_FOLLOWUP + SEX + batch,
              family = binomial, data = hosp1)

# ICU vs population controls
g_icu2 <- glm(formula = case ~ score_std + PC1 + PC2 + PC3 + PC4 + AGE_AT_DEATH_OR_END_OF_FOLLOWUP + SEX + batch,
             family = binomial, data = icu2)

# Hospitalized vs population controls
g_hosp2 <- glm(formula = case ~ score_std + PC1 + PC2 + PC3 + PC4 + AGE_AT_DEATH_OR_END_OF_FOLLOWUP + SEX + batch,
              family = binomial, data = hosp2)

# ICU vs population controls, BMI adjusted
g_icu_bmi2 <- glm(formula = case ~ score_std + BMI + PC1 + PC2 + PC3 + PC4 + AGE_AT_DEATH_OR_END_OF_FOLLOWUP + SEX + batch,
              family = binomial, data = icu2)

# Hospitalized vs population controls, BMI adjusted
g_hosp_bmi2 <- glm(formula = case ~ score_std + BMI + PC1 + PC2 + PC3 + PC4 + AGE_AT_DEATH_OR_END_OF_FOLLOWUP + SEX + batch,
               family = binomial, data = hosp2)


# # #  Summarize and write results
res <- data.frame(Analysis = c("ICU vs covid19+", "ICU vs covid19+",
                               "Hosp. vs covid19+", "Hosp. vs covid19+",
                               "ICU vs pop controls", "ICU vs pop controls",
                               "Hosp. vs pop controls", "Hosp. vs pop controls"),
                  Adjust = rep(c("-", "BMI"), 4),
                  Ncases = c(table(icu1$case)["1"], table(icu1$case, !is.na(icu1$BMI))["1", "TRUE"],
                             table(hosp1$case)["1"], table(hosp1$case, !is.na(hosp1$BMI))["1", "TRUE"],
                             table(icu2$case)["1"], table(icu2$case, !is.na(icu2$BMI))["1", "TRUE"],
                             table(hosp2$case)["1"], table(hosp2$case, !is.na(hosp2$BMI))["1", "TRUE"]),
                  Ncontrols = c(table(icu1$case)["0"], table(icu1$case, !is.na(icu1$BMI))["0", "TRUE"],
                                table(hosp1$case)["0"], table(hosp1$case, !is.na(hosp1$BMI))["0", "TRUE"],
                                table(icu2$case)["0"], table(icu2$case, !is.na(icu2$BMI))["0", "TRUE"],
                                table(hosp2$case)["0"], table(hosp2$case, !is.na(hosp2$BMI))["0", "TRUE"]),
                  beta_cyc = c(coef(summary(g_icu))["score_std",1], coef(summary(g_icu_bmi))["score_std",1],
                               coef(summary(g_hosp))["score_std",1], coef(summary(g_hosp_bmi))["score_std",1],
                               coef(summary(g_icu2))["score_std",1], coef(summary(g_icu_bmi2))["score_std",1],
                               coef(summary(g_hosp2))["score_std",1], coef(summary(g_hosp_bmi2))["score_std",1]),
                  se_beta_cyc = c(coef(summary(g_icu))["score_std",2], coef(summary(g_icu_bmi))["score_std",2],
                                  coef(summary(g_hosp))["score_std",2], coef(summary(g_hosp_bmi))["score_std",2],
                                  coef(summary(g_icu2))["score_std",2], coef(summary(g_icu_bmi2))["score_std",2],
                                  coef(summary(g_hosp2))["score_std",2], coef(summary(g_hosp_bmi2))["score_std",2]),
                  pval_cyc = c(coef(summary(g_icu))["score_std",4], coef(summary(g_icu_bmi))["score_std",4],
                               coef(summary(g_hosp))["score_std",4], coef(summary(g_hosp_bmi))["score_std",4],
                               coef(summary(g_icu2))["score_std",4], coef(summary(g_icu_bmi2))["score_std",4],
                               coef(summary(g_hosp2))["score_std",4], coef(summary(g_hosp_bmi2))["score_std",4]),
                  beta_bmi = c("-", coef(summary(g_icu_bmi))["BMI",1],
                               "-", coef(summary(g_hosp_bmi))["BMI",1],
                               "-", coef(summary(g_icu_bmi2))["BMI",1],
                               "-", coef(summary(g_hosp_bmi2))["BMI",1]),
                  se_beta_bmi = c("-", coef(summary(g_icu_bmi))["BMI",2],
                                  "-", coef(summary(g_hosp_bmi))["BMI",2],
                                  "-", coef(summary(g_icu_bmi2))["BMI",2],
                                  "-", coef(summary(g_hosp_bmi2))["BMI",2]),
                  pval_bmi = c("-", coef(summary(g_icu_bmi))["BMI",4],
                               "-", coef(summary(g_hosp_bmi))["BMI",4],
                               "-", coef(summary(g_icu_bmi2))["BMI",4],
                               "-", coef(summary(g_hosp_bmi2))["BMI",4])
                  )

fwrite(res, "/home/ivm/cyc_prs/cyc_prs_results.tsv", sep = "\t", quote=F)
