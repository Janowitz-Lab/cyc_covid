#!/usr/bin/env Rscript

#--------------------------------------------------------------------------------------------------------------------
# Setup environment
#----

options(width=200)
options(digits=6)
source("../scripts/print_script_name.R")

# Packages
library(data.table)
library(dplyr)

# Static
data_dir <- "/exports/igmm/eddie/UK-BioBank-proj19655/phenotypes/2021_08_18/dataportal/data/"



#--------------------------------------------------------------------------------------------------------------------
# START
#----

sessionInfo()

#----
# COVID tests
#----

covid19_result_england <- fread(paste0(data_dir,"covid19_result_england.txt"))
pos_england <- subset(covid19_result_england, result==1)$eid

covid19_result_scotland <- fread(paste0(data_dir,"covid19_result_scotland.txt"))
pos_scotland <- subset(covid19_result_scotland, result==1)$eid

covid19_result_wales <- fread(paste0(data_dir,"covid19_result_wales.txt"))
pos_wales <- subset(covid19_result_wales, result==1)$eid

pos_all <- unique(c(pos_england, pos_scotland, pos_wales))

#----
# Hospital records
#----

# In-hospital diagnoses
hesin_diag <- fread(paste0(data_dir,"hesin_diag.txt"))
                         
metabolic <- c("U071")
hesin_diag <- subset(hesin_diag, diag_icd10 %in% metabolic & level==1)
hosp_all <- unique(hesin_diag$eid)

# In-hospital mortality
hesin <- fread(paste0(data_dir,"hesin.txt"))
hesin <- inner_join(hesin, subset(hesin_diag, level==1), by=c('eid','ins_index'))
hesin <- subset(hesin, disdest_uni=="11001")

death_in_hospital <- unique(hesin$eid)

#----
# Critical care records
#----

hesin_critical <- fread(paste0(data_dir,"hesin_critical.txt"))
hesin_critical <- inner_join(hesin_critical, subset(hesin_diag, level==1), by=c('eid','ins_index'))
critical_all <- unique(hesin_critical$eid)


#----
# Deaths
#----

death_cause <- fread(paste0(data_dir,"death_cause.txt"))
death_cause <- subset(death_cause, cause_icd10 %in% metabolic & level==1)
death_all <- unique(death_cause$eid)
death_all <- unique(c(subset(hesin_critical, ccdisdest==6)$eid, death_all, death_in_hospital))

death_and_hospital <- intersect(death_all, hosp_all)


#----
# Total cases
#----

total <- unique(c(pos_all, hosp_all, critical_all, death_all))
covid <- data.frame(eid=total)
covid$hospital <- ifelse(covid$eid %in% hosp_all,1,0)
covid$critical <- ifelse(covid$eid %in% critical_all,1,0)
covid$death <- ifelse(covid$eid %in% death_all,1,0)
covid$death_hospital <- ifelse(covid$eid %in% c(death_and_hospital,death_in_hospital),1,0)

#----
# Write to file
#----

fwrite(covid, "st01_01_covid_cases.tsv", sep="\t", na="NA", quote=FALSE)