# ------------------------------------------------------------
# 
# Script name: bmi_pca.R
#
# Purpose: plot genetic PC1-PC2 distribution coloured by proportion of having/not having BMI measured 
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
library(ggplot2)
library(viridis)


covs <- fread('/finngen/library-red/finngen_R8/analysis_covariates/finngen_R8_cov_1.0.txt.gz')

covs$IS_BMI <- ifelse(is.na(covs$BMI),0,1)

fwrite(covs, "/home/ivm/cyc_prs/IS_BMI_PHENO.tsv", sep = "\t", na = "NA", quote = F)
png("cyc_prs/BMI_PCs.png", width = 8, height = 8, unit = "in", res = 300)
ggplot(covs, aes(x = PC1, y = PC2, z = IS_BMI)) +
  stat_summary_hex(fun = function(z) sum(z)/length(z), alpha = 0.8, bins = c(150,150)) +
  coord_equal() +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_fill_viridis(option = "D")
dev.off()