#!/usr/bin/env Rscript

#----
# Setup environment
#----

options(width=200)
options(digits=6)
library(data.table)

# Variables
hwe_p_thresh <- 1e-6
maf_thresh <- 0.01
info_thresh <- 0.8
missingness_thresh <- 0.01


# Standard
max_maf_diff <- 0.1
max_info_diff <- 1


# Strict
max_maf_diff_strict <- 0.01
max_info_diff_strict <- 0.02


#----
# START
#----

dat <- fread("st02_01_snpstats.tsv")
n_ukbb <- as.numeric(system("tail -n+2 ../01_ukbb/st02_01_ukbb_stats.tsv | wc -l ", intern=T))
n_genomicc <- as.numeric(system("tail -n+2 ../02_genomicc/st02_01_genomicc_stats.tsv | wc -l ", intern=T))
n <- data.table(study=c("ukbb","genomicc"),N=c(n_ukbb,n_genomicc))

#----
# Comparison
#----

# Apply filters
dat1 <- dat[!(hwe_exact_p < hwe_p_thresh | hwe_lrt_p < hwe_p_thresh | maf < maf_thresh | info < info_thresh | impute_info < info_thresh | missingness > missingness_thresh)]

n_filter <- dat1[,.N,by="study"]
n_snps <- merge(n, n_filter,by="study", suffixes=c("_orig","_filter"))
print(n_snps[,.(study, N_orig, N_filter, `%`=100*N_filter/N_orig)])


# Format for comparison
mdat1 <- dcast(dat1, snpid_ukbb+snpid_genomicc ~ study, value.var=c("maf","info"))

# Remove SNPs that are not shared between studies
mdat1 <- mdat1[complete.cases(mdat1)] 
print(mdat1[,.N])

# Plot MAF and INFO differences
pdf("st02_02_ukbb_genomicc_differences_preqc.pdf", width=8, height=6)
par(mfrow=c(1,2))
with(mdat1, hist(maf_ukbb-maf_genomicc))
with(mdat1, hist(info_ukbb-info_genomicc))
dev.off()

# Subset to highly similar SNPs
mdat1[,c("maf_diff", "info_diff"):=.(maf_ukbb-maf_genomicc, info_ukbb-info_genomicc)]
mdat2 <- mdat1[abs(maf_diff) < max_maf_diff & abs(info_diff) < max_info_diff] # Standard
mdat3 <- mdat1[abs(maf_diff) < max_maf_diff_strict & abs(info_diff) < max_info_diff_strict] # Strict


# Plot MAF and INFO differences post-qc

pdf("st02_02_ukbb_genomicc_differences_postqc.pdf", width=8, height=6)
par(mfrow=c(1,2))
with(mdat2, hist(maf_ukbb-maf_genomicc))
with(mdat2, hist(info_ukbb-info_genomicc))
dev.off()

pdf("st02_02_ukbb_genomicc_differences_postqc_strict.pdf", width=8, height=6)
par(mfrow=c(1,2))
with(mdat3, hist(maf_ukbb-maf_genomicc))
with(mdat3, hist(info_ukbb-info_genomicc))
dev.off()



t1 <- mdat2[,.(delta_maf=max_maf_diff, delta_info=max_info_diff, n_snps=.N, maf_ukbb=mean(maf_ukbb), maf_genomicc=mean(maf_genomicc), info_ukbb=mean(info_ukbb), info_genomicc=mean(info_genomicc))]
t2 <- mdat3[,.(delta_maf=max_maf_diff_strict, delta_info=max_info_diff_strict, n_snps=.N, maf_ukbb=mean(maf_ukbb), maf_genomicc=mean(maf_genomicc), info_ukbb=mean(info_ukbb), info_genomicc=mean(info_genomicc))]

print(rbind(t1,t2))

#----
# Write to file
#----

snp_dat <- copy(mdat2)
snp_dat[,quality:=ifelse(snpid_ukbb %in% mdat3$snpid_ukbb, "strict", "standard")]

fwrite(data.table(format(snp_dat,digits=4)), "st02_03_shared_snps.tsv", sep="\t", quote=FALSE, na="NA")
