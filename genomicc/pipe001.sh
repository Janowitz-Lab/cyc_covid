#!/bin/bash

#----
# Setup environment
#----

MAIN=/exports/igmm/eddie/wilson-lab/projects/prj_187_genomicc_cystatin_prs
cd ${MAIN}

mkdir -p ${MAIN}/data
cd ${MAIN}/data

# Postcode to Easting/Northing
# https://osdatahub.os.uk/downloads/open/CodePointOpen
mkdir -p ${MAIN}/data/d000_postcode_coordinates
wget 'https://api.os.uk/downloads/v1/products/CodePointOpen/downloads?area=GB&format=CSV&redirect' -O d000_postcode_coordinates/code_point_open.zip
unzip d000_postcode_coordinates/code_point_open.zip && rm d000_postcode_coordinates/code_point_open.zip
awk -v OFS="\t" 'NR > 1' d000_postcode_coordinates/Doc/Code-Point_Open_Column_Headers.csv d000_postcode_coordinates/Data/CSV/*.csv > d000_postcode_coordinates.csv && rm -r d000_postcode_coordinates


# Polygenic risk score weights
## Downloaded from Google Drive as 'd001_pgs_weights.tsv'

# Base phenotypes
cp /exports/igmm/datastore/ISARIC4C/wp5-gwas/pheno/r7/genomicc.r7.pheno.csv d002_genomicc_phenotypes.csv
ln -s /exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/phenotypes/ukbb_19655_base_phenotypes.tsv d002_ukbb_phenotypes.tsv

# Withdrawals
awk -v FS="," -v OFS="\t" '$5 == "Full Withdrawal" {print $1}' d002_genomicc_phenotypes.csv | sort -u > d003_genomicc_withdrawals.txt
ln -s /exports/igmm/eddie/wilson-lab/data/processing/ukbb_19655/phenotypes/p013_pt_withdrawals_aug_2021/data/withdrawn_20210809.txt d003_ukbb_withdrawals.txt


# Cystatin GWAS
## Downloaded from Google Drive as 'd004_pgs_gwas.tsv.gz'


#~~~
#~~~~~~
# Step 1: Classify UKBB Covid cases
#~~~~~~
#~~~

mkdir -p ${MAIN}/p01_ukbb_covid
cd ${MAIN}/p01_ukbb_covid

Rscript ../scripts/a01_01_classify_ukbb.R |& tee log_classify_ukbb.log


#~~~
#~~~~~~
# Step 2: QC SNPs
#~~~~~~
#~~~

mkdir -p ${MAIN}/p02_snpstats
cd ${MAIN}/p02_snpstats

#----
# UK Biobank
#----

mkdir -p ${MAIN}/p02_snpstats/01_ukbb
cd ${MAIN}/p02_snpstats/01_ukbb

# Make list of SNPs
cut -f5 ../../data/d001_pgs_weights.tsv > st02_00_ukbb_snps.txt

# Get list of samples
zcat /exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/phenotypes/ukbb_19655_genetic_phenotypes.tsv.gz | awk -v OFS="\t" '$3=="TRUE" && $6==1 {print $1}' > st02_00_ukbb_samples.txt

# Exclude withdrawn individuals
awk 'NR == FNR {withdrawn[$1]++; next} withdrawn[$1] == 0' \
../../data/d003_ukbb_withdrawals.txt st02_00_ukbb_samples.txt \
> st02_00_ukbb_samples.tmp && mv st02_00_ukbb_samples.tmp st02_00_ukbb_samples.txt  

# Calculate SNP stats
qsub -N stats.ukbb -l h_vmem=8G -l h_rt=01:00:00 -t 1-522 -j y -o log_ukbb_stats.\$TASK_ID.log -cwd -V <<"QSUB"
truncate -s 0 log_ukbb_stats.${SGE_TASK_ID}.log
bgen=/exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/genotypes/hrc/bgen/chunks/ukbb_19655_chunk_${SGE_TASK_ID}.bgen
sample=/exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/genotypes/hrc/bgen/ukbb_19655.sample
qctool -g ${bgen} -s ${sample} -incl-rsids st02_00_ukbb_snps.txt -incl-samples st02_00_ukbb_samples.txt -snp-stats -osnp st02_01_ukbb_stats.${SGE_TASK_ID}
QSUB

# Merge when done
awk '$0 !~ /#/ {print; exit}' st02_01_ukbb_stats.1 > st02_01_ukbb_stats.tsv
for i in {1..522}
do

if [[ `wc -l < st02_01_ukbb_stats.${i}` -lt 15 ]]; then
    rm -f log_ukbb_stats.${i}.log st02_01_ukbb_stats.${i} # Remove empty
else
    awk '$0 !~ /#/' st02_01_ukbb_stats.${i}  | awk 'NR > 1' >> st02_01_ukbb_stats.tsv # Merge non-empty
fi
done


# Remove intermediate files after checking log files
rm -f st02_01_ukbb_stats.{1..522} log_ukbb_stats.{1..522}.log


#----
# GenOMICC
#----

mkdir -p ${MAIN}/p02_snpstats/02_genomicc
cd ${MAIN}/p02_snpstats/02_genomicc

# Get list of SNPs 
## Since we don't have rsids, we'll allow for any combination of 
## chr:pos:a1:a0, chr:pos:a0:a1, or complementary strand alleles)
awk -v OFS="\t" 'BEGIN{comp["A"]="T"; comp["T"]="A"; comp["C"]="G"; comp["G"]="C"}
NR > 1 {
    print "chr"$1":"$21":"$4":"$3; 
    print "chr"$1":"$21":"$3":"$4;
    print "chr"$1":"$21":"comp[$4]":"comp[$3];
    print "chr"$1":"$21":"comp[$3]":"comp[$4]
}' \
../../data/d001_pgs_weights.tsv > st02_00_genomicc_snps.txt

# Get list of samples
awk -v OFS="\t" 'NR == FNR {if($3=="EUR"){iid[$1]++}; next} $1 in iid {print $1}' \
/exports/igmm/datastore/ISARIC4C/wp5-gwas/QC/wp5-gwas-r7/wp5-r7-covid19_estimatedancestry.txt \
/exports/igmm/datastore/ISARIC4C/wp5-gwas/QC/wp5-gwas-r7/wp5-r7-unrelatedindividuals.txt  > st02_00_genomicc_samples.txt


# Exclude withdrawn individuals
awk 'NR == FNR {withdrawn[$1]++; next} withdrawn[$1] == 0' \
../../data/d003_genomicc_withdrawals.txt st02_00_genomicc_samples.txt \
> st02_00_genomicc_samples.tmp && mv st02_00_genomicc_samples.tmp st02_00_genomicc_samples.txt  

# Calculate SNPs stats
qsub -N stage.genomicc -q staging -l h_vmem=1G -j y -o log_genomicc_stats_stage.log -cwd -V <<"QSUB"
truncate -s 0 log_genomicc_stats_stage.log
for SGE_TASK_ID in {1..22}
do
bgen=/exports/igmm/datastore/ISARIC4C/wp5-gwas/imputation/wp5-r7/filtered/bgen/GenOMICC_chr${SGE_TASK_ID}_TOPMedFreeze5_nomono_I4.bgen
sample=/exports/igmm/datastore/ISARIC4C/wp5-gwas/imputation/wp5-r7/filtered/bgen/GenOMICC_chr${SGE_TASK_ID}_TOPMedFreeze5_nomono_I4.sample
cp ${bgen} GenOMICC_chr${SGE_TASK_ID}_TOPMedFreeze5_nomono_I4.bgen
cp ${sample} GenOMICC_chr${SGE_TASK_ID}_TOPMedFreeze5_nomono_I4.sample
done
[ -s log_genomicc_stats_stage.log ] || rm log_genomicc_stats_stage.log
QSUB

qsub -N stats.genomicc -hold_jid stage.genomicc -l h_rt=04:00:00 -t 1-22 -j y -o log_genomicc_stats.\$TASK_ID.log -cwd -V <<"QSUB"
truncate -s 0 log_genomicc_stats.${SGE_TASK_ID}.log
bgen=GenOMICC_chr${SGE_TASK_ID}_TOPMedFreeze5_nomono_I4.bgen
sample=GenOMICC_chr${SGE_TASK_ID}_TOPMedFreeze5_nomono_I4.sample
qctool -g ${bgen} -s ${sample} -incl-rsids st02_00_genomicc_snps.txt -incl-samples st02_00_genomicc_samples.txt -snp-stats -osnp st02_01_genomicc_stats.${SGE_TASK_ID} && rm ${bgen} ${sample}
QSUB


# Merge when done
awk '$0 !~ /#/ {print; exit}' st02_01_genomicc_stats.1 > st02_01_genomicc_stats.tsv
for i in {1..22}
do
if [[ `wc -l < st02_01_genomicc_stats.${i}` -lt 15 ]]; then
    rm -f log_genomicc_stats.${i}.log st02_01_genomicc_stats.${i} # Remove empty
else
    awk '$0 !~ /#/' st02_01_genomicc_stats.${i}  | awk 'NR > 1' >> st02_01_genomicc_stats.tsv # Merge non-empty
fi
done


# Remove intermediate files after checking log files
rm -f st02_01_genomicc_stats.{1..22} log_genomicc_stats.{1..22}.log



#----
# Merge & Compare
#----

mkdir -p ${MAIN}/p02_snpstats/03_merge
cd ${MAIN}/p02_snpstats/03_merge

# Create 1-to-1 map of SNPs
awk -v OFS="\t" 'BEGIN{comp["A"]="T"; comp["T"]="A"; comp["C"]="G"; comp["G"]="C"; print "ukbb","genomicc"} 
FNR == 1 {next}
ARGIND==1 {snpid[$3"_"$4]=$2; next}
ARGIND==2 && "chr"$1"_"$21 in snpid {print $5,snpid["chr"$1"_"$21]}' \
../02_genomicc/st02_01_genomicc_stats.tsv ../../data/d001_pgs_weights.tsv > st02_00_study_snps.tsv


# Extract data of interest 
cat st02_00_study_snps.tsv | awk -v OFS="\t" \
'BEGIN{print "snpid_ukbb","snpid_genomicc","a1","a0","hwe_exact_p","hwe_lrt_p","maf","info","impute_info","missingness","study"}
FNR == 1 {next} NR == FNR {rsid[$2]=$1; snpid[$1]=$2; next}
$2 in snpid {print $2,($2 in snpid ? snpid[$2] : "NA"),$5,$6,$8,$9,$14,$17,$18,$19,"ukbb"; next}
$2 in rsid {print ($2 in rsid ? rsid[$2] : "NA"),$2,$5,$6,$8,$9,$14,$17,$18,$19,"genomicc"}' \
- ../01_ukbb/st02_01_ukbb_stats.tsv ../02_genomicc/st02_01_genomicc_stats.tsv > st02_01_snpstats.tsv 


# Compare studies
Rscript ../../scripts/a02_03_compare_snpstats.R



# Print number of strand flips (these are fixed automatically in PRSice)
awk -v FS=":|\t" 'NR ==1 {next} $4 != $6 && $5 != $6' st02_01_snpstats.tsv  | wc -l



#~~~
#~~~~~~
# Step 3: Scores UKBB
#~~~~~~
#~~~


mkdir -p ${MAIN}/p03_scores_ukbb
cd ${MAIN}/p03_scores_ukbb


# Format weights
mkdir -p ${MAIN}/p03_scores_ukbb/p01_weights
cd ${MAIN}/p03_scores_ukbb/p01_weights


ln -sf ../../p02_snpstats/01_ukbb/st02_00_ukbb_samples.txt st03_01_ukbb_samples.txt
awk -v OFS="\t" 'NR > 1 {print $1}' ../../p02_snpstats/03_merge/st02_03_shared_snps.tsv > st03_01_ukbb_snps.txt

awk -v OFS="\t" 'NR==1 {print "rsid","chr","pos","a1","a0","beta1","p"; next}
NR == FNR {rsid[$1]++; quality[$1]=$NF; next}
$5 in rsid {print $5,$1,$2,$4,$3,$(NF-1),(quality[$5] == "strict" ? 0 : 0.5)}' \
../../p02_snpstats/03_merge/st02_03_shared_snps.tsv ../../data/d001_pgs_weights.tsv > st03_01_weights.ukbb.tsv



# Calculate scores
mkdir -p ${MAIN}/p03_scores_ukbb/p02_scores
cd ${MAIN}/p03_scores_ukbb/p02_scores

qsub -N prs.ukbb -l h_vmem=8G -l h_rt=02:30:00 -t 1-522 -j y -o log_calculate_score.\$TASK_ID.log -cwd -V <<"QSUB"
truncate -s 0 log_calculate_score.${SGE_TASK_ID}.log
/gpfs/igmmfs01/eddie/wilson-lab/apps_by_us_full_stuff/prs_pipeline/scripts/prsice/PRSice_linux \
    --base ../p01_weights/st03_01_weights.ukbb.tsv \
    --extract ../p01_weights/st03_01_ukbb_snps.txt \
    --keep ../p01_weights/st03_01_ukbb_samples.txt \
    --ignore-fid \
    --snp rsid \
    --chr chr \
    --bp pos \
    --a1 a1 \
    --a2 a0 \
    --stat beta1 \
    --pvalue p \
    --bar-levels 1,0 \
    --beta  \
    --binary-target F \
    --fastscore  \
    --lower 1 \
    --no-clump  \
    --no-regress  \
    --num-auto 22 \
    --out st03_02_ukbb_score.${SGE_TASK_ID} \
    --score sum \
    --seed 1680002709 \
    --target /exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/genotypes/hrc/bgen/chunks/ukbb_19655_chunk_${SGE_TASK_ID},/exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/genotypes/hrc/bgen/ukbb_19655.sample \
    --print-snp \
    --thread 1 \
    --type bgen \
    --upper 1
QSUB



# Combine when done
awk -v OFS="\t" 'BEGIN{print "iid","score_strict","score"} FNR > 1 {score_strict[$1]+=$3; score[$1]+=$4} END{for(iid in score){print iid,score_strict[iid],score[iid]}}' st03_02_ukbb_score.*.all_score > st03_02_ukbb_score.tsv
awk -v OFS="\t" 'NR == 1 {print $0} FNR > 1 {pheno[$3]=$1; set[$3]=$2; nsnp[$3]+=$4} END{for (thresh in nsnp) {print pheno[thresh],set[thresh],thresh,nsnp[thresh]}}' st03_02_ukbb_score.*.prsice > st03_02_ukbb_score.prsice
awk -v OFS="\t" 'NR ==1 || FNR > 1' st03_02_ukbb_score.*.snp | sort -k1g,1 -k3g,3 > st03_02_ukbb_score.snp

for log in log_calculate_score.{1..522}.log; do grep -q "No vairant remained" $log && rm $log; done


rm st03_02_ukbb_score.{1..522}.log st03_02_ukbb_score.*.all_score st03_02_ukbb_score.*.snp st03_02_ukbb_score.*.prsice

#~~~
#~~~~~~
# Step 4: Scores GenOMICC 
#~~~~~~
#~~~


mkdir -p ${MAIN}/p04_scores_genomicc
cd ${MAIN}/p04_scores_genomicc


# Format weights
mkdir -p ${MAIN}/p04_scores_genomicc/p01_weights
cd ${MAIN}/p04_scores_genomicc/p01_weights


ln -sf ../../p02_snpstats/02_genomicc/st02_00_genomicc_samples.txt st04_01_genomicc_samples.txt
awk -v OFS="\t" 'NR > 1 {print $2}' ../../p02_snpstats/03_merge/st02_03_shared_snps.tsv > st04_01_genomicc_snps.txt


awk -v OFS="\t" 'NR==1  {print "rsid","chr","pos","a1","a0","beta1","p"; next}
NR == FNR {snpid[$1]=$2; quality[$1]=$NF; next}
$5 in snpid {print snpid[$5],"chr"$1,$21,$4,$3,$(NF-1),(quality[$5] == "strict" ? 0 : 0.5)}' \
../../p02_snpstats/03_merge/st02_03_shared_snps.tsv ../../data/d001_pgs_weights.tsv > st04_01_weights.genomicc.tsv

 

# Calculate scores

mkdir -p ${MAIN}/p04_scores_genomicc/p02_scores
cd ${MAIN}/p04_scores_genomicc/p02_scores

qsub -N stage.genomicc -q staging -l h_vmem=1G -j y -o log_genomicc_stats_stage.log -cwd -V <<"QSUB"
truncate -s 0 log_genomicc_stats_stage.log
for SGE_TASK_ID in {1..22}
do
bgen=/exports/igmm/datastore/ISARIC4C/wp5-gwas/imputation/wp5-r7/filtered/bgen/GenOMICC_chr${SGE_TASK_ID}_TOPMedFreeze5_nomono_I4.bgen
sample=/exports/igmm/datastore/ISARIC4C/wp5-gwas/imputation/wp5-r7/filtered/bgen/GenOMICC_chr${SGE_TASK_ID}_TOPMedFreeze5_nomono_I4.sample
cp ${bgen} GenOMICC_chr${SGE_TASK_ID}_TOPMedFreeze5_nomono_I4.bgen
cp ${sample} GenOMICC_chr${SGE_TASK_ID}_TOPMedFreeze5_nomono_I4.sample
done
[ -s log_genomicc_stats_stage.log ] || rm log_genomicc_stats_stage.log
QSUB


qsub -N prs.genomicc -hold_jid stage.genomicc -l h_vmem=8G -l h_rt=00:30:00 -j y -o log_calculate_score.log -cwd -V <<"QSUB"
truncate -s 0 log_calculate_score.log

/gpfs/igmmfs01/eddie/wilson-lab/apps_by_us_full_stuff/prs_pipeline/scripts/prsice/PRSice_linux \
    --base ../p01_weights/st04_01_weights.genomicc.tsv \
    --extract ../p01_weights/st04_01_genomicc_snps.txt \
    --keep ../p01_weights/st04_01_genomicc_samples.txt \
    --ignore-fid \
    --snp rsid \
    --chr chr \
    --bp pos \
    --a1 a1 \
    --a2 a0 \
    --stat beta1 \
    --pvalue p \
    --bar-levels 1,0 \
    --beta  \
    --binary-target F \
    --fastscore  \
    --lower 1 \
    --no-clump  \
    --no-regress  \
    --num-auto 22 \
    --out st04_02_genomicc_score \
    --score sum \
    --seed 1680002709 \
    --target GenOMICC_chr#_TOPMedFreeze5_nomono_I4,GenOMICC_chr1_TOPMedFreeze5_nomono_I4.sample \
    --print-snp \
    --thread 1 \
    --type bgen \
    --upper 1 \
&& rm GenOMICC_chr{1..22}_TOPMedFreeze5_nomono_I4.bgen GenOMICC_chr{1..22}_TOPMedFreeze5_nomono_I4.sample 

awk -v OFS="\t" 'BEGIN{print "iid","score_strict","score"} FNR > 1 {print $1, $3, $4}' st04_02_genomicc_score.all_score > st04_02_genomicc_score.tsv
QSUB

# Combine when done




#~~~
#~~~~~~
# Step 5: Calculate principal components
#~~~~~~
#~~~


mkdir -p ${MAIN}/p05_pca
cd ${MAIN}/p05_pca


#----
# Format GenOMICC array files
#----


# Get GenOMICC array genotype data

qsub -N stage.genomicc -q staging -l h_vmem=1G -j y -o log_genomicc_stage.log -cwd -V <<"QSUB"
truncate -s 0 log_genomicc_stage.log
bfile=/exports/igmm/datastore/ISARIC4C/wp5-gwas/QC/wp5-gwas-r7/wp5-gwas-r7-20210630-miss95-geno99-maf01-miss97-nodup
for ext in bim bed fam
do
cp ${bfile}.${ext} genomicc.${ext}
done
[ -s log_genomicc_stage.log ] || rm log_genomicc_stage.log
QSUB


# Subset to unrelated EUR ancestry
plink \
--bfile genomicc \
--keep <(awk '{print $1,$1}' ../p04_scores_genomicc/p01_weights/st04_01_genomicc_samples.txt) \
--make-bed \
--out st05_01_genomicc && rm genomicc.bim genomicc.bed genomicc.fam


wc -l st05_01_genomicc.bim st05_01_genomicc.fam # Print number of SNPs and samples


# Fix SNPs missing chromosome or position
rsids=`mktemp`
awk -v OFS="\t" '$1 == 0 || $4 == 0 {print $2}' st05_01_genomicc.bim > ${rsids}
snp2gene.sh ${rsids} > ${rsids}.dat 
awk -v OFS="\t" 'NR == FNR {chr[$1]=$2; pos[$1]=$3; next} $2 in chr {$1=chr[$2]; $4=pos[$2]} 1' \
${rsids}.dat st05_01_genomicc.bim > st05_01_genomicc.bim.tmp && \
mv st05_01_genomicc.bim.tmp st05_01_genomicc.bim && rm ${rsids} ${rsids}.dat


# Rename rsids to match UK Biobank
ukbb_bim=`mktemp`; echo /exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/genotypes/affy/plink/ukbb_19655_chr_{1..22}.bim  | tr -s ' ' '\n' > ${ukbb_bim}
wc -l `cat ${ukbb_bim}`

awk -v OFS="\t" '(ARGIND+1) < ARGC {rsid[$1"_"$4]=$2; next} $1"_"$4 in rsid {$2=rsid[$1"_"$4]} 1' \
`cat ${ukbb_bim}` st05_01_genomicc.bim > st05_01_genomicc.bim.tmp && \
mv st05_01_genomicc.bim.tmp st05_01_genomicc.bim && rm ${ukbb_bim}


# Find multi-allelic SNPs
awk -v OFS="\t" '{id=$1"_"$4; snpid[id]++} snpid[id] > 1 {print $2}' st05_01_genomicc.bim | sort -u > st05_01_genomicc_dups.txt

# Get SNPs shared with UK Biobank + same allele coding
ukbb_bim=`mktemp`; echo /exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/genotypes/affy/plink/ukbb_19655_chr_{1..22}.bim  | tr -s ' ' '\n' > ${ukbb_bim}
awk -v OFS="\t" '(ARGIND+1) < ARGC {a1[$2]=$5; a0[$2]=$6; next} (a1[$2]==$5 && a0[$2]==$6) || (a0[$2]==$5 && a1[$2]==$6) {print $2}' \
`cat ${ukbb_bim}` st05_01_genomicc.bim | sort -u > st05_02_shared_snps.txt

# Remove multi-allelic SNPs from shared list
grep -v -wFf st05_01_genomicc_dups.txt st05_02_shared_snps.txt > st05_02_shared_snps.tmp && mv st05_02_shared_snps.tmp st05_02_shared_snps.txt && rm st05_01_genomicc_dups.txt
wc -l st05_02_shared_snps.txt

# Subset GenOMICC SNPs to bi-allelic SNPs shared with UKBB
for i in {1..22}
do
plink \
--bfile st05_01_genomicc \
--extract st05_02_shared_snps.txt \
--chr ${i} \
--make-bed \
--out st05_01_genomicc_chr_${i}
done

rm -f st05_01_genomicc.* *.nosex

#----
# Merge with UK Biobank
#----

# Get list of samples
cat ../p02_snpstats/02_genomicc/st02_00_genomicc_samples.txt > st05_02_samples.txt
cat ../p02_snpstats/01_ukbb/st02_00_ukbb_samples.txt >> st05_02_samples.txt
paste st05_02_samples.txt st05_02_samples.txt > st05_02_samples.tmp && mv st05_02_samples.tmp st05_02_samples.txt


# Merge with UK Biobank data
qsub -N merge.ukbb -l h_vmem=8G -l h_rt=00:10:00 -t 1-22 -j y -o log_merge_ukbb_chr_\$TASK_ID.log -cwd -V <<"QSUB"
i=${SGE_TASK_ID}
truncate -s 0 log_merge_ukbb_chr_${i}.log

echo -e "\n============\nChromosome ${i}\n\n========================\n\n"

echo -e ">> Subsetting UKBB\n"

/exports/igmm/software/pkg/el7/apps/plink/1.90b4/plink \
--bfile /exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/genotypes/affy/plink/ukbb_19655_chr_${i} \
--extract st05_02_shared_snps.txt \
--keep st05_02_samples.txt \
--make-bed \
--out st05_03_ukbb_chr_${i}

echo -e ">> Merging GenOMICC + UKBB\n"

/exports/igmm/software/pkg/el7/apps/plink/1.90b4/plink \
--bfile st05_01_genomicc_chr_${i} \
--bmerge st05_03_ukbb_chr_${i} \
--make-bed \
--out st05_04_merged_chr_${i} && rm st05_01_genomicc_chr_${i}.* st05_03_ukbb_chr_${i}.*
QSUB



#----
# Calculate PCs
#----

# https://sites.google.com/a/broadinstitute.org/ricopili/pca -- QC recommendations
# https://github.com/gabraham/flashpca -- Regions to exclude + LD pruning parameters

# Get list of ambiguous AT CG SNPs, MHC, and chromosome 8 inversion
awk -v OFS="\t" 'BEGIN{AT["A"]++; AT["T"]++; CG["C"]++; CG["G"]++}
$1 == 5 && $4 >= 44000000 && $4 <= 51500000 {print $2,"r1"; next}
$1 == 6 && $4 >= 25000000 && $4 <= 33500000 {print $2,"r2"; next}
$1 == 8 && $4 >= 8000000 && $4 <= 12000000 {print $2,"r3"; next}
$1 == 11 && $4 >= 45000000 && $4 <= 57000000 {print $2,"r4"; next}
($5 in AT) && ($6 in AT) {print $2, "AT"; next}
($5 in CG) && ($6 in CG) {print $2, "CG"; next}' \
st05_04_merged_chr_*.bim > st05_05_pca.exclude


# Exclude the above SNPs, perform initial QC, and start first LD pruning pass

qsub -N pca.sample -l h_vmem=8G -l h_rt=00:30:00 -t 1-22 -j y -o log_pca_sample_chr_\$TASK_ID.log -cwd -V <<"QSUB"
i=${SGE_TASK_ID}
truncate -s 0 log_pca_sample_chr_${i}.log

/exports/igmm/software/pkg/el7/apps/plink/1.90b4/plink \
--bfile st05_04_merged_chr_${i} \
--maf 0.01 \
--hwe 1e-3 \
--geno 0.02 \
--exclude st05_05_pca.exclude \
--indep-pairwise 1000 50 0.05 \
--threads `nproc` \
--out st05_05_pca_sample_chr_${i}

# Second round of LD pruning
/exports/igmm/software/pkg/el7/apps/plink/1.90b4/plink \
--bfile st05_04_merged_chr_${i} \
--extract st05_05_pca_sample_chr_${i}.prune.in \
--indep-pairwise 1000 50 0.05 \
--threads `nproc` \
--out st05_05_pca_sample

# Create PLINK files with the robust, pruned SNP list
/exports/igmm/software/pkg/el7/apps/plink/1.90b4/plink \
--bfile st05_04_merged_chr_${i} \
--extract st05_05_pca_sample_chr_${i}.prune.in \
--threads `nproc` \
--make-bed \
--out st05_05_pca_sample_chr_${i} && rm st05_04_merged_chr_${i}.*
QSUB


# Merge chromosomes
pca_sample_list=`mktemp`
echo st05_05_pca_sample_chr_{1..22} | tr -s ' ' '\n' > ${pca_sample_list}

plink \
--merge-list ${pca_sample_list} \
--make-bed \
--out st05_05_pca_sample && rm st05_05_pca_sample_chr_*


# Calculate PCs (recommend using 8+ cores for this)
# qlogin -l h_vmem=2G -pe interactivemem 16
plink2 \
--bfile st05_05_pca_sample \
--threads `nproc` \
--pca approx \
--out st05_06_pca_sample


#----
# Identify duplicates across cohorts
#----


# Calculate duplicates (can be done with 1 core & 16G memory; ~10 minutes)

king \
-b st05_05_pca_sample.bed \
--prefix st05_07_duplicates \
--duplicate \
--cpus `nproc` &> st05_07_duplicates.log




#~~~
#~~~~~~
# Step 6: Analyse scores
#~~~~~~
#~~~


mkdir -p ${MAIN}/p06_analyse
cd ${MAIN}/p06_analyse


Rscript ../scripts/a06_01_analyse.R



#~~~
#~~~~~~
# Step 7: Perform MR
#~~~~~~
#~~~


mkdir -p ${MAIN}/p07_cystatin_mr
cd ${MAIN}/p07_cystatin_mr

# Get BMI and Townsend Deprivation Index GWAS
wget 'https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz' -O 21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget 'https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/189_irnt.gwas.imputed_v3.both_sexes.tsv.bgz' -O 189_irnt.gwas.imputed_v3.both_sexes.tsv.bgz

#----
# Format sumstats
#----

# Cystatin C
awk -v OFS="\t" 'NR == 1 {print "SNP","chr","pos","A1","A2","N","freq1","b","se"; next}
$6 > 0.001 {print $3,$1,$2,$4,$5,$10,$11,$7,$8}'  ../data/d004_prs_gwas.tsv > st07_01_cystatin.tsv

# BMI
zcat 21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz | \
awk -v OFS="\t" 'NR == 1 {print "SNP","chr","pos","A1","A2","N","freq1","b","se"; next}
NR == FNR {rsid[$2"_"$3]=$1; next}
$3 > 0.001 {split($1,chrpos,":"); if(chrpos[1]"_"chrpos[2] in rsid) {print rsid[chrpos[1]"_"chrpos[2]],chrpos[1],chrpos[2],chrpos[4],chrpos[3],$5,($2 == chrpos[4] ? $3 : 1-$3),$8,$9}}' \
st07_01_cystatin.tsv - > st07_01_bmi.tsv && rm 21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz

# Townsend
zcat 189_irnt.gwas.imputed_v3.both_sexes.tsv.bgz | \
awk -v OFS="\t" 'NR == 1 {print "SNP","chr","pos","A1","A2","N","freq1","b","se"; next}
NR == FNR {rsid[$2"_"$3]=$1; next}
$3 > 0.001 {split($1,chrpos,":"); if(chrpos[1]"_"chrpos[2] in rsid) {print rsid[chrpos[1]"_"chrpos[2]],chrpos[1],chrpos[2],chrpos[4],chrpos[3],$5,($2 == chrpos[4] ? $3 : 1-$3),$8,$9}}' \
st07_01_cystatin.tsv - > st07_01_townsend.tsv && rm 189_irnt.gwas.imputed_v3.both_sexes.tsv.bgz



# Calculate genetic correlations with Cystatin PRS
for trait1 in cystatin bmi
do
    for trait2 in bmi townsend
    do
        if [[ "${trait1}" == "${trait2}" ]]; then
            continue
        fi

qsub -N rg.${trait1:0:3}.${trait2:0:3} -l h_vmem=8G -l h_rt=00:30:00 -j y -o st07_02_rg_${trait1}_${trait2}.log -cwd -V <<QSUB
truncate -s 0 st07_02_rg_${trait1}_${trait2}.log
Rscript ../scripts/a07_02_calculate_correlations.R ${trait1} ${trait2}
QSUB

    done
done


# Combine results
awk -v OFS="\t" 'NR == 1 || FNR > 1' st07_02_rg_*_*.tsv > st07_02_rg_cystatin.tsv

