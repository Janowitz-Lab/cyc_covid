#!/usr/bin/env Rscript

#--------------------------------------------------------------------------------------------------------
# Setup environment
#----

options(width=200)
options(digits=6)

# Load stats packages
library(data.table)
library(speedglm)
library(rgdal)
library(sp)

# Load plotting packages
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggthemes)


# UKBB Fields
traits <- c("f.eid","f.21003.0.0","f.31.0.0", "f.20075.0.0", "f.20074.0.0", "f.21001.0.0")
trait_names <- c("iid", "age", "sex", "northing", "easting", "bmi")


#--------------------------------------------------------------------------------------------------------
# START
#----

sessionInfo()


#----
# Load data
#----

# Samples
samples <- fread("../p05_pca/st05_02_samples.txt", header=FALSE, col.names=c("fid","iid"))[,iid]
duplicates <- fread("../p05_pca/st05_07_duplicates.con")
duplicates[,ID1:=as.character(ID1)]


# Postcodes 
postcodes  <- fread("../data/d000_postcode_coordinates.csv")
postcodes[,postcode_partial:=gsub(" ", "",substr(Postcode, start=0, stop=4))]
postcodes_list <- postcodes[,.(northing=median(Northings), easting=median(Eastings)), by="postcode_partial"]


# UK Biobank Covid Classifier
classifier <- fread("../p01_ukbb_covid/st01_01_covid_cases.tsv")
classifier <- classifier[eid %in% samples,]

# Scores
ukbb <- fread("../p03_scores_ukbb/p02_scores/st03_02_ukbb_score.tsv")
ukbb[,study:="ukbb"]
genomicc <- fread("../p04_scores_genomicc/p02_scores/st04_02_genomicc_score.tsv")
genomicc[,study:="genomicc"]

scores <- rbind(ukbb, genomicc)
scores_desc <- fread("../p03_scores_ukbb/p02_scores/st03_02_ukbb_score.prsice")[,.(score=c("score_strict","score"),n_snps=Num_SNP)]


# Principal components
genetic_pcs <- fread("../p05_pca/st05_06_pca_sample.eigenvec", col.names=c("fid","iid",paste0("pc",1:10)))


# Phenotypes
ph_ukbb <- fread("../data/d002_ukbb_phenotypes.tsv", select=traits, col.names=trait_names)
ph_ukbb[,study:="ukbb"]
ph_ukbb[,covid_status:=ifelse(iid %in% classifier$eid, 1, 0)]


ph_genomicc <- fread("../data/d002_genomicc_phenotypes.csv", col.names=c("iid", "age", "sex", "postcode_partial","withdrawn","covid_status"), na.strings=c("","NA"))
ph_genomicc[,study:="genomicc"]


#----
# Standardise score
#----
 
# Remove duplicates and calculate population characteristics

pop_mean <- scores[!iid %in% duplicates[,c(ID1,ID2)],mean(score)]
pop_sd <- scores[!iid %in% duplicates[,c(ID1,ID2)],sd(score)]
scores[,score:=.((score-pop_mean)/pop_sd)]



pop_mean <- scores[!iid %in% duplicates[,c(ID1,ID2)],mean(score_strict)]
pop_sd <- scores[!iid %in% duplicates[,c(ID1,ID2)],sd(score_strict)]
scores[,score_strict:=.((score_strict-pop_mean)/pop_sd)]




#----
# Compare duplicates
#----

d1 <- scores[duplicates,,on=c(iid="ID1")][,.(iid=ID2,score,score_strict,study)]
d2 <- scores[duplicates,,on=c(iid="ID2")][,.(iid, score, score_strict,study)]
duplicate_scores <- dcast(rbind(d1,d2), iid~study, value.var=c("score", "score_strict"))
duplicate_scores[,c("score_diff","score_strict_diff"):=.(score_genomicc - score_ukbb, score_strict_genomicc - score_strict_ukbb)]

pdf("st06_01_duplicates_score_diff.pdf", height=6, width=8)
par(mfrow=c(1,2))
with(duplicate_scores, plot(score_ukbb, score_genomicc))
abline(0,1)
with(duplicate_scores, plot(score_strict_ukbb, score_strict_genomicc))
abline(0,1)
dev.off()


t_diff <- with(duplicate_scores, t.test(x=score_ukbb, y=score_genomicc, paired=TRUE))
t_diff_strict <- with(duplicate_scores, t.test(x=score_strict_ukbb, y=score_strict_genomicc, paired=TRUE))

print(t_diff)
print(t_diff_strict)


# Remove UKBB duplicates from rest of analysis

samples <- samples[!samples %in% duplicates$ID1]
scores <- scores[iid %in% samples,]
ph_ukbb <- ph_ukbb[iid %in% samples,]
ph_genomicc <- ph_genomicc[iid %in% samples,]



#----
# Spatial coordinates and deprivation index
#----

ph_genomicc1 <- postcodes_list[ph_genomicc,.(iid, age, sex, northing, easting, covid_status, study),on="postcode_partial"]
ph <- rbind(ph_genomicc1, ph_ukbb, fill=TRUE)

geo_dir <- "/exports/igmm/eddie/wilson-lab/data/processing/ukbb_19655/phenotypes/p012_pt_geography_2021/data/"
sp_lsoa <- readOGR(geo_dir, paste0("infuse_lsoa_lyr_2011_clipped"))
sp_dist <- readOGR(geo_dir, paste0("infuse_dist_lyr_2011_clipped"))
townsend_mapper <- fread(paste0(geo_dir,"townsend_lsoa_2011.csv"), col.names=c("id","geo_code","geo_label","town_di","town_di_quintile"))

# Select complete coordinate data
ph1 <- ph[,.(iid, northing, easting)] 
coord_dat <- ph1[!is.na(easting) & !is.na(northing), ]
coord_dat_iid <- coord_dat$iid

# Convert raw easting and northing data to spacial coordinates
coordinates(coord_dat) <- c("easting","northing")
proj4string(coord_dat) <- CRS("+init=epsg:27700")

# Reproject easting and northing data to match lsoa/msoa coordinate system
coord_sp <- spTransform(coord_dat, CRS(proj4string(sp_lsoa)))

# Map easting northing points onto lsoa polygons
area_dat <- coord_sp %over% sp_lsoa
area_dat <- townsend_mapper[area_dat,,on=c("geo_code","geo_label")]

# Add results onto original data frame
ph[match(coord_dat_iid,iid),"lsoa"] <- area_dat$geo_code
ph[match(coord_dat_iid,iid),"lsoa_label"] <- area_dat$geo_label
ph[match(coord_dat_iid,iid),"townsend"] <- area_dat$town_di

# Map easting northing points onto msoa polygons
area_dat <- coord_sp %over% sp_dist

# Add results onto original data frame
ph[match(coord_dat_iid,iid),"dist"] <- area_dat$geo_code
ph[match(coord_dat_iid,iid),"dist_label"] <- area_dat$geo_label




#----
# Classify scores
#----

mild_positive <- classifier[hospital==0 & critical==0 & death==0 & death_hospital==0, eid]
scores[,mild_positive:=iid %in% mild_positive]
scores[,negative:=!iid %in% classifier$eid]


#----
# Plot Severe vs. Mild
#----

density1 <- scores[study=="genomicc",density(score, n=1000)]
density2 <- scores[mild_positive==TRUE, density(score, n=1000)]
density_df <- rbind(data.table(x=density1$x, y=density1$y, study="GenOMICC"), data.table(x=density2$x, y=density2$y, study="UK Biobank"))

text_df <- density_df[,.SD[which.max(y),],by=study]
text_df$label <- c(density1$n, density2$n)

rf <- density_df[,.(x=pretty(x,4), y=pretty(y,4))]


p1 <- ggplot(density_df, aes(x=x, y=y)) + 
geom_line(aes(colour=study), size=1.2) +
geom_polygon(aes(fill=study), alpha=0.25) +
geom_text_repel(data=text_df, aes(colour=study, label=label), nudge_y=0.01) +
geom_rangeframe(data=rf) +
scale_colour_manual(values=c("#b2182b","#053061"), name="", guide="none") + 
scale_fill_manual(values=c("#b2182b","#053061"), labels=c("Severe","Mild"), name="") + 
guides(fill = guide_legend(override.aes = list(alpha=1))) +
scale_x_continuous(breaks=unique(rf$x)) +
scale_y_continuous(breaks=unique(rf$y)) +
coord_cartesian(xlim=range(rf$x), ylim=range(rf$y)) +
labs(x="Cystatin-C PRS\n(Z score)", y="Density") +
theme_pubclean() + 
theme(legend.position="right")



#----
# Plot Severe vs. Negative
#----

density1 <- scores[study=="genomicc",density(score, n=1000)]
density2 <- scores[study=="ukbb" & negative==TRUE, density(score, n=1000)]
density_df <- rbind(data.table(x=density1$x, y=density1$y, study="GenOMICC"), data.table(x=density2$x, y=density2$y, study="UK Biobank"))

text_df <- density_df[,.SD[which.max(y),],by=study]
text_df$label <- c(density1$n, density2$n)

rf <- density_df[,.(x=pretty(x,4), y=pretty(y,4))]


p2 <- ggplot(density_df, aes(x=x, y=y)) + 
geom_line(aes(colour=study), size=1.2) +
geom_polygon(aes(fill=study), alpha=0.25) +
geom_text_repel(data=text_df, aes(colour=study, label=label), nudge_y=0.01) +
geom_rangeframe(data=rf) +
scale_colour_manual(values=c("#b2182b","#1a1a1a"), name="", guide="none") + 
scale_fill_manual(values=c("#b2182b","#1a1a1a"), labels=c("Severe","Negative"), name="") + 
guides(fill = guide_legend(override.aes = list(alpha=1))) +
scale_x_continuous(breaks=unique(rf$x)) +
scale_y_continuous(breaks=unique(rf$y)) +
coord_cartesian(xlim=range(rf$x), ylim=range(rf$y)) +
labs(x="Cystatin-C PRS\n(Z score)", y="Density") +
theme_pubclean() + 
theme(legend.position="right")


p12 <- ggarrange(p1, p2, nrow=2, align="hv")

pdf("st06_01_score_density.pdf", width=6, height=8)
plot(p12)
dev.off()




#----
# Model
#----

dat <- scores[ph,,on="iid"]
dat <- genetic_pcs[dat,,on="iid"]

# Define variables
dat[,age1 := ifelse(age > 40, age, NA)]
dat[,icu := study=="genomicc" & covid_status == 1]
dat[,is_positive_only := ifelse(covid_status == 1, TRUE, NA)]
dat[,is_icu_or_negative := ifelse(icu == 1 | covid_status == 0, TRUE, NA)]

# Define local authorities
dist_size <- dat[,.N,by="dist"][order(-N)]
dist_common <- dist_size[N > 50, as.character(dist)]
dat[,local_authority := ifelse(dist %in% dist_common, as.character(dist_label), NA)]

writeLines(paste0("Fitting ",length(dist_common)," common districts"))
writeLines(paste0("Reference district: ",dat[,sort(local_authority)[1]],"\n"))


# Define score deciles
deciles <- quantile(dat$score, probs=seq(0.1,0.9,0.1), na.rm=T)
dat[,score_dec := factor(as.numeric(cut(score, breaks=c(-Inf,deciles,Inf))))]


# ICU vs. mild-positive
#----------------------

res <- data.table()
res_dec <- data.table()

for (s in c("score")) { # for (s in c("score", "score_strict")) {

# All

f1 <- paste0("icu ~ is_positive_only + age + sex + local_authority + townsend + ",paste0(paste0("pc",1:10), collapse=" + "), " + ", s)
m1 <- try(speedglm(formula=as.formula(f1), family=binomial(link="logit"), data=dat, model=TRUE))
if("try-error" %in% class(m1)) m1 <- glm(formula=as.formula(f1), family=binomial(link="logit"), data=dat, model=TRUE)
s1 <- summary(m1); s1$call[[2]] <- as.formula(f1); print(s1)

x1 <- data.table(m1$model)
res1 <- data.table(analysis="severe_vs_mild", n_cases=sum(x1$icu), n_controls=sum(!x1$icu), n_snps=scores_desc[score==s,n_snps], s1$coefficients[s,])
names(res1)[5:8] <- c("beta","se","z","p")

model1 <- s1$coefficients[complete.cases(s1$coefficients),]
model1 <- data.table(variable=gsub("local_authority","",rownames(model1)), model1)
colnames(model1) <- c("variable", "beta", "se", "z", "p")
res <- rbind(res, res1)

# 40+

f1p <- paste0("icu ~ is_positive_only + age1 + sex + local_authority + townsend + ",paste0(paste0("pc",1:10), collapse=" + "), " + ", s)
m1p <- try(speedglm(formula=as.formula(f1p), family=binomial(link="logit"), data=dat, model=TRUE))
if("try-error" %in% class(m1p)) m1p <- glm(formula=as.formula(f1p), family=binomial(link="logit"), data=dat, model=TRUE)
s1p <- summary(m1p); s1p$call[[2]] <- as.formula(f1p)

x1p <- data.table(m1p$model)
res1p <- data.table(analysis="severe_vs_mild_40_plus", n_cases=sum(x1p$icu), n_controls=sum(!x1p$icu), n_snps=scores_desc[score==s,n_snps], s1p$coefficients[s,,drop=FALSE])
names(res1p)[5:8] <- c("beta","se","z","p")

model1p <- s1p$coefficients[complete.cases(s1p$coefficients),]
model1p <- data.table(variable=gsub("local_authority","",rownames(model1p)), model1p)
colnames(model1p) <- c("variable", "beta", "se", "z", "p")
res <- rbind(res, res1p)



# ICU vs. negative
#-----------------

# All

f2 <- paste0("icu ~ is_icu_or_negative + age + sex + local_authority + townsend + ",paste0(paste0("pc",1:10), collapse=" + "), " + ", s)
m2 <- try(speedglm(formula=as.formula(f2), family=binomial(link="logit"), data=dat, model=TRUE))
if("try-error" %in% class(m2)) m2 <- glm(formula=as.formula(f2), family=binomial(link="logit"), data=dat, model=TRUE)
s2 <- summary(m2); s2$call[[2]] <- as.formula(f2); print(s2)

x2 <- data.table(m2$model)
res2 <- data.table(analysis="severe_vs_negative", n_cases=sum(x2$icu), n_controls=sum(!x2$icu), n_snps=scores_desc[score==s,n_snps], s2$coefficients[s,,drop=FALSE])
names(res2)[5:8] <- c("beta","se","z","p")

model2 <- s2$coefficients[complete.cases(s2$coefficients),]
model2 <- data.table(variable=gsub("local_authority","",rownames(model2)), model2)
colnames(model2) <- c("variable", "beta", "se", "z", "p")
res <- rbind(res, res2)

# 40+

f2p <- paste0("icu ~ is_icu_or_negative + age1 + sex + local_authority + townsend + ",paste0(paste0("pc",1:10), collapse=" + "), " + ", s)
m2p <- try(speedglm(formula=as.formula(f2p), family=binomial(link="logit"), data=dat, model=TRUE))
if("try-error" %in% class(m2p)) m2p <- glm(formula=as.formula(f2p), family=binomial(link="logit"), data=dat, model=TRUE)
s2p <- summary(m2p); s2p$call[[2]] <- as.formula(f2p); print(s2p)

x2p <- data.table(m2p$model)
res2p <- data.table(analysis="severe_vs_negative_40_plus", n_cases=sum(x2p$icu), n_controls=sum(!x2p$icu), n_snps=scores_desc[score==s,n_snps], s2p$coefficients[s,,drop=FALSE])
names(res2p)[5:8] <- c("beta","se","z","p")

model2p <- s2p$coefficients[complete.cases(s2p$coefficients),]
model2p <- data.table(variable=gsub("local_authority","",rownames(model2p)), model2p)
colnames(model2p) <- c("variable", "beta", "se", "z", "p")
res <- rbind(res, res2p)


#----
# Deciles
#----

f3 <- paste0("icu ~ is_icu_or_negative + age + sex + local_authority + townsend + ",paste0(paste0("pc",1:10), collapse=" + "), " + ", s, "_dec")
m3 <- try(speedglm(formula=as.formula(f3), family=binomial(link="logit"), data=dat, model=TRUE))
if("try-error" %in% class(m3)) m3 <- glm(formula=as.formula(f3), family=binomial(link="logit"), data=dat, model=TRUE)
s3 <- summary(m3); s3$call[[2]] <- as.formula(f3); print(s3)


x3 <- data.table(m3$model)
res3 <- data.table(analysis=factor(paste0(s,"_dec",2:10),levels=paste0(s,"_dec",1:10)), n_cases=sum(x2p$icu), n_controls=sum(!x2p$icu), n_snps=scores_desc[score==s,n_snps], s3$coefficients[paste0(s,"_dec",2:10),,drop=FALSE])
names(res3)[5:8] <- c("beta","se","z","p")

model3 <- s3$coefficients[complete.cases(s3$coefficients),]
model3 <- data.table(variable=gsub("local_authority","",rownames(model3)), model3)
colnames(model3) <- c("variable", "beta", "se", "z", "p")

ggplot(res3, aes(x=beta, y=analysis)) + 
geom_errorbarh(aes(xmin=beta - qnorm(0.975)*se, xmax=beta + qnorm(0.975)*se)) + 
geom_point()


res_dec <- rbind(res_dec, res3)

#----
# Create tables
#----


fwrite(model1, paste0("st06_02_model1_",s,"_coef.csv"), sep=",", na="NA", quote=TRUE)
fwrite(model1p, paste0("st06_02_model1_",s,"_coef_40_plus.csv"), sep=",", na="NA", quote=TRUE)
fwrite(model2, paste0("st06_02_model2_",s,"_coef.csv"), sep=",", na="NA", quote=TRUE)
fwrite(model2p, paste0("st06_02_model2_",s,"_coef_40_plus.csv"), sep=",", na="NA", quote=TRUE)

fwrite(res, paste0("st06_03_score_coef.csv"), sep=",", na="NA", quote=FALSE)

#----
# Get descriptives
#----


names(x1)[ncol(x1)] <- "score"
d1 <- rbind(
	x1[,.(icu="ALL", sex="ALL", n_total=.N, n_female=sum(sex=="Female"), min_age=min(age), max_age=max(age), mean_age=mean(age), sd_age=sd(age), mean_score=mean(score), sd_score=sd(score))],
	x1[,.(icu="ALL", n_total=.N, n_female=c(Female=.N, Male=0)[sex], min_age=min(age), max_age=max(age), mean_age=mean(age), sd_age=sd(age), mean_score=mean(score), sd_score=sd(score)), by="sex"],
	x1[,.(sex="ALL", n_total=.N, n_female=sum(sex=="Female"), min_age=min(age), max_age=max(age), mean_age=mean(age), sd_age=sd(age), mean_score=mean(score), sd_score=sd(score)), by="icu"],
	x1[,.(n_total=.N, n_female=c(Female=.N, Male=0)[sex], min_age=min(age), max_age=max(age), mean_age=mean(age), sd_age=sd(age), mean_score=mean(score), sd_score=sd(score)),by=c("icu","sex")]
)

names(x1p)[ncol(x1p)] <- "score"
d1p <- rbind(
	x1p[,.(icu="ALL", sex="ALL", n_total=.N, n_female=sum(sex=="Female"), min_age=min(age1), max_age=max(age1), mean_age=mean(age1), sd_age=sd(age1), mean_score=mean(score), sd_score=sd(score))],
	x1p[,.(icu="ALL", n_total=.N, n_female=c(Female=.N, Male=0)[sex], min_age=min(age1), max_age=max(age1), mean_age=mean(age1), sd_age=sd(age1), mean_score=mean(score), sd_score=sd(score)), by="sex"],
	x1p[,.(sex="ALL", n_total=.N, n_female=sum(sex=="Female"), min_age=min(age1), max_age=max(age1), mean_age=mean(age1), sd_age=sd(age1), mean_score=mean(score), sd_score=sd(score)), by="icu"],
	x1p[,.(n_total=.N, n_female=c(Female=.N, Male=0)[sex], min_age=min(age1), max_age=max(age1), mean_age=mean(age1), sd_age=sd(age1), mean_score=mean(score), sd_score=sd(score)),by=c("icu","sex")]
)

names(x2)[ncol(x2)] <- "score"
d2 <- rbind(
	x2[,.(icu="ALL", sex="ALL", n_total=.N, n_female=sum(sex=="Female"), min_age=min(age), max_age=max(age), mean_age=mean(age), sd_age=sd(age), mean_score=mean(score), sd_score=sd(score))],
	x2[,.(icu="ALL", n_total=.N, n_female=c(Female=.N, Male=0)[sex], min_age=min(age), max_age=max(age), mean_age=mean(age), sd_age=sd(age), mean_score=mean(score), sd_score=sd(score)), by="sex"],
	x2[,.(sex="ALL", n_total=.N, n_female=sum(sex=="Female"), min_age=min(age), max_age=max(age), mean_age=mean(age), sd_age=sd(age), mean_score=mean(score), sd_score=sd(score)), by="icu"],
	x2[,.(n_total=.N, n_female=c(Female=.N, Male=0)[sex], min_age=min(age), max_age=max(age), mean_age=mean(age), sd_age=sd(age), mean_score=mean(score), sd_score=sd(score)),by=c("icu","sex")]
)

names(x2p)[ncol(x2p)] <- "score"
d2p <- rbind(
	x2p[,.(icu="ALL", sex="ALL", n_total=.N, n_female=sum(sex=="Female"), min_age=min(age1), max_age=max(age1), mean_age=mean(age1), sd_age=sd(age1), mean_score=mean(score), sd_score=sd(score))],
	x2p[,.(icu="ALL", n_total=.N, n_female=c(Female=.N, Male=0)[sex], min_age=min(age1), max_age=max(age1), mean_age=mean(age1), sd_age=sd(age1), mean_score=mean(score), sd_score=sd(score)), by="sex"],
	x2p[,.(sex="ALL", n_total=.N, n_female=sum(sex=="Female"), min_age=min(age1), max_age=max(age1), mean_age=mean(age1), sd_age=sd(age1), mean_score=mean(score), sd_score=sd(score)), by="icu"],
	x2p[,.(n_total=.N, n_female=c(Female=.N, Male=0)[sex], min_age=min(age1), max_age=max(age1), mean_age=mean(age1), sd_age=sd(age1), mean_score=mean(score), sd_score=sd(score)),by=c("icu","sex")]
)

d1$analysis <- "severe_vs_mild"
d1p$analysis <- "severe_vs_mild_40_plus"
d2$analysis <- "severe_vs_negative"
d2p$analysis <- "severe_vs_negative_40_plus"

desc <- rbind(d1, d1p, d2, d2p)

fwrite(desc, paste0("st06_00_model_",s,"_descriptives.csv"), sep=",", na="NA", quote=FALSE)
}