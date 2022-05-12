#######################

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  Sys.time()

f1 <- read.table(args[1], header=T) #### PRS
f2 <- read.delim(args[2],header=T) ####  PCA
f3 <- read.delim(args[3],header=T) ####  Pheno


#########
f33 <- merge(f1,f2,by.x="IID", by.y="IID")

f4 <- read.delim(args[4],header=T) ####  Sex and Age



f3 <- merge(f33,f3,by.x="IID", by.y="IID")

f3 <- merge(f4,f3,by.x="IID", by.y="IID")


print(dim(f3))

f3$SCORESUM <- scale(f3$SCORESUM)
f3 <- f3[which(f3$Sex!="Unknown"),]

#print(head(f3))

########
f3$PC1 <- f3$PC2
f3$PC2 <- f3$PC3
f3$PC3 <- f3$PC4
f3$PC4 <- f3$PC5
f3$PC5 <- f3$PC6
f3$PC6 <- f3$PC7
f3$PC7 <- f3$PC8
f3$PC8 <- f3$PC9





#print(summary(glm(Pheno ~ SCORESUM + PC1 + PC2 + PC3 + as.factor(Sex) + Age, data=f3, family = "binomial")))

print(summary(glm(Pheno ~ SCORESUM +Age + as.factor(Sex) + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8, data=f3, family = "binomial")))
#print(summary(glm(Pheno ~ SCORE + PC1 + PC2 + SEX ,data=f3, family = "binomial")))

########
#print(summary(glm(Pheno~ SCORE,data=f3, family = "binomial")))


#print(head(f3))

write.table(f3, "PRS_PCS_Pheno_v1", col.names=TRUE,row.names=FALSE, append=FALSE, quote=FALSE, sep="\t")

case <-f3[which(f3$Pheno==1),]
ctl <- f3[which(f3$Pheno==0),]

print(t.test(case$SCORESUM,ctl$SCORESUM))
print(table(f3$Pheno))
}
