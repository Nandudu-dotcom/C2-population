
#creating deregressed_BLUPS
rm(list=ls())
setwd("~/Desktop/phd_Defense/2021/2022_readings/Publications/GWAS_1_5_scores/Data_sets_and_code/BLUPS/deregressed_BLUPs")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(corrplot)
library(lme4)
library(minque)
library(tidyverse)
library(rrBLUP)
library(dplyr)
library(corrplot)
library(qqman)
library(tidyr)
library(lme4)
library(magrittr)
library(dplyr)
library(ggplot2)
library(EMMREML)
CETS <-read.csv("~/Desktop/phd_Defense/2021/2022_readings/Publications/GWAS_1_5_scores/Data_sets_and_code/Variance_and_convariance/Data/CETS.csv")
dim(CETS) # 1662   36
CMD_phenotypes <- read.csv("~/Desktop/phd_Defense/2021/2022_readings/Publications/Data_sets/CMD_phenotypes.csv")
dim(CMD_phenotypes) #1659    9

CET <- merge(CETS,CMD_phenotypes , by = c("studyYear", "locationName","germplasmName","blockNumber",'plotNumber'))
write.csv(CET , file ="CETS.csv" )
CETS <- CET 
dim(CETS) #1658   42


CETS$Env <- paste(CETS$locationName,CETS$studyYear)
CETS$Block <- paste(CETS$locationName,CETS$studyYear,CETS$blockNumber)
CETS[is.na(CETS)] = 0

closeAllConnections()
print(xtabs(~blockNumber +locationName +studyYear+ germplasmName, CETS))

#making the clone random but also separating the check and clone effects 
#WE will also fix block effects and calculate the habitability
# Adding a new column that will help us divide data set into two groups with one group being new entries and one group being checks
closeAllConnections()
data2 <- within(CETS, {
  new <- ifelse(CETS$entryType!="check",0,1)
});
tail(data2)

#check is 1 while the test clone is 0 

data2 <- within(data2, {
  entryc <- ifelse(data2$new > 0, 999,data2$germplasmName)})
str(data2)
tail(data2)

data2$blockNumber <- as.factor(data2$blockNumber)
data2$Env <- as.factor(data2$Env)
data2$germplasmName<- as.factor(data2$germplasmName)
data2$new <- as.factor(data2$new)
data2$entryc <- as.factor(data2$entryc)
str(data2)  # Display the structure of an arbitrary R Object compactly

#randomize clone 
#fix entryc that accouts for the checks vs clone effects,ENV, grand mean




# Deregressed_BLUPS -------------------------------------------------------


model1 <- lmer(CBSDi3 ~ entryc  +(1|germplasmName:new) +(1|Env) +(1|germplasmName:Env) , data=data2, REML = T)
modelSum <- summary(model1)
print (modelSum)

BLUP <- ranef(model1, condVar=TRUE)$`germplasmName:new`
write.csv(BLUP, file = "Combined_BLUPs_CBSDi3.csv")
PEV <- c(attr(BLUP, "postVar"))
VARCOMPS <- as.data.frame(VarCorr(model1),comp="Variance") 
CLONEbyENV <- VARCOMPS[1,4]
CLONE <- VARCOMPS[2,4]
ENV <- VARCOMPS[3,4]
ResidVar <- (VARCOMPS[4,4])
#deregressing OF blupS
out <- BLUP/(1-(PEV/CLONE))

write.csv(out, file = "Combined_deregressedBLUPs_CBSDi3.csv")

model2 <- lmer(CBSDi6 ~ entryc  +(1|germplasmName:new) +(1|Env) +(1|germplasmName:Env) , data=data2, REML = T)
modelSum <- summary(model2)
print (modelSum)

BLUP <- ranef(model2, condVar=TRUE)$`germplasmName:new`
write.csv(BLUP, file = "Combined_BLUPs_CBSDi6.csv")
PEV <- c(attr(BLUP, "postVar"))
VARCOMPS <- as.data.frame(VarCorr(model2),comp="Variance") 
CLONEbyENV <- VARCOMPS[1,4]
CLONE <- VARCOMPS[2,4]
ENV <- VARCOMPS[3,4]
ResidVar <- (VARCOMPS[4,4])
#deregressing OF blupS
out <- BLUP/(1-(PEV/CLONE))

write.csv(out, file = "Combined_deregressedBLUPs_CBSDi6.csv")
model3 <- lmer(CBSDi12 ~ entryc  +(1|germplasmName:new) +(1|Env) +(1|germplasmName:Env) , data=data2, REML = T)
modelSum <- summary(model3)
print (modelSum)

BLUP <- ranef(model3, condVar=TRUE)$`germplasmName:new`
write.csv(BLUP, file = "Combined_BLUPs_CBSDi12.csv")
PEV <- c(attr(BLUP, "postVar"))
VARCOMPS <- as.data.frame(VarCorr(model3),comp="Variance") 
CLONEbyENV <- VARCOMPS[1,4]
CLONE <- VARCOMPS[2,4]
ENV <- VARCOMPS[3,4]
ResidVar <- (VARCOMPS[4,4])
#deregressing OF blupS
out <- BLUP/(1-(PEV/CLONE))

write.csv(out, file = "Combined_deregressedBLUPs_CBSDi12.csv")

model4 <- lmer(CBSDs3 ~ entryc  +(1|germplasmName:new) +(1|Env) +(1|germplasmName:Env) , data=data2, REML = T)
modelSum <- summary(model4)
print (modelSum)

BLUP <- ranef(model4, condVar=TRUE)$`germplasmName:new`
write.csv(BLUP, file = "Combined_BLUPs_CBSDs3.csv")
PEV <- c(attr(BLUP, "postVar"))
VARCOMPS <- as.data.frame(VarCorr(model4),comp="Variance") 
CLONEbyENV <- VARCOMPS[1,4]
CLONE <- VARCOMPS[2,4]
ENV <- VARCOMPS[3,4]
ResidVar <- (VARCOMPS[4,4])
#deregressing OF blupS
out <- BLUP/(1-(PEV/CLONE))

write.csv(out, file = "Combined_deregressedBLUPs_CBSDs3.csv")

model5 <- lmer(CBSDs6 ~ entryc  +(1|germplasmName:new) +(1|Env) +(1|germplasmName:Env) , data=data2, REML = T)
modelSum <- summary(model5)
print (modelSum)

BLUP <- ranef(model5, condVar=TRUE)$`germplasmName:new`
write.csv(BLUP, file = "Combined_BLUPs_CBSDs6.csv")
PEV <- c(attr(BLUP, "postVar"))
VARCOMPS <- as.data.frame(VarCorr(model5),comp="Variance") 
CLONEbyENV <- VARCOMPS[1,4]
CLONE <- VARCOMPS[2,4]
ENV <- VARCOMPS[3,4]
ResidVar <- (VARCOMPS[4,4])
#deregressing OF blupS
out <- BLUP/(1-(PEV/CLONE))

write.csv(out, file = "Combined_deregressedBLUPs_CBSDs6.csv")

model6 <- lmer(CBSDs12 ~ entryc  +(1|germplasmName:new) +(1|Env) +(1|germplasmName:Env) , data=data2, REML = T)
modelSum <- summary(model6)
print (modelSum)

BLUP <- ranef(model6, condVar=TRUE)$`germplasmName:new`
write.csv(BLUP, file = "Combined_BLUPs_CBSDs12.csv")
PEV <- c(attr(BLUP, "postVar"))
VARCOMPS <- as.data.frame(VarCorr(model6),comp="Variance") 
CLONEbyENV <- VARCOMPS[1,4]
CLONE <- VARCOMPS[2,4]
ENV <- VARCOMPS[3,4]
ResidVar <- (VARCOMPS[4,4])
#deregressing OF blupS
out <- BLUP/(1-(PEV/CLONE))

write.csv(out, file = "Combined_deregressedBLUPs_CBSDs12.csv")

model7 <- lmer(CMDi3 ~ entryc  +(1|germplasmName:new) +(1|Env) +(1|germplasmName:Env) , data=data2, REML = T)
modelSum <- summary(model7)
print (modelSum)
BLUP <- ranef(model7, condVar=TRUE)$`germplasmName:new`
write.csv(BLUP, file = "Combined_BLUPs_CMDi3.csv")
PEV <- c(attr(BLUP, "postVar"))
VARCOMPS <- as.data.frame(VarCorr(model7),comp="Variance") 
CLONEbyENV <- VARCOMPS[1,4]
CLONE <- VARCOMPS[2,4]
ENV <- VARCOMPS[3,4]
ResidVar <- (VARCOMPS[4,4])
#deregressing OF blupS
out <- BLUP/(1-(PEV/CLONE))

write.csv(out, file = "Combined_deregressedBLUPs_CMDi3.csv")


model8 <- lmer(CMDi6 ~ entryc  +(1|germplasmName:new) +(1|Env) +(1|germplasmName:Env) , data=data2, REML = T)
modelSum <- summary(model8)
print (modelSum)



BLUP <- ranef(model8, condVar=TRUE)$`germplasmName:new`
write.csv(BLUP, file = "Combined_BLUPs_CMDi6.csv")
PEV <- c(attr(BLUP, "postVar"))
VARCOMPS <- as.data.frame(VarCorr(model8),comp="Variance") 
CLONEbyENV <- VARCOMPS[1,4]
CLONE <- VARCOMPS[2,4]
ENV <- VARCOMPS[3,4]
ResidVar <- (VARCOMPS[4,4])
#deregressing OF blupS
out <- BLUP/(1-(PEV/CLONE))

write.csv(out, file = "Combined_deregressedBLUPs_CMDi6.csv")



model9 <- lmer(CMDs3 ~ entryc  +(1|germplasmName:new) +(1|Env) +(1|germplasmName:Env) , data=data2, REML = T)
modelSum <- summary(model9)
print (modelSum)



BLUP <- ranef(model9, condVar=TRUE)$`germplasmName:new`
write.csv(BLUP, file = "Combined_BLUPs_CMDs3.csv")
PEV <- c(attr(BLUP, "postVar"))
VARCOMPS <- as.data.frame(VarCorr(model9),comp="Variance") 
CLONEbyENV <- VARCOMPS[1,4]
CLONE <- VARCOMPS[2,4]
ENV <- VARCOMPS[3,4]
ResidVar <- (VARCOMPS[4,4])
#deregressing OF blupS
out <- BLUP/(1-(PEV/CLONE))

write.csv(out, file = "Combined_deregressedBLUPs_CMDs3.csv")



model10 <- lmer(CMDs6 ~ entryc  +(1|germplasmName:new) +(1|Env) +(1|germplasmName:Env) , data=data2, REML = T)
modelSum <- summary(model10)
print (modelSum)

BLUP <- ranef(model10, condVar=TRUE)$`germplasmName:new`
write.csv(BLUP, file = "Combined_BLUPs_CMDs6.csv")
PEV <- c(attr(BLUP, "postVar"))
VARCOMPS <- as.data.frame(VarCorr(model10),comp="Variance") 
CLONEbyENV <- VARCOMPS[1,4]
CLONE <- VARCOMPS[2,4]
ENV <- VARCOMPS[3,4]
ResidVar <- (VARCOMPS[4,4])
#deregressing OF blupS
out <- BLUP/(1-(PEV/CLONE))

write.csv(out, file = "Combined_deregressedBLUPs_CMDs6.csv")



# Narrow sense heritability -----------------------------------------------


#narrow sense heritability 

Marker_data <- readRDS("~/Desktop/phd_Defense/2021/marker_data/New_combined_marker_data/New_marker_data_51862_SNPs.rds")
dim(Marker_data) #51862  1598
head(Marker_data[1:100,1:100])

#if i transpose directly, the marker column name appears as a row name creating NAs 
####changing marker data to be coded as -1,0,1
#Eariler the error was due to R considering the marker column and everything in them as characters 
#plan , remove the marker row and then remove the first column and then transpose
#removing SNP names 
SNP <- Marker_data$MARKER
####tracking filtering
rownames(Marker_data) <- SNP
#removing the first column 
markerdataNoSNP <- Marker_data[,-1]
#transposing and changing the code to 1,0,-1
Markerdata1= t(markerdataNoSNP -1)
dim(Markerdata1) #1597 51862

#works well 

# Data filtering  ---------------------------------------------------------

#### filtering for missing individuals and missing markers 
Data_filtered <- function(Markerdata1, NIM, NMG){
  
  missing_indids <-apply(Markerdata1,1, function(x){
    return(length(which(is.na(x)))/ncol(Markerdata1))
  })
  
  missing_genos <- apply(Markerdata1,2,function(x) {
    return(length(which(is.na(x)))/nrow(Markerdata1))
  })
  
  filtering <- Markerdata1[which( missing_indids< NIM),which(missing_genos < NMG)]
  
  return(filtering)
}


#remove individuals with >10% data and makers with >5% missing data
filtered <-  Data_filtered (Markerdata1,0.10,0.05)
dim(filtered) #  1597 51862
dim(Markerdata1)# 1597 51862  

###from the filtering no genotype is filtered


# Filtering for minor allele frequency ------------------------------------


#Calculating allele freq and filtering for minor allele frequency using the filtered marker data 

MAf_filtered <-filtered + 1

alleleFreq <- colMeans(MAf_filtered )/2

hist(alleleFreq)

summary(alleleFreq)

MAF <- apply(MAf_filtered, 2, function(x) sum(x) / (length(x) * 2))

#filtering for MAF

Genotypes_filtered <- filtered[,which(MAF > 0.05 & MAF < 0.95)]

dim(Genotypes_filtered)  #1597 41165 Reduction after filtering for  MAF 
write_rds(Genotypes_filtered,file = "Marker_data.rds")
##########checking if there are no Nas 
sum(is.na(Genotypes_filtered)) # =0 

length(which(is.na(Genotypes_filtered)))  # =0 



# Phenotypic data ---------------------------------------------------------

CETS <-read.csv("~/Desktop/phd_Defense/2021/2022_readings/Publications/Data_sets/CETs_data.csv")
dim(CETS) # 1658   44

common_code <- read.csv("~/Desktop/phd_Defense/2021/2022_readings/Publications/ELISA/Data_sets/Set_1/common_code.csv")
dim(common_code) #6570    2

names(common_code)[names(common_code) == "accession_name"] <- "germplasmName"
CETS<-merge(CETS,common_code, by="germplasmName")
dim(CETS) #1500   45


commonInd <- intersect(rownames(Genotypes_filtered),CETS$observationUnitName.y) ####commonindividuals 
length(commonInd)  #302


C2_Geno_1_subset <- Genotypes_filtered[(rownames(Genotypes_filtered)%in%commonInd),] 
dim(C2_Geno_1_subset) # 302 41165



Phenotypes_subset <- CETS[CETS$observationUnitName.y %in% commonInd,]
dim(Phenotypes_subset) #  1044   45


# Additive relationship matrix --------------------------------------------
#when working with sommer the variables have to be the same length
CETS$blo <-paste(CETS$locationName,CETS$studyYear)
CETS$ENV <-paste(CETS$locationName,CETS$studyYear,CETS$blockNumber)
CETS$blockNumber <- as.factor(CETS$blockNumber)
#TO ENSURE THAT NAs don't affect the length of the traits, i will replace the NAs with 0


#cbsdi3
CETS_CBSDi3 <-CETS[,c(45,47,31)]
dim(CETS_CBSDi3 ) #1503    3
CETS_CBSDi3 [1:5,]
# remove NAs
CETS_CBSDi3_na <- na.omit(CETS_CBSDi3)
dim(CETS_CBSDi3_na) #1476    3
## dosage file with common clones
CETS_dosageC2_123format_CBSDi3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%CETS_CBSDi3_na$observationUnitName.y,]
CETS_dosageC2_123format_CBSDi3[1:5,1:5]
dim(CETS_dosageC2_123format_CBSDi3) #297 41165

CETS_CBSDi3_sub <- CETS_CBSDi3_na[CETS_CBSDi3_na$observationUnitName.y%in%rownames(CETS_dosageC2_123format_CBSDi3),]
dim(CETS_CBSDi3_sub) #1039    3

CETS_CBSDi3k <-A.mat(CETS_dosageC2_123format_CBSDi3-1)
dim(CETS_CBSDi3k) #  297 297
Z_model_CETS_CBSDi3<-model.matrix(~observationUnitName.y-1, data = CETS_CBSDi3_sub)
dim(Z_model_CETS_CBSDi3) # 1039  297
#making block as a factor
CETS_CBSDi3_sub$ENV <- as.factor(CETS_CBSDi3_sub$ENV)

X_model_CETS_CBSDi3 <- model.matrix(~ENV,data = CETS_CBSDi3_sub)

#Narrow sense heritability 
Model_CETS_CBSDi3 <-emmreml(y=CETS_CBSDi3_sub$CBSDi3, X=X_model_CETS_CBSDi3, Z=Z_model_CETS_CBSDi3, K=CETS_CBSDi3k)

Model_CETS_CBSDi3$Vu/(Model_CETS_CBSDi3$Vu+Model_CETS_CBSDi3$Ve) 
# Narrowsense 0.3347411
#Extract blups
CETS_CBSDi3_blups <- Model_CETS_CBSDi3$uhat
write.csv(CETS_CBSDi3_blups, file = "CETS_CBSDi3emmreml_blups.csv")


#CBSDi6
CETS_CBSDi6 <-CETS[,c(45,47,32)]
dim(CETS_CBSDi6 ) #1500    3
CETS_CBSDi6 [1:5,]
# remove NAs
CETS_CBSDi6_na <- na.omit(CETS_CBSDi6)
dim(CETS_CBSDi6_na) #1488    3
## dosage file with common clones
CETS_dosageC2_123format_CBSDi6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%CETS_CBSDi6_na$observationUnitName.y,]
CETS_dosageC2_123format_CBSDi6[1:5,1:5]
dim(CETS_dosageC2_123format_CBSDi6) # 302 41165

CETS_CBSDi6_sub <- CETS_CBSDi6_na[CETS_CBSDi6_na$observationUnitName.y%in%rownames(CETS_dosageC2_123format_CBSDi6),]
dim(CETS_CBSDi6_sub) # 1037    3

CETS_CBSDi6k <-A.mat(CETS_dosageC2_123format_CBSDi6-1)
dim(CETS_CBSDi6k) #302 302
Z_model_CETS_CBSDi6<-model.matrix(~observationUnitName.y-1, data = CETS_CBSDi6_sub)
dim(Z_model_CETS_CBSDi6) # 1037  302
#making block as a factor
CETS_CBSDi6_sub$ENV <- as.factor(CETS_CBSDi6_sub$ENV)

X_model_CETS_CBSDi6 <- model.matrix(~ENV,data = CETS_CBSDi6_sub)

#Narrow sense heritability 
Model_CETS_CBSDi6 <-emmreml(y=CETS_CBSDi6_sub$CBSDi6, X=X_model_CETS_CBSDi6, Z=Z_model_CETS_CBSDi6, K=CETS_CBSDi6k)

Model_CETS_CBSDi6$Vu/(Model_CETS_CBSDi6$Vu+Model_CETS_CBSDi6$Ve) 
# Narrowsense  0.3808686
#Extract blups
CETS_CBSDi6_blups <- Model_CETS_CBSDi6$uhat
write.csv(CETS_CBSDi6_blups, file = "CETS_CBSDi6emmreml_blups.csv")



#CBSDi12
CETS_CBSDi12 <-CETS[,c(45,47,35)]
dim(CETS_CBSDi12 ) # 1500    3
CETS_CBSDi12 [1:5,]
# remove NAs
CETS_CBSDi12_na <- na.omit(CETS_CBSDi12)
dim(CETS_CBSDi12_na) #  866   3
## dosage file with common clones
CETS_dosageC2_123format_CBSDi12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%CETS_CBSDi12_na$observationUnitName.y,]
CETS_dosageC2_123format_CBSDi12[1:5,1:5]
dim(CETS_dosageC2_123format_CBSDi12) #270 41165

CETS_CBSDi12_sub <- CETS_CBSDi12_na[CETS_CBSDi12_na$observationUnitName.y%in%rownames(CETS_dosageC2_123format_CBSDi12),]
dim(CETS_CBSDi12_sub) #647   3

CETS_CBSDi12k <-A.mat(CETS_dosageC2_123format_CBSDi12-1)
dim(CETS_CBSDi12k) # 270 270
Z_model_CETS_CBSDi12<-model.matrix(~observationUnitName.y-1, data = CETS_CBSDi12_sub)
dim(Z_model_CETS_CBSDi12) #  647 270
#making block as a factor
CETS_CBSDi12_sub$ENV <- as.factor(CETS_CBSDi12_sub$ENV)

X_model_CETS_CBSDi12 <- model.matrix(~ENV,data = CETS_CBSDi12_sub)

#Narrow sense heritability 
Model_CETS_CBSDi12 <-emmreml(y=CETS_CBSDi12_sub$CBSDi12, X=X_model_CETS_CBSDi12, Z=Z_model_CETS_CBSDi12, K=CETS_CBSDi12k)

Model_CETS_CBSDi12$Vu/(Model_CETS_CBSDi12$Vu+Model_CETS_CBSDi12$Ve) 
# Narrowsense  0.1952025
#Extract blups
CETS_CBSDi12_blups <- Model_CETS_CBSDi12$uhat
write.csv(CETS_CBSDi12_blups, file = "CETS_CBSDi12emmreml_blups.csv")



#cbsds3
CETS_CBSDs3 <-CETS[,c(45,47,33)]
dim(CETS_CBSDs3 ) #1500    3
CETS_CBSDs3 [1:5,]
# remove NAs
CETS_CBSDs3_na <- na.omit(CETS_CBSDs3)
dim(CETS_CBSDs3_na) #1474    3
## dosage file with common clones
CETS_dosageC2_123format_CBSDs3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%CETS_CBSDs3_na$observationUnitName.y,]
CETS_dosageC2_123format_CBSDs3[1:5,1:5]
dim(CETS_dosageC2_123format_CBSDs3) #297 41165

CETS_CBSDs3_sub <- CETS_CBSDs3_na[CETS_CBSDs3_na$observationUnitName.y%in%rownames(CETS_dosageC2_123format_CBSDs3),]
dim(CETS_CBSDs3_sub) #1038    3

CETS_CBSDs3k <-A.mat(CETS_dosageC2_123format_CBSDs3-1)
dim(CETS_CBSDs3k) # 297 297
Z_model_CETS_CBSDs3<-model.matrix(~observationUnitName.y-1, data = CETS_CBSDs3_sub)
dim(Z_model_CETS_CBSDs3) #1038  297
#making block as a factor
CETS_CBSDs3_sub$ENV <- as.factor(CETS_CBSDs3_sub$ENV)

X_model_CETS_CBSDs3 <- model.matrix(~ENV,data = CETS_CBSDs3_sub)

#Narrow sense heritability 
Model_CETS_CBSDs3 <-emmreml(y=CETS_CBSDs3_sub$CBSDs3, X=X_model_CETS_CBSDs3, Z=Z_model_CETS_CBSDs3, K=CETS_CBSDs3k)

Model_CETS_CBSDs3$Vu/(Model_CETS_CBSDs3$Vu+Model_CETS_CBSDs3$Ve) 
# Narrowsense  0.3772576
#Extract blups
CETS_CBSDs3_blups <- Model_CETS_CBSDs3$uhat
write.csv(CETS_CBSDs3_blups, file = "CETS_CBSDs3emmreml_blups.csv")


#CBSDi6
CETS_CBSDs6 <-CETS[,c(45,47,34)]
dim(CETS_CBSDs6 ) #1500    3
CETS_CBSDs6 [1:5,]
# remove NAs
CETS_CBSDs6_na <- na.omit(CETS_CBSDs6)
dim(CETS_CBSDs6_na) #1453    3
## dosage file with common clones
CETS_dosageC2_123format_CBSDs6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%CETS_CBSDs6_na$observationUnitName.y,]
CETS_dosageC2_123format_CBSDs6[1:5,1:5]
dim(CETS_dosageC2_123format_CBSDs6) #297 41165

CETS_CBSDs6_sub <- CETS_CBSDs6_na[CETS_CBSDs6_na$observationUnitName.y%in%rownames(CETS_dosageC2_123format_CBSDs6),]
dim(CETS_CBSDs6_sub) #1030    3

CETS_CBSDs6k <-A.mat(CETS_dosageC2_123format_CBSDs6-1)
dim(CETS_CBSDs6k) #297 297
Z_model_CETS_CBSDs6<-model.matrix(~observationUnitName.y-1, data = CETS_CBSDs6_sub)
dim(Z_model_CETS_CBSDs6) #  1030  297
#making block as a factor
CETS_CBSDs6_sub$ENV <- as.factor(CETS_CBSDs6_sub$ENV)

X_model_CETS_CBSDs6 <- model.matrix(~ENV,data = CETS_CBSDs6_sub)

#Narrow sense heritability 
Model_CETS_CBSDs6 <-emmreml(y=CETS_CBSDs6_sub$CBSDs6, X=X_model_CETS_CBSDs6, Z=Z_model_CETS_CBSDs6, K=CETS_CBSDs6k)

Model_CETS_CBSDs6$Vu/(Model_CETS_CBSDs6$Vu+Model_CETS_CBSDs6$Ve) 
# Narrowsense 0.402012
#Extract blups
CETS_CBSDs6_blups <- Model_CETS_CBSDs6$uhat
write.csv(CETS_CBSDs6_blups, file = "CETS_CBSDs6emmreml_blups.csv")



#CBSDs12
CETS_CBSDs12 <-CETS[,c(45,47,36)]
dim(CETS_CBSDs12 ) #1500    3
CETS_CBSDs12 [1:5,]
# remove NAs
CETS_CBSDs12_na <- na.omit(CETS_CBSDs12)
dim(CETS_CBSDs12_na) # 1099    3
## dosage file with common clones
CETS_dosageC2_123format_CBSDs12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%CETS_CBSDs12_na$observationUnitName.y,]
CETS_dosageC2_123format_CBSDs12[1:5,1:5]
dim(CETS_dosageC2_123format_CBSDs12) #274 41165

CETS_CBSDs12_sub <- CETS_CBSDs12_na[CETS_CBSDs12_na$observationUnitName.y%in%rownames(CETS_dosageC2_123format_CBSDs12),]
dim(CETS_CBSDs12_sub) #828   3

CETS_CBSDs12k <-A.mat(CETS_dosageC2_123format_CBSDs12-1)
dim(CETS_CBSDs12k) #275 275
Z_model_CETS_CBSDs12<-model.matrix(~observationUnitName.y-1, data = CETS_CBSDs12_sub)
dim(Z_model_CETS_CBSDs12) # 828 274

#making block as a factor
CETS_CBSDs12_sub$ENV <- as.factor(CETS_CBSDs12_sub$ENV)

X_model_CETS_CBSDs12 <- model.matrix(~ENV,data = CETS_CBSDs12_sub)

#Narrow sense heritability 
Model_CETS_CBSDs12 <-emmreml(y=CETS_CBSDs12_sub$CBSDs12, X=X_model_CETS_CBSDs12, Z=Z_model_CETS_CBSDs12, K=CETS_CBSDs12k)

Model_CETS_CBSDs12$Vu/(Model_CETS_CBSDs12$Vu+Model_CETS_CBSDs12$Ve) 
# Narrowsense 0.3752487
#Extract blups
CETS_CBSDs12_blups <- Model_CETS_CBSDs12$uhat
write.csv(CETS_CBSDs12_blups, file = "CETS_CBSDs12emmreml_blups.csv")


#zz3
CETS_CMDi3 <-CETS[,c(45,47,39)]
dim(CETS_CMDi3 ) #1500    3
CETS_CMDi3 [1:5,]
# remove NAs
CETS_CMDi3_na <- na.omit(CETS_CMDi3)
dim(CETS_CMDi3_na) # 1476    3
## dosage file with common clones
CETS_dosageC2_123format_CMDi3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%CETS_CMDi3_na$observationUnitName.y,]
CETS_dosageC2_123format_CMDi3[1:5,1:5]
dim(CETS_dosageC2_123format_CMDi3) #274 41165

CETS_CMDi3_sub <- CETS_CMDi3_na[CETS_CMDi3_na$observationUnitName.y%in%rownames(CETS_dosageC2_123format_CMDi3),]
dim(CETS_CMDi3_sub) #1039    3

CETS_CMDi3k <-A.mat(CETS_dosageC2_123format_CMDi3-1)
dim(CETS_CMDi3k) #297 297
Z_model_CETS_CMDi3<-model.matrix(~observationUnitName.y-1, data = CETS_CMDi3_sub)
dim(Z_model_CETS_CMDi3) # 1039  297

#making block as a factor
CETS_CMDi3_sub$ENV <- as.factor(CETS_CMDi3_sub$ENV)

X_model_CETS_CMDi3 <- model.matrix(~ENV,data = CETS_CMDi3_sub)

#Narrow sense heritability 

Model_CETS_CMDi3 <-emmreml(y=CETS_CMDi3_sub$CMDi3, X=X_model_CETS_CMDi3, Z=Z_model_CETS_CMDi3, K=CETS_CMDi3k)

Model_CETS_CMDi3$Vu/(Model_CETS_CMDi3$Vu+Model_CETS_CMDi3$Ve) 
# Narrowsense 0.1718519
#Extract blups
CETS_CMDi3_blups <- Model_CETS_CMDi3$uhat
write.csv(CETS_CMDi3_blups, file = "CETS_CMDi3emmreml_blups.csv")




#CMDi6
CETS_CMDi6 <-CETS[,c(45,47,40)]
dim(CETS_CMDi6 ) #1500    3
CETS_CMDi6 [1:5,]
# remove NAs
CETS_CMDi6_na <- na.omit(CETS_CMDi6)
dim(CETS_CMDi6_na) # 1485    3
## dosage file with common clones
CETS_dosageC2_123format_CMDi6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%CETS_CMDi6_na$observationUnitName.y,]
CETS_dosageC2_123format_CMDi6[1:5,1:5]
dim(CETS_dosageC2_123format_CMDi6) #302 41165

CETS_CMDi6_sub <- CETS_CMDi6_na[CETS_CMDi6_na$observationUnitName.y%in%rownames(CETS_dosageC2_123format_CMDi6),]
dim(CETS_CMDi6_sub) #1037    3

CETS_CMDi6k <-A.mat(CETS_dosageC2_123format_CMDi6-1)
dim(CETS_CMDi6k) #302 302
Z_model_CETS_CMDi6<-model.matrix(~observationUnitName.y-1, data = CETS_CMDi6_sub)
dim(Z_model_CETS_CMDi6) # 1037  302

#making block as a factor
CETS_CMDi6_sub$ENV <- as.factor(CETS_CMDi6_sub$ENV)

X_model_CETS_CMDi6 <- model.matrix(~ENV,data = CETS_CMDi6_sub)

#Narrow sense heritability 

Model_CETS_CMDi6 <-emmreml(y=CETS_CMDi6_sub$CMDi6, X=X_model_CETS_CMDi6, Z=Z_model_CETS_CMDi6, K=CETS_CMDi6k)

Model_CETS_CMDi6$Vu/(Model_CETS_CMDi6$Vu+Model_CETS_CMDi6$Ve) 
# Narrowsense 0.1706189
#Extract blups
CETS_CMDi6_blups <- Model_CETS_CMDi6$uhat
write.csv(CETS_CMDi6_blups, file = "CETS_CMDi6emmreml_blups.csv")


#CMDs3
CETS_CMDs3 <-CETS[,c(45,47,41)]
dim(CETS_CMDs3 ) #1500    3
CETS_CMDs3 [1:5,]
# remove NAs
CETS_CMDs3_na <- na.omit(CETS_CMDs3)
dim(CETS_CMDs3_na) # 1474    3
## dosage file with common clones
CETS_dosageC2_123format_CMDs3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%CETS_CMDs3_na$observationUnitName.y,]
CETS_dosageC2_123format_CMDs3[1:5,1:5]
dim(CETS_dosageC2_123format_CMDs3) #297 41165

CETS_CMDs3_sub <- CETS_CMDs3_na[CETS_CMDs3_na$observationUnitName.y%in%rownames(CETS_dosageC2_123format_CMDs3),]
dim(CETS_CMDs3_sub) #1038    3

CETS_CMDs3k <-A.mat(CETS_dosageC2_123format_CMDs3-1)
dim(CETS_CMDs3k) #297 297
Z_model_CETS_CMDs3<-model.matrix(~observationUnitName.y-1, data = CETS_CMDs3_sub)
dim(Z_model_CETS_CMDs3) # 1038  297

#making block as a factor
CETS_CMDs3_sub$ENV <- as.factor(CETS_CMDs3_sub$ENV)

X_model_CETS_CMDs3 <- model.matrix(~ENV,data = CETS_CMDs3_sub)

#Narrow sense heritability 

Model_CETS_CMDs3 <-emmreml(y=CETS_CMDs3_sub$CMDs3, X=X_model_CETS_CMDs3, Z=Z_model_CETS_CMDs3, K=CETS_CMDs3k)

Model_CETS_CMDs3$Vu/(Model_CETS_CMDs3$Vu+Model_CETS_CMDs3$Ve) 
# Narrowsense 0.7403347
#Extract blups
CETS_CMDs3_blups <- Model_CETS_CMDs3$uhat
write.csv(CETS_CMDs3_blups, file = "CETS_CMDs3emmreml_blups.csv")



#CMDs6
CETS_CMDs6 <-CETS[,c(45,47,42)]
dim(CETS_CMDs6 ) #1500    3
CETS_CMDs6 [1:5,]
# remove NAs
CETS_CMDs6_na <- na.omit(CETS_CMDs6)
dim(CETS_CMDs6_na) # 1133    3
## dosage file with common clones
CETS_dosageC2_123format_CMDs6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%CETS_CMDs6_na$observationUnitName.y,]
CETS_dosageC2_123format_CMDs6[1:5,1:5]
dim(CETS_dosageC2_123format_CMDs6) #291 41165

CETS_CMDs6_sub <- CETS_CMDs6_na[CETS_CMDs6_na$observationUnitName.y%in%rownames(CETS_dosageC2_123format_CMDs6),]
dim(CETS_CMDs6_sub) #819   3

CETS_CMDs6k <-A.mat(CETS_dosageC2_123format_CMDs6-1)
dim(CETS_CMDs6k) #291 291
Z_model_CETS_CMDs6<-model.matrix(~observationUnitName.y-1, data = CETS_CMDs6_sub)
dim(Z_model_CETS_CMDs6) # 819 291

#making block as a factor
CETS_CMDs6_sub$ENV <- as.factor(CETS_CMDs6_sub$ENV)

X_model_CETS_CMDs6 <- model.matrix(~ENV,data = CETS_CMDs6_sub)

#Narrow sense heritability 

Model_CETS_CMDs6 <-emmreml(y=CETS_CMDs6_sub$CMDs6, X=X_model_CETS_CMDs6, Z=Z_model_CETS_CMDs6, K=CETS_CMDs6k)

Model_CETS_CMDs6$Vu/(Model_CETS_CMDs6$Vu+Model_CETS_CMDs6$Ve) 
# Narrowsense 0.6560799
#Extract blups
CETS_CMDs6_blups <- Model_CETS_CMDs6$uhat
write.csv(CETS_CMDs6_blups, file = "CETS_CMDs6emmreml_blups.csv")
















#namulonge 2019
Namulonge_2019 <-CETS %>% filter(locationName == "Namulonge", studyYear == "2019")
dim(Namulonge_2019) #444  40



#cbsdi3
Namulonge_CBSDi3 <-Namulonge_2019[,c(45,47,31)]
dim(Namulonge_CBSDi3 ) #444   3
Namulonge_CBSDi3 [1:5,]
# remove NAs
Namulonge_CBSDi3_na <- na.omit(Namulonge_CBSDi3)
dim(Namulonge_CBSDi3_na) #428   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDi3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDi3_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDi3[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDi3) # 288 41165

Namulonge_CBSDi3_sub <- Namulonge_CBSDi3_na[Namulonge_CBSDi3_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDi3),]
dim(Namulonge_CBSDi3_sub) #288   3

Namulonge_CBSDi3k <-A.mat(Namulonge_dosageC2_123format_CBSDi3-1)
dim(Namulonge_CBSDi3k) #  288 288
Z_model_Namulonge_CBSDi3 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDi3_sub)
dim(Z_model_Namulonge_CBSDi3) # 288 288
#making block as a factor
Namulonge_CBSDi3_sub$ENV <- as.factor(Namulonge_CBSDi3_sub$ENV)

X_model_Namulonge_CBSDi3 <- model.matrix(~ENV,data = Namulonge_CBSDi3_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDi3 <-emmreml(y=Namulonge_CBSDi3_sub$CBSDi3, X=X_model_Namulonge_CBSDi3, Z=Z_model_Namulonge_CBSDi3, K=Namulonge_CBSDi3k)

Model_Namulonge_CBSDi3$Vu/(Model_Namulonge_CBSDi3$Vu+Model_Namulonge_CBSDi3$Ve) 
# Narrowsense 0.2356727
#Extract blups
Namulonge_CBSDi3_blups <- Model_Namulonge_CBSDi3$uhat
write.csv(Namulonge_CBSDi3_blups, file = "Namulonge2019_CBSDi3emmreml_blups.csv")


#CBSDi6
Namulonge_CBSDi6 <-Namulonge_2019[,c(45,47,32)]
dim(Namulonge_CBSDi6 ) #444   3
Namulonge_CBSDi6 [1:5,]
# remove NAs
Namulonge_CBSDi6_na <- na.omit(Namulonge_CBSDi6)
dim(Namulonge_CBSDi6_na) #444   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDi6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDi6_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDi6[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDi6) # 293 41165

Namulonge_CBSDi6_sub <- Namulonge_CBSDi6_na[Namulonge_CBSDi6_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDi6),]
dim(Namulonge_CBSDi6_sub) #293   3

Namulonge_CBSDi6k <-A.mat(Namulonge_dosageC2_123format_CBSDi6-1)
dim(Namulonge_CBSDi6k) #  293 293
Z_model_Namulonge_CBSDi6 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDi6_sub)
dim(Z_model_Namulonge_CBSDi6) # 293 293
#making block as a factor
Namulonge_CBSDi6_sub$ENV <- as.factor(Namulonge_CBSDi6_sub$ENV)

X_model_Namulonge_CBSDi6 <- model.matrix(~ENV,data = Namulonge_CBSDi6_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDi6 <-emmreml(y=Namulonge_CBSDi6_sub$CBSDi6, X=X_model_Namulonge_CBSDi6, Z=Z_model_Namulonge_CBSDi6, K=Namulonge_CBSDi6k)

Model_Namulonge_CBSDi6$Vu/(Model_Namulonge_CBSDi6$Vu+Model_Namulonge_CBSDi6$Ve) 
# Narrowsense 0.2643503
#Extract blups
Namulonge_CBSDi6_blups <- Model_Namulonge_CBSDi6$uhat
write.csv(Namulonge_CBSDi6_blups, file = "Namulonge2019_CBSDi6emmreml_blups.csv")



#CBSDi12
Namulonge_CBSDi12 <-Namulonge_2019[,c(45,47,35)]
dim(Namulonge_CBSDi12 ) #444   3
Namulonge_CBSDi12 [1:5,]
# remove NAs
Namulonge_CBSDi12_na <- na.omit(Namulonge_CBSDi12)
dim(Namulonge_CBSDi12_na) #315   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDi12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDi12_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDi12[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDi12) # 244 41165

Namulonge_CBSDi12_sub <- Namulonge_CBSDi12_na[Namulonge_CBSDi12_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDi12),]
dim(Namulonge_CBSDi12_sub) #244   3

Namulonge_CBSDi12k <-A.mat(Namulonge_dosageC2_123format_CBSDi12-1)
dim(Namulonge_CBSDi12k) #  244 244
Z_model_Namulonge_CBSDi12 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDi12_sub)
dim(Z_model_Namulonge_CBSDi12) # 244 244
#making block as a factor
Namulonge_CBSDi12_sub$ENV <- as.factor(Namulonge_CBSDi12_sub$ENV)

X_model_Namulonge_CBSDi12 <- model.matrix(~ENV,data = Namulonge_CBSDi12_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDi12 <-emmreml(y=Namulonge_CBSDi12_sub$CBSDi12, X=X_model_Namulonge_CBSDi12, Z=Z_model_Namulonge_CBSDi12, K=Namulonge_CBSDi12k)

Model_Namulonge_CBSDi12$Vu/(Model_Namulonge_CBSDi12$Vu+Model_Namulonge_CBSDi12$Ve) 
# Narrowsense 0.4211937
#Extract blups
Namulonge_CBSDi12_blups <- Model_Namulonge_CBSDi12$uhat
write.csv(Namulonge_CBSDi12_blups, file = "Namulonge2019_CBSDi12emmreml_blups.csv")




#CBSDi12
Namulonge_CBSDi12 <-Namulonge_2019[,c(45,47,35)]
dim(Namulonge_CBSDi12 ) #444   3
Namulonge_CBSDi12 [1:5,]
# remove NAs
Namulonge_CBSDi12_na <- na.omit(Namulonge_CBSDi12)
dim(Namulonge_CBSDi12_na) #315   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDi12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDi12_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDi12[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDi12) # 244 41165

Namulonge_CBSDi12_sub <- Namulonge_CBSDi12_na[Namulonge_CBSDi12_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDi12),]
dim(Namulonge_CBSDi12_sub) #244   3

Namulonge_CBSDi12k <-A.mat(Namulonge_dosageC2_123format_CBSDi12-1)
dim(Namulonge_CBSDi12k) #  244 244
Z_model_Namulonge_CBSDi12 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDi12_sub)
dim(Z_model_Namulonge_CBSDi12) # 244 244
#making block as a factor
Namulonge_CBSDi12_sub$ENV <- as.factor(Namulonge_CBSDi12_sub$ENV)

X_model_Namulonge_CBSDi12 <- model.matrix(~ENV,data = Namulonge_CBSDi12_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDi12 <-emmreml(y=Namulonge_CBSDi12_sub$CBSDi12, X=X_model_Namulonge_CBSDi12, Z=Z_model_Namulonge_CBSDi12, K=Namulonge_CBSDi12k)

Model_Namulonge_CBSDi12$Vu/(Model_Namulonge_CBSDi12$Vu+Model_Namulonge_CBSDi12$Ve) 
# Narrowsense 0.4211937
#Extract blups
Namulonge_CBSDi12_blups <- Model_Namulonge_CBSDi12$uhat
write.csv(Namulonge_CBSDi12_blups, file = "Namulonge2019_CBSDi12emmreml_blups.csv")






#CBSDs3
Namulonge_CBSDs3 <-Namulonge_2019[,c(45,47,33)]
dim(Namulonge_CBSDs3 ) #444   3
Namulonge_CBSDs3 [1:5,]
# remove NAs
Namulonge_CBSDs3_na <- na.omit(Namulonge_CBSDs3)
dim(Namulonge_CBSDs3_na) #315   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDs3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDs3_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDs3[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDs3) # 288 41165

Namulonge_CBSDs3_sub <- Namulonge_CBSDs3_na[Namulonge_CBSDs3_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDs3),]
dim(Namulonge_CBSDs3_sub) #288   3

Namulonge_CBSDs3k <-A.mat(Namulonge_dosageC2_123format_CBSDs3-1)
dim(Namulonge_CBSDs3k) #  288 288
Z_model_Namulonge_CBSDs3 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDs3_sub)
dim(Z_model_Namulonge_CBSDs3) # 288 288
#making block as a factor
Namulonge_CBSDs3_sub$ENV <- as.factor(Namulonge_CBSDs3_sub$ENV)

X_model_Namulonge_CBSDs3 <- model.matrix(~ENV,data = Namulonge_CBSDs3_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDs3 <-emmreml(y=Namulonge_CBSDs3_sub$CBSDs3, X=X_model_Namulonge_CBSDs3, Z=Z_model_Namulonge_CBSDs3, K=Namulonge_CBSDs3k)

Model_Namulonge_CBSDs3$Vu/(Model_Namulonge_CBSDs3$Vu+Model_Namulonge_CBSDs3$Ve) 
# Narrowsense 0.2695323
#Extract blups
Namulonge_CBSDs3_blups <- Model_Namulonge_CBSDs3$uhat
write.csv(Namulonge_CBSDs3_blups, file = "Namulonge2019_CBSDs3emmreml_blups.csv")



#CBSDs12
Namulonge_CBSDs12 <-Namulonge_2019[,c(45,47,36)]
dim(Namulonge_CBSDs12 ) #444   3
Namulonge_CBSDs12 [1:5,]
# remove NAs
Namulonge_CBSDs12_na <- na.omit(Namulonge_CBSDs12)
dim(Namulonge_CBSDs12_na) # 315   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDs12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDs12_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDs12[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDs12) # 244 41165

Namulonge_CBSDs12_sub <- Namulonge_CBSDs12_na[Namulonge_CBSDs12_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDs12),]
dim(Namulonge_CBSDs12_sub) #244   3

Namulonge_CBSDs12k <-A.mat(Namulonge_dosageC2_123format_CBSDs12-1)
dim(Namulonge_CBSDs12k) #   244 244
Z_model_Namulonge_CBSDs12 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDs12_sub)
dim(Z_model_Namulonge_CBSDs12) #  244 244
#making block as a factor
Namulonge_CBSDs12_sub$ENV <- as.factor(Namulonge_CBSDs12_sub$ENV)

X_model_Namulonge_CBSDs12 <- model.matrix(~ENV,data = Namulonge_CBSDs12_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDs12 <-emmreml(y=Namulonge_CBSDs12_sub$CBSDs12, X=X_model_Namulonge_CBSDs12, Z=Z_model_Namulonge_CBSDs12, K=Namulonge_CBSDs12k)

Model_Namulonge_CBSDs12$Vu/(Model_Namulonge_CBSDs12$Vu+Model_Namulonge_CBSDs12$Ve) 
# Narrowsense 0.6062639
#Extract blups
Namulonge_CBSDs12_blups <- Model_Namulonge_CBSDs12$uhat
write.csv(Namulonge_CBSDs12_blups, file = "Namulonge2019_CBSDs12emmreml_blups.csv")




#CMDi3
Namulonge_CMDi3 <-Namulonge_2019[,c(45,47,39)]
dim(Namulonge_CMDi3 ) #444   3
Namulonge_CMDi3 [1:5,]
# remove NAs
Namulonge_CMDi3_na <- na.omit(Namulonge_CMDi3)
dim(Namulonge_CMDi3_na) # 428   3
## dosage file with common clones
Namulonge_dosageC2_123format_CMDi3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CMDi3_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CMDi3[1:5,1:5]
dim(Namulonge_dosageC2_123format_CMDi3) # 288 41165

Namulonge_CMDi3_sub <- Namulonge_CMDi3_na[Namulonge_CMDi3_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CMDi3),]
dim(Namulonge_CMDi3_sub) #288   3

Namulonge_CMDi3k <-A.mat(Namulonge_dosageC2_123format_CMDi3-1)
dim(Namulonge_CMDi3k) #   288 288
Z_model_Namulonge_CMDi3 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CMDi3_sub)
dim(Z_model_Namulonge_CMDi3) #  288 288
#making block as a factor
Namulonge_CMDi3_sub$ENV <- as.factor(Namulonge_CMDi3_sub$ENV)

X_model_Namulonge_CMDi3 <- model.matrix(~ENV,data = Namulonge_CMDi3_sub)

#Narrow sense heritability 
Model_Namulonge_CMDi3 <-emmreml(y=Namulonge_CMDi3_sub$CMDi3, X=X_model_Namulonge_CMDi3, Z=Z_model_Namulonge_CMDi3, K=Namulonge_CMDi3k)

Model_Namulonge_CMDi3$Vu/(Model_Namulonge_CMDi3$Vu+Model_Namulonge_CMDi3$Ve) 
# Narrowsense 0.3727123
#Extract blups
Namulonge_CMDi3_blups <- Model_Namulonge_CMDi3$uhat
write.csv(Namulonge_CMDi3_blups, file = "Namulonge2019_CMDi3emmreml_blups.csv")




#CMDi6
Namulonge_CMDi6 <-Namulonge_2019[,c(45,47,40)]
dim(Namulonge_CMDi6 ) #444   3
Namulonge_CMDi6 [1:5,]
# remove NAs
Namulonge_CMDi6_na <- na.omit(Namulonge_CMDi6)
dim(Namulonge_CMDi6_na) # 428   3
## dosage file with common clones
Namulonge_dosageC2_123format_CMDi6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CMDi6_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CMDi6[1:5,1:5]
dim(Namulonge_dosageC2_123format_CMDi6) # 288 41165

Namulonge_CMDi6_sub <- Namulonge_CMDi6_na[Namulonge_CMDi6_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CMDi6),]
dim(Namulonge_CMDi6_sub) #288   3

Namulonge_CMDi6k <-A.mat(Namulonge_dosageC2_123format_CMDi6-1)
dim(Namulonge_CMDi6k) #   288 288
Z_model_Namulonge_CMDi6 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CMDi6_sub)
dim(Z_model_Namulonge_CMDi6) #  288 288
#making block as a factor
Namulonge_CMDi6_sub$ENV <- as.factor(Namulonge_CMDi6_sub$ENV)

X_model_Namulonge_CMDi6 <- model.matrix(~ENV,data = Namulonge_CMDi6_sub)

#Narrow sense heritability 
Model_Namulonge_CMDi6 <-emmreml(y=Namulonge_CMDi6_sub$CMDi6, X=X_model_Namulonge_CMDi6, Z=Z_model_Namulonge_CMDi6, K=Namulonge_CMDi6k)

Model_Namulonge_CMDi6$Vu/(Model_Namulonge_CMDi6$Vu+Model_Namulonge_CMDi6$Ve) 
# Narrowsense  0.3226885
#Extract blups
Namulonge_CMDi6_blups <- Model_Namulonge_CMDi6$uhat
write.csv(Namulonge_CMDi6_blups, file = "Namulonge2019_CMDi6emmreml_blups.csv")



#CMDs3
Namulonge_CMDs3 <-Namulonge_2019[,c(45,47,41)]
dim(Namulonge_CMDs3 ) #444   3
Namulonge_CMDs3 [1:5,]
# remove NAs
Namulonge_CMDs3_na <- na.omit(Namulonge_CMDs3)
dim(Namulonge_CMDs3_na) # 428   3
## dosage file with common clones
Namulonge_dosageC2_123format_CMDs3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CMDs3_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CMDs3[1:5,1:5]
dim(Namulonge_dosageC2_123format_CMDs3) # 288 41165

Namulonge_CMDs3_sub <- Namulonge_CMDs3_na[Namulonge_CMDs3_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CMDs3),]
dim(Namulonge_CMDs3_sub) #288   3

Namulonge_CMDs3k <-A.mat(Namulonge_dosageC2_123format_CMDs3-1)
dim(Namulonge_CMDs3k) #   288 288
Z_model_Namulonge_CMDs3 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CMDs3_sub)
dim(Z_model_Namulonge_CMDs3) #  288 288
#making block as a factor
Namulonge_CMDs3_sub$ENV <- as.factor(Namulonge_CMDs3_sub$ENV)

X_model_Namulonge_CMDs3 <- model.matrix(~ENV,data = Namulonge_CMDs3_sub)

#Narrow sense heritability 
Model_Namulonge_CMDs3 <-emmreml(y=Namulonge_CMDs3_sub$CMDs3, X=X_model_Namulonge_CMDs3, Z=Z_model_Namulonge_CMDs3, K=Namulonge_CMDs3k)

Model_Namulonge_CMDs3$Vu/(Model_Namulonge_CMDs3$Vu+Model_Namulonge_CMDs3$Ve) 
# Narrowsense 0.3497928
#Extract blups
Namulonge_CMDs3_blups <- Model_Namulonge_CMDs3$uhat
write.csv(Namulonge_CMDs3_blups, file = "Namulonge2019_CMDs3emmreml_blups.csv")




#CMDs6
Namulonge_CMDs6 <-Namulonge_2019[,c(45,47,42)]
dim(Namulonge_CMDs6 ) #444   3
Namulonge_CMDs6 [1:5,]
# remove NAs
Namulonge_CMDs6_na <- na.omit(Namulonge_CMDs6)
dim(Namulonge_CMDs6_na) # 428   3
## dosage file with common clones
Namulonge_dosageC2_123format_CMDs6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CMDs6_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CMDs6[1:5,1:5]
dim(Namulonge_dosageC2_123format_CMDs6) # 288 41165

Namulonge_CMDs6_sub <- Namulonge_CMDs6_na[Namulonge_CMDs6_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CMDs6),]
dim(Namulonge_CMDs6_sub) #288   3

Namulonge_CMDs6k <-A.mat(Namulonge_dosageC2_123format_CMDs6-1)
dim(Namulonge_CMDs6k) #   288 288
Z_model_Namulonge_CMDs6 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CMDs6_sub)
dim(Z_model_Namulonge_CMDs6) #  288 288
#making block as a factor
Namulonge_CMDs6_sub$ENV <- as.factor(Namulonge_CMDs6_sub$ENV)

X_model_Namulonge_CMDs6 <- model.matrix(~ENV,data = Namulonge_CMDs6_sub)

#Narrow sense heritability 
Model_Namulonge_CMDs6 <-emmreml(y=Namulonge_CMDs6_sub$CMDs6, X=X_model_Namulonge_CMDs6, Z=Z_model_Namulonge_CMDs6, K=Namulonge_CMDs6k)

Model_Namulonge_CMDs6$Vu/(Model_Namulonge_CMDs6$Vu+Model_Namulonge_CMDs6$Ve) 
# Narrowsense  0.2093126
#Extract blups
Namulonge_CMDs6_blups <- Model_Namulonge_CMDs6$uhat
write.csv(Namulonge_CMDs6_blups, file = "Namulonge2019_CMDs6emmreml_blups.csv")


















#namulonge 2020
Namulonge_2020 <-CETS %>% filter(locationName == "Namulonge", studyYear == "2020")
dim(Namulonge_2020) #365  41



#cbsdi3
Namulonge_CBSDi3 <-Namulonge_2020[,c(45,47,31)]
dim(Namulonge_CBSDi3 ) #365   3
Namulonge_CBSDi3 [1:5,]
# remove NAs
Namulonge_CBSDi3_na <- na.omit(Namulonge_CBSDi3)
dim(Namulonge_CBSDi3_na) #365   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDi3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDi3_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDi3[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDi3) # 267 41165

Namulonge_CBSDi3_sub <- Namulonge_CBSDi3_na[Namulonge_CBSDi3_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDi3),]
dim(Namulonge_CBSDi3_sub) #269   3

Namulonge_CBSDi3k <-A.mat(Namulonge_dosageC2_123format_CBSDi3-1)
dim(Namulonge_CBSDi3k) #  267 267
Z_model_Namulonge_CBSDi3 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDi3_sub)
dim(Z_model_Namulonge_CBSDi3) # 267 267
#making block as a factor
Namulonge_CBSDi3_sub$ENV <- as.factor(Namulonge_CBSDi3_sub$ENV)

X_model_Namulonge_CBSDi3 <- model.matrix(~ENV,data = Namulonge_CBSDi3_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDi3 <-emmreml(y=Namulonge_CBSDi3_sub$CBSDi3, X=X_model_Namulonge_CBSDi3, Z=Z_model_Namulonge_CBSDi3, K=Namulonge_CBSDi3k)

Model_Namulonge_CBSDi3$Vu/(Model_Namulonge_CBSDi3$Vu+Model_Namulonge_CBSDi3$Ve) 
# Narrowsense 0.1818876
#Extract blups
Namulonge_CBSDi3_blups <- Model_Namulonge_CBSDi3$uhat
write.csv(Namulonge_CBSDi3_blups, file = "Namulonge2020_CBSDi3emmreml_blups.csv")


#CBSDi6
Namulonge_CBSDi6 <-Namulonge_2020[,c(45,47,32)]
dim(Namulonge_CBSDi6 ) #444   3
Namulonge_CBSDi6 [1:5,]
# remove NAs
Namulonge_CBSDi6_na <- na.omit(Namulonge_CBSDi6)
dim(Namulonge_CBSDi6_na) #444   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDi6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDi6_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDi6[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDi6) # 293 41165

Namulonge_CBSDi6_sub <- Namulonge_CBSDi6_na[Namulonge_CBSDi6_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDi6),]
dim(Namulonge_CBSDi6_sub) #269   3

Namulonge_CBSDi6k <-A.mat(Namulonge_dosageC2_123format_CBSDi6-1)
dim(Namulonge_CBSDi6k) #  267 267
Z_model_Namulonge_CBSDi6 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDi6_sub)
dim(Z_model_Namulonge_CBSDi6) # 267 267
#making block as a factor
Namulonge_CBSDi6_sub$ENV <- as.factor(Namulonge_CBSDi6_sub$ENV)

X_model_Namulonge_CBSDi6 <- model.matrix(~ENV,data = Namulonge_CBSDi6_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDi6 <-emmreml(y=Namulonge_CBSDi6_sub$CBSDi6, X=X_model_Namulonge_CBSDi6, Z=Z_model_Namulonge_CBSDi6, K=Namulonge_CBSDi6k)

Model_Namulonge_CBSDi6$Vu/(Model_Namulonge_CBSDi6$Vu+Model_Namulonge_CBSDi6$Ve) 
# Narrowsense 0.2294016
#Extract blups
Namulonge_CBSDi6_blups <- Model_Namulonge_CBSDi6$uhat
write.csv(Namulonge_CBSDi6_blups, file = "Namulonge2020_CBSDi6emmreml_blups.csv")



#CBSDi12
Namulonge_CBSDi12 <-Namulonge_2020[,c(45,47,35)]
dim(Namulonge_CBSDi12 ) #444   3
Namulonge_CBSDi12 [1:5,]
# remove NAs
Namulonge_CBSDi12_na <- na.omit(Namulonge_CBSDi12)
dim(Namulonge_CBSDi12_na) #315   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDi12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDi12_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDi12[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDi12) # 244 41165

Namulonge_CBSDi12_sub <- Namulonge_CBSDi12_na[Namulonge_CBSDi12_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDi12),]
dim(Namulonge_CBSDi12_sub) #244   3

Namulonge_CBSDi12k <-A.mat(Namulonge_dosageC2_123format_CBSDi12-1)
dim(Namulonge_CBSDi12k) #  244 244
Z_model_Namulonge_CBSDi12 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDi12_sub)
dim(Z_model_Namulonge_CBSDi12) # 244 244
#making block as a factor
Namulonge_CBSDi12_sub$ENV <- as.factor(Namulonge_CBSDi12_sub$ENV)

X_model_Namulonge_CBSDi12 <- model.matrix(~ENV,data = Namulonge_CBSDi12_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDi12 <-emmreml(y=Namulonge_CBSDi12_sub$CBSDi12, X=X_model_Namulonge_CBSDi12, Z=Z_model_Namulonge_CBSDi12, K=Namulonge_CBSDi12k)

Model_Namulonge_CBSDi12$Vu/(Model_Namulonge_CBSDi12$Vu+Model_Namulonge_CBSDi12$Ve) 
# Narrowsense 0.4211937
#Extract blups
Namulonge_CBSDi12_blups <- Model_Namulonge_CBSDi12$uhat
write.csv(Namulonge_CBSDi12_blups, file = "Namulonge2020_CBSDi12emmreml_blups.csv")




#CBSDi12
Namulonge_CBSDi12 <-Namulonge_2020[,c(45,47,35)]
dim(Namulonge_CBSDi12 ) #444   3
Namulonge_CBSDi12 [1:5,]
# remove NAs
Namulonge_CBSDi12_na <- na.omit(Namulonge_CBSDi12)
dim(Namulonge_CBSDi12_na) #315   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDi12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDi12_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDi12[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDi12) # 244 41165

Namulonge_CBSDi12_sub <- Namulonge_CBSDi12_na[Namulonge_CBSDi12_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDi12),]
dim(Namulonge_CBSDi12_sub) #244   3

Namulonge_CBSDi12k <-A.mat(Namulonge_dosageC2_123format_CBSDi12-1)
dim(Namulonge_CBSDi12k) #  244 244
Z_model_Namulonge_CBSDi12 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDi12_sub)
dim(Z_model_Namulonge_CBSDi12) # 244 244
#making block as a factor
Namulonge_CBSDi12_sub$ENV <- as.factor(Namulonge_CBSDi12_sub$ENV)

X_model_Namulonge_CBSDi12 <- model.matrix(~ENV,data = Namulonge_CBSDi12_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDi12 <-emmreml(y=Namulonge_CBSDi12_sub$CBSDi12, X=X_model_Namulonge_CBSDi12, Z=Z_model_Namulonge_CBSDi12, K=Namulonge_CBSDi12k)

Model_Namulonge_CBSDi12$Vu/(Model_Namulonge_CBSDi12$Vu+Model_Namulonge_CBSDi12$Ve) 
# Narrowsense 0.4211937
#Extract blups
Namulonge_CBSDi12_blups <- Model_Namulonge_CBSDi12$uhat
write.csv(Namulonge_CBSDi12_blups, file = "Namulonge2020_CBSDi12emmreml_blups.csv")






#CBSDs3
Namulonge_CBSDs3 <-Namulonge_2020[,c(45,47,33)]
dim(Namulonge_CBSDs3 ) #444   3
Namulonge_CBSDs3 [1:5,]
# remove NAs
Namulonge_CBSDs3_na <- na.omit(Namulonge_CBSDs3)
dim(Namulonge_CBSDs3_na) #364   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDs3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDs3_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDs3[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDs3) # 266 41165

Namulonge_CBSDs3_sub <- Namulonge_CBSDs3_na[Namulonge_CBSDs3_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDs3),]
dim(Namulonge_CBSDs3_sub) #268   3

Namulonge_CBSDs3k <-A.mat(Namulonge_dosageC2_123format_CBSDs3-1)
dim(Namulonge_CBSDs3k) #  266 266
Z_model_Namulonge_CBSDs3 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDs3_sub)
dim(Z_model_Namulonge_CBSDs3) #  268 266
#making block as a factor
Namulonge_CBSDs3_sub$ENV <- as.factor(Namulonge_CBSDs3_sub$ENV)

X_model_Namulonge_CBSDs3 <- model.matrix(~ENV,data = Namulonge_CBSDs3_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDs3 <-emmreml(y=Namulonge_CBSDs3_sub$CBSDs3, X=X_model_Namulonge_CBSDs3, Z=Z_model_Namulonge_CBSDs3, K=Namulonge_CBSDs3k)

Model_Namulonge_CBSDs3$Vu/(Model_Namulonge_CBSDs3$Vu+Model_Namulonge_CBSDs3$Ve) 
# Narrowsense 0.1641091
#Extract blups
Namulonge_CBSDs3_blups <- Model_Namulonge_CBSDs3$uhat
write.csv(Namulonge_CBSDs3_blups, file = "Namulonge2020_CBSDs3emmreml_blups.csv")


#CBSDs6
Namulonge_CBSDs6 <-Namulonge_2020[,c(45,47,34)]
dim(Namulonge_CBSDs6 ) #365   3
Namulonge_CBSDs6 [1:5,]
# remove NAs
Namulonge_CBSDs6_na <- na.omit(Namulonge_CBSDs6)
dim(Namulonge_CBSDs6_na) #364   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDs6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDs6_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDs6[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDs6) # 266 41165

Namulonge_CBSDs6_sub <- Namulonge_CBSDs6_na[Namulonge_CBSDs6_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDs6),]
dim(Namulonge_CBSDs6_sub) #268   3

Namulonge_CBSDs6k <-A.mat(Namulonge_dosageC2_123format_CBSDs6-1)
dim(Namulonge_CBSDs6k) #  266 266
Z_model_Namulonge_CBSDs6 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDs6_sub)
dim(Z_model_Namulonge_CBSDs6) #  268 266
#making block as a factor
Namulonge_CBSDs6_sub$ENV <- as.factor(Namulonge_CBSDs6_sub$ENV)

X_model_Namulonge_CBSDs6 <- model.matrix(~ENV,data = Namulonge_CBSDs6_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDs6 <-emmreml(y=Namulonge_CBSDs6_sub$CBSDs6, X=X_model_Namulonge_CBSDs6, Z=Z_model_Namulonge_CBSDs6, K=Namulonge_CBSDs6k)

Model_Namulonge_CBSDs6$Vu/(Model_Namulonge_CBSDs6$Vu+Model_Namulonge_CBSDs6$Ve) 
# Narrowsense 0.2689294
#Extract blups
Namulonge_CBSDs6_blups <- Model_Namulonge_CBSDs6$uhat
write.csv(Namulonge_CBSDs6_blups, file = "Namulonge2020_CBSDs6emmreml_blups.csv")







#CBSDs12
Namulonge_CBSDs12 <-Namulonge_2020[,c(45,47,36)]
dim(Namulonge_CBSDs12 ) #365   3
Namulonge_CBSDs12 [1:5,]
# remove NAs
Namulonge_CBSDs12_na <- na.omit(Namulonge_CBSDs12)
dim(Namulonge_CBSDs12_na) # 232   3
## dosage file with common clones
Namulonge_dosageC2_123format_CBSDs12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CBSDs12_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CBSDs12[1:5,1:5]
dim(Namulonge_dosageC2_123format_CBSDs12) # 179 41165

Namulonge_CBSDs12_sub <- Namulonge_CBSDs12_na[Namulonge_CBSDs12_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CBSDs12),]
dim(Namulonge_CBSDs12_sub) #180   3

Namulonge_CBSDs12k <-A.mat(Namulonge_dosageC2_123format_CBSDs12-1)
dim(Namulonge_CBSDs12k) #   179 179
Z_model_Namulonge_CBSDs12 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CBSDs12_sub)
dim(Z_model_Namulonge_CBSDs12) #  180 179
#making block as a factor
Namulonge_CBSDs12_sub$ENV <- as.factor(Namulonge_CBSDs12_sub$ENV)

X_model_Namulonge_CBSDs12 <- model.matrix(~ENV,data = Namulonge_CBSDs12_sub)

#Narrow sense heritability 
Model_Namulonge_CBSDs12 <-emmreml(y=Namulonge_CBSDs12_sub$CBSDs12, X=X_model_Namulonge_CBSDs12, Z=Z_model_Namulonge_CBSDs12, K=Namulonge_CBSDs12k)

Model_Namulonge_CBSDs12$Vu/(Model_Namulonge_CBSDs12$Vu+Model_Namulonge_CBSDs12$Ve) 
# Narrowsense 0.3872835
#Extract blups
Namulonge_CBSDs12_blups <- Model_Namulonge_CBSDs12$uhat
write.csv(Namulonge_CBSDs12_blups, file = "Namulonge2020_CBSDs12emmreml_blups.csv")


#CMDi3
Namulonge_CMDi3 <-Namulonge_2020[,c(45,47,39)]
dim(Namulonge_CMDi3 ) #365   3
Namulonge_CMDi3 [1:5,]
# remove NAs
Namulonge_CMDi3_na <- na.omit(Namulonge_CMDi3)
dim(Namulonge_CMDi3_na) # 232   3
## dosage file with common clones
Namulonge_dosageC2_123format_CMDi3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CMDi3_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CMDi3[1:5,1:5]
dim(Namulonge_dosageC2_123format_CMDi3) # 179 41165

Namulonge_CMDi3_sub <- Namulonge_CMDi3_na[Namulonge_CMDi3_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CMDi3),]
dim(Namulonge_CMDi3_sub) #180   3

Namulonge_CMDi3k <-A.mat(Namulonge_dosageC2_123format_CMDi3-1)
dim(Namulonge_CMDi3k) #   179 179
Z_model_Namulonge_CMDi3 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CMDi3_sub)
dim(Z_model_Namulonge_CMDi3) #  180 179
#making block as a factor
Namulonge_CMDi3_sub$ENV <- as.factor(Namulonge_CMDi3_sub$ENV)

X_model_Namulonge_CMDi3 <- model.matrix(~ENV,data = Namulonge_CMDi3_sub)

#Narrow sense heritability 
Model_Namulonge_CMDi3 <-emmreml(y=Namulonge_CMDi3_sub$CMDi3, X=X_model_Namulonge_CMDi3, Z=Z_model_Namulonge_CMDi3, K=Namulonge_CMDi3k)

Model_Namulonge_CMDi3$Vu/(Model_Namulonge_CMDi3$Vu+Model_Namulonge_CMDi3$Ve) 
# Narrowsense 0.1908764
#Extract blups
Namulonge_CMDi3_blups <- Model_Namulonge_CMDi3$uhat
write.csv(Namulonge_CMDi3_blups, file = "Namulonge2020_CMDi3emmreml_blups.csv")


#CMDi6
Namulonge_CMDi6 <-Namulonge_2020[,c(45,47,40)]
dim(Namulonge_CMDi6 ) #365   3
Namulonge_CMDi6 [1:5,]
# remove NAs
Namulonge_CMDi6_na <- na.omit(Namulonge_CMDi6)
dim(Namulonge_CMDi6_na) # 232   3
## dosage file with common clones
Namulonge_dosageC2_123format_CMDi6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CMDi6_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CMDi6[1:5,1:5]
dim(Namulonge_dosageC2_123format_CMDi6) # 179 41165

Namulonge_CMDi6_sub <- Namulonge_CMDi6_na[Namulonge_CMDi6_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CMDi6),]
dim(Namulonge_CMDi6_sub) #180   3

Namulonge_CMDi6k <-A.mat(Namulonge_dosageC2_123format_CMDi6-1)
dim(Namulonge_CMDi6k) #   179 179
Z_model_Namulonge_CMDi6 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CMDi6_sub)
dim(Z_model_Namulonge_CMDi6) #  180 179
#making block as a factor
Namulonge_CMDi6_sub$ENV <- as.factor(Namulonge_CMDi6_sub$ENV)

X_model_Namulonge_CMDi6 <- model.matrix(~ENV,data = Namulonge_CMDi6_sub)

#Narrow sense heritability 
Model_Namulonge_CMDi6 <-emmreml(y=Namulonge_CMDi6_sub$CMDi6, X=X_model_Namulonge_CMDi6, Z=Z_model_Namulonge_CMDi6, K=Namulonge_CMDi6k)

Model_Namulonge_CMDi6$Vu/(Model_Namulonge_CMDi6$Vu+Model_Namulonge_CMDi6$Ve) 
# Narrowsense 0.1597638
#Extract blups
Namulonge_CMDi6_blups <- Model_Namulonge_CMDi6$uhat
write.csv(Namulonge_CMDi6_blups, file = "Namulonge2020_CMDi6emmreml_blups.csv")




#CMDs3
Namulonge_CMDs3 <-Namulonge_2020[,c(45,47,41)]
dim(Namulonge_CMDs3 ) #365   3
Namulonge_CMDs3 [1:5,]
# remove NAs
Namulonge_CMDs3_na <- na.omit(Namulonge_CMDs3)
dim(Namulonge_CMDs3_na) # 232   3
## dosage file with common clones
Namulonge_dosageC2_123format_CMDs3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CMDs3_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CMDs3[1:5,1:5]
dim(Namulonge_dosageC2_123format_CMDs3) # 179 41165

Namulonge_CMDs3_sub <- Namulonge_CMDs3_na[Namulonge_CMDs3_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CMDs3),]
dim(Namulonge_CMDs3_sub) #180   3

Namulonge_CMDs3k <-A.mat(Namulonge_dosageC2_123format_CMDs3-1)
dim(Namulonge_CMDs3k) #   179 179
Z_model_Namulonge_CMDs3 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CMDs3_sub)
dim(Z_model_Namulonge_CMDs3) #  180 179
#making block as a factor
Namulonge_CMDs3_sub$ENV <- as.factor(Namulonge_CMDs3_sub$ENV)

X_model_Namulonge_CMDs3 <- model.matrix(~ENV,data = Namulonge_CMDs3_sub)

#Narrow sense heritability 
Model_Namulonge_CMDs3 <-emmreml(y=Namulonge_CMDs3_sub$CMDs3, X=X_model_Namulonge_CMDs3, Z=Z_model_Namulonge_CMDs3, K=Namulonge_CMDs3k)

Model_Namulonge_CMDs3$Vu/(Model_Namulonge_CMDs3$Vu+Model_Namulonge_CMDs3$Ve) 
# Narrowsense 0.1261503
#Extract blups
Namulonge_CMDs3_blups <- Model_Namulonge_CMDs3$uhat
write.csv(Namulonge_CMDs3_blups, file = "Namulonge2020_CMDs3emmreml_blups.csv")


#CMDs6
Namulonge_CMDs6 <-Namulonge_2020[,c(45,47,42)]
dim(Namulonge_CMDs6 ) #365   3
Namulonge_CMDs6 [1:5,]
# remove NAs
Namulonge_CMDs6_na <- na.omit(Namulonge_CMDs6)
dim(Namulonge_CMDs6_na) # 232   3
## dosage file with common clones
Namulonge_dosageC2_123format_CMDs6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Namulonge_CMDs6_na$observationUnitName.y,]
Namulonge_dosageC2_123format_CMDs6[1:5,1:5]
dim(Namulonge_dosageC2_123format_CMDs6) # 179 41165

Namulonge_CMDs6_sub <- Namulonge_CMDs6_na[Namulonge_CMDs6_na$observationUnitName.y%in%rownames(Namulonge_dosageC2_123format_CMDs6),]
dim(Namulonge_CMDs6_sub) #180   3

Namulonge_CMDs6k <-A.mat(Namulonge_dosageC2_123format_CMDs6-1)
dim(Namulonge_CMDs6k) #   179 179
Z_model_Namulonge_CMDs6 <-model.matrix(~observationUnitName.y-1, data = Namulonge_CMDs6_sub)
dim(Z_model_Namulonge_CMDs6) #  180 179
#making block as a factor
Namulonge_CMDs6_sub$ENV <- as.factor(Namulonge_CMDs6_sub$ENV)

X_model_Namulonge_CMDs6 <- model.matrix(~ENV,data = Namulonge_CMDs6_sub)

#Narrow sense heritability 
Model_Namulonge_CMDs6 <-emmreml(y=Namulonge_CMDs6_sub$CMDs6, X=X_model_Namulonge_CMDs6, Z=Z_model_Namulonge_CMDs6, K=Namulonge_CMDs6k)

Model_Namulonge_CMDs6$Vu/(Model_Namulonge_CMDs6$Vu+Model_Namulonge_CMDs6$Ve) 
# Narrowsense 0.1168979
#Extract blups
Namulonge_CMDs6_blups <- Model_Namulonge_CMDs6$uhat
write.csv(Namulonge_CMDs6_blups, file = "Namulonge2020_CMDs6emmreml_blups.csv")









# serere ------------------------------------------------------------------

#Serere 2019
Serere_2019 <-CETS %>% filter(locationName == "Serere", studyYear == "2019")
dim(Serere_2019) #444  40



#cbsdi3
Serere_CBSDi3 <-Serere_2019[,c(45,47,31)]
dim(Serere_CBSDi3 ) #444   3
Serere_CBSDi3 [1:5,]
# remove NAs
Serere_CBSDi3_na <- na.omit(Serere_CBSDi3)
dim(Serere_CBSDi3_na) #428   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDi3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDi3_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDi3[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDi3) # 288 41165

Serere_CBSDi3_sub <- Serere_CBSDi3_na[Serere_CBSDi3_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDi3),]
dim(Serere_CBSDi3_sub) #288   3

Serere_CBSDi3k <-A.mat(Serere_dosageC2_123format_CBSDi3-1)
dim(Serere_CBSDi3k) #  288 288
Z_model_Serere_CBSDi3 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDi3_sub)
dim(Z_model_Serere_CBSDi3) # 288 288
#making block as a factor
Serere_CBSDi3_sub$ENV <- as.factor(Serere_CBSDi3_sub$ENV)

X_model_Serere_CBSDi3 <- model.matrix(~ENV,data = Serere_CBSDi3_sub)

#Narrow sense heritability 
Model_Serere_CBSDi3 <-emmreml(y=Serere_CBSDi3_sub$CBSDi3, X=X_model_Serere_CBSDi3, Z=Z_model_Serere_CBSDi3, K=Serere_CBSDi3k)

Model_Serere_CBSDi3$Vu/(Model_Serere_CBSDi3$Vu+Model_Serere_CBSDi3$Ve) 
# Narrowsense 0.08957445
#Extract blups
Serere_CBSDi3_blups <- Model_Serere_CBSDi3$uhat
write.csv(Serere_CBSDi3_blups, file = "Serere2019_CBSDi3emmreml_blups.csv")


#CBSDi6
Serere_CBSDi6 <-Serere_2019[,c(45,47,32)]
dim(Serere_CBSDi6 ) #444   3
Serere_CBSDi6 [1:5,]
# remove NAs
Serere_CBSDi6_na <- na.omit(Serere_CBSDi6)
dim(Serere_CBSDi6_na) #444   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDi6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDi6_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDi6[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDi6) # 293 41165

Serere_CBSDi6_sub <- Serere_CBSDi6_na[Serere_CBSDi6_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDi6),]
dim(Serere_CBSDi6_sub) #293   3

Serere_CBSDi6k <-A.mat(Serere_dosageC2_123format_CBSDi6-1)
dim(Serere_CBSDi6k) #  293 293
Z_model_Serere_CBSDi6 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDi6_sub)
dim(Z_model_Serere_CBSDi6) # 293 293
#making block as a factor
Serere_CBSDi6_sub$ENV <- as.factor(Serere_CBSDi6_sub$ENV)

X_model_Serere_CBSDi6 <- model.matrix(~ENV,data = Serere_CBSDi6_sub)

#Narrow sense heritability 
Model_Serere_CBSDi6 <-emmreml(y=Serere_CBSDi6_sub$CBSDi6, X=X_model_Serere_CBSDi6, Z=Z_model_Serere_CBSDi6, K=Serere_CBSDi6k)

Model_Serere_CBSDi6$Vu/(Model_Serere_CBSDi6$Vu+Model_Serere_CBSDi6$Ve) 
# Narrowsense 0.1649301
#Extract blups
Serere_CBSDi6_blups <- Model_Serere_CBSDi6$uhat
write.csv(Serere_CBSDi6_blups, file = "Serere2019_CBSDi6emmreml_blups.csv")



#CBSDi12
Serere_CBSDi12 <-Serere_2019[,c(45,47,35)]
dim(Serere_CBSDi12 ) #444   3
Serere_CBSDi12 [1:5,]
# remove NAs
Serere_CBSDi12_na <- na.omit(Serere_CBSDi12)
dim(Serere_CBSDi12_na) #315   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDi12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDi12_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDi12[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDi12) # 244 41165

Serere_CBSDi12_sub <- Serere_CBSDi12_na[Serere_CBSDi12_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDi12),]
dim(Serere_CBSDi12_sub) #244   3

Serere_CBSDi12k <-A.mat(Serere_dosageC2_123format_CBSDi12-1)
dim(Serere_CBSDi12k) #  244 244
Z_model_Serere_CBSDi12 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDi12_sub)
dim(Z_model_Serere_CBSDi12) # 244 244
#making block as a factor
Serere_CBSDi12_sub$ENV <- as.factor(Serere_CBSDi12_sub$ENV)

X_model_Serere_CBSDi12 <- model.matrix(~ENV,data = Serere_CBSDi12_sub)

#Narrow sense heritability 
Model_Serere_CBSDi12 <-emmreml(y=Serere_CBSDi12_sub$CBSDi12, X=X_model_Serere_CBSDi12, Z=Z_model_Serere_CBSDi12, K=Serere_CBSDi12k)

Model_Serere_CBSDi12$Vu/(Model_Serere_CBSDi12$Vu+Model_Serere_CBSDi12$Ve) 
# Narrowsense 0.03022201
#Extract blups
Serere_CBSDi12_blups <- Model_Serere_CBSDi12$uhat
write.csv(Serere_CBSDi12_blups, file = "Serere2019_CBSDi12emmreml_blups.csv")




#CBSDi12
Serere_CBSDi12 <-Serere_2019[,c(45,47,35)]
dim(Serere_CBSDi12 ) #444   3
Serere_CBSDi12 [1:5,]
# remove NAs
Serere_CBSDi12_na <- na.omit(Serere_CBSDi12)
dim(Serere_CBSDi12_na) #315   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDi12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDi12_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDi12[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDi12) # 244 41165

Serere_CBSDi12_sub <- Serere_CBSDi12_na[Serere_CBSDi12_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDi12),]
dim(Serere_CBSDi12_sub) #244   3

Serere_CBSDi12k <-A.mat(Serere_dosageC2_123format_CBSDi12-1)
dim(Serere_CBSDi12k) #  244 244
Z_model_Serere_CBSDi12 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDi12_sub)
dim(Z_model_Serere_CBSDi12) # 244 244
#making block as a factor
Serere_CBSDi12_sub$ENV <- as.factor(Serere_CBSDi12_sub$ENV)

X_model_Serere_CBSDi12 <- model.matrix(~ENV,data = Serere_CBSDi12_sub)

#Narrow sense heritability 
Model_Serere_CBSDi12 <-emmreml(y=Serere_CBSDi12_sub$CBSDi12, X=X_model_Serere_CBSDi12, Z=Z_model_Serere_CBSDi12, K=Serere_CBSDi12k)

Model_Serere_CBSDi12$Vu/(Model_Serere_CBSDi12$Vu+Model_Serere_CBSDi12$Ve) 
# Narrowsense 0.03003499
#Extract blups
Serere_CBSDi12_blups <- Model_Serere_CBSDi12$uhat
write.csv(Serere_CBSDi12_blups, file = "Serere2019_CBSDi12emmreml_blups.csv")






#CBSDs3
Serere_CBSDs3 <-Serere_2019[,c(45,47,33)]
dim(Serere_CBSDs3 ) #444   3
Serere_CBSDs3 [1:5,]
# remove NAs
Serere_CBSDs3_na <- na.omit(Serere_CBSDs3)
dim(Serere_CBSDs3_na) #315   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDs3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDs3_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDs3[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDs3) # 288 41165

Serere_CBSDs3_sub <- Serere_CBSDs3_na[Serere_CBSDs3_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDs3),]
dim(Serere_CBSDs3_sub) #288   3

Serere_CBSDs3k <-A.mat(Serere_dosageC2_123format_CBSDs3-1)
dim(Serere_CBSDs3k) #  288 288
Z_model_Serere_CBSDs3 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDs3_sub)
dim(Z_model_Serere_CBSDs3) # 288 288
#making block as a factor
Serere_CBSDs3_sub$ENV <- as.factor(Serere_CBSDs3_sub$ENV)

X_model_Serere_CBSDs3 <- model.matrix(~ENV,data = Serere_CBSDs3_sub)

#Narrow sense heritability 
Model_Serere_CBSDs3 <-emmreml(y=Serere_CBSDs3_sub$CBSDs3, X=X_model_Serere_CBSDs3, Z=Z_model_Serere_CBSDs3, K=Serere_CBSDs3k)

Model_Serere_CBSDs3$Vu/(Model_Serere_CBSDs3$Vu+Model_Serere_CBSDs3$Ve) 
# Narrowsense 0.08072132
#Extract blups
Serere_CBSDs3_blups <- Model_Serere_CBSDs3$uhat
write.csv(Serere_CBSDs3_blups, file = "Serere2019_CBSDs3emmreml_blups.csv")



#CBSDs6
Serere_CBSDs6 <-Serere_2019[,c(45,47,34)]
dim(Serere_CBSDs6 ) #444   3
Serere_CBSDs6 [1:5,]
# remove NAs
Serere_CBSDs6_na <- na.omit(Serere_CBSDs6)
dim(Serere_CBSDs6_na) #315   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDs6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDs6_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDs6[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDs6) # 288 41165

Serere_CBSDs6_sub <- Serere_CBSDs6_na[Serere_CBSDs6_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDs6),]
dim(Serere_CBSDs6_sub) #288   3

Serere_CBSDs6k <-A.mat(Serere_dosageC2_123format_CBSDs6-1)
dim(Serere_CBSDs6k) #  288 288
Z_model_Serere_CBSDs6 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDs6_sub)
dim(Z_model_Serere_CBSDs6) # 288 288
#making block as a factor
Serere_CBSDs6_sub$ENV <- as.factor(Serere_CBSDs6_sub$ENV)

X_model_Serere_CBSDs6 <- model.matrix(~ENV,data = Serere_CBSDs6_sub)

#Narrow sense heritability 
Model_Serere_CBSDs6 <-emmreml(y=Serere_CBSDs6_sub$CBSDs6, X=X_model_Serere_CBSDs6, Z=Z_model_Serere_CBSDs6, K=Serere_CBSDs6k)

Model_Serere_CBSDs6$Vu/(Model_Serere_CBSDs6$Vu+Model_Serere_CBSDs6$Ve) 
# Narrowsense  0.168187
#Extract blups
Serere_CBSDs6_blups <- Model_Serere_CBSDs6$uhat
write.csv(Serere_CBSDs6_blups, file = "Serere2019_CBSDs6emmreml_blups.csv")









#CBSDs12
Serere_CBSDs12 <-Serere_2019[,c(45,47,36)]
dim(Serere_CBSDs12 ) #444   3
Serere_CBSDs12 [1:5,]
# remove NAs
Serere_CBSDs12_na <- na.omit(Serere_CBSDs12)
dim(Serere_CBSDs12_na) # 315   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDs12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDs12_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDs12[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDs12) # 244 41165

Serere_CBSDs12_sub <- Serere_CBSDs12_na[Serere_CBSDs12_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDs12),]
dim(Serere_CBSDs12_sub) #244   3

Serere_CBSDs12k <-A.mat(Serere_dosageC2_123format_CBSDs12-1)
dim(Serere_CBSDs12k) #   244 244
Z_model_Serere_CBSDs12 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDs12_sub)
dim(Z_model_Serere_CBSDs12) #  244 244
#making block as a factor
Serere_CBSDs12_sub$ENV <- as.factor(Serere_CBSDs12_sub$ENV)

X_model_Serere_CBSDs12 <- model.matrix(~ENV,data = Serere_CBSDs12_sub)

#Narrow sense heritability 
Model_Serere_CBSDs12 <-emmreml(y=Serere_CBSDs12_sub$CBSDs12, X=X_model_Serere_CBSDs12, Z=Z_model_Serere_CBSDs12, K=Serere_CBSDs12k)

Model_Serere_CBSDs12$Vu/(Model_Serere_CBSDs12$Vu+Model_Serere_CBSDs12$Ve) 
# Narrowsense 0.07819979
#Extract blups
Serere_CBSDs12_blups <- Model_Serere_CBSDs12$uhat
write.csv(Serere_CBSDs12_blups, file = "Serere2019_CBSDs12emmreml_blups.csv")


#CMDi3
Serere_CMDi3 <-Serere_2019[,c(45,47,39)]
dim(Serere_CMDi3 ) #444   3
Serere_CMDi3 [1:5,]
# remove NAs
Serere_CMDi3_na <- na.omit(Serere_CMDi3)
dim(Serere_CMDi3_na) # 315   3
## dosage file with common clones
Serere_dosageC2_123format_CMDi3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CMDi3_na$observationUnitName.y,]
Serere_dosageC2_123format_CMDi3[1:5,1:5]
dim(Serere_dosageC2_123format_CMDi3) # 244 41165

Serere_CMDi3_sub <- Serere_CMDi3_na[Serere_CMDi3_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CMDi3),]
dim(Serere_CMDi3_sub) #244   3

Serere_CMDi3k <-A.mat(Serere_dosageC2_123format_CMDi3-1)
dim(Serere_CMDi3k) #   244 244
Z_model_Serere_CMDi3 <-model.matrix(~observationUnitName.y-1, data = Serere_CMDi3_sub)
dim(Z_model_Serere_CMDi3) #  244 244
#making block as a factor
Serere_CMDi3_sub$ENV <- as.factor(Serere_CMDi3_sub$ENV)

X_model_Serere_CMDi3 <- model.matrix(~ENV,data = Serere_CMDi3_sub)

#Narrow sense heritability 
Model_Serere_CMDi3 <-emmreml(y=Serere_CMDi3_sub$CMDi3, X=X_model_Serere_CMDi3, Z=Z_model_Serere_CMDi3, K=Serere_CMDi3k)

Model_Serere_CMDi3$Vu/(Model_Serere_CMDi3$Vu+Model_Serere_CMDi3$Ve) 
# Narrowsense 0.208784
#Extract blups
Serere_CMDi3_blups <- Model_Serere_CMDi3$uhat
write.csv(Serere_CMDi3_blups, file = "Serere2019_CMDi3emmreml_blups.csv")


#CMDi6
Serere_CMDi6 <-Serere_2019[,c(45,47,40)]
dim(Serere_CMDi6 ) #444   3
Serere_CMDi6 [1:5,]
# remove NAs
Serere_CMDi6_na <- na.omit(Serere_CMDi6)
dim(Serere_CMDi6_na) # 315   3
## dosage file with common clones
Serere_dosageC2_123format_CMDi6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CMDi6_na$observationUnitName.y,]
Serere_dosageC2_123format_CMDi6[1:5,1:5]
dim(Serere_dosageC2_123format_CMDi6) # 244 41165

Serere_CMDi6_sub <- Serere_CMDi6_na[Serere_CMDi6_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CMDi6),]
dim(Serere_CMDi6_sub) #244   3

Serere_CMDi6k <-A.mat(Serere_dosageC2_123format_CMDi6-1)
dim(Serere_CMDi6k) #   244 244
Z_model_Serere_CMDi6 <-model.matrix(~observationUnitName.y-1, data = Serere_CMDi6_sub)
dim(Z_model_Serere_CMDi6) #  244 244
#making block as a factor
Serere_CMDi6_sub$ENV <- as.factor(Serere_CMDi6_sub$ENV)

X_model_Serere_CMDi6 <- model.matrix(~ENV,data = Serere_CMDi6_sub)

#Narrow sense heritability 
Model_Serere_CMDi6 <-emmreml(y=Serere_CMDi6_sub$CMDi6, X=X_model_Serere_CMDi6, Z=Z_model_Serere_CMDi6, K=Serere_CMDi6k)

Model_Serere_CMDi6$Vu/(Model_Serere_CMDi6$Vu+Model_Serere_CMDi6$Ve) 
# Narrowsense  0.1803754
#Extract blups
Serere_CMDi6_blups <- Model_Serere_CMDi6$uhat
write.csv(Serere_CMDi6_blups, file = "Serere2019_CMDi6emmreml_blups.csv")



#CMDs3
Serere_CMDs3 <-Serere_2019[,c(45,47,41)]
dim(Serere_CMDs3 ) #444   3
Serere_CMDs3 [1:5,]
# remove NAs
Serere_CMDs3_na <- na.omit(Serere_CMDs3)
dim(Serere_CMDs3_na) # 315   3
## dosage file with common clones
Serere_dosageC2_123format_CMDs3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CMDs3_na$observationUnitName.y,]
Serere_dosageC2_123format_CMDs3[1:5,1:5]
dim(Serere_dosageC2_123format_CMDs3) # 244 41165

Serere_CMDs3_sub <- Serere_CMDs3_na[Serere_CMDs3_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CMDs3),]
dim(Serere_CMDs3_sub) #244   3

Serere_CMDs3k <-A.mat(Serere_dosageC2_123format_CMDs3-1)
dim(Serere_CMDs3k) #   244 244
Z_model_Serere_CMDs3 <-model.matrix(~observationUnitName.y-1, data = Serere_CMDs3_sub)
dim(Z_model_Serere_CMDs3) #  244 244
#making block as a factor
Serere_CMDs3_sub$ENV <- as.factor(Serere_CMDs3_sub$ENV)

X_model_Serere_CMDs3 <- model.matrix(~ENV,data = Serere_CMDs3_sub)

#Narrow sense heritability 
Model_Serere_CMDs3 <-emmreml(y=Serere_CMDs3_sub$CMDs3, X=X_model_Serere_CMDs3, Z=Z_model_Serere_CMDs3, K=Serere_CMDs3k)

Model_Serere_CMDs3$Vu/(Model_Serere_CMDs3$Vu+Model_Serere_CMDs3$Ve) 
# Narrowsense 0.07753807
#Extract blups
Serere_CMDs3_blups <- Model_Serere_CMDs3$uhat
write.csv(Serere_CMDs3_blups, file = "Serere2019_CMDs3emmreml_blups.csv")


#CMDs6
Serere_CMDs6 <-Serere_2019[,c(45,47,42)]
dim(Serere_CMDs6 ) #444   3
Serere_CMDs6 [1:5,]
# remove NAs
Serere_CMDs6_na <- na.omit(Serere_CMDs6)
dim(Serere_CMDs6_na) # 315   3
## dosage file with common clones
Serere_dosageC2_123format_CMDs6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CMDs6_na$observationUnitName.y,]
Serere_dosageC2_123format_CMDs6[1:5,1:5]
dim(Serere_dosageC2_123format_CMDs6) # 244 41165

Serere_CMDs6_sub <- Serere_CMDs6_na[Serere_CMDs6_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CMDs6),]
dim(Serere_CMDs6_sub) #244   3

Serere_CMDs6k <-A.mat(Serere_dosageC2_123format_CMDs6-1)
dim(Serere_CMDs6k) #   244 244
Z_model_Serere_CMDs6 <-model.matrix(~observationUnitName.y-1, data = Serere_CMDs6_sub)
dim(Z_model_Serere_CMDs6) #  244 244
#making block as a factor
Serere_CMDs6_sub$ENV <- as.factor(Serere_CMDs6_sub$ENV)

X_model_Serere_CMDs6 <- model.matrix(~ENV,data = Serere_CMDs6_sub)

#Narrow sense heritability 
Model_Serere_CMDs6 <-emmreml(y=Serere_CMDs6_sub$CMDs6, X=X_model_Serere_CMDs6, Z=Z_model_Serere_CMDs6, K=Serere_CMDs6k)

Model_Serere_CMDs6$Vu/(Model_Serere_CMDs6$Vu+Model_Serere_CMDs6$Ve) 
# Narrowsense  0.1803754
#Extract blups
Serere_CMDs6_blups <- Model_Serere_CMDs6$uhat
write.csv(Serere_CMDs6_blups, file = "Serere2019_CMDs6emmreml_blups.csv")


















#Serere 2020
Serere_2020 <-CETS %>% filter(locationName == "Serere", studyYear == "2020")
dim(Serere_2020) #365  41



#cbsdi3
Serere_CBSDi3 <-Serere_2020[,c(45,47,31)]
dim(Serere_CBSDi3 ) #365   3
Serere_CBSDi3 [1:5,]
# remove NAs
Serere_CBSDi3_na <- na.omit(Serere_CBSDi3)
dim(Serere_CBSDi3_na) #365   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDi3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDi3_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDi3[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDi3) # 267 41165

Serere_CBSDi3_sub <- Serere_CBSDi3_na[Serere_CBSDi3_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDi3),]
dim(Serere_CBSDi3_sub) #269   3

Serere_CBSDi3k <-A.mat(Serere_dosageC2_123format_CBSDi3-1)
dim(Serere_CBSDi3k) #  267 267
Z_model_Serere_CBSDi3 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDi3_sub)
dim(Z_model_Serere_CBSDi3) # 267 267
#making block as a factor
Serere_CBSDi3_sub$ENV <- as.factor(Serere_CBSDi3_sub$ENV)

X_model_Serere_CBSDi3 <- model.matrix(~ENV,data = Serere_CBSDi3_sub)

#Narrow sense heritability 
Model_Serere_CBSDi3 <-emmreml(y=Serere_CBSDi3_sub$CBSDi3, X=X_model_Serere_CBSDi3, Z=Z_model_Serere_CBSDi3, K=Serere_CBSDi3k)

Model_Serere_CBSDi3$Vu/(Model_Serere_CBSDi3$Vu+Model_Serere_CBSDi3$Ve) 
# Narrowsense 0.4441812
#Extract blups
Serere_CBSDi3_blups <- Model_Serere_CBSDi3$uhat
write.csv(Serere_CBSDi3_blups, file = "Serere2020_CBSDi3emmreml_blups.csv")


#CBSDi6
Serere_CBSDi6 <-Serere_2020[,c(45,47,32)]
dim(Serere_CBSDi6 ) #444   3
Serere_CBSDi6 [1:5,]
# remove NAs
Serere_CBSDi6_na <- na.omit(Serere_CBSDi6)
dim(Serere_CBSDi6_na) #444   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDi6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDi6_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDi6[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDi6) # 293 41165

Serere_CBSDi6_sub <- Serere_CBSDi6_na[Serere_CBSDi6_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDi6),]
dim(Serere_CBSDi6_sub) #269   3

Serere_CBSDi6k <-A.mat(Serere_dosageC2_123format_CBSDi6-1)
dim(Serere_CBSDi6k) #  267 267
Z_model_Serere_CBSDi6 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDi6_sub)
dim(Z_model_Serere_CBSDi6) # 267 267
#making block as a factor
Serere_CBSDi6_sub$ENV <- as.factor(Serere_CBSDi6_sub$ENV)

X_model_Serere_CBSDi6 <- model.matrix(~ENV,data = Serere_CBSDi6_sub)

#Narrow sense heritability 
Model_Serere_CBSDi6 <-emmreml(y=Serere_CBSDi6_sub$CBSDi6, X=X_model_Serere_CBSDi6, Z=Z_model_Serere_CBSDi6, K=Serere_CBSDi6k)

Model_Serere_CBSDi6$Vu/(Model_Serere_CBSDi6$Vu+Model_Serere_CBSDi6$Ve) 
# Narrowsense 0.2758978
#Extract blups
Serere_CBSDi6_blups <- Model_Serere_CBSDi6$uhat
write.csv(Serere_CBSDi6_blups, file = "Serere2020_CBSDi6emmreml_blups.csv")



#CBSDi12
Serere_CBSDi12 <-Serere_2020[,c(45,47,35)]
dim(Serere_CBSDi12 ) #444   3
Serere_CBSDi12 [1:5,]
# remove NAs
Serere_CBSDi12_na <- na.omit(Serere_CBSDi12)
dim(Serere_CBSDi12_na) #315   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDi12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDi12_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDi12[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDi12) # 244 41165

Serere_CBSDi12_sub <- Serere_CBSDi12_na[Serere_CBSDi12_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDi12),]
dim(Serere_CBSDi12_sub) #244   3

Serere_CBSDi12k <-A.mat(Serere_dosageC2_123format_CBSDi12-1)
dim(Serere_CBSDi12k) #  244 244
Z_model_Serere_CBSDi12 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDi12_sub)
dim(Z_model_Serere_CBSDi12) # 244 244
#making block as a factor
Serere_CBSDi12_sub$ENV <- as.factor(Serere_CBSDi12_sub$ENV)

X_model_Serere_CBSDi12 <- model.matrix(~ENV,data = Serere_CBSDi12_sub)

#Narrow sense heritability 
Model_Serere_CBSDi12 <-emmreml(y=Serere_CBSDi12_sub$CBSDi12, X=X_model_Serere_CBSDi12, Z=Z_model_Serere_CBSDi12, K=Serere_CBSDi12k)

Model_Serere_CBSDi12$Vu/(Model_Serere_CBSDi12$Vu+Model_Serere_CBSDi12$Ve) 
# Narrowsense 0.1325479
#Extract blups
Serere_CBSDi12_blups <- Model_Serere_CBSDi12$uhat
write.csv(Serere_CBSDi12_blups, file = "Serere2020_CBSDi12emmreml_blups.csv")




#CBSDi12
Serere_CBSDi12 <-Serere_2020[,c(45,47,35)]
dim(Serere_CBSDi12 ) #444   3
Serere_CBSDi12 [1:5,]
# remove NAs
Serere_CBSDi12_na <- na.omit(Serere_CBSDi12)
dim(Serere_CBSDi12_na) #315   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDi12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDi12_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDi12[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDi12) # 244 41165

Serere_CBSDi12_sub <- Serere_CBSDi12_na[Serere_CBSDi12_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDi12),]
dim(Serere_CBSDi12_sub) #244   3

Serere_CBSDi12k <-A.mat(Serere_dosageC2_123format_CBSDi12-1)
dim(Serere_CBSDi12k) #  244 244
Z_model_Serere_CBSDi12 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDi12_sub)
dim(Z_model_Serere_CBSDi12) # 244 244
#making block as a factor
Serere_CBSDi12_sub$ENV <- as.factor(Serere_CBSDi12_sub$ENV)

X_model_Serere_CBSDi12 <- model.matrix(~ENV,data = Serere_CBSDi12_sub)

#Narrow sense heritability 
Model_Serere_CBSDi12 <-emmreml(y=Serere_CBSDi12_sub$CBSDi12, X=X_model_Serere_CBSDi12, Z=Z_model_Serere_CBSDi12, K=Serere_CBSDi12k)

Model_Serere_CBSDi12$Vu/(Model_Serere_CBSDi12$Vu+Model_Serere_CBSDi12$Ve) 
# Narrowsense 0.4211937
#Extract blups
Serere_CBSDi12_blups <- Model_Serere_CBSDi12$uhat
write.csv(Serere_CBSDi12_blups, file = "Serere2020_CBSDi12emmreml_blups.csv")






#CBSDs3
Serere_CBSDs3 <-Serere_2020[,c(45,47,33)]
dim(Serere_CBSDs3 ) #444   3
Serere_CBSDs3 [1:5,]
# remove NAs
Serere_CBSDs3_na <- na.omit(Serere_CBSDs3)
dim(Serere_CBSDs3_na) #364   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDs3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDs3_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDs3[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDs3) # 266 41165

Serere_CBSDs3_sub <- Serere_CBSDs3_na[Serere_CBSDs3_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDs3),]
dim(Serere_CBSDs3_sub) #268   3

Serere_CBSDs3k <-A.mat(Serere_dosageC2_123format_CBSDs3-1)
dim(Serere_CBSDs3k) #  266 266
Z_model_Serere_CBSDs3 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDs3_sub)
dim(Z_model_Serere_CBSDs3) #  268 266
#making block as a factor
Serere_CBSDs3_sub$ENV <- as.factor(Serere_CBSDs3_sub$ENV)

X_model_Serere_CBSDs3 <- model.matrix(~ENV,data = Serere_CBSDs3_sub)

#Narrow sense heritability 
Model_Serere_CBSDs3 <-emmreml(y=Serere_CBSDs3_sub$CBSDs3, X=X_model_Serere_CBSDs3, Z=Z_model_Serere_CBSDs3, K=Serere_CBSDs3k)

Model_Serere_CBSDs3$Vu/(Model_Serere_CBSDs3$Vu+Model_Serere_CBSDs3$Ve) 
# Narrowsense 0.2785055
#Extract blups
Serere_CBSDs3_blups <- Model_Serere_CBSDs3$uhat
write.csv(Serere_CBSDs3_blups, file = "Serere2020_CBSDs3emmreml_blups.csv")


#CBSDs6
Serere_CBSDs6 <-Serere_2020[,c(45,47,34)]
dim(Serere_CBSDs6 ) #365   3
Serere_CBSDs6 [1:5,]
# remove NAs
Serere_CBSDs6_na <- na.omit(Serere_CBSDs6)
dim(Serere_CBSDs6_na) #364   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDs6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDs6_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDs6[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDs6) # 266 41165

Serere_CBSDs6_sub <- Serere_CBSDs6_na[Serere_CBSDs6_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDs6),]
dim(Serere_CBSDs6_sub) #268   3

Serere_CBSDs6k <-A.mat(Serere_dosageC2_123format_CBSDs6-1)
dim(Serere_CBSDs6k) #  266 266
Z_model_Serere_CBSDs6 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDs6_sub)
dim(Z_model_Serere_CBSDs6) #  268 266
#making block as a factor
Serere_CBSDs6_sub$ENV <- as.factor(Serere_CBSDs6_sub$ENV)

X_model_Serere_CBSDs6 <- model.matrix(~ENV,data = Serere_CBSDs6_sub)

#Narrow sense heritability 
Model_Serere_CBSDs6 <-emmreml(y=Serere_CBSDs6_sub$CBSDs6, X=X_model_Serere_CBSDs6, Z=Z_model_Serere_CBSDs6, K=Serere_CBSDs6k)

Model_Serere_CBSDs6$Vu/(Model_Serere_CBSDs6$Vu+Model_Serere_CBSDs6$Ve) 
# Narrowsense 0.1962771
#Extract blups
Serere_CBSDs6_blups <- Model_Serere_CBSDs6$uhat
write.csv(Serere_CBSDs6_blups, file = "Serere2020_CBSDs6emmreml_blups.csv")







#CBSDs12
Serere_CBSDs12 <-Serere_2020[,c(45,47,36)]
dim(Serere_CBSDs12 ) #365   3
Serere_CBSDs12 [1:5,]
# remove NAs
Serere_CBSDs12_na <- na.omit(Serere_CBSDs12)
dim(Serere_CBSDs12_na) # 232   3
## dosage file with common clones
Serere_dosageC2_123format_CBSDs12 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CBSDs12_na$observationUnitName.y,]
Serere_dosageC2_123format_CBSDs12[1:5,1:5]
dim(Serere_dosageC2_123format_CBSDs12) # 179 41165

Serere_CBSDs12_sub <- Serere_CBSDs12_na[Serere_CBSDs12_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CBSDs12),]
dim(Serere_CBSDs12_sub) #180   3

Serere_CBSDs12k <-A.mat(Serere_dosageC2_123format_CBSDs12-1)
dim(Serere_CBSDs12k) #   179 179
Z_model_Serere_CBSDs12 <-model.matrix(~observationUnitName.y-1, data = Serere_CBSDs12_sub)
dim(Z_model_Serere_CBSDs12) #  180 179
#making block as a factor
Serere_CBSDs12_sub$ENV <- as.factor(Serere_CBSDs12_sub$ENV)

X_model_Serere_CBSDs12 <- model.matrix(~ENV,data = Serere_CBSDs12_sub)

#Narrow sense heritability 
Model_Serere_CBSDs12 <-emmreml(y=Serere_CBSDs12_sub$CBSDs12, X=X_model_Serere_CBSDs12, Z=Z_model_Serere_CBSDs12, K=Serere_CBSDs12k)

Model_Serere_CBSDs12$Vu/(Model_Serere_CBSDs12$Vu+Model_Serere_CBSDs12$Ve) 
# Narrowsense 0.2704852
#Extract blups
Serere_CBSDs12_blups <- Model_Serere_CBSDs12$uhat
write.csv(Serere_CBSDs12_blups, file = "Serere2020_CBSDs12emmreml_blups.csv")



#CMDi3
Serere_CMDi3 <-Serere_2020[,c(45,47,39)]
dim(Serere_CMDi3 ) #365   3
Serere_CMDi3 [1:5,]
# remove NAs
Serere_CMDi3_na <- na.omit(Serere_CMDi3)
dim(Serere_CMDi3_na) # 232   3
## dosage file with common clones
Serere_dosageC2_123format_CMDi3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CMDi3_na$observationUnitName.y,]
Serere_dosageC2_123format_CMDi3[1:5,1:5]
dim(Serere_dosageC2_123format_CMDi3) # 179 41165

Serere_CMDi3_sub <- Serere_CMDi3_na[Serere_CMDi3_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CMDi3),]
dim(Serere_CMDi3_sub) #180   3

Serere_CMDi3k <-A.mat(Serere_dosageC2_123format_CMDi3-1)
dim(Serere_CMDi3k) #   179 179
Z_model_Serere_CMDi3 <-model.matrix(~observationUnitName.y-1, data = Serere_CMDi3_sub)
dim(Z_model_Serere_CMDi3) #  180 179
#making block as a factor
Serere_CMDi3_sub$ENV <- as.factor(Serere_CMDi3_sub$ENV)

X_model_Serere_CMDi3 <- model.matrix(~ENV,data = Serere_CMDi3_sub)

#Narrow sense heritability 
Model_Serere_CMDi3 <-emmreml(y=Serere_CMDi3_sub$CMDi3, X=X_model_Serere_CMDi3, Z=Z_model_Serere_CMDi3, K=Serere_CMDi3k)

Model_Serere_CMDi3$Vu/(Model_Serere_CMDi3$Vu+Model_Serere_CMDi3$Ve) 
# Narrowsense 0.2452161
#Extract blups
Serere_CMDi3_blups <- Model_Serere_CMDi3$uhat
write.csv(Serere_CMDi3_blups, file = "Serere2020_CMDi3emmreml_blups.csv")


#CMDi6
Serere_CMDi6 <-Serere_2020[,c(45,47,40)]
dim(Serere_CMDi6 ) #365   3
Serere_CMDi6 [1:5,]
# remove NAs
Serere_CMDi6_na <- na.omit(Serere_CMDi6)
dim(Serere_CMDi6_na) # 232   3
## dosage file with common clones
Serere_dosageC2_123format_CMDi6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CMDi6_na$observationUnitName.y,]
Serere_dosageC2_123format_CMDi6[1:5,1:5]
dim(Serere_dosageC2_123format_CMDi6) # 179 41165

Serere_CMDi6_sub <- Serere_CMDi6_na[Serere_CMDi6_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CMDi6),]
dim(Serere_CMDi6_sub) #180   3

Serere_CMDi6k <-A.mat(Serere_dosageC2_123format_CMDi6-1)
dim(Serere_CMDi6k) #   179 179
Z_model_Serere_CMDi6 <-model.matrix(~observationUnitName.y-1, data = Serere_CMDi6_sub)
dim(Z_model_Serere_CMDi6) #  180 179
#making block as a factor
Serere_CMDi6_sub$ENV <- as.factor(Serere_CMDi6_sub$ENV)

X_model_Serere_CMDi6 <- model.matrix(~ENV,data = Serere_CMDi6_sub)

#Narrow sense heritability 
Model_Serere_CMDi6 <-emmreml(y=Serere_CMDi6_sub$CMDi6, X=X_model_Serere_CMDi6, Z=Z_model_Serere_CMDi6, K=Serere_CMDi6k)

Model_Serere_CMDi6$Vu/(Model_Serere_CMDi6$Vu+Model_Serere_CMDi6$Ve) 
# Narrowsense 0.04910497
#Extract blups
Serere_CMDi6_blups <- Model_Serere_CMDi6$uhat
write.csv(Serere_CMDi6_blups, file = "Serere2020_CMDi6emmreml_blups.csv")


#CMDs3
Serere_CMDs3 <-Serere_2020[,c(45,47,41)]
dim(Serere_CMDs3 ) #365   3
Serere_CMDs3 [1:5,]
# remove NAs
Serere_CMDs3_na <- na.omit(Serere_CMDs3)
dim(Serere_CMDs3_na) # 232   3
## dosage file with common clones
Serere_dosageC2_123format_CMDs3 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CMDs3_na$observationUnitName.y,]
Serere_dosageC2_123format_CMDs3[1:5,1:5]
dim(Serere_dosageC2_123format_CMDs3) # 179 41165

Serere_CMDs3_sub <- Serere_CMDs3_na[Serere_CMDs3_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CMDs3),]
dim(Serere_CMDs3_sub) #180   3

Serere_CMDs3k <-A.mat(Serere_dosageC2_123format_CMDs3-1)
dim(Serere_CMDs3k) #   179 179
Z_model_Serere_CMDs3 <-model.matrix(~observationUnitName.y-1, data = Serere_CMDs3_sub)
dim(Z_model_Serere_CMDs3) #  180 179
#making block as a factor
Serere_CMDs3_sub$ENV <- as.factor(Serere_CMDs3_sub$ENV)

X_model_Serere_CMDs3 <- model.matrix(~ENV,data = Serere_CMDs3_sub)

#Narrow sense heritability 
Model_Serere_CMDs3 <-emmreml(y=Serere_CMDs3_sub$CMDs3, X=X_model_Serere_CMDs3, Z=Z_model_Serere_CMDs3, K=Serere_CMDs3k)

Model_Serere_CMDs3$Vu/(Model_Serere_CMDs3$Vu+Model_Serere_CMDs3$Ve) 
# Narrowsense 0.1472728
#Extract blups
Serere_CMDs3_blups <- Model_Serere_CMDs3$uhat
write.csv(Serere_CMDs3_blups, file = "Serere2020_CMDs3emmreml_blups.csv")


#CMDs6
Serere_CMDs6 <-Serere_2020[,c(45,47,42)]
dim(Serere_CMDs6 ) #365   3
Serere_CMDs6 [1:5,]
# remove NAs
Serere_CMDs6_na <- na.omit(Serere_CMDs6)
dim(Serere_CMDs6_na) # 232   3
## dosage file with common clones
Serere_dosageC2_123format_CMDs6 <- C2_Geno_1_subset[rownames(C2_Geno_1_subset)%in%Serere_CMDs6_na$observationUnitName.y,]
Serere_dosageC2_123format_CMDs6[1:5,1:5]
dim(Serere_dosageC2_123format_CMDs6) # 179 41165

Serere_CMDs6_sub <- Serere_CMDs6_na[Serere_CMDs6_na$observationUnitName.y%in%rownames(Serere_dosageC2_123format_CMDs6),]
dim(Serere_CMDs6_sub) #180   3

Serere_CMDs6k <-A.mat(Serere_dosageC2_123format_CMDs6-1)
dim(Serere_CMDs6k) #   179 179
Z_model_Serere_CMDs6 <-model.matrix(~observationUnitName.y-1, data = Serere_CMDs6_sub)
dim(Z_model_Serere_CMDs6) #  180 179
#making block as a factor
Serere_CMDs6_sub$ENV <- as.factor(Serere_CMDs6_sub$ENV)

X_model_Serere_CMDs6 <- model.matrix(~ENV,data = Serere_CMDs6_sub)

#Narrow sense heritability 
Model_Serere_CMDs6 <-emmreml(y=Serere_CMDs6_sub$CMDs6, X=X_model_Serere_CMDs6, Z=Z_model_Serere_CMDs6, K=Serere_CMDs6k)

Model_Serere_CMDs6$Vu/(Model_Serere_CMDs6$Vu+Model_Serere_CMDs6$Ve) 
# Narrowsense 0.02693887
#Extract blups
Serere_CMDs6_blups <- Model_Serere_CMDs6$uhat
write.csv(Serere_CMDs6_blups, file = "Serere2020_CMDs6emmreml_blups.csv")

####end 



# one_step ----------------------------------------------------------------


common_code <- read.csv("~/Desktop/phd_Defense/2021/2022_readings/Publications/ELISA/Data_sets/Set_1/common_code.csv")
dim(common_code) #6570    2

names(common_code)[names(common_code) == "accession_name"] <- "germplasmName"
CETS<-merge(CETS,common_code, by="germplasmName")
dim(CETS) #1500   45


commonInd <- intersect(rownames(Genotypes_filtered),CETS$observationUnitName.y) ####commonindividuals 
length(commonInd)  #302 


C2_Geno_1_subset <- Genotypes_filtered[(rownames(Genotypes_filtered)%in%commonInd),] 
dim(C2_Geno_1_subset) # 303 41165



Phenotypes_subset <- CETS[CETS$observationUnitName.y %in% commonInd,]
dim(Phenotypes_subset) #  302 41165

write.csv(Phenotypes_subset, file = "C2_CET_PHENOTYPES.csv")




# Additive relationship matrix --------------------------------------------
relation_matrix <- A.mat(C2_Geno_1_subset,impute.method="EM",min.MAF = 0.05)

pca=prcomp(relation_matrix, center = TRUE, scale. = TRUE, rank. = 10)#conducting pca analysis 
loadings=pca$rotation

eigen1 <- (pca$sdev)^2 

# Variances in percentage
variance1 <- eigen1*100/sum(eigen1) 

plot(variance1 , main = 'Eigen values ') #### most variation explained by the first 7 PCs 


# cumulative variance
cumvar1 <- cumsum(variance1)

eigenvals <- data.frame(eigen = eigen1, variance = variance1, cumvariance = cumvar1)

loadings_=as.data.frame(loadings[,1:10])

write.csv(loadings_,'PCA_LOADINGS1.csv')

data=read.csv('PCA_LOADINGS1.csv')



#PCAS 
summary(pca)$importance[,1:10] 

pc_scores<-pca$x
pc_scores[1:5,1:5]
A=pc_scores %>% 
  as.data.frame %>% 
  ggplot(.,aes(x=PC1,y=PC2)) + 
  geom_point() + 
  theme_classic() +
  labs(x="PC1 (27.6%)",y="PC2 (18%)")

B=pc_scores %>% 
  as.data.frame %>% 
  ggplot(.,aes(x=PC3,y=PC4)) + 
  geom_point() + 
  theme_classic() +
  labs(x="PC3 (10.2%)",y="PC4 (7.2%)")




#Marker data object 
SNP <- colnames(C2_Geno_1_subset)
data<- data.frame(SNP)

head(data)



#seperating the data points \

data2<-data %>%separate(SNP, c("chrom", "pos"), "_")

data2$chrom <- as.numeric(sub("S", "",data2$chrom, perl=TRUE))

data2$pos <- as.integer(data2$pos)

data2$chrom <- as.character(data2$chrom)

head(data2)
###creaating the number of rows that will be combined to mark an object of the marker data for GWA analysis 
n <- nrow(data2)
n


##creating the marker objective that will contain the number of rows for the data object with the chromosome and position plus the actual marker ranking

#was getting an error because of different lengths Error in data.frame(..., check.names = FALSE) : arguments imply differing number of rows: 24040, 1596
C2_Geno_1_subset1 = t(C2_Geno_1_subset)

#marker object
markerObj <- cbind(1:n,data2,C2_Geno_1_subset1 )

head(markerObj)

dim(markerObj) #41165   306

class(markerObj) #"data.frame"





# GWAS --------------------------------------------------------------------

#object for the trait of evaluation
#REMOVING THE BLOCK since its accounted for under location
Phenotypes_subset$ENV <-paste(Phenotypes_subset$locationName,Phenotypes_subset$studyYear)


phenotype1 <- Phenotypes_subset[,c(45,46,31)]

phenotype1.1 <- na.omit(phenotype1)

phenotype1.1 = data.frame(phenotype1.1)

dim(phenotype1.1) # 1039    5





GWASMODEL_CBSDi3_4_1 <- GWAS(pheno = phenotype1.1, geno = markerObj , P3D=TRUE, K=relation_matrix, fixed = "ENV" , n.core=10, min.MAF = 0.05, plot=TRUE, n.PC = 4)
#relation_matrix
#SAVING THE MODEL

#4 Pcs
save(GWASMODEL_CBSDi3_4_1, file = "GWASMODEL_CBSDi3_4_1.RData")

#CBSDI6
phenotype1 <- Phenotypes_subset[,c(45,46,32)]

phenotype1.1 <- na.omit(phenotype1)

phenotype1.1 = data.frame(phenotype1.1)

dim(phenotype1.1) # 1039    5
#Gwas model which is a mixed model 


GWASMODEL_CBSDi6_4_1 <- GWAS(pheno = phenotype1.1, geno = markerObj , P3D=TRUE, K=relation_matrix, fixed = "ENV" , n.core=10, min.MAF = 0.05, plot=TRUE, n.PC = 4)
#relation_matrix
#SAVING THE MODEL
#4 Pcs
save(GWASMODEL_CBSDi6_4_1, file = "GWASMODEL_CBSDi6_4_1.RData")


#CBSDi12
phenotype1 <- Phenotypes_subset[,c(45,46,35)]

phenotype1.1 <- na.omit(phenotype1)

phenotype1.1 = data.frame(phenotype1.1)

dim(phenotype1.1) # 1039    5
#Gwas model which is a mixed model 



GWASMODEL_CBSDi12_4_1 <- GWAS(pheno = phenotype1.1, geno = markerObj , P3D=TRUE, K=relation_matrix, fixed = "ENV" , n.core=10, min.MAF = 0.05, plot=TRUE, n.PC = 4)
#relation_matrix
#SAVING THE MODEL
#4 Pcs
save(GWASMODEL_CBSDi12_4_1, file = "GWASMODEL_CBSDi12_4_1.RData")



#####################################SEVERITIES 
phenotype1 <- Phenotypes_subset[,c(45,46,33)]

phenotype1.1 <- na.omit(phenotype1)

phenotype1.1 = data.frame(phenotype1.1)

dim(phenotype1.1) # 1039    5





GWASMODEL_CBSDs3_4_1 <- GWAS(pheno = phenotype1.1, geno = markerObj , P3D=TRUE, K=relation_matrix, fixed = "ENV" , n.core=10, min.MAF = 0.05, plot=TRUE, n.PC = 4)
#relation_matrix
#SAVING THE MODEL
#4 Pcs
save(GWASMODEL_CBSDs3_4_1, file = "GWASMODEL_CBSDs3_4_1.RData")


#CBSDs6
phenotype1 <- Phenotypes_subset[,c(45,46,34)]

phenotype1.1 <- na.omit(phenotype1)

phenotype1.1 = data.frame(phenotype1.1)

dim(phenotype1.1) # 1039    5
#Gwas model which is a mixed model 



GWASMODEL_CBSDs6_4_1 <- GWAS(pheno = phenotype1.1, geno = markerObj , P3D=TRUE, K=relation_matrix, fixed = "ENV" , n.core=10, min.MAF = 0.05, plot=TRUE, n.PC = 4)
#relation_matrix
#SAVING THE MODEL
#4 Pcs
save(GWASMODEL_CBSDs6_4_1, file = "GWASMODEL_CBSDs6_4_1.RData")


#CBSDs12
phenotype1 <- Phenotypes_subset[,c(45,46,36)]

phenotype1.1 <- na.omit(phenotype1)


GWASMODEL_CBSDs12_4_1 <- GWAS(pheno = phenotype1.1, geno = markerObj , P3D=TRUE, K=relation_matrix, fixed = "ENV" , n.core=10, min.MAF = 0.05, plot=TRUE, n.PC = 4)
#relation_matrix
#SAVING THE MODEL
#4 Pcs
save(GWASMODEL_CBSDs12_4_1, file = "GWASMODEL_CBSDs12_4_1.RData")






#CMDi 
#
phenotype1 <- Phenotypes_subset[,c(45,46,39)]

phenotype1.1 <- na.omit(phenotype1)

phenotype1.1 = data.frame(phenotype1.1)

dim(phenotype1.1) # 1039    5
#Gwas model which is a mixed model 


GWASMODEL_CMDi3_4_1 <- GWAS(pheno = phenotype1.1, geno = markerObj , P3D=TRUE, K=relation_matrix, fixed = "ENV" , n.core=10, min.MAF = 0.05, plot=TRUE, n.PC = 4)
#relation_matrix
#SAVING THE MODEL
#4 Pcs
save(GWASMODEL_CMDi3_4_1, file = "GWASMODEL_CMDi3_4_1.RData")

#CMDi6
phenotype1 <- Phenotypes_subset[,c(45,46,40)]

phenotype1.1 <- na.omit(phenotype1)

phenotype1.1 = data.frame(phenotype1.1)

dim(phenotype1.1) # 1039    5


GWASMODEL_CMDi6_4_1 <- GWAS(pheno = phenotype1.1, geno = markerObj , P3D=TRUE, K=relation_matrix, fixed = "ENV" , n.core=10, min.MAF = 0.05, plot=TRUE, n.PC = 4)
#relation_matrix
#SAVING THE MODEL
#4 Pcs
save(GWASMODEL_CMDi6_4_1, file = "GWASMODEL_CMDi6_4_1.RData")






#CMDs 
#
phenotype1 <- Phenotypes_subset[,c(45,46,41)]

phenotype1.1 <- na.omit(phenotype1)

phenotype1.1 = data.frame(phenotype1.1)

dim(phenotype1.1) # 1039    5
#Gwas model which is a mixed model 


GWASMODEL_CMDs3_4_1 <- GWAS(pheno = phenotype1.1, geno = markerObj , P3D=TRUE, K=relation_matrix, fixed = "ENV" , n.core=10, min.MAF = 0.05, plot=TRUE, n.PC = 4)
#relation_matrix
#SAVING THE MODEL
#4 Pcs
save(GWASMODEL_CMDs3_4_1, file = "GWASMODEL_CMDs3_4_1.RData")

#CMDs6
phenotype1 <- Phenotypes_subset[,c(45,46,42)]

phenotype1.1 <- na.omit(phenotype1)

phenotype1.1 = data.frame(phenotype1.1)

dim(phenotype1.1) # 1039    5
#Gwas model which is a mixed model 


GWASMODEL_CMDs6_4_1 <- GWAS(pheno = phenotype1.1, geno = markerObj , P3D=TRUE, K=relation_matrix, fixed = "ENV" , n.core=10, min.MAF = 0.05, plot=TRUE, n.PC = 4)
#relation_matrix
#SAVING THE MODEL
#4 Pcs
save(GWASMODEL_CMDs6_4_1, file = "GWASMODEL_CMDs6_4_1.RData")



#TWO STEP GWAS

# Two step GWAS - Univariate and Multivariate in GEMMA --------------------

#Used default settings in Gemma 
