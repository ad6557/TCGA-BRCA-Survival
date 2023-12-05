library(skimr)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)

load("preprocessing+genefinding.RData")

######### prepare dataset for cox regression #########
# RNAseq only
d_mat = t(limma_res$voomObj$E) # 1224*27289
coxdata = as.data.frame(d_mat)

# firsth merge more clinical
# https://docs.gdc.cancer.gov/Data_Dictionary/viewer/

coxdata$race=as.factor(limma_res$voomObj$targets$race)
coxdata$age=as.numeric(limma_res$voomObj$targets$age_at_index)

coxdata$stage=as.factor(limma_res$voomObj$targets$paper_pathologic_stage)
coxdata$stage[which(is.na(coxdata$stage))] <- "NA"
# levels(coxdata$stage)[levels(coxdata$stage)=="NA"] <- NA
# too many NAs in this variable -- can be problematic if just toss them

# then merge survival
coxdata$vital_status=limma_res$voomObj$targets$vital_status
coxdata$days_to_last_follow_up=limma_res$voomObj$targets$days_to_last_follow_up
coxdata$days_to_death=limma_res$voomObj$targets$days_to_death

coxdata = coxdata[limma_res$voomObj$targets$definition == "Primary solid Tumor",] # 1111
coxdata$overall_survival <- ifelse(coxdata$vital_status == "Alive",
                                   coxdata$days_to_last_follow_up,
                                   coxdata$days_to_death)
table(coxdata$vital_status) # 956+154=1110
coxdata$vital_status =ifelse(coxdata$vital_status=="Dead",1,0)
table(coxdata$vital_status)

######### addressing abnormal values #########
naos = which(is.na(coxdata$overall_survival)) # 401 758
negativeost = which(coxdata$overall_survival<=0) # 3
coxdata = coxdata[-c(naos,negativeost),] # 1088

# map from gene_id to gene_name 
colnames(coxdata)[1:27289] == limma_res$voomObj$genes$gene_id
colnames(coxdata)[1:27289] <- limma_res$voomObj$genes$gene_name

######### split train test #########

set.seed(2023)
id = sample(1:1088, 800)
train = coxdata[id,]
test = coxdata[-id,]

######### Multivariate Cox #########

# race not significant, remove it!
train.cox <- train[,c(select.final,"overall_survival","vital_status",
                      "age","stage")]
test.cox <- test[,c(select.final,"overall_survival","vital_status",
                    "age","stage")]

mdrcox <- coxph(Surv(overall_survival, vital_status) ~ .,data=train.cox)
summary(mdrcox)
cox.zph(mdrcox) # stage must be taken out and the model must be stratified by stage. 

mdrcox <- coxph(Surv(overall_survival, vital_status) ~ strata(stage)+., data=train.cox)
summary(mdrcox)
cox.zph(mdrcox)

######### calculate risk score #########

beta = coef(mdrcox)[1:16]
train.score = as.matrix(train.cox[,1:16]) %*% beta
test.score = as.matrix(test.cox[,1:16]) %*% beta

train.cox$risk = as.vector(train.score)
train.cox$group = ifelse(train.cox$risk>median(train.score),"high","low")
test.cox$risk = as.vector(test.score)
test.cox$group = ifelse(test.cox$risk>median(test.score),"high","low")


save(limma_res, train.cox, test.cox,
     beta, select.final, file="signature_building.RData")
