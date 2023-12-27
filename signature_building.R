library(skimr)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)

load("preprocessing+genefinding.RData")

######### Multivariate Cox #########

# round 1
train.cox <- train[,c(select.final,"overall_survival","vital_status",
                      "age","stage")]
mdrcox <- coxph(Surv(overall_survival, vital_status) ~ .,data=train.cox)
summary(mdrcox)

coef.names = names(summary(mdrcox)$coefficients[,1])
select.driver = coef.names[summary(mdrcox)$coefficients[ ,5] < 0.05]
select.driver = setdiff(select.driver, c("age","stageStage_I","stageStage_II","stageStage_III","stageStage_IV"))
select = replace(select.driver, select.driver=="`SPACA6P-AS`", "SPACA6P-AS")
select

# round 2
train.cox <- train[,c(select,"overall_survival","vital_status",
                      "age","stage")]
mdrcox <- coxph(Surv(overall_survival, vital_status) ~ .,data=train.cox)
summary(mdrcox)
cox.zph(mdrcox)

######### calculate risk score #########
test.cox <- test[,c(select,"overall_survival","vital_status",
                    "age","stage")]

beta = coef(mdrcox)[1:length(select)]
train.score = as.matrix(train.cox[,1:length(select)]) %*% beta
test.score = as.matrix(test.cox[,1:length(select)]) %*% beta

train.cox$risk = as.vector(train.score)
train.cox$group = ifelse(train.cox$risk>median(train.score),"high","low")
test.cox$risk = as.vector(test.score)
test.cox$group = ifelse(test.cox$risk>median(test.score),"high","low")


save(limma_res, select.final, train.cox, test.cox,
     beta, select, file="signature_building.RData")
