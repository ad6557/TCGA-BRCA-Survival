library(survival)
library(coxme)
library(tibble)
load("signature_building.RData")


train.cox <- train[,c(select,"overall_survival","vital_status",
                      "age","stage")]
train.f <- tibble::rownames_to_column(train.cox, "ID")
cox.f <- coxme(Surv(overall_survival, vital_status) ~ 
                 SORBS1+AC104211.1+AC007686.3+AC016583.1+WDR17+
                 AC007014.2+AC016582.3+`SPACA6P-AS`+QPRT+ACSM4+
                 age+stage+(1|ID),data=train.f)

anova(cox.f, mdrcox) # 0.9106

AIC(cox.f) # 1052
AIC(mdrcox) # 1052
