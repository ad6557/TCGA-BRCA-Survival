library(survival)
library(ggplot2)
library(survminer)
library(survivalROC)

load("signature_building.RData")

######### KM plot #########

# fitting survival curve
train.KM <- survfit(Surv(overall_survival, vital_status) ~ group, data = train.cox)
ggsurvplot(train.KM,
           data = train.cox,
           pval = T,
           risk.table = T)

test.KM <- survfit(Surv(overall_survival, vital_status) ~ group, data = test.cox)
ggsurvplot(test.KM,
           data = test.cox,
           pval = T,
           risk.table = T)

######### time-dependent ROC and AUC #########

train.ROC <-survivalROC(Stime=train.cox$overall_survival,
                        status=train.cox$vital_status,
                        marker=train.cox$risk,
                        predict.time = 365*10,
                        method = 'KM')
train.AUC = round(train.ROC[["AUC"]],3)

## Plot
train_ROC_data <- data.frame(TP=train.ROC$TP,FP=train.ROC$FP)
train_ROC_data %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank()) + 
  geom_label(label=train.AUC, x=0.5,y=0.5,
             label.padding = unit(0.55, "lines"))


test.ROC <-survivalROC(Stime=test.cox$overall_survival,
                       status=test.cox$vital_status,
                       marker=test.cox$risk,
                       predict.time = 365*10,
                       method = 'KM')
test.AUC = round(test.ROC[["AUC"]],3)

## Plot
test_ROC_data <- data.frame(TP=test.ROC$TP,FP=test.ROC$FP)
test_ROC_data %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank()) + 
  geom_label(label=test.AUC, x=0.5,y=0.5,
    label.padding = unit(0.55, "lines"))

######### correlation analysis #########

# age
train.cox$age_group = ifelse(train.cox$age>60,"age>60","age<=60")
ggplot(train.cox, aes(x = age_group, y = risk, fill = age_group)) +
  geom_jitter(aes(color = age_group),width = 0.2, size = 1) +
  labs(title = "Distribution of risk scores across age groups", x = "Stage", y = "Score") +
  theme_minimal()

# stage

ggplot(train.cox, aes(x = stage, y = risk, fill = stage)) +
  geom_jitter(aes(color = stage),width = 0.2, size = 1) +
  labs(title = "Distribution of risk scores across stages", x = "Stage", y = "Score") +
  theme_minimal()

# age
test.cox$age_group = ifelse(test.cox$age>60,"age>60","age<=60")
ggplot(test.cox, aes(x = age_group, y = risk, fill = age_group)) +
  geom_jitter(aes(color = age_group),width = 0.2, size = 1) +
  labs(title = "Distribution of risk scores across age groups", x = "Stage", y = "Score") +
  theme_minimal()

# stage

ggplot(test.cox, aes(x = stage, y = risk, fill = stage)) +
  geom_jitter(aes(color = stage),width = 0.2, size = 1) +
  labs(title = "Distribution of risk scores across stages", x = "Stage", y = "Score") +
  theme_minimal()

######### PCA #########

PCA <- prcomp(train.cox[,1:16])
pca_data <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = train.cox$group)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point() +
  labs(title = "PCA of training set", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

PCA <- prcomp(test.cox[,1:16])
pca_data <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = test.cox$group)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point() +
  labs(title = "PCA of testing set", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

