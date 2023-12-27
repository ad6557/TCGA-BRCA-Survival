library(WGCNA)
library(DESeq2)
library(tidyverse)
library(gridExtra)

allowWGCNAThreads()
load("preprocessing+genefinding.RData")

# RNAseq only
d_mat = t(limma_res$voomObj$E) # 1224*27289
coxdata = as.data.frame(d_mat)

# merge demographic
coxdata$race=as.factor(limma_res$voomObj$targets$race)
coxdata$age=as.numeric(limma_res$voomObj$targets$age_at_index)

coxdata$stage=as.factor(limma_res$voomObj$targets$paper_pathologic_stage)
coxdata$stage[which(is.na(coxdata$stage))] <- "NA"
# levels(coxdata$stage)[levels(coxdata$stage)=="NA"] <- NA
# too many NAs in this variable -- can be problematic if just toss them

# merge clinical
coxdata$vital_status=limma_res$voomObj$targets$vital_status
coxdata$days_to_last_follow_up=limma_res$voomObj$targets$days_to_last_follow_up
coxdata$days_to_death=limma_res$voomObj$targets$days_to_death
coxdata$gender=as.factor(limma_res$voomObj$targets$gender)

coxdata$definition = limma_res$voomObj$targets$definition

coxdata = coxdata[coxdata$gender == "female",]
coxdata$overall_survival <- ifelse(coxdata$vital_status == "Alive",
                                   coxdata$days_to_last_follow_up,
                                   coxdata$days_to_death)
table(coxdata$vital_status) # 1013+197
coxdata$vital_status =ifelse(coxdata$vital_status=="Dead",1,0)
table(coxdata$vital_status)

######### addressing abnormal values #########
naos = which(is.na(coxdata$overall_survival)) 
negativeost = which(coxdata$overall_survival<=0) 
coxdata = coxdata[-c(naos,negativeost),] # 1187

# map from gene_id to gene_name 
# colnames(coxdata)[1:27289] == limma_res$voomObj$genes$gene_id
colnames(coxdata)[1:27289] <- limma_res$voomObj$genes$gene_name

# Get all differentially expressed genes into a dataframe

EDGdata = coxdata[,EDGs$gene_name]

# detect outlier genes 
# (nonono I should have done it before normalization!)

gsg <- goodSamplesGenes(EDGdata)
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(EDGdata,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)


# convert matrix to numeric
EDGdata[] <- sapply(EDGdata, as.numeric)

soft_power <- 10 # high R square low connectivity
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(EDGdata,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)




# grey module = all genes that doesn't fall into other modules were assigned to the grey module


# 6A. Relate modules to traits --------------------------------------------------
# module trait associations


traits <- data.frame(status = coxdata$vital_status)


# 6B. Intramodular analysis: Identifying driver genes ---------------


# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

nSamples = 1187
module.membership.measure <- cor(module_eigengenes, EDGdata, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:10,1:10]


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(EDGdata, traits$status, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


head30 = gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(30)

save(gene.signf.corr.pvals, file = "WGCNA.RData")


