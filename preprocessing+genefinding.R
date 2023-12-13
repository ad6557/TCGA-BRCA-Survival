library(SummarizedExperiment)
library(TCGAbiolinks)
library(edgeR)
library(skimr)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(purrr)
library(glmnet)
library(EnhancedVolcano)
library(randomForestSRC)


######### get gene expression data #########
query_brca = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

#GDCdownload(query = query_brca)

tcga_data_BRCA = GDCprepare(query_brca)

table(tcga_data_BRCA@colData$definition)

######### RNAseq preprocessing #########

limma_pipeline = function(
    tcga_data,
    condition_variable,
    reference_group=NULL){
  
  design_factor = colData(tcga_data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(tcga_data),
                samples=colData(tcga_data),
                genes=as.data.frame(rowData(tcga_data)))
  
  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), n = Inf, sort.by="p")
  
  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the 100 most differentially expressed genes
    )
  )
}

limma_res = limma_pipeline(
  tcga_data=tcga_data_BRCA,
  condition_variable="definition",
  reference_group = "Solid Tissue Normal"
)

######### get DEGS #########

toptable <- limma_res$topGenes
pCutoff = 0.05/27289
EnhancedVolcano(toptable,
                lab = toptable$gene_name,
                x = 'logFC',
                y = 'P.Value',
                pCutoff = pCutoff,
                FCcutoff = 1,
                title = 'Differentially expressed genes from TCGA-BTCA',
                subtitle = 'cutoff: Fold Change = 2; p value = 0.05 (Bonferroni correction)')
EDGs = toptable %>% filter(abs(logFC)>1 & P.Value<=pCutoff) # 6640

######### prepare dataset for cox regression #########
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

coxdata = coxdata[limma_res$voomObj$targets$definition == "Primary solid Tumor",] # 1111
coxdata = coxdata[coxdata$gender == "female",] # 1099 female 12 male
coxdata$overall_survival <- ifelse(coxdata$vital_status == "Alive",
                                   coxdata$days_to_last_follow_up,
                                   coxdata$days_to_death)
table(coxdata$vital_status) # 945+153=1098
coxdata$vital_status =ifelse(coxdata$vital_status=="Dead",1,0)
table(coxdata$vital_status)

######### addressing abnormal values #########
naos = which(is.na(coxdata$overall_survival)) 
negativeost = which(coxdata$overall_survival<=0) 
coxdata = coxdata[-c(naos,negativeost),] # 1076

# map from gene_id to gene_name 
colnames(coxdata)[1:27289] == limma_res$voomObj$genes$gene_id
colnames(coxdata)[1:27289] <- limma_res$voomObj$genes$gene_name

######### split train test #########

set.seed(2023)
id = sample(1:1076, 800)
train = coxdata[id,]
test = coxdata[-id,]

######### univariate Cox model #########
uni_cox_p = data.frame(gene_id = EDGs$gene_id,
                       gene_name = EDGs$gene_name,
                       p=NA)

unicoxpvalue = function(x){
  data = train[,c(x,"overall_survival","vital_status")]
  fit <- coxph(Surv(overall_survival, vital_status) ~ ., data = data)
  return(summary(fit)$waldtest[[3]])
}

for(x in uni_cox_p$gene_name){
  uni_cox_p[uni_cox_p$gene_name==x,3]=unicoxpvalue(x)
}

select1 = uni_cox_p[uni_cox_p$p<0.05,2] # 736

dta <- train[,c(select1,"overall_survival","vital_status",
                "race","age","stage")]

######### Random survival forest #########

srf <- rfsrc(Surv(overall_survival,vital_status)~., dta,
             ntree = 1000, nodesize = 5, importance = TRUE, seed=2023)
print(srf)

length(srf$importance[srf$importance>0]) # 618
length(srf$importance[srf$importance>0.005]) # 9
select.srf = names(srf$importance[srf$importance>0])

# visualize it
srf.imp = sort(srf$importance,decreasing = TRUE)

impvarplot = data.frame(Gene=names(srf.imp)[1:50],
                        Importance=srf.imp[1:50])

rsf.var.plot = ggplot(impvarplot,aes(Importance, reorder(Gene,Importance))) +
  geom_bar(stat ="identity") +
  ylab("Gene Name") +
  ggtitle("Top Important Variables")
rsf.var.plot

######### Cox LASSO #########
Y<-cbind(time=dta$overall_survival, status=dta$vital_status)
Z<-as.matrix(dta[,select1])
cv.coxlasso <- cv.glmnet(Z, Y, alpha = 1, family = "cox", nfolds = 5)
plot(cv.coxlasso)
text(x=-5, y=200, "lambda.min = 0.028", cex=2, font=3)
text(x=-5, y=190, "log lambda.min = -3.58", cex=2, font=3)

cv.coxlasso$lambda.min # 0.028
bestlam = which(cv.coxlasso$lambda==cv.coxlasso$lambda.min)

coxlasso <- glmnet(Z, Y, alpha = 1, family = "cox")

# cv.coxlasso$lambda == coxlasso$lambda
sum(abs(coxlasso$beta[,bestlam])>1e-20) # 49
select.coxlasso = coxlasso[["beta"]]@Dimnames[[1]][abs(coxlasso$beta[,bestlam])>1e-20]


select.final = intersect(select.coxlasso,select.srf)


######### StepCox(both) analyses #########
data.step = train[,c(select.final,"overall_survival","vital_status",
                     "race","age","stage")]
StepCox <- step(coxph(Surv(overall_survival,vital_status)~.,data.step),
                direction = "both")
StepCox$coefficients
select.step = names(StepCox$coefficients)






save(limma_res, EDGs, train, test, dta, select.final,
     srf, cv.coxlasso, coxlasso, # StepCox, 
     select1, select.srf, select.coxlasso, # select.step, 
     file="preprocessing+genefinding.RData")

