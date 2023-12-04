library(SummarizedExperiment)
library(TCGAbiolinks)
library(edgeR)
library(glmnet)
library(skimr)
library(survival)
library(survminer)
library(survivalROC)
library(dplyr)
library(tidyr)
library(purrr)
library(EnhancedVolcano)
library(survival)
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
                title = 'limma results',
                subtitle = 'Differential expression')
EDGs = toptable %>% filter(abs(logFC)>1 & P.Value<=pCutoff) # 6640
# FC=1.5
######### prepare dataset for cox regression #########
# RNAseq only
d_mat = t(limma_res$voomObj$E) # 1224*27289
coxdata = as.data.frame(d_mat)
# merge clinical
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

select1 = uni_cox_p[uni_cox_p$p<0.005,2] # 97

dta <- train[,c(select1,"overall_survival","vital_status")]

######### Random survival forest #########

srf <- rfsrc(Surv(overall_survival,vital_status)~., dta,
             ntree = 1000, nodesize = 5, importance = TRUE, seed=123)
print(srf)

length(srf$importance[srf$importance>0.005]) # 55
select.srf = names(srf$importance[srf$importance>0.005])

# visualize it
srf.imp = sort(srf$importance,decreasing = TRUE)

impvarplot = data.frame(Gene=names(srf.imp)[1:50],
                        Importance=srf.imp[1:50])

rsf.var.plot = ggplot(impvarplot,aes(Importance, reorder(Gene,Importance))) +
  geom_bar(stat ="identity") +
  ylab("Gene Name") +
  ggtitle("Top Important Variables")
rsf.var.plot

######### StepCox(both) analyses #########
StepCox <- step(coxph(Surv(overall_survival,vital_status)~.,dta),direction = "both")
StepCox$coefficients
select.step = names(StepCox$coefficients) # 50

######### Multivariate Cox #########

select.final = intersect(select.srf,select.step)


save(limma_res, EDGs, train, test, dta, srf, StepCox, 
     select.srf, select.step, select.final, file="12042023.RData")

