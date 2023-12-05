# TCGA-BRCA-Survival
Analysis of RNAseq and clinical data of TCGA-BRCA female patients

## `preprocessing+genefinding.R` 
- download and preprocess RNAseq data
  - remove NAs in survival times
  - remove 12 males
- find differentially expressed genes (DEGs)
- **Univariate Cox proportional hazard regression analysis** over DEGs
  - **random survival forest (RSF)**
  - **stepCox (both)**
- **marker genes**: intersect of genes selected by RSF and stepCox

## `signature_building.R`
- **Multivariate Cox proportional hazard regression analysis**
  - adjusted by ~~race,~~ age, stage and marker genes
- **risk scores** = sum of [multicox coef(RNA)*count(RNA)]
- divide into high risk and low risk groups cut off by median of risk scores

## `validating.R`
- Kaplan‒Meier (K‒M) survival curves & log-rank test
- time-dependent ROC and AUC
- Principal Component Analysis (PCA)
- correlation analysis
  - risk score distribution by age
  - risk score distribution by stages

## `enrichment_analysis.R`
- gene ontology (GO)
- KEGG
