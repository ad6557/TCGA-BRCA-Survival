# TCGA-BRCA-Survival
Analysis of RNAseq and clinical data of TCGA-BRCA female patients

## Course: Advanced Survival Analysis-Fall 2023
### Instructor: Dr. Yongzhao Shao

## `preprocessing+genefinding.R` 
- download and preprocess RNAseq data
  - filtering
  - normalization
  - remove NAs in survival times
  - remove 12 males
- find **Differentially Expressed Cancer Driver Genes** (DEGs)
  - both male and female
- **Univariate Cox proportional hazard regression analysis** over DEGs
  - **random survival forest (RSF)**
    - adjusted by race, age, stage
  - ~~stepCox (both)~~
  - **Cox LASSO**
- **marker genes**: intersect of genes selected by RSF and Cox LASSO

## `signature_building.R`
- **Multivariate Cox proportional hazard regression analysis**
  - toss insignificant marker genes and race
  - in the end adjusted by ~~race,~~ age, stage and 10 of marker genes
- **risk scores** = sum of [multicox coef(RNA)*count(RNA)] + coef(age)*age + coef(stage)*factor(stage)
- divide into high risk and low risk groups cut off by median of risk scores

## `signature_building_frailty.R`
- add intercept frailty term and test/compare AIC

## `validating.R`
- Kaplan‒Meier (K‒M) survival curves & log-rank test
- time-dependent ROC and AUC
  - 5 years
  - 10 years
- Principal Component Analysis (PCA)
- correlation analysis
  - risk score distribution by age
  - risk score distribution by stages
  - kruskal-wallis test

## `WCGNA.R`
- determined that a soft-thresholding power
- dendrogram
- significant genes

## `enrichment_analysis.R`
- gene ontology (GO)
- KEGG
