# R-Based-Projects

## 1. Biomarker Identification
**Goal:** Segment patients into groups by biomarkers to promote precision medicine <br/>
**Problem:** Identify prognosis biomarkers for cervical cancer by survival analysis <br/>
**Methods:** <br/>
  * preliminarily select significant covariates by correlation heatmap, univariate cox regression, and log-rank test
  * subset the dataset by different cell types and focus on squamous cell carcinoma and adenocarcinoma
  * find 2~3 optimal cut points for each biomarker using `surv_cutpoint()` or `rhier()`
  * conduct full model selection by stepwise cox regression

