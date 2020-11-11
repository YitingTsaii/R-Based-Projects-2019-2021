# R-Based-Projects

## 1. Biomarker Identification (Survival Analysis)
* **Goal:** Segment patients into groups by biomarkers to promote precision medicine <br/>
* **Problem:** Identify prognosis biomarkers for cervical cancer by survival analysis <br/>
* **Methods:**  <br/>
  * preliminarily select covariates by correlation heatmap, univariate cox regression, and log-rank test 
  * subset the dataset by different cell types and focus on squamous cell carcinoma and adenocarcinoma 
  * find 2~3 optimal cut points for each biomarker using `surv_cutpoint()` or `rhier()` 
  * conduct full model selection by stepwise cox regression

## 2. Material Price Prediction (Hidden Markov Model)
* **Goal:** Grasp the trend of future material price to improve inventory control plan <br/>
* **Problem:** Predict future material price by Hidden Markov Model <br/>
* **Methods:**  <br/>
  * assign 3 hidden states for the model, each following a normal distribution 
  * estimate all the parameters, including mean and variance of normal distributions and elements in the transition matrix, by forward, backward, and Viterbi algorithm (by Stan)
  
## 3. Bond/Stock Return Simulation (Hull-White Model & Copula)
