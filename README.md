# R-Based-Projects

## 1. Biomarker Identification (Survival Analysis)
* **Goal:** Segment cervical cancer patients into groups by biomarkers to promote precision medicine <br/>
* **Problem:** Identify prognosis biomarkers for cervical cancer by survival analysis <br/>
* **Methods:**  <br/>
  * Select covariates preliminarily: correlation heatmap, univariate cox regression, log-rank test 
  * Subset the data by cell types: focus on squamous cell carcinoma and adenocarcinoma 
  * Find 2~3 optimal cut points for each biomarker: `surv_cutpoint()` and `rhier()` 
  * Select biomarkers and covariates: stepwise cox regression

## 2. Material Price Prediction (Hidden Markov Model)
* **Goal:** Grasp the trend of future material price to improve inventory control plan <br/>
* **Problem:** Predict future material price by Hidden Markov Model <br/>
* **Methods:**  <br/>
  * Assign 3 hidden states: each follows a normal distribution, representing the low, medium, and high status
  * 
  * estimate all the parameters, including mean and variance of normal distributions and elements in the transition matrix, by forward, backward, and Viterbi algorithm (by Stan)
  
## 3. Bond/Stock Return Simulation (Hull-White Model & Copula)
