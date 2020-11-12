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
* **Goal:** Grasp the trend of future material prices to improve inventory control plan <br/>
* **Problem:** Predict future material prices by Hidden Markov Model <br/>
* **Methods:**  <br/>
  * Assign 3 hidden states: each follows a normal distribution, representing the low, medium, and high status
  * Construct the transition matrix: a 3 by 3 matrix storing the transition probabilities 
  * Estimate parameters: forward-backward algorithm, Viterbi algorithm (using Bayesian approach by `RStan`)
  
## 3. Bond/Stock Return Simulation (Stocahstic Models & Copula)
* **Goal:** Generate economic scenarios to help determine the optimal declared interest rate <br/>
* **Problem:** Simulate future bond return and stock return by stochastic models and Copula <br/>
* **Methods:**  <br/>
  * Simulate future bond return: Hull-White Model (short rate -> bond price -> bond return)
  * Simulate future stock return: Geometric Brownian Motion
  * Capture the correlation between bond and stock return: Gaussian copula or Archimedean copula 

## 4. Bayesian Variable Selection & GDP Forecast (EMVS & Regression)
* **Goal:** Construct proper linear regression models to forecast GDP <br/>
* **Problem:** Select variables by Expectation Maximization Variable Selection (EMVS) <br/>
* **Methods:**  <br/>
  * 
  * 
  * 
 
