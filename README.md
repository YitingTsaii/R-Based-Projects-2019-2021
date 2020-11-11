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
  * assign 3 hidden states for the model, each following a normal distribution ![normal dist](https://latex.codecogs.com/gif.latex?N%28%5Cmu_i%2C%20%5Csigma_i%29)
  * estimate all the parameters by forward, backward, and Viterbi algorithm (by Stan)
  
  
  ![transition matrix](https://latex.codecogs.com/gif.latex?A%20%3D%20%5Cbegin%7Bbmatrix%7D%20a_%7B11%7D%20%26%20a_%7B12%7D%20%26%20a_%7B13%7D%5C%5C%20a_%7B21%7D%20%26%20a_%7B22%7D%20%26%20a_%7B23%7D%5C%5C%20a_%7B31%7D%20%26%20a_%7B32%7D%20%26%20a_%7B33%7D%20%5Cend%7Bbmatrix%7D)
  ![noraml with reg](https://latex.codecogs.com/gif.latex?N%28%5Calpha_i&plus;%5Cbeta_%7Bi1%7Dx_1&plus;%5Cbeta_%7Bi2%7Dx_2%2C%5C%3B%20%5Csigma_i%29)
