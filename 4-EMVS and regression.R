#import packages
library(EMVS)
library(bayesm)
library(Matrix)
library(ggplot2)

#import data
data_TEJ_orig <- read.csv("/Users/data/TEJ_and_original_data_to_r.csv") 
data_1858 <- read.csv("/Users/data/EMVS_1858_data_to_r.csv") 

#transform data (use 1858 inputs)
X1858 <- data_TEJ_orig[1:123, 4:1861] #123 periods
Y <- data_TEJ_orig[2:124, 3]


################################
### select variables by EMVS ###
################################

#settings (change v0 and v1 according to the setting of prior 1~5)
v0 <- seq(0.1, 2, length.out = 20)
v1 <- 1000
beta_init <- rep(1,ncol(X1858))  #(1,p) where p=the num of input variables 
sigma_init <- 1 
a <- 1
b <- 1
epsilon <- 10^{-5}

#variable selection with EMVS (5 results)
beta_init <- rep(1,ncol(X1858))
result1858 <- EMVS(Y, X1858, v0 = v0, v1 = v1, type = "betabinomial",
                  independent = FALSE, beta_init = beta_init, sigma_init = sigma_init,
                  epsilon = epsilon, a = a, b = b)
EMVSplot(result1858, "both", FALSE)
EMVSbest(result1858)


####################################################
### take the result of EMVS as regression inputs ###
####################################################

#use the variables selected by EMVS (prior 1, prior 4, prior 5)
Y <- data_1858[2:124, 3]
prior1 <- data_1858[1:130, c(4,6,7,8,11,16,20,21,22,25,26,27,29,31,32,35,36,39,43,46,48,49,51)]
prior4 <- data_1858[1:130, -c(1,2,3,18,22,36,40,47)]
prior5 <- data_1858[1:130, c(6,7,8,11,16,18,20,26,29,40,43,47,51)]
X_prior1 <- data_1858[1:123, c(4,6,7,8,11,16,20,21,22,25,26,27,29,31,32,35,36,39,43,46,48,49,51)]
X_prior4 <- data_1858[1:123, -c(1,2,3,18,22,36,40,47)]
X_prior5 <- data_1858[1:123, c(6,7,8,11,16,18,20,26,29,40,43,47,51)]

#take log for variables with large magnitude
prior1_log <- prior1
prior1_log[,c(1,7,8,9,20,22,23)] <- log(prior1[,c(1,7,8,9,20,22,23)])
prior4_log <- prior4
prior4_log[,c(1,2,7,16,17,19,38,39,41,42,43)] <- log(prior4[,c(1,2,7,16,17,19,38,39,41,42,43)])
prior5_log <- prior5
prior5_log[,c(7,12,13)] <- log(prior5[,c(7,12,13)])
X_prior1_log <- prior1_log[1:123,]
X_prior4_log <- prior4_log[1:123,]
X_prior5_log <- prior5_log[1:123,]

#data for prediction
pred_prior1_log <- prior1_log[124:130,]
pred_prior4_log <- prior4_log[124:130,]
pred_prior5_log <- prior5_log[124:130,]


####################
### GDP forecast ###
####################

#substitute prior 1 for prior 4 or 5
#regression --OLS 
Z1 <- as.data.frame(cbind(Y,X_prior1_log))
result1 <- lm(Z1)
summary(result1)

#regression --bayes 
ddata1<-list(y=Y, X=model.matrix(~. , data=X_prior1_log))
mmcmc1 <- list(R = 50000, nprint = 10000)
bayesreg1 <- runireg(Data = ddata1, Mcmc = mmcmc1)
summary(bayesreg1$betadraw)

#prediction --OLS
pred_prior1_OLS <- predict(lm(Z1), pred_prior1_log)  

#prediction --bayes
s1 <- as.data.frame(summary(bayesreg1$betadraw))
coeff_bayes1 <- s1$mean
pred_data1 <- cbind(1,pred_prior1_log)
pred_prior1_bayes <- rep(1,7)
for(i in 1:7){
  pred_prior1_bayes[i] <- as.numeric(pred_data1[i,]) %*% coeff_bayes1
}

#plot the result 
actual_y <- data_1858[125:131,3]
plot(pred_prior1_OLS, type="l")
plot(pred_prior1_bayes, type="l")
plot(actual_y, type="l")

