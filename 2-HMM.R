#import packages
library(rstan)
library(tidyverse)

#import data
AN_data0717 <- read.csv("/Users/Desktop/0717_AN_data.csv") 

#linear regression (entire data)
AN_0717_reg <- cbind(AN_data0717[2:119,2],AN_data0717[1:118,3:11])
AN_result <- lm(AN_0717_reg)
summary(AN_result)

#linear regression (now = 17)
AN_reg17 <- AN_0717_reg[1:17,]
AN_result17 <- lm(AN_reg17)
summary(AN_result17)

#use the real regression input data (now=17)
pre <- predict(lm(AN_reg17), AN_0717_reg[18:47,2:10])
plott17 <- as.data.frame(cbind(1:30, pre, AN_0717_reg[18:47, 1] ))

ggplot(plott17, aes(x=V1)) + 
  geom_line(aes(y = pre, colour = "predict"),linetype = "longdash") + 
  geom_line(aes(y = V3, colour = "actual")) +
  labs(y="AN price", x="time") +ggtitle("AN - Once")+
  scale_colour_manual("", values=c("black","red"))+
  theme(plot.title = element_text(hjust = 0.5))

#Stan code for HMM
hmm_model = "
functions {
  vector normalize(vector x) {
    return x / sum(x);
  }
}
data {
  int<lower=1> T;
  int<lower=1> K;
  real y[T];
}
parameters {
  // Discrete state model
  simplex[K] pi1;
  simplex[K] A[K];
  
  // Continuous observation model
  ordered[K] mu;
  real<lower=0> sigma[K];
}

transformed parameters {
  vector[K] logalpha[T];
  { // Forward algorithm log p(z_t = j | y_{1:t})
    real accumulator[K];
    logalpha[1] = log(pi1) + normal_lpdf(y[1] | mu, sigma);  //this is a vector
    for (t in 2:T) {
      for (j in 1:K) { // j = current (t)
      for (i in 1:K) { // i = previous (t-1)
                      // Murphy (2012) p. 609 eq. 17.48
                      // belief state      + transition prob + local evidence at t
        accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(y[t] | mu[j], sigma[j]);
        }
      logalpha[t, j] = log_sum_exp(accumulator);  //log_sum_exp: take exp, sum the exps, then take log
      }
    }
  } // Forward
}

model {   
  target += log_sum_exp(logalpha[T]); // Note: update based only on last logalpha
                                      //sum the K elts of logalpha[T], this is the P(O|lamda)
  mu ~ normal(mean(y), sd(y));
  sigma ~ exponential(1);
}

generated quantities {
  vector[K] alpha[T];
  
  vector[K] logbeta[T];
  vector[K] loggamma[T];
  
  vector[K] beta[T];
  vector[K] gamma[T];
  
  int<lower=1, upper=K> zstar[T];
  real logp_zstar;
  
  { // Forward algortihm
    for (t in 1:T)
      alpha[t] = softmax(logalpha[t]);
  } // Forward
  
  { // Backward algorithm log p(y_{t+1:T} | z_t = j)
    real accumulator[K];
    
    for (j in 1:K)
      logbeta[T, j] = 1;
    
    for (tforward in 0:(T-2)) {
      int t;
      t = T - tforward;
      
      for (j in 1:K) {    // j = previous (t-1)
      for (i in 1:K) {  // i = next (t)
                        // Murphy (2012) Eq. 17.58
                        // backwards t    + transition prob + local evidence at t
        accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(y[t] | mu[i], sigma[i]);
        }
      logbeta[t-1, j] = log_sum_exp(accumulator);
      }
    }
    
    for (t in 1:T)
      beta[t] = softmax(logbeta[t]);
  } // Backward
  
  { // forward-backward algorithm log p(z_t = j | y_{1:T})
    for(t in 1:T) {
      loggamma[t] = alpha[t] .* beta[t];
    }
    
    for(t in 1:T)
      gamma[t] = normalize(loggamma[t]);
  } // forward-backward
  
  
  
  { // Viterbi algorithm
    int bpointer[T, K]; // backpointer to the most likely previous state on the most probable path
    real delta[T, K];   // max prob for the sequence up to t
    // that ends with an emission from state k
    
    for (j in 1:K)
      delta[1, K] = normal_lpdf(y[1] | mu[j], sigma[j]);
    
    for (t in 2:T) {
      for (j in 1:K) { // j = current (t)
        delta[t, j] = negative_infinity();
        for (i in 1:K) { // i = previous (t-1)
          real logp;
          logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(y[t] | mu[j], sigma[j]);
          if (logp > delta[t, j]) {   //save the max till now
             bpointer[t, j] = i;
             delta[t, j] = logp;
          }
        }
      }
    }
    
    logp_zstar = max(delta[T]);
    
    for (j in 1:K)
      if (delta[T, j] == logp_zstar)
        zstar[T] = j;   //zstar[T] is the most probable state in time T
    
    for (t in 1:(T - 1)) {
      zstar[T - t] = bpointer[T - t + 1, zstar[T - t + 1]];
    }
  }
}"

#HMM functions
hmm_init <- function(K, y) {  # K: number of states, y: price training data
  clasif <- kmeans(y, K)
  init.mu <- by(y, clasif$cluster, mean) 
  init.sigma <- by(y, clasif$cluster, sd) 
  init.order <- order(init.mu)
  
  list(
    mu = init.mu[init.order], 
    sigma = init.sigma[init.order]
  ) 
}

hmm_VI_fit <- function(K, y) {   # K: number of states, y: price training data
  rstan_options(auto_write = TRUE) 
  options(mc.cores = parallel::detectCores())
  
  stan.model = stan_model(model_code = hmm_model) 
  stan.data = list(
    T = length(y), 
    K = K,
    y = y
  )
  vb(stan.model,
     data = stan.data,
     iter = 20000,
     #init = function(){hmm_init(K,y)}
     init = "random"
     )
}

hmm_predict <- function(data, T, K) {  # data: put in  x_test, T: predict next T periods, K: number of states
  # 1. Parameters
  A = matrix(data[1:K^2], nrow = K)
  mu = as_vector(data[(K^2+1):(K^2+K)])
  sigma = as_vector(data[(K^2+K+1):(K^2+2*K)])
  Z = unlist(data[(K^2+2*K+1)])
  
  # 2. Hidden path
  z = vector("numeric", T+1)
  z[1] = Z
  for (t in 2:(T+1)){
    z[t] = sample(1:K, size = 1, prob = A[z[t - 1], ])
  }
  
  # 3. prediction
  y = vector("numeric", T)
  for (t in 1:T){
    y[t] = rnorm(1, mu[z[t+1]], sigma[z[t+1]])
  }
  list(y = y, z = z[2:length(z)])
}

# Try K = 3, T = 10, train = 10, now = 17
#predict
input_AN17 <- as.data.frame(matrix(1:270, nrow=30))
K <- 3
TT <- 30

for(i in 2:10){
  fit <- hmm_VI_fit(K, AN_0717_reg[8:17, i])
  result <- as.data.frame(fit)
  x_test <- result[1:1000,c((K+1):(K^2+3*K), ncol(result)-2)]
  pred =  x_test %>%
    apply(., MARGIN=1, FUN=hmm_predict, T=TT, K=K) %>%
    lapply(., function (x) x[c('y')]) %>% 
    unlist() %>% 
    matrix(, nrow = TT)   
  pred_t <- t(pred) 

  for(j in 1:30){
    input_AN17[j, (i-1)] <- mean(pred_t[,j])
  }
}

# use the predict regression input data (now=17)
colnames(input_AN17) <- c("ABS_demand", "AF_demand", "AN_Supply", "C3_FOB_KOR", "C3_CFR_China", "AN_CNY", "ABS_OP", "AF_OP", "AN_OP")
pre <- predict(lm(AN_reg17), input_AN17)
plott17 <- as.data.frame(cbind(1:30, pre, AN_0717_reg[18:47, 1] ))

ggplot(plott17, aes(x=V1)) + 
  geom_line(aes(y = pre, colour = "predict"),linetype = "longdash") + 
  geom_line(aes(y = V3, colour = "actual")) +
  labs(y="AN price", x="time") +ggtitle("AN - Once")+
  scale_colour_manual("", values=c("black","red"))+
  theme(plot.title = element_text(hjust = 0.5))

#scale
plott17_scale <- as.data.frame(cbind(1:30, scale(pre), scale(AN_0717_reg[18:47, 1] )))

ggplot(plott17_scale, aes(x=V1)) + 
  geom_line(aes(y = V2, colour = "predict"),linetype = "longdash") + 
  geom_line(aes(y = V3, colour = "actual")) +
  labs(y="AN price", x="time") +ggtitle("AN - Once")+
  scale_colour_manual("", values=c("black","red"))+
  theme(plot.title = element_text(hjust = 0.5))

#write the data frame to csv
AN17 <- as.data.frame(cbind(AN_0717_reg[18:47, 1],pre,scale(AN_0717_reg[18:47, 1]),scale(pre)))
colnames(AN17) = c( "AN","AN_predict","AN_scale", "AN_predict_scale")
write.csv(AN17,"/Users/Desktop/0717_predict_price/AN_reg17.csv", row.names = FALSE)


