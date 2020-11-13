#import packages
install.packages("copula")
install.packages("QRM")
library(QRM)
library(copula)

#import data
data_Gaussian <- read.csv('/Users/Desktop/excel/gaussian copula data_r.csv')


##################################################
### marginal distribution parameter estimation ###
##################################################

#marginal: stock return (GBM), short rate r (Hull White)
#GBM parameter estimation (sigma) and calculate v 
stock_price <- data_Gaussian$stock.price..month.[1:144] #2008/01-2019/12
s_month <- length(stock_price)
stock_return <- log(stock_price[2:s_month]/stock_price[1:(s_month-1)])
GBM_sigma <- (var(stock_return)*12)^0.5
stock_hist_return <- data_Gaussian$stock.return.year.[1:11]
GBM_hist_mean <- rep(mean(stock_hist_return) , 11) #mu=mean+sigma^2/2
v <- pnorm(stock_hist_return, mean=GBM_hist_mean, sd=GBM_sigma) #change (1)
inv_v <- qnorm(v) #change (2) #turn into standard normal

#Hull White parameter estimation (a, sigma) and calculate u
short_rate <- data_Gaussian$short.rate[2:3047] #use up to 2019/12/31
r_t <- short_rate[3:3047]
r_t_minus1 <- short_rate[2:3046]
reg <- lm(r_t ~ r_t_minus1) 
summary <- summary(reg)
a_hat <- as.numeric(reg$coefficients[2])
sigma_hat <- summary$sigma
h <- 1/360
HW_a <- -log(a_hat)/h
HW_sigma <- sigma_hat*(-2*log(a_hat)/((1-a_hat^2)*h))^0.5 #HW_sigma means the parameter sigma
HW_hist_r <- data_Gaussian$short.rate.year.[2:12] #correspond to GBM

#calculate HW_hist_mean
hist_term <- cbind(data_Gaussian$term_1.year[2:12],data_Gaussian$term_2.year[2:12]) #2008-2018
hist_PM <- cbind(1/(1+hist_term[,1]), 1/(1+hist_term[,2])^2)
hist_fM01 <- -(log(hist_PM[,2])-log(hist_PM[,1]))
r_s <- data_Gaussian$short.rate.year.[1:11]
alpha_t <- hist_fM01 + 0.5*(HW_sigma^2/HW_a^2)*(1-exp(-HW_a))^2
alpha_s <- -log(hist_PM[,1])
HW_hist_mean <- r_s*exp(-HW_a) + alpha_t - alpha_s*exp(-HW_a)
HW_distribution_sigma <- (0.5*(HW_sigma^2/HW_a)*(1-exp(-2*HW_a)))^0.5
u <- pnorm(HW_hist_r, mean=HW_hist_mean, sd=HW_distribution_sigma ) #change (1)
inv_u <- qnorm(u) #change (2)


###################################
### copula parameter estimation ###
###################################

#Gaussian copula: calculate correlation matrix directly
corr <- cor(cbind(inv_v,inv_u))
BVN_rand <- rmvnorm(1000*120, mean=c(0,0), sigma=corr) #sampling from the BVN #mvtnorm package
Rho <- corr[1,2]
rand_vu <- pnorm(BVN_rand)  #reverse change (2) #change back to u and v (cumulative prob)

#Gaussian copula with copula package
data <- cbind(v,u)
copula <- summary(fitCopula(normalCopula(dim=2, dispstr="un"), data, method="ml")) 
Rho <- copula$coefficients[1]
norm.cop <- normalCopula(param=Rho,dim = 2, dispstr = "ex")
rand_vu <- rCopula(120*1000, norm.cop)

#Gaussian copula with QRM package
data <- cbind(v,u)
copula <- fit.gausscopula(Udata=data) #data is the area below
sigma <- copula$P
rand_vu <- rcopula.gauss(n=120*1000, Sigma=sigma)

#t copula with QRM package
data <- cbind(v,u)
copula <- fit.tcopula(Udata=data, method='Kendall')
sigma <- copula$P
df <- copula$nu
rand_vu <- rcopula.t(n=120*1000, Sigma=sigma, df=df)

#Archimedean copula (Frank, Clayton, AMH, Gumbel, Joe) with copula package
data <- cbind(v,u)
fit_copu <- fitCopula(copula=frankCopula(), data=data, method="ml") #or claytonCopula() and so on
summary <- summary(fit_copu)
para <- summary$coefficients[1]
copu <- archmCopula(family="frank", param=para) #or family="clayton" ...
rand_vu <- rCopula(copula=copu, n=120*1000)


###########################################
### derive stock return and bond return ###
###########################################

rand_v <- as.data.frame(matrix(0,1000,120))
rand_u <- as.data.frame(matrix(0,1000,120))
for(i in 1:120){
  rand_v[1:1000, i] <- rand_vu[(i*1000-999):(i*1000),1] #1 is the v part
  rand_u[1:1000, i] <- rand_vu[(i*1000-999):(i*1000),2] #2 is the u part
}

#use the random sample to derive bond/stock return simulation
#simulate stock return (GBM) 
term_struc_forward <- data_Gaussian$term_forward_2019.12.31[1:120] #term structure right now
GBM_future_mean <- term_struc_forward - GBM_sigma^2/2
simu_stock_return <- as.data.frame(matrix(0,1000,120))
for(i in 1:120){
  simu_stock_return[1:1000, i] <- qnorm(rand_v[1:1000,i], mean=GBM_future_mean[i], sd=GBM_sigma)
}
simu_stock_Rmean <- colMeans(simu_stock_return)
matplot(time(simu_stock_Rmean),simu_stock_Rmean, type = 'l',xlab="year",ylab="stock return", main="simulated stock return")

#simulate short rate (Hull-White)
t <- 1:120
term_struc_spot <- data_Gaussian$term_spot_2019.12.31[1:120]
price <- 1/(1+term_struc_spot)^t
log_price <- log(price)
forward <- -diff(log_price)/diff(t) #diff(t)=1 for all time
alpha <- forward[1:119] + (HW_sigma^2/(2*HW_a^2))*(1-exp(-HW_a*t[1:119]))^2
alpha_0 <- -log_price[1]
HW_future_mean <- as.data.frame(matrix(0,1000,119))
HW_future_mean[1:1000,1] <- rep(data_Gaussian$short.rate.year.[13],1000) #the short rate at 2019/12
for(i in 2:119){
  HW_future_mean[1:1000,i] <- HW_future_mean[1:1000,i-1]*exp(-HW_a)+alpha[i]-alpha[i-1]*exp(-HW_a)
}
matplot(time(colMeans(HW_future_mean)),colMeans(HW_future_mean), type = 'l',
        xlab="year",ylab="short rate", main="short rate distribution mean")
simu_r <- as.data.frame(matrix(0,1000,119))
for(i in 1:119){
  simu_r[1:1000, i] <- qnorm(rand_u[1:1000,i], mean=HW_future_mean[1:1000,i], sd=HW_distribution_sigma)
}
simu_r_mean <- colMeans(simu_r)
matplot(time(simu_r_mean),simu_r_mean, type = 'l',xlab="year",ylab="short rate", main="simulated r")

#turn short rate back to bond return
#calculate bond price
r0 <- data_Gaussian$short.rate.year.[13]
B <- (1/HW_a)*(1-exp(-HW_a)) #B is fixed
A <- rep(0,119)
A[1] <- price[1]*exp(-B*log_price[1]) #A[1]: do it separately
for(i in 2:119){
  A[i] <- (price[i]/price[i-1])*exp(B*forward[i-1]-(HW_sigma^2/(4*HW_a))*(1-exp(-2*HW_a*(i-1)))*B^2)
}
P <- as.data.frame(matrix(0,1000,119)) 
P[1:1000, 1] <- rep(A[1]*exp(-B*r0),1000) #the first column of P
for(i in 2:119){
  P[1:1000, i] <- A[i]*exp(-B*simu_r[1:1000,i-1])  #use r only to 118
}

#bond return (matrix)
R <- as.data.frame(matrix(0,1000,119)) 
R <- -log(P)
R_mean_copula <- colMeans(R)
matplot(time(R_mean_copula), R_mean_copula, type = 'l',xlab="year",ylab="bond return", main="bond return by H-W and copula")


#######################
### plot the result ###
#######################

matplot(time(GBM_return_ave), cbind(GBM_return_ave, simu_stock_Rmean), type = 'l',lty=c(1,1),
        col=c("blue","red"),xlab="year",ylab="stock return",main="stock return (t)")
legend(x=65,y=0,legend=c("marginal","copula"),lty=c(1,1),col=c("blue","red"))

matplot(time(R_ave), cbind(R_ave, R_mean_copula), type = 'l',lty=c(1,1),col=c("blue","red"),
        xlab="year",ylab="bond return",main="bond return (t)")
legend(x=65,y=0.02,legend=c("marginal","copula"),lty=c(1,1),col=c("blue","red"))

