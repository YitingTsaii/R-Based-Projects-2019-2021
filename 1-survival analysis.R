#import packages
library(tidyverse)  
library(survminer)
library(survival)
library(reshape2)  #correlation heatmap, melt
library(MASS)  #stepwise regression, stepAIC
library(rolr)  #rhier()

##########################
### Data Preprocessing ###
##########################

## < import data/change column names/merge data >
data <- read.csv(file="/Users/maggie/Desktop/vgh/excel/cervical_cancer_data.csv")
data_page4 <- data[1:41]
data_page3 <- data[42:48]
names(data_page4)[names(data_page4)=="id_page4"] <- "id"
names(data_page3)[names(data_page3)=="id_page3"] <- "id"
data_merge <- merge(data_page4, data_page3, by="id", all=T)
data_merge$name_page3 <- NULL
data_merge$name_page4 <- NULL

## < data wrangling >
#age, height, weight, BMI (class = numeric)
#summary(data_merge[2:5])  #we can discover some problems, such as the min of weight = 5.05, max = 999
data_merge[,3:4][sapply(data_merge[,3:4], function(x) x == 999)] <- NA
data_merge$BMI[is.na(data_merge$height)] <- NA
idx <- which(data_merge$weight == 5.05)
data_merge$weight[idx] <- 50.5
data_merge$BMI[idx] <- data_merge$weight[idx]/ (data_merge$height[idx]/100)^2

#marriage, menopause, high bloos pressure, diabetes (0/1, class = factor)
#G, P, A, grading (categorical data, class = factor)
#summary(data_merge[,6:13])
data_merge$Marriage[data_merge$Marriage == 999] <- NA
data_merge[,7:9][sapply(data_merge[,7:9], function(x) x == "" | x == "無資料" )] <- NA
data_merge[,7:9][sapply(data_merge[,7:9], function(x) x == "2-3" )] <- "2"
data_merge$grading[data_merge$grading == 9] <- NA
data_merge[,6:13] <- as.data.frame(lapply(data_merge[,6:13], as.factor))
data_merge$Marriage <- as.factor(data_merge$Marriage)

#c.stage, p.stage (categorical data, class = factor)
data_merge[,14:15][sapply(data_merge[,14:15], function(x) x == "BBB" | x == "888" | x == "999" | x == "1")] <- NA
data_merge[,14:15][sapply(data_merge[,14:15], function(x) x == "2b")] <- "2B"
data_merge[,14:15] <- as.data.frame(lapply(data_merge[,14:15], as.factor))

#FIGO.Stage.B and FIGO.Stage.A (use B) (categorical data, class = factor)
data_merge[c(16,17)] <- sapply(data_merge[c(16,17)], trimws)  #trim white space
data_merge$FIGO.Stage.B <- apply(data_merge[c(16,17)], 1, function(x) {ifelse(x[1] == "", x[2], x[1])}) #first use B, use A if B is missing
#which(data_merge$FIGO.Stage.B == "") #at row 153, 233, both A, B = ""
#not_same <- which(data_merge$FIGO.Stage.B != data_merge$FIGO.Stage.A & data_merge$FIGO.Stage.B != "" & data_merge$FIGO.Stage.A != "") #both are not "", and the two are not the same：114 126 131 136 137 141 168 174 183 217 235 243
#data_merge[not_same, c(16,17)]  #see the rows that A, B are not the same. Use A in row 126, otherwise, use B
data_merge$FIGO.Stage.B[126] <- data_merge$FIGO.Stage.A[126]
data_merge$FIGO.Stage.B[data_merge$FIGO.Stage.B == ""] <- NA
data_merge$FIGO.Stage.B <- as.factor(data_merge$FIGO.Stage.B)
data_merge$FIGO.Stage.A <- NULL
names(data_merge)[names(data_merge)=="FIGO.Stage.B"] <- "FIGO.Stage"

#tumor size (class = numeric)
data_merge$tumor.size[data_merge$tumor.size==999 | data_merge$tumor.size==990 | data_merge$tumor.size==994 | data_merge$tumor.size==998 | data_merge$tumor.size==992] <- NA

#chemotherapy, radiotherapy (delete these two columns, use RT, CT instead)
data_merge$Radiotherapy <- NULL
data_merge$Chemotherapy <- NULL

#recurrence (0/1, class = factor)
data_merge$Recurrence[data_merge$Recurrence=="" | data_merge$Recurrence=="nil"] <- NA
data_merge$Recurrence[data_merge$Recurrence=="1\n(meta)" | data_merge$Recurrence=="1 (meta)"] <- "1"
data_merge$Recurrence <- as.factor(data_merge$Recurrence)

#death status (0/1, class = factor)
#unique(data_merge$Death.date.1.death.0.Live.)
names(data_merge)[names(data_merge)=="Death.date.1.death.0.Live."] <- "death.status"
data_merge$death.status[data_merge$death.status==11] <- 1
#data_merge$death.status <- as.factor(data_merge$death.status)

#disease free interval
data_merge$Disease.free.interval..month.[data_merge$Disease.free.interval..month. == 0] <- NA

#recurrence date, first hospitalization
data_merge[c(19,20)] <- sapply(data_merge[c(19,20)], function(x) as.POSIXlt(x, format="%Y/%m/%d")) #convert to time variable, incorrect inputs will turn into NA automatically
diff <- as.double(difftime(data_merge$Recurrence.Date, data_merge$first.hospitalization, units="days"))
#diff[!is.na(diff)] #check that recurrence is no sooner than hospitaliztion
diff <- diff/30  #change the interval from day to month
#check whether there are rows with recurrence=1, but no interval data
#data_merge$Disease.free.interval..month.[data_merge$Recurrence == 1]  #recurrence = 1, but no disease free interval
#sum(is.na(data_merge$Disease.free.interval..month.[data_merge$Recurrence == 1]))  #38 of them

#disease free interval, disease free status
#calculate disease free interval on my own (recurrence data - first hospitalization)
data_merge$Disease.free.interval..month. <- as.double(difftime(data_merge$Recurrence.Date, data_merge$first.hospitalization, units="days"))/30
data_merge$Disease.free.interval..month. <-  apply(data_merge[c(21,24)], 1, function(x) {ifelse(is.na(x[1]), x[2], x[1])}) #first use disease free interval, if NA then use overall survival time
names(data_merge)[names(data_merge) == "Recurrence"] <- "disease.free.status"
#if disease free status = 0, death status = 1, then change disease free status to 1
data_merge$disease.free.status <- apply(data_merge[c(18,22)], 1, function(x) {ifelse(x[1]==0 & x[2]==1, x[2], x[1])})
data_merge$disease.free.status <- apply(data_merge[c(18,22)], 1, function(x) {ifelse(is.na(x[1]), x[2], x[1])})
data_merge$disease.free.status <- as.numeric(data_merge$disease.free.status)

#delete unnecessary columns
data_merge[c(19,20,23,25)] <- NULL

#CA125, CA153, CA199, CEA, SCC (class = numeric)
data_merge$CA125 <- as.numeric(data_merge$CA125) #nil is turned to NA automatically
data_merge[,22:26] <- sapply(data_merge[,22:26], as.numeric)
data_merge$CA153 <- NULL

#Biomarkers (class = numeric)
data_merge[,26:31] <- sapply(data_merge[,26:31], as.numeric)

#RT, CT
data_merge[,c(33,34)][sapply(data_merge[,c(33,34)], function(x) x == "N")] <- 0
data_merge[,c(33,34)][sapply(data_merge[,c(33,34)], function(x) x == "Y")] <- 1
data_merge[,33:34] <- as.data.frame(lapply(data_merge[,33:34], as.factor))

#smoke, betel nut, alcohol
data_merge[,35:37][sapply(data_merge[,35:37], function(x) x == "no")] <- 0
data_merge[,35:37][sapply(data_merge[,35:37], function(x) x == "yes")] <- 1
data_merge[,35:37][sapply(data_merge[,35:37], function(x) x == "null")] <- NA
data_merge[,35:37] <- as.data.frame(lapply(data_merge[,35:37], as.factor))

## < create a new column: cell type (squamous/adeno/others) >
squamous <- grepl("squamous",data_merge$Final.pathology, ignore.case = TRUE)
adeno <- grepl("adeno",data_merge$Final.pathology, ignore.case = TRUE)
data_merge$cell.type <- rep(0, nrow(data_merge))  #create a new column "cell.type"
data_merge$cell.type[squamous == TRUE & adeno == FALSE] <- "squamous"
data_merge$cell.type[squamous == FALSE & adeno == TRUE] <- "adeno"
data_merge$cell.type[squamous == FALSE & adeno == FALSE] <- "others"
data_merge$cell.type[squamous == TRUE & adeno == TRUE] <- "others"
data_merge$Final.pathology <- NULL 

# <delete CA125, CA199 due to too many NAs>
data_merge$CA125 <- NULL
data_merge$CA199 <- NULL


###################################
### Screen Out Noisy Covariates ###
###################################

## < correlation heatmap for continuous covariates (Pearson) >
#get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

#plot correlation heatmap (Pearson)
cts_input <- data.frame(data_merge[c(2:5,17,22:29)])
names(cts_input) <- c("age", "height", "weight", "BMI", "tumor_size", "CEA", "SCC", "UBE2C", "PDL1", "HPV18", "HPV58", "HPV16", "ASCC2")
cor_mar <- round(cor(cts_input, method="pearson", use="complete.obs"), 2)
upper_cor_mat <-  get_upper_tri(cor_mar)
melted_cor_mat <- melt(upper_cor_mat, na.rm = TRUE)
ggplot(data = melted_cor_mat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
  ggtitle("Heatmap for continuous inputs") +
  theme(plot.title = element_text(hjust = 0.5))

#test whether the corr between weight and BMI is significantly different from 0
cor.test(cts_input$weight, cts_input$BMI, method="pearson", use="complete.obs")

## < univariate cox regression >
inputs <- names(cts_input)
y <- data_merge[c(21,20)]
names(y) <- c("time", "status")
uni_cox_data <- cbind(y, cts_input)
for(input in inputs){
  fit <- coxph(formula = as.formula(paste("Surv(time, status) ~", input)), data = uni_cox_data)
  s <- summary(fit)
  #print(s)
  print(sprintf("%s: p-value=%f", input, s$coefficients[,5]))
}

## < Kaplan-Meier curve for categorical covariates >
categorical_input <- data.frame(data_merge[c(6:16,30:34)])
inputs <-  names(categorical_input)
categorical_input <- cbind(y, categorical_input)
splots <- list()
i <- 1
for(input in inputs){
  print(input)
  fit <- survfit(as.formula(paste("Surv(time, status) ~", input)), data = categorical_input) 
  print(fit)
  splots[[i]]= ggsurvplot(fit, pval = TRUE)
  i <- i + 1
}
arrange_ggsurvplots(splots,print = TRUE,ncol = 4,nrow = 4)
#only consider the categorical covariates with the p-value of log-rank test < 0.05

## < correlation heatmap for ordinal categorical covariates (Spearman) >
ordinal_input <- data.frame(data_merge[c(14:16)])
stages <- sapply(ordinal_input, as.numeric)
cor_mar <- round(cor(stages, method="spearman",use="complete.obs"), 4)
upper_cor_mat <-  get_upper_tri(cor_mar)
melted_cor_mat <- melt(upper_cor_mat, na.rm = TRUE)
ggplot(data = melted_cor_mat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
  ggtitle("Heatmap for 3 stage measurements (spearman)") +
  theme(plot.title = element_text(hjust = 0.5))

#test whether the correlations are significantly different from 0
cor.test(stages[,1], stages[,3], method="spearman",use="complete.obs") #c stage vs FIGO stage
cor.test(stages[,1], stages[,2], method="spearman",use="complete.obs") #c stage vs p stage
cor.test(stages[,2], stages[,3], method="spearman",use="complete.obs") #p stage vs FIGO stage

## < merging groups for FIGO stage >
#try FIGO.stage with level: 1.I 2.II 3.III 4.IV
FIGO.stage.level <- rep(NA, nrow(data_merge))
FIGO.stage.level[data_merge$FIGO.Stage=="IA" | data_merge$FIGO.Stage=="IA1" | data_merge$FIGO.Stage=="IA2" | data_merge$FIGO.Stage=="IB" | data_merge$FIGO.Stage=="IB1" | data_merge$FIGO.Stage=="IB2"] <- "I"
FIGO.stage.level[data_merge$FIGO.Stage=="IIA" | data_merge$FIGO.Stage=="IIA1" | data_merge$FIGO.Stage=="IIA2" | data_merge$FIGO.Stage=="IIB"] <- "II"
FIGO.stage.level[data_merge$FIGO.Stage=="IIIA" | data_merge$FIGO.Stage=="IIIB"] <- "III"
FIGO.stage.level[data_merge$FIGO.Stage=="IVA" | data_merge$FIGO.Stage=="IVB"] <- "IV"
data_merge <- cbind(data_merge, FIGO.stage.level)
data_merge$FIGO.stage.level <- factor(data_merge$FIGO.stage.level)
#plot survival curve from this level
fit <- survfit(Surv(Overall.survival.time..month., death.status) ~ FIGO.stage.level, data = data_merge) 
ggsurvplot(fit, pval = TRUE)
##pairwise log-rank test for this level
fit_pairwise <- pairwise_survdiff(Surv(Overall.survival.time..month., death.status) ~ FIGO.stage.level, data = data_merge)
round(fit_pairwise[['p.value']],4) #II and III are not significantly apart
data_merge$FIGO.stage.level <- NULL

#try FIGO.stage with level: 1.I 2.II 3.III 4.IV
FIGO.stage.level <- rep(NA, nrow(data_merge))
FIGO.stage.level[data_merge$FIGO.Stage=="IA" | data_merge$FIGO.Stage=="IA1" | data_merge$FIGO.Stage=="IA2" | data_merge$FIGO.Stage=="IB" | data_merge$FIGO.Stage=="IB1" | data_merge$FIGO.Stage=="IB2"] <- "I"
FIGO.stage.level[data_merge$FIGO.Stage=="IIA" | data_merge$FIGO.Stage=="IIA1" | data_merge$FIGO.Stage=="IIA2" | data_merge$FIGO.Stage=="IIB"] <- "II/III"
FIGO.stage.level[data_merge$FIGO.Stage=="IIIA" | data_merge$FIGO.Stage=="IIIB"] <- "II/III"
FIGO.stage.level[data_merge$FIGO.Stage=="IVA" | data_merge$FIGO.Stage=="IVB"] <- "IV"
data_merge <- cbind(data_merge, FIGO.stage.level)
data_merge$FIGO.stage.level <- factor(data_merge$FIGO.stage.level)
##plot survival curve from this level
fit <- survfit(Surv(Overall.survival.time..month., death.status) ~ FIGO.stage.level, data = data_merge) 
ggsurvplot(fit, pval = TRUE, pval.coord = c(0, 0.1))
##pairwise log-rank test for this level
fit_pairwise <- pairwise_survdiff(Surv(Overall.survival.time..month., death.status) ~ FIGO.stage.level, data = data_merge)
round(fit_pairwise[['p.value']],4) 

## < the variables remained after preliminary variable selection >
#biomarkers: UBE2C, PDL1, HPV18, HPV58, HPV16, ASCC2 
#continuous covariates: age, weight, tumor size, SCC
#categorical covariates: marriage, RT, CT, FIGO stage


###############################################################
### Find Cutpoints for Biomarkers & Stepwise Cox Regression ###
###############################################################

## < choosing cell type and endpoint >
#subset data_merge with different cell types
data_squamous <- subset(data_merge, cell.type == "squamous")
data_adeno <- subset(data_merge, cell.type == "adeno")
data_others <- subset(data_merge, cell.type == "others")

#choose cell type (squamous/adeno/others)
data_now <- data_squamous
#data_now <- data_adeno
#data_now <- data_others

#choose endpoint (overall survival/disease free survival)
y <- data.frame(data_now[c(21,20)])  #overall survival
#y <- data.frame(data_now[c(19,18)])  #disease free survival
names(y) <- c("time", "status")

#extract the variables 
cts_input <- data.frame(data_now[c(2,4,17,23)])
names(cts_input) <- c("age", "weight", "tumor_size", "SCC")
cat_input <- data.frame(data_now[c(6,30,31,36)])
names(cat_input) <- c("marriage", "RT", "CT", "FIGO_stage")
biomarkers <- data.frame(data_now[24:29])
names(biomarkers) <- c("UBE2C", "PDL1", "HPV18", "HPV58", "HPV16", "ASCC2")

## < one optimal cutpoint by survminer >
y_and_biomarkers <- cbind(y, biomarkers)
cutpoints <- surv_cutpoint(y_and_biomarkers, time = "time", event = "status", variables = c("UBE2C", "PDL1", "HPV18", "HPV58", "HPV16", "ASCC2"))
summary(cutpoints)
plot(cutpoints) #draw the cutpoint of all the variables
cut_category <- surv_categorize(cutpoints) #a dataframe with time, status, and all the categorized biomarkers
cut_category[,3:8] <- as.data.frame(lapply(cut_category[,3:8], as.factor))
cut_category[,3:7] <- as.data.frame(lapply(cut_category[,3:7], relevel, ref = "low")) #adjust the reference to "low" except ASCC2

#plot KM curves
inputs <-  c("UBE2C", "PDL1", "HPV18", "HPV58", "HPV16", "ASCC2")
splots <- list()
i <- 1
for(input in inputs){
  print(input)
  fit <- survfit(as.formula(paste("Surv(time, status) ~", input)), data = cut_category) 
  print(fit)
  splots[[i]]= ggsurvplot(fit, pval = TRUE)
  i <- i + 1
}
arrange_ggsurvplots(splots,print = TRUE,ncol = 3,nrow = 2)

#cox regression
cox_reg_data <- cbind(cut_category, cts_input, cat_input) 
fit <- coxph(formula = Surv(time, status) ~ ., data = cox_reg_data)
summary(fit)
#stepwise cox regression
cox_reg_data <- cox_reg_data[rowSums(is.na(cox_reg_data)) == 0,]
fit <- coxph(formula = Surv(time, status) ~ ., data = cox_reg_data) %>%
  stepAIC(trace = TRUE, direction = "both") 
summary(fit)


## < two optimal cutpoints by rolr >
y_and_biomarkers <- cbind(y, biomarkers)
y_and_biomarkers_no_NA <- y_and_biomarkers[rowSums(is.na(y_and_biomarkers)) == 0,]
cutpoints_2 <- data.frame(matrix(0, ncol = 6, nrow = 2))
names(cutpoints_2) <- inputs
row.names(cutpoints_2) <- c("left_cutpoint", "right_cutpoint")
for(i in 3:8){
  print(names(y_and_biomarkers_no_NA)[i])
  res=rhier(times=y_and_biomarkers_no_NA$time, status=y_and_biomarkers_no_NA$status, x=y_and_biomarkers_no_NA[,i], ns=15, alt='decrease')
  j <- i - 2
  cutpoints_2[1, j]<- res$best.splits.hier[1]
  cutpoints_2[2, j]<- res$best.splits.hier[2]
}

#transform the cts inputs to category inputs 
cut_category_2 <- data.frame(matrix(0, ncol = 6, nrow = nrow(data_now)))
names(cut_category_2) <- inputs
for(i in 1:6){
  cut_category_2[,i] <- ifelse(biomarkers[,i]> cutpoints_2[2,i], "high", ifelse(biomarkers[,i] < cutpoints_2[1,i], "low", "medium"))
}

#combine overall survival/disease free survival
cut_category_2 <- cbind(y, cut_category_2)
cut_category_2[,3:8] <- as.data.frame(lapply(cut_category_2[,3:8], as.factor))
cut_category_2[,3:7] <- as.data.frame(lapply(cut_category_2[,3:7], relevel, ref = "low")) #除了ASCC2之外，把reference都改成low

#KM curve
splots <- list()
i <- 1
for(input in inputs){
  print(input)
  fit <- survfit(as.formula(paste("Surv(time, status) ~", input)), data = cut_category_2) 
  print(fit)
  splots[[i]]= ggsurvplot(fit, pval = TRUE)
  i <- i + 1
}
arrange_ggsurvplots(splots,print = TRUE,ncol = 3,nrow = 2)

#cox regression
cox_reg_data <- cbind(cut_category_2, cts_input, cat_input) 
fit <- coxph(formula = Surv(time, status) ~ ., data = cox_reg_data)
summary(fit)
#stepwise cox regression
cox_reg_data <- cox_reg_data[rowSums(is.na(cox_reg_data)) == 0,]
fit <- coxph(formula = Surv(time, status) ~ ., data = cox_reg_data) %>%
  stepAIC(trace = TRUE, direction = "both") 
summary(fit)




