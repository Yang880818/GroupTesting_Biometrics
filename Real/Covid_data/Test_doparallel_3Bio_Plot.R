rm(list = ls())
library(foreach)
library(doParallel)
library(combinat)
library(bbmle)
library(Rfast)
library(MASS) 
library(moments)
library(latex2exp)

CurrentData <- Sys.Date()

SEED.NUM      = 10000
#KK            = 1
Pi0           = 0.96
Pi1           = 0.952
SEED.NUM.BOOT = 1000

time1 <- Sys.time()
seed.num = SEED.NUM



#setwd("C:/Users/yangj21/Desktop/Real/Covid_Data")
Covid_Data <- read.csv("./Data/Data_final_new_2023.csv", header = T)
dim(Covid_Data)

Covid_Data$Self_Antibody[which(Covid_Data$Self_method == 2)] <- NA

## Diseased status of self antibody test
sum(is.na(Covid_Data$Self_Antibody))
length(which(Covid_Data$Self_Antibody == 0))
sum(Covid_Data$Self_Antibody[which(Covid_Data$Self_Antibody == 1)])
length(which(Covid_Data$Self_Antibody == 2))


#### We will use self antibody test result as corresponding variables
Covid_Data$Self_Antibody[which(Covid_Data$Self_Antibody == 2)] <- NA ### Ignore Ig G positive
length(which(Covid_Data$Self_Antibody == 0))
sum(Covid_Data$Self_Antibody[which(Covid_Data$Self_Antibody == 1)])
length(which(Covid_Data$Self_Antibody == 2))



#### We will use self antibody test result as corresponding variables
Covid_Data$Self_Antibody[which(Covid_Data$Self_Antibody == 2)] <- NA  ### Ignore Ig G positive

Data_resample <- na.omit(Covid_Data[c("Self_Antibody", "Total.bilirubin", "Alanine.aminotransferase", "Aspartate.aminotransferase")])
dim(Data_resample)
N = 11830

# Data_resample <- na.omit(Covid_Data[c("Self_Antibody", "Gamma.glutamyltransferase", "Triglycerides", "LDL.direct")])
# dim(Data_resample)
# N = 11000

colnames(Data_resample) <- c("Disease_status", "Biomarker_1", "Biomarker_2", "Biomarker_3")
set.seed(seed.num)
Data_resample_final <- Data_resample[sample(nrow(Data_resample), N, replace = FALSE), ]
dim(Data_resample_final)

Ori_Biomarker_1 <- Data_resample_final$Biomarker_1
Ori_Biomarker_2 <- Data_resample_final$Biomarker_2
Ori_Biomarker_3 <- Data_resample_final$Biomarker_3


Biomarker_1 <- scale(Ori_Biomarker_1, center = FALSE, scale = TRUE)
Biomarker_2 <- scale(Ori_Biomarker_2, center = FALSE, scale = TRUE)
Biomarker_3 <- scale(Ori_Biomarker_3, center = FALSE, scale = TRUE)

# Test_Biomarker_1_new <- log(Biomarker_1)
# Data_resample_final$Biomarker_1 <- Test_Biomarker_1_new
# Test_Biomarker_2_new <- log(Biomarker_2)
# Data_resample_final$Biomarker_2 <- Test_Biomarker_2_new
# Test_Biomarker_3_new <- log(Biomarker_3)
# Data_resample_final$Biomarker_3 <- Test_Biomarker_3_new

# r1 <- lm(Test_Biomarker_1_new~1)
# hist(r1$residuals)
# skewness(r1$residuals)
# 
# r2 <- lm(Test_Biomarker_2_new~1)
# hist(r2$residuals)
# skewness(r2$residuals)
# 
# r3 <- lm(Test_Biomarker_3_new~1)
# hist(r3$residuals)
# skewness(r3$residuals)
# 


## Test Biomarker_1
Test_Biomarker_1 <- lm(Biomarker_1~1)
hist(Test_Biomarker_1$residuals)
skewness(Test_Biomarker_1$residuals)

b <- boxcox(Biomarker_1~1)
lambda <- b$x
lik <- b$y
bc <- cbind(lambda, lik)
Test_lambda_Biomarker_1 <- bc[order(-lik),][lik][1]
Test_Biomarker_1_new = (Biomarker_1^Test_lambda_Biomarker_1 - 1)/Test_lambda_Biomarker_1
hist(Test_Biomarker_1_new)
r0 <- lm(Test_Biomarker_1_new~1)
hist(r0$residuals)
skewness(r0$residuals)
Data_resample_final$Biomarker_1 <- Test_Biomarker_1_new



### Test Biomarker_2
Test_Biomarker_2 <- lm(Biomarker_2~1)
hist(Test_Biomarker_2$residuals)
skewness(Test_Biomarker_2$residuals)

b <- boxcox(Biomarker_2~1)
lambda <- b$x
lik <- b$y
bc <- cbind(lambda, lik)
Test_lambda_Biomarker_2 <- bc[order(-lik),][lik][1]
Test_Biomarker_2_new = (Biomarker_2^Test_lambda_Biomarker_2 - 1)/Test_lambda_Biomarker_2
hist(Test_Biomarker_2_new)
r0 <- lm(Test_Biomarker_2_new~1)
hist(r0$residuals)
skewness(r0$residuals)
Data_resample_final$Biomarker_2 <- Test_Biomarker_2_new


### Test Biomarker_3
Test_Biomarker_3 <- lm(Biomarker_3~1)
hist(Test_Biomarker_3$residuals)
skewness(Test_Biomarker_3$residuals)

b <- boxcox(Biomarker_3~1)
lambda <- b$x
lik <- b$y
bc <- cbind(lambda, lik)
Test_lambda_Biomarker_3 <- bc[order(-lik),][lik][1]
Test_Biomarker_3_new = (Biomarker_3^Test_lambda_Biomarker_3 - 1)/Test_lambda_Biomarker_3
hist(Test_Biomarker_3_new)
r0 <- lm(Test_Biomarker_3_new~1)
hist(r0$residuals)
skewness(r0$residuals)
Data_resample_final$Biomarker_3 <- Test_Biomarker_3_new


### Non-Dieseasd Biomarkers
Test_Biomarker_1_0_new <- Data_resample_final$Biomarker_1[which(Data_resample_final$Disease_status == 0)]
Test_Biomarker_2_0_new <- Data_resample_final$Biomarker_2[which(Data_resample_final$Disease_status == 0)]
Test_Biomarker_3_0_new <- Data_resample_final$Biomarker_3[which(Data_resample_final$Disease_status == 0)]
### Dieseasd Biomarkers
Test_Biomarker_1_1_new <- Data_resample_final$Biomarker_1[which(Data_resample_final$Disease_status == 1)]
Test_Biomarker_2_1_new <- Data_resample_final$Biomarker_2[which(Data_resample_final$Disease_status == 1)]
Test_Biomarker_3_1_new <- Data_resample_final$Biomarker_3[which(Data_resample_final$Disease_status == 1)]




## Mean of Non-diseased
mu01 = mean(Test_Biomarker_1_0_new)
mu01
mu02 = mean(Test_Biomarker_2_0_new)
mu02
mu03 = mean(Test_Biomarker_3_0_new)
mu03

## Mean of Diseased
mu11 = mean(Test_Biomarker_1_1_new)
mu11
mu12 = mean(Test_Biomarker_2_1_new)
mu12
mu13 = mean(Test_Biomarker_3_1_new)
mu13



## Sd of Non-diseased
sd01 = sd(Test_Biomarker_1_0_new)
sd01
sd02 = sd(Test_Biomarker_2_0_new)
sd02
sd03 = sd(Test_Biomarker_3_0_new)
sd03

## Sd of Diseased
sd11 = sd(Test_Biomarker_1_1_new)
sd11
sd12 = sd(Test_Biomarker_2_1_new)
sd12
sd13 = sd(Test_Biomarker_3_1_new)
sd13



## Correlation Coefficients in Non-diseased
rho0_12 = cor(Test_Biomarker_1_0_new, Test_Biomarker_2_0_new)
rho0_12
rho0_13 = cor(Test_Biomarker_1_0_new, Test_Biomarker_3_0_new)
rho0_13
rho0_23 = cor(Test_Biomarker_2_0_new, Test_Biomarker_3_0_new)
rho0_23


## Correlation Coefficients in Diseased
rho1_12 = cor(Test_Biomarker_1_1_new, Test_Biomarker_2_1_new)
rho1_12
rho1_13 = cor(Test_Biomarker_1_1_new, Test_Biomarker_3_1_new)
rho1_13
rho1_23 = cor(Test_Biomarker_2_1_new, Test_Biomarker_3_1_new)
rho1_23


#####################################
Bio0 <- cbind(Test_Biomarker_1_0_new, Test_Biomarker_2_0_new, Test_Biomarker_3_0_new)
Bio1 <- cbind(Test_Biomarker_1_1_new, Test_Biomarker_2_1_new, Test_Biomarker_3_1_new)

mu0 = c(mu01, mu02, mu03)
mu0
mu1 = c(mu11, mu12, mu13)
mu1

sd0 = c(sd01, sd02, sd03)
sd0
sd1 = c(sd11, sd12, sd13)
sd1

rho0 = c(rho0_12, rho0_13, rho0_23)
rho0
rho1 = c(rho1_12, rho1_13, rho1_23)
rho1

c(mu0, mu1, sd0, sd1, rho0, rho1)

###################################################### AUC True ######################################################

### 3 Biomarkers (1/2/3)
mu_true <- as.matrix(mu1 - mu0)
mu_true
sigma_0_true = matrix(c(sd01^2, rho0_12 * sd01 * sd02, rho0_13 * sd01 * sd03, 
                        rho0_12 * sd02 * sd01, sd02^2, rho0_23 * sd02 * sd03,
                        rho0_13 * sd03 * sd01, rho0_23 * sd03 * sd02, sd03^2), 3, 3)
sigma_0_true
sigma_1_true = matrix(c(sd11^2, rho1_12 * sd11 * sd12, rho1_13 * sd11 * sd13, 
                        rho1_12 * sd12 * sd11, sd12^2, rho1_23 * sd12 * sd13, 
                        rho1_13 * sd13 * sd11, rho1_23 * sd13 * sd12, sd13^2), 3, 3)
sigma_1_true
AUC_true <- pnorm(sqrt(t(mu_true) %*% solve(sigma_1_true + sigma_0_true) %*% mu_true))
AUC_true

AUC_1 <- pnorm(sqrt(t(mu11-mu01) %*% solve(sd01^2 + sd11^2) %*% (mu11-mu01)))
AUC_1
AUC_2 <- pnorm(sqrt(t(mu12-mu02) %*% solve(sd02^2 + sd12^2) %*% (mu12-mu02)))
AUC_2
AUC_3 <- pnorm(sqrt(t(mu13-mu03) %*% solve(sd03^2 + sd13^2) %*% (mu13-mu03)))
AUC_3
###########################################################################
nodeNum = 100
CDF_true <-  seq(from = 1/nodeNum, to = 1 - 1/nodeNum, by = 1/nodeNum)
CDF_true <- c(0, CDF_true, 1)
knotNum <- length(CDF_true)
aa <- CDF_true

## ROC_1
ROC.Fun_1 <- function(u)
{
  return(pnorm((mu11 - mu01 + sd01 * qnorm(u))/sd11))
}
ROC.true_1 <- ROC.Fun_1(aa)
plot(aa, ROC.true_1, type = "l", xlab = "u", ylab = TeX(r"($ROC(u)$)"), mgp = c(1.5,0.5,0), cex.lab = 2, lwd = 2, cex.axis = 1.5)
abline(0, 1, col = "red", lwd = 1, lty = 1)

## ROC_2
ROC.Fun_2 <- function(u)
{
  return(pnorm((mu12 - mu02 + sd02 * qnorm(u))/sd12))
}
ROC.true_2 <- ROC.Fun_2(aa)
plot(aa, sort(1 - ROC.true_2), type = "l", xlab = "u", ylab = TeX(r"($ROC(u)$)"), mgp = c(1.5,0.5,0), cex.lab = 2, lwd = 2, cex.axis = 1.5)
abline(0, 1, col = "red", lwd = 1, lty = 1)

## ROC_3
ROC.Fun_3 <- function(u)
{
  return(pnorm((mu13 - mu03 + sd03 * qnorm(u))/sd13))
}
ROC.true_3 <- ROC.Fun_3(aa)
plot(aa, ROC.true_3, type = "l", xlab = "u", ylab = TeX(r"($ROC(u)$)"), mgp = c(1.5,0.5,0), cex.lab = 2, lwd = 2, cex.axis = 1.5)
abline(0, 1, col = "red", lwd = 1, lty = 1)


###### AUC Combination
c <- t(mu_true) %*% solve(sigma_1_true + sigma_0_true)
c

linear_alpha <- (1/c[,1]) * t(mu_true) %*% solve(sigma_1_true + sigma_0_true)
linear_alpha

new_mu <- c(linear_alpha %*% mu_true)
new_mu
new_sigma0 <- c(linear_alpha %*% sigma_0_true %*% t(linear_alpha))
new_sigma0
new_sigma1 <- c(linear_alpha %*% sigma_1_true %*% t(linear_alpha))
new_sigma1

AUC_try <- pnorm(sqrt(t(new_mu) %*% solve(new_sigma0 + new_sigma1) %*% new_mu))
AUC_try
ROC.Fun <- function(u)
{
  return(pnorm((new_mu + sqrt(new_sigma0) * qnorm(u))/sqrt(new_sigma1)))
}

ROC_true <- ROC.Fun(aa)
plot(aa, ROC_true, type = "l", xlab = "u", ylab = TeX(r"($ROC(u)$)"), mgp = c(1.5,0.5,0), cex.lab = 2, lwd = 2, cex.axis = 1.5)
abline(0, 1, col = "red", lwd = 1, lty = 1)
######################################################################################################################

###################################################### AUC Est. ######################################################
CDF_est <-  seq(from = 1/nodeNum, to = 1 - 1/nodeNum, by = 1/nodeNum)
CDF_est <- c(0, CDF_est, 1)
knotNum <- length(CDF_est)
aa <- CDF_est


### For K = 1
Result_1 <- read.csv("./results/Boot/TTB_scaleSep_Seed.boot=1000_K=1_0.96_0.952_2023-06-27.csv", header = T)
Result_1
## Mean of Diseased and Non-diseased
mu01_est_1 = Result_1[1,2]
mu01_est_1
mu02_est_1 = Result_1[2,2]
mu02_est_1
mu03_est_1 = Result_1[3,2]
mu03_est_1
mu11_est_1 = Result_1[1,3]
mu11_est_1
mu12_est_1 = Result_1[2,3]
mu12_est_1
mu13_est_1 = Result_1[3,3]
mu13_est_1
## Sd of Diseased and Non-diseased
sd01_est_1 = Result_1[1,4]
sd01_est_1
sd02_est_1 = Result_1[2,4]
sd02_est_1
sd03_est_1 = Result_1[3,4]
sd03_est_1
sd11_est_1 = Result_1[1,5]
sd11_est_1
sd12_est_1 = Result_1[2,5]
sd12_est_1
sd13_est_1 = Result_1[3,5]
sd13_est_1

## Correlation Coefficients in Diseased and Non-diseased
rho0_12_est_1 = Result_1[1,6]
rho0_12_est_1
rho0_13_est_1 = Result_1[2,6]
rho0_13_est_1
rho0_23_est_1 = Result_1[3,6]
rho0_23_est_1
rho1_12_est_1 = Result_1[1,7]
rho1_12_est_1
rho1_13_est_1 = Result_1[2,7]
rho1_13_est_1
rho1_23_est_1 = Result_1[3,7]
rho1_23_est_1

### ROC_1_Est
ROC.Fun_1_1 <- function(u)
{
  return(pnorm((mu11_est_1 - mu01_est_1 + sd01_est_1 * qnorm(u))/sd11_est_1))
}
ROC.est_1_1 <- ROC.Fun_1_1(aa)
### ROC_2_Est
ROC.Fun_2_1 <- function(u)
{
  return(pnorm((mu12_est_1 - mu02_est_1 + sd02_est_1 * qnorm(u))/sd12_est_1))
}
ROC.est_2_1 <- ROC.Fun_2_1(aa)
### ROC_3_Est
ROC.Fun_3_1 <- function(u)
{
  return(pnorm((mu13_est_1 - mu03_est_1 + sd03_est_1 * qnorm(u))/sd13_est_1))
}
ROC.est_3_1 <- ROC.Fun_3_1(aa)

#########################################################################################################################################
### For K = 2
Result_2 <- read.csv("./results/Boot/TTB_scaleSep_Seed.boot=1000_K=2_0.96_0.952_2023-06-28.csv", header = T)
Result_2
## Mean of Diseased and Non-diseased
mu01_est_2 = Result_2[1,2]
mu02_est_2 = Result_2[2,2]
mu03_est_2 = Result_2[3,2]

mu11_est_2 = Result_2[1,3]
mu12_est_2 = Result_2[2,3]
mu13_est_2 = Result_2[3,3]
## Sd of Diseased and Non-diseased
sd01_est_2 = Result_2[1,4]
sd02_est_2 = Result_2[2,4]
sd03_est_2 = Result_2[3,4]

sd11_est_2 = Result_2[1,5]
sd12_est_2 = Result_2[2,5]
sd13_est_2 = Result_2[3,5]

## Correlation Coefficients in Diseased and Non-diseased
rho0_12_est_2 = Result_2[1,6]
rho0_13_est_2 = Result_2[2,6]
rho0_23_est_2 = Result_2[3,6]

rho1_12_est_2 = Result_2[1,7]
rho1_13_est_2 = Result_2[2,7]
rho1_23_est_2 = Result_2[3,7]

### ROC_1_Est
ROC.Fun_1_2 <- function(u)
{
  return(pnorm((mu11_est_2 - mu01_est_2 + sd01_est_2 * qnorm(u))/sd11_est_2))
}
ROC.est_1_2 <- ROC.Fun_1_2(aa)
### ROC_2_Est
ROC.Fun_2_2 <- function(u)
{
  return(pnorm((mu12_est_2 - mu02_est_2 + sd02_est_2 * qnorm(u))/sd12_est_2))
}
ROC.est_2_2 <- ROC.Fun_2_2(aa)
### ROC_3_Est
ROC.Fun_3_2 <- function(u)
{
  return(pnorm((mu13_est_2 - mu03_est_2 + sd03_est_2 * qnorm(u))/sd13_est_2))
}
ROC.est_3_2 <- ROC.Fun_3_2(aa)

#########################################################################################################################################
## For K = 5
Result_5 <- read.csv("./results/Boot/TTB_scaleSep_Seed.boot=1000_K=5_0.96_0.952_2023-06-28.csv", header = T)
Result_5
## Mean of Diseased and Non-diseased
mu01_est_5 = Result_5[1,2]
mu02_est_5 = Result_5[2,2]
mu03_est_5 = Result_5[3,2]

mu11_est_5 = Result_5[1,3]
mu12_est_5 = Result_5[2,3]
mu13_est_5 = Result_5[3,3]
## Sd of Diseased and Non-diseased
sd01_est_5 = Result_5[1,4]
sd02_est_5 = Result_5[2,4]
sd03_est_5 = Result_5[3,4]

sd11_est_5 = Result_5[1,5]
sd12_est_5 = Result_5[2,5]
sd13_est_5 = Result_5[3,5]

## Correlation Coefficients in Diseased and Non-diseased
rho0_12_est_5 = Result_5[1,6]
rho0_13_est_5 = Result_5[2,6]
rho0_23_est_5 = Result_5[3,6]

rho1_12_est_5 = Result_5[1,7]
rho1_13_est_5 = Result_5[2,7]
rho1_23_est_5 = Result_5[3,7]

### ROC_1_Est
ROC.Fun_1_5 <- function(u)
{
  return(pnorm((mu11_est_5 - mu01_est_5 + sd01_est_5 * qnorm(u))/sd11_est_5))
}
ROC.est_1_5 <- ROC.Fun_1_5(aa)
### ROC_2_Est
ROC.Fun_2_5 <- function(u)
{
  return(pnorm((mu12_est_5 - mu02_est_5 + sd02_est_5 * qnorm(u))/sd12_est_5))
}
ROC.est_2_5 <- ROC.Fun_2_5(aa)
### ROC_3_Est
ROC.Fun_3_5 <- function(u)
{
  return(pnorm((mu13_est_5 - mu03_est_5 + sd03_est_5 * qnorm(u))/sd13_est_5))
}
ROC.est_3_5 <- ROC.Fun_3_5(aa)

########################## Plot #################################
### Biomarker 1
plot(aa, ROC.true_1, type = "l", col = "red",  lty = "solid", lwd = 2, xlab = "u", ylab = TeX(r"($\widehat{ROC}(u)$)"), mgp = c(1.5,0.5,0), cex.lab = 2, cex.axis = 1.5)
lines(aa, sort(1 - ROC.est_1_1), type = "l", col = "purple",  lty = "dotted", lwd = 2)
abline(0, 1, col = "red", lwd = 1, lty = 1)
lines(aa, sort(1 - ROC.est_1_2), type = "l", col = "blue",  lty = "longdash", lwd = 2)
lines(aa, sort(1 - ROC.est_1_5), type = "l", col = "deeppink",  lty = "dotdash", lwd = 2)
legend(0.55, 0.65, legend=c("K = 1", "K = 2", "K = 5"), col=c("purple", "blue", "deeppink"), lty=c(3,5,4), bty = "n", lwd = 3, cex=1.8, y.intersp = 0.35)
### Biomarker 2
plot(aa, ROC.true_2, type = "l", col = "red",  lty = "solid", lwd = 2, xlab = "u", ylab = TeX(r"($\widehat{ROC}(u)$)"), mgp = c(1.5,0.5,0), cex.lab = 2, cex.axis = 1.5)
lines(aa, ROC.est_2_1, type = "l", col = "purple",  lty = "dotted", lwd = 2)
abline(0, 1, col = "red", lwd = 1, lty = 1)
lines(aa, ROC.est_2_2, type = "l", col = "blue",  lty = "longdash", lwd = 2)
lines(aa, ROC.est_2_5, type = "l", col = "deeppink",  lty = "dotdash", lwd = 2)
legend(0.55, 0.65, legend=c("K = 1", "K = 2", "K = 5"), col=c("purple", "blue", "deeppink"), lty=c(3,5,4), bty = "n", lwd = 3, cex=1.8, y.intersp = 0.35)
### Biomarker 3
plot(aa, ROC.true_3, type = "l", col = "red",  lty = "solid", lwd = 2, xlab = "u", ylab = TeX(r"($\widehat{ROC}(u)$)"), mgp = c(1.5,0.5,0), cex.lab = 2, cex.axis = 1.5)
lines(aa, ROC.est_3_1, type = "l", col = "purple",  lty = "dotted", lwd = 2)
abline(0, 1, col = "red", lwd = 1, lty = 1)
lines(aa, ROC.est_3_2, type = "l", col = "blue",  lty = "longdash", lwd = 2)
lines(aa, ROC.est_3_5, type = "l", col = "deeppink",  lty = "dotdash", lwd = 2)
legend(0.55, 0.65, legend=c("K = 1", "K = 2", "K = 5"), col=c("purple", "blue", "deeppink"), lty=c(3,5,4), bty = "n", lwd = 3, cex=1.8, y.intersp = 0.35)

###### AUC Combination Estimation
mu1_est_1 = c(mu11_est_1, mu12_est_1, mu13_est_1)
mu0_est_1 = c(mu01_est_1, mu02_est_1, mu03_est_1)
mu_est_1 <- as.matrix(mu1_est_1 - mu0_est_1)
mu_est_1
sigma_0_est_1 = matrix(c(sd01_est_1^2, rho0_12_est_1 * sd01_est_1 * sd02_est_1, rho0_13_est_1 * sd01_est_1 * sd03_est_1,
                         rho0_12_est_1 * sd02_est_1 * sd01_est_1, sd02_est_1^2, rho0_23_est_1 * sd02_est_1 * sd03_est_1,
                         rho0_13_est_1 * sd01_est_1 * sd03_est_1, rho0_23_est_1 * sd02_est_1 * sd03_est_1, sd03_est_1^2), 3, 3)
sigma_0_est_1
sigma_1_est_1 = matrix(c(sd11_est_1^2, rho1_12_est_1 * sd11_est_1 * sd12_est_1, rho1_13_est_1 * sd11_est_1 * sd13_est_1,
                         rho1_12_est_1 * sd12_est_1 * sd11_est_1, sd12_est_1^2, rho1_23_est_1 * sd12_est_1 * sd13_est_1,
                         rho1_13_est_1 * sd11_est_1 * sd13_est_1, rho1_23_est_1 * sd12_est_1 * sd13_est_1, sd13_est_1^2), 3, 3)
sigma_1_est_1

c_1 <- t(mu_est_1) %*% solve(sigma_1_est_1 + sigma_0_est_1)
c_1

linear_alpha_1 <- (1/c_1[,1]) * t(mu_est_1) %*% solve(sigma_1_est_1 + sigma_0_est_1)
linear_alpha_1

new_mu_1 <- c(linear_alpha_1 %*% mu_est_1)
new_mu_1
new_sigma0_1 <- c(linear_alpha_1 %*% sigma_0_est_1 %*% t(linear_alpha_1))
new_sigma0_1
new_sigma1_1 <- c(linear_alpha_1 %*% sigma_1_est_1 %*% t(linear_alpha_1))
new_sigma1_1

AUC_try_1 <- pnorm(sqrt(t(new_mu_1) %*% solve(new_sigma0_1 + new_sigma1_1) %*% new_mu_1))
AUC_try_1
ROC.Curve.fun_1 <- function(u)
{
  return(pnorm((new_mu_1 + sqrt(new_sigma0_1) * qnorm(u))/sqrt(new_sigma1_1)))
}

ROC.curve_1 <- ROC.Curve.fun_1(aa)
#########################################################################
mu1_est_2 = c(mu11_est_2, mu12_est_2, mu13_est_2)
mu0_est_2 = c(mu01_est_2, mu02_est_2, mu03_est_2)
mu_est_2 <- as.matrix(mu1_est_2 - mu0_est_2)
mu_est_2
sigma_0_est_2 = matrix(c(sd01_est_2^2, rho0_12_est_2 * sd01_est_2 * sd02_est_2, rho0_13_est_2 * sd01_est_2 * sd03_est_2,
                         rho0_12_est_2 * sd02_est_2 * sd01_est_2, sd02_est_2^2, rho0_23_est_2 * sd02_est_2 * sd03_est_2,
                         rho0_13_est_2 * sd01_est_2 * sd03_est_2, rho0_23_est_2 * sd02_est_2 * sd03_est_2, sd03_est_2^2), 3, 3)
sigma_0_est_2
sigma_1_est_2 = matrix(c(sd11_est_2^2, rho1_12_est_2 * sd11_est_2 * sd12_est_2, rho1_13_est_2 * sd11_est_2 * sd13_est_2,
                         rho1_12_est_2 * sd12_est_2 * sd11_est_2, sd12_est_2^2, rho1_23_est_2 * sd12_est_2 * sd13_est_2,
                         rho1_13_est_2 * sd11_est_2 * sd13_est_2, rho1_23_est_2 * sd12_est_2 * sd13_est_2, sd13_est_2^2), 3, 3)
sigma_1_est_2

c_2 <- t(mu_est_2) %*% solve(sigma_1_est_2 + sigma_0_est_2)
c_2

linear_alpha_2 <- (1/c_2[,1]) * t(mu_est_2) %*% solve(sigma_1_est_2 + sigma_0_est_2)
linear_alpha_2

new_mu_2 <- c(linear_alpha_2 %*% mu_est_2)
new_mu_2
new_sigma0_2 <- c(linear_alpha_2 %*% sigma_0_est_2 %*% t(linear_alpha_2))
new_sigma0_2
new_sigma1_2 <- c(linear_alpha_2 %*% sigma_1_est_2 %*% t(linear_alpha_2))
new_sigma1_2

AUC_try_2 <- pnorm(sqrt(t(new_mu_2) %*% solve(new_sigma0_2 + new_sigma1_2) %*% new_mu_2))
AUC_try_2
ROC.Curve.fun_2 <- function(u)
{
  return(pnorm((new_mu_2 + sqrt(new_sigma0_2) * qnorm(u))/sqrt(new_sigma1_2)))
}

ROC.curve_2 <- ROC.Curve.fun_2(aa)
#########################################################################
mu1_est_5 = c(mu11_est_5, mu12_est_5, mu13_est_5)
mu0_est_5 = c(mu01_est_5, mu02_est_5, mu03_est_5)
mu_est_5 <- as.matrix(mu1_est_5 - mu0_est_5)
mu_est_5
sigma_0_est_5 = matrix(c(sd01_est_5^2, rho0_12_est_5 * sd01_est_5 * sd02_est_5, rho0_13_est_5 * sd01_est_5 * sd03_est_5,
                         rho0_12_est_5 * sd02_est_5 * sd01_est_5, sd02_est_5^2, rho0_23_est_5 * sd02_est_5 * sd03_est_5,
                         rho0_13_est_5 * sd01_est_5 * sd03_est_5, rho0_23_est_5 * sd02_est_5 * sd03_est_5, sd03_est_5^2), 3, 3)
sigma_0_est_5
sigma_1_est_5 = matrix(c(sd11_est_5^2, rho1_12_est_5 * sd11_est_5 * sd12_est_5, rho1_13_est_5 * sd11_est_5 * sd13_est_5,
                         rho1_12_est_5 * sd12_est_5 * sd11_est_5, sd12_est_5^2, rho1_23_est_5 * sd12_est_5 * sd13_est_5,
                         rho1_13_est_5 * sd11_est_5 * sd13_est_5, rho1_23_est_5 * sd12_est_5 * sd13_est_5, sd13_est_5^2), 3, 3)
sigma_1_est_5

c_5 <- t(mu_est_5) %*% solve(sigma_1_est_5 + sigma_0_est_5)
c_5

linear_alpha_5 <- (1/c_5[,1]) * t(mu_est_5) %*% solve(sigma_1_est_5 + sigma_0_est_5)
linear_alpha_5

new_mu_5 <- c(linear_alpha_5 %*% mu_est_5)
new_mu_5
new_sigma0_5 <- c(linear_alpha_5 %*% sigma_0_est_5 %*% t(linear_alpha_5))
new_sigma0_5
new_sigma1_5 <- c(linear_alpha_5 %*% sigma_1_est_5 %*% t(linear_alpha_5))
new_sigma1_5

AUC_try_5 <- pnorm(sqrt(t(new_mu_5) %*% solve(new_sigma0_5 + new_sigma1_5) %*% new_mu_5))
AUC_try_5
ROC.Curve.fun_5 <- function(u)
{
  return(pnorm((new_mu_5 + sqrt(new_sigma0_5) * qnorm(u))/sqrt(new_sigma1_5)))
}

ROC.curve_5 <- ROC.Curve.fun_5(aa)

# plot(aa, ROC.curve_1, type = "l", col = "purple",  lty = "dotted", lwd = 2, xlab = "u", ylab = TeX(r"($\widehat{ROC}(u)$)"), mgp = c(1.5,0.5,0), cex.lab = 2, cex.axis = 1.5)
plot(aa, ROC.curve_1, type = "l", col = "purple",  lty = "dotted", lwd = 2, xlab = "", ylab = "", mgp = c(1.5,0.5,0), cex.lab = 2, cex.axis = 1.5)
lines(aa, ROC.curve_2, type = "l", col = "blue",  lty = "longdash", lwd = 2)
lines(aa, ROC.curve_5, type = "l", col = "deeppink",  lty = "dotdash", lwd = 2)
#abline(0, 1, col = "red", lwd = 1, lty = 1)
#lines(aa, ROC_true, type = "l", col = "red",  lty = "solid", lwd = 2)
legend(0.68, 0.45, legend=c("J = 1", "J = 2", "J = 5"), col=c("purple", "blue", "deeppink"), lty=c(3,5,4), bty = "n", lwd = 3, cex = 1.8, y.intersp = 1.8)


# c(AUC_true, AUC_try_1, AUC_try_2, AUC_try_5)
# ##########################################################
# ### Biomarker 1
# AUC_1_1 <- pnorm(sqrt(t(mu11_est_1-mu01_est_1) %*% solve(sd01_est_1^2 + sd11_est_1^2) %*% (mu11_est_1-mu01_est_1)))
# AUC_1_1
# AUC_1_2 <- pnorm(sqrt(t(mu11_est_2-mu01_est_2) %*% solve(sd01_est_2^2 + sd11_est_2^2) %*% (mu11_est_2-mu01_est_2)))
# AUC_1_2
# AUC_1_5 <- pnorm(sqrt(t(mu11_est_5-mu01_est_5) %*% solve(sd01_est_5^2 + sd11_est_5^2) %*% (mu11_est_5-mu01_est_5)))
# AUC_1_5
# c(AUC_1, AUC_1_1, AUC_1_2, AUC_1_5)
# ### Biomarker 2
# AUC_2_1 <- pnorm(sqrt(t(mu12_est_1-mu02_est_1) %*% solve(sd02_est_1^2 + sd12_est_1^2) %*% (mu12_est_1-mu02_est_1)))
# AUC_2_1
# AUC_2_2 <- pnorm(sqrt(t(mu12_est_2-mu02_est_2) %*% solve(sd02_est_2^2 + sd12_est_2^2) %*% (mu12_est_2-mu02_est_2)))
# AUC_2_2
# AUC_2_5 <- pnorm(sqrt(t(mu12_est_5-mu02_est_5) %*% solve(sd02_est_5^2 + sd12_est_5^2) %*% (mu12_est_5-mu02_est_5)))
# AUC_2_5
# c(AUC_2, AUC_2_1, AUC_2_2, AUC_2_5)
# ### Biomarker 3
# AUC_3_1 <- pnorm(sqrt(t(mu13_est_1-mu03_est_1) %*% solve(sd03_est_1^2 + sd13_est_1^2) %*% (mu13_est_1-mu03_est_1)))
# AUC_3_1
# AUC_3_2 <- pnorm(sqrt(t(mu13_est_2-mu03_est_2) %*% solve(sd03_est_2^2 + sd13_est_2^2) %*% (mu13_est_2-mu03_est_2)))
# AUC_3_2
# AUC_3_5 <- pnorm(sqrt(t(mu13_est_5-mu03_est_5) %*% solve(sd03_est_5^2 + sd13_est_5^2) %*% (mu13_est_5-mu03_est_5)))
# AUC_3_5
# c(AUC_3, AUC_3_1, AUC_3_2, AUC_3_5)
