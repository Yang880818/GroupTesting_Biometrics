rm(list = ls())
library(foreach)
library(doParallel)
library(combinat)
library(bbmle)
library(Rfast)
library(MASS) 
library(moments)
#library(latex2exp)

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
plot(aa, ROC.true_1, type = "l", xlab = "u", ylab = "ROC", mgp = c(1.5,0.5,0), cex.lab = 2, lwd = 2, cex.axis = 1.5)
abline(0, 1, col = "red", lwd = 1, lty = 1)

## ROC_2
ROC.Fun_2 <- function(u)
{
  return(pnorm((mu12 - mu02 + sd02 * qnorm(u))/sd12))
}
ROC.true_2 <- ROC.Fun_2(aa)
plot(aa, sort(1 - ROC.true_2), type = "l", xlab = "u", ylab = "ROC", mgp = c(1.5,0.5,0), cex.lab = 2, lwd = 2, cex.axis = 1.5)
abline(0, 1, col = "red", lwd = 1, lty = 1)

## ROC_3
ROC.Fun_3 <- function(u)
{
  return(pnorm((mu13 - mu03 + sd03 * qnorm(u))/sd13))
}
ROC.true_3 <- ROC.Fun_3(aa)
plot(aa, sort(1 - ROC.true_3), type = "l", xlab = "u", ylab = "ROC", mgp = c(1.5,0.5,0), cex.lab = 2, lwd = 2, cex.axis = 1.5)
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
plot(aa, ROC_true, type = "l", xlab = "u", ylab = "ROC", mgp = c(1.5,0.5,0), cex.lab = 2, lwd = 2, cex.axis = 1.5)
abline(0, 1, col = "red", lwd = 1, lty = 1)


######################################################################################################################
Bio0_comb = comb_n(Bio0[1, ], 2)
Bio0_local <- as.matrix(apply(comb_n(Bio0[1, ], 2), 2, function(x){return(which(Bio0[1, ] %in% x))})) 
Bio0_local

Bio1_comb = comb_n(Bio1[1, ], 2)
Bio1_local <- as.matrix(apply(comb_n(Bio1[1, ], 2), 2, function(x){return(which(Bio1[1, ] %in% x))})) 
Bio1_local

mu0_comb = comb_n(mu0, 2)
mu0_local <- as.matrix(apply(comb_n(mu0, 2), 2, function(x){return(which(mu0 %in% x))})) 
mu0_local

mu1_comb = comb_n(mu1, 2)
mu1_local <- as.matrix(apply(comb_n(mu1, 2), 2, function(x){return(which(mu1 %in% x))}))
mu1_local

sd0_comb = comb_n(sd0, 2)
sd0_local <- as.matrix(apply(comb_n(sd0, 2), 2, function(x){return(which(sd0 %in% x))})) 
sd0_local

sd1_comb = comb_n(sd1, 2)
sd1_local <- as.matrix(apply(comb_n(sd1, 2), 2, function(x){return(which(sd1 %in% x))})) 
sd1_local

rho_leg = length(rho0)





p = mean(Data_resample_final$Disease_status == 1)
p
#N = dim(Data_resample_final)[1]
K = 1
n = N/K
lambda = 0.02
pi0 = Pi0
pi1 = Pi1
alpha = 0.05


EstPMeanSD  <- rep(NA, 11)
p_est_final <- rep(NA, 1)
AUC.value   <- rep(NA, 1)
mu0.final  <- matrix(NA, 1, length(mu0))
mu1.final  <- matrix(NA, 1, length(mu1))
sd0.final  <- matrix(NA, 1, length(sd0))
sd1.final  <- matrix(NA, 1, length(sd1))
rho0.final <- matrix(NA, 1, length(rho0))
rho1.final <- matrix(NA, 1, length(rho1))

mu0.mat <- matrix(0, rho_leg, length(mu0))
mu0.mat
mu1.mat <- matrix(0, rho_leg, length(mu1))
mu1.mat
sd0.mat <- matrix(0, rho_leg, length(sd0))
sd0.mat
sd1.mat <- matrix(0, rho_leg, length(sd1))
sd1.mat
rho0.mat <- rep(0, rho_leg)
rho0.mat
rho1.mat <- rep(0, rho_leg)
rho1.mat
p.vec    <- rep(0, rho_leg)
p.vec

for(l in 1 : rho_leg)
{
  Dmat <- t(matrix(Data_resample_final$Disease_status, nrow = K, ncol = n))
  XYvec <- matrix(NA, N, 2)
  
  mu_1 <- c(mu1_comb[ ,l])
  sigma_1 <- matrix(c(sd1_comb[1, l]^2, rho1[l] * sd1_comb[1, l] * sd1_comb[2, l], rho1[l] * sd1_comb[1, l] * sd1_comb[2, l], sd1_comb[2, l]^2), 2, 2)
  XYvec[Dmat==1, ] <- Bio1[, Bio1_local[,l]]
  
  mu_0 <- c(mu0_comb[ ,l])
  sigma_0 <- matrix(c(sd0_comb[1, l]^2, rho0[l] * sd0_comb[1, l] * sd0_comb[2, l], rho0[l] * sd0_comb[1, l] * sd0_comb[2, l], sd0_comb[2, l]^2), 2, 2)
  XYvec[Dmat==0, ] <- Bio0[, Bio0_local[,l]]
  
  
  # ########### GroupM ############
  dVal <- rowSums(Dmat) ###GroupM_true !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ## The case
  groupD <- ifelse(dVal == 0, yes = 0, no = 1)
  print(sum(groupD == 1))
  
  ### For the negative group
  D0Loc <- which(groupD == 0)
  D0M1.num <- sum(rbinom(sum(groupD == 0), size = 1, prob = 1 - pi0) == 1)
  #D0M1.num
  groupM <- groupD
  D0M1Loc <- sample(D0Loc, D0M1.num, replace = F)
  groupM[D0M1Loc] <- 1
  #groupM
  
  SenFun <- function(pi1, K, dnum, lambda)
  {
    res = pi1 * dnum /(dnum + lambda * (K - dnum))
    return(res)
  }
  
  ## For the positive group
  pi1.vec <- apply(matrix(1:K, ncol = 1), 1, SenFun, pi1 = pi1, K = K, lambda = lambda)
  for(dd in 1:K)
  {
    tempLoc <- which(dVal == dd)
    if(length(tempLoc) > 0)
    {
      tempNum <- sum(dVal == dd)
      D1M0.num <- sum(rbinom(tempNum, size = 1, prob = pi1.vec[dd]) == 0)
      if(D1M0.num > 0)
      {
        D1M0Loc <- sample(tempLoc, D1M0.num, replace = F)
        groupM[D1M0Loc] <- 0
      }
    }
  }
  #groupM
  print(sum(groupM == 1))
  
  pi1_lik = apply(matrix(1:K, ncol = 1), 1, SenFun, pi1 = pi1, K = K, lambda = lambda)
  pi1_lik.mat = as.matrix(pi1_lik)
  pi1_lik.mat
  
  ############ binary matrix ############
  a.mat <- rbind()
  pi1_lik.matrix <- rbind()
  for(d in 1:K)
  {
    for(i in 1 : choose(K, d))
    {
      a <- rep(0, K)
      a[as.matrix(comb_n(K, d))[ ,i]] <- 1
      pi1_lik.matrix <- rbind(pi1_lik.matrix, apply(matrix(d, ncol = 1), 1, SenFun, pi1 = pi1, K = K, lambda = lambda))
      a.mat <- rbind(a.mat,a)
    }
  }
  aa.mat <- matrix(a.mat, dim(a.mat)[1], dim(a.mat)[2])
  pi1_lik.matrix
  ############ Key Matrix #############
  aa.mat
  
  ####################################
  
  LoglikPMeanSD_est <- function(theta, pi0, pi1, N, n, K, XYvec, lambda, groupM, aa.mat, pi1_lik.matrix)
  {
    mu01_est <- theta[1]
    mu02_est <- theta[2]
    mu11_est <- theta[3]
    mu12_est <- theta[4]
    sd01_est <- theta[5]
    sd02_est <- theta[6]
    sd11_est <- theta[7]
    sd12_est <- theta[8]
    rho0_est <- theta[9]
    rho1_est <- theta[10]
    p_est    <- theta[11]
    
    LoglikPMeanSD = rep(NA, n)
    
    f0.vec  <- rep(NA, N)
    f0.vec  <- (1/(2 * pi * sd01_est * sd02_est * sqrt(1 - rho0_est^2))) * exp(-1/(2 * (1 - rho0_est^2)) * (((XYvec[ ,1] - mu01_est)/sd01_est)^2 - 2 * rho0_est * ((XYvec[ ,1] - mu01_est)/sd01_est) * ((XYvec[ ,2] - mu02_est)/sd02_est) + ((XYvec[ ,2] - mu02_est)/sd02_est)^2))
    f0.mat  <- t(matrix(f0.vec, nrow = K, ncol = n))
    
    f1.vec  <- rep(NA, N)
    f1.vec  <- (1/(2 * pi * sd11_est * sd12_est * sqrt(1 - rho1_est^2))) * exp(-1/(2 * (1 - rho1_est^2)) * (((XYvec[ ,1] - mu11_est)/sd11_est)^2 - 2 * rho1_est * ((XYvec[ ,1] - mu11_est)/sd11_est) * ((XYvec[ ,2] - mu12_est)/sd12_est) + ((XYvec[ ,2] - mu12_est)/sd12_est)^2))
    f1.mat  <- t(matrix(f1.vec, nrow = K, ncol = n))
    
    mat.time <- dim(aa.mat)[1]
    f0.mat_final = matrix(rep(f0.mat, mat.time), n, K * mat.time)
    f1.mat_final = matrix(rep(f1.mat, mat.time), n, K * mat.time)
    
    aa.mat_new = t(replicate(n, as.vector(t(aa.mat))))
    
    if(K==1)
    {
      f0.matrix = f0.mat_final * t(1 - aa.mat_new)
      f1.matrix = f1.mat_final * t(aa.mat_new)
      
      f.matrix = f0.matrix + f1.matrix
    }else
    {
      f0.matrix = f0.mat_final * (1 - aa.mat_new)
      f1.matrix = f1.mat_final * aa.mat_new
      
      f.matrix = f0.matrix + f1.matrix
    }
    
    
    f.matrix_final <- matrix(NA, n, mat.time)
    for(j in 1 : mat.time)
    {
      for(i in 1 : n)
        f.matrix_final[i, j] <- prod(f.matrix[i , (1 + K * (j - 1)) : (K + K * (j - 1))])
    }
    
    
    ############### pi1 * p^d * (1 - p)^(K - d) ###############
    pd <- p_est^(apply(aa.mat, 1, sum)) * (1 - p_est)^(K-apply(aa.mat, 1, sum))
    pd.mat <- as.matrix(pd)
    pd.mat
    
    LoglikPMeanSD_final <- 0
    for(j in 1 : n)
    {
      
      if(groupM[j] == 0)
      {
        tempA <- 0
        
        tempA = tempA + f.matrix_final[j, ] %*% ((1 - pi1_lik.matrix) * pd.mat)
        
        LoglikPMeanSD[j] <- log(pi0 * (1 - p_est)^K * prod(f0.mat[j, ]) + tempA)
        
      }
      else
      {
        tempE <- 0
        
        tempE = tempE + f.matrix_final[j, ] %*% (pi1_lik.matrix * pd.mat)
        
        LoglikPMeanSD[j] <- log((1 - pi0) * (1 - p_est)^K * prod(f0.mat[j, ]) + tempE)
        
      }
      
      LoglikPMeanSD_final =  LoglikPMeanSD_final + LoglikPMeanSD[j]
    }
    
    res = -LoglikPMeanSD_final
    return(res)
  }
  
  
  EstPMeanSD  <- as.vector(optim(c(mu0_comb[1,l], mu0_comb[2,l], mu1_comb[1,l], mu1_comb[2,l], sd0_comb[1,l], sd0_comb[2,l], sd1_comb[1,l], sd1_comb[2,l], rho0[l], rho1[l], p),
                                 LoglikPMeanSD_est, pi0 = pi0, pi1 = pi1, N = N, n = n, K = K, XYvec = XYvec, lambda = lambda,
                                 groupM = groupM, aa.mat = aa.mat, pi1_lik.matrix = pi1_lik.matrix)$par)
  
  # 
  mu01_final <- EstPMeanSD[1]
  mu02_final <- EstPMeanSD[2]
  mu11_final <- EstPMeanSD[3]
  mu12_final <- EstPMeanSD[4]
  sd01_final <- EstPMeanSD[5]
  sd02_final <- EstPMeanSD[6]
  sd11_final <- EstPMeanSD[7]
  sd12_final <- EstPMeanSD[8]
  rho0_final <- EstPMeanSD[9]
  rho1_final <- EstPMeanSD[10]
  p_est_final <- EstPMeanSD[11]
  
  print(l)
  
  mu0.mat[l, (mu0_local[ ,l])] <- mu0.mat[l, (mu0_local[ ,l])] + c(mu01_final, mu02_final)
  mu1.mat[l, (mu1_local[ ,l])] <- mu1.mat[l, (mu1_local[ ,l])] + c(mu11_final, mu12_final)
  sd0.mat[l, (sd0_local[ ,l])] <- sd0.mat[l, (sd0_local[ ,l])] + c(sd01_final, sd02_final)
  sd1.mat[l, (sd1_local[ ,l])] <- sd1.mat[l, (sd1_local[ ,l])] + c(sd11_final, sd12_final)
  rho0.mat[l] <- rho0.mat[l] + rho0_final
  rho1.mat[l] <- rho1.mat[l] + rho1_final
  p.vec[l] <- p.vec[l] + p_est_final
}

mu0.final <- apply(mu0.mat, 2, sum)/(length(mu0) - 1)
mu1.final <- apply(mu1.mat, 2, sum)/(length(mu1) - 1)
sd0.final <- apply(sd0.mat, 2, sum)/(length(sd0) - 1)
sd1.final <- apply(sd1.mat, 2, sum)/(length(sd1) - 1)
rho0.final <- rho0.mat
rho1.final <- rho1.mat


mu_final <- mu1.final - mu0.final
mu_final
sigma_case_final <- sigma_control_final <- matrix(NA, 3, 3)
sigma_case_final <- matrix(c(sd1.final[1]^2, rho1.final[1] * sd1.final[1] * sd1.final[2], rho1.final[2] * sd1.final[1] * sd1.final[3], 
                             rho1.final[1] * sd1.final[2] * sd1.final[1], sd1.final[2]^2, rho1.final[3] * sd1.final[2] * sd1.final[3],
                             rho1.final[2] * sd1.final[3] * sd1.final[1], rho1.final[3] * sd1.final[3] * sd1.final[2], sd1.final[3]^2), 3, 3)
sigma_case_final
sigma_control_final <- matrix(c(sd0.final[1]^2, rho0.final[1] * sd0.final[1] * sd0.final[2], rho0.final[2] * sd0.final[1] * sd0.final[3], 
                                rho0.final[1] * sd0.final[2] * sd0.final[1], sd0.final[2]^2, rho0.final[3] * sd0.final[2] * sd0.final[3],
                                rho0.final[2] * sd0.final[3] * sd0.final[1], rho0.final[3] * sd0.final[3] * sd0.final[2], sd0.final[3]^2), 3, 3)
sigma_control_final
AUC.value <- pnorm(sqrt(t(mu_final) %*% solve(sigma_control_final + sigma_case_final) %*% mu_final))
AUC.value

mu_final
sigma_control_final
sigma_case_final

mu_true
sigma_0_true
sigma_1_true
AUC_true
mu0.final - mu0
mu1.final - mu1
sd0.final - sd0
sd1.final - sd1
rho0.final - rho0
rho1.final - rho1

EstPMeanSD


p.final <- mean(na.omit(p.vec))
p.bias <- p.final - p
p.MSE  <- mean((na.omit(p.vec) - p)^2)

final_results =
  as.vector(
    list("mu0.final" = mu0.final, "mu1.final" = mu1.final, "sd0.final" = sd0.final, "sd1.final" = sd1.final, 
         "rho0.final" = rho0.final, "rho1.final" = rho1.final, "AUC_true" = AUC_true, "AUC.value" = AUC.value, "p" = p, "p.final" = p.final
    )
  )

write.csv(final_results, file = paste0("./results/TTB_scaleSep_K=", K, "_", pi0, "_", pi1, "_", CurrentData, ".csv"))





seed.num.boot = SEED.NUM.BOOT

#Data_resample_final
UseCores <- 100 # 
detectCores() - 1
#UseCores
cl <- makeCluster(UseCores)
registerDoParallel(cl) 

EstPMeanSD.boot  <- matrix(NA, seed.num.boot, 11)

p_est_final.boot <- rep(NA, seed.num.boot)
AUC.value.boot   <- rep(NA, seed.num.boot)
AUC_1.value.boot <- rep(NA, seed.num.boot)
AUC_2.value.boot <- rep(NA, seed.num.boot)
AUC_3.value.boot <- rep(NA, seed.num.boot)
mu0.final.boot  <- matrix(NA, seed.num.boot, length(mu0))
mu1.final.boot  <- matrix(NA, seed.num.boot, length(mu1))
sd0.final.boot  <- matrix(NA, seed.num.boot, length(sd0))
sd1.final.boot  <- matrix(NA, seed.num.boot, length(sd1))
rho0.final.boot <- matrix(NA, seed.num.boot, length(rho0))
rho1.final.boot <- matrix(NA, seed.num.boot, length(rho1))
p.final.boot <- rep(NA, seed.num.boot)

### Bio_1
Est_Bio1.boot <- matrix(NA, seed.num.boot, 5)
mu0_Bio1.boot <- rep(NA, seed.num.boot)
mu1_Bio1.boot <- rep(NA, seed.num.boot)
sd0_Bio1.boot <- rep(NA, seed.num.boot)
sd1_Bio1.boot <- rep(NA, seed.num.boot)
p_Bio1.boot   <- rep(NA, seed.num.boot)
AUC_1.boot <- rep(NA, seed.num.boot)
### Bio_2
Est_Bio2.boot <- matrix(NA, seed.num.boot, 5)
mu0_Bio2.boot <- rep(NA, seed.num.boot)
mu1_Bio2.boot <- rep(NA, seed.num.boot)
sd0_Bio2.boot <- rep(NA, seed.num.boot)
sd1_Bio2.boot <- rep(NA, seed.num.boot)
p_Bio2.boot   <- rep(NA, seed.num.boot)
AUC_2.boot <- rep(NA, seed.num.boot)
### Bio_3
Est_Bio3.boot <- matrix(NA, seed.num.boot, 5)
mu0_Bio3.boot <- rep(NA, seed.num.boot)
mu1_Bio3.boot <- rep(NA, seed.num.boot)
sd0_Bio3.boot <- rep(NA, seed.num.boot)
sd1_Bio3.boot <- rep(NA, seed.num.boot)
p_Bio3.boot   <- rep(NA, seed.num.boot)
AUC_3.boot <- rep(NA, seed.num.boot)


Main_Fun <- foreach(j = 1 : seed.num.boot, .packages= c('MASS', 'Rfast', 'bbmle', 'combinat')) %dopar%
  {
    set.seed(j)
    
    Xgroup.boot.Loc <- sample(c(1:N), replace = TRUE)
    Xgroup.boot.Loc
    
    Data_resample_new <- Data_resample_final[Xgroup.boot.Loc, ]
    
    
    Xgroup.boot <- Data_resample_new$Disease_status
    
    n0 <- sum(Xgroup.boot == 0)
    n0
    n1 <- sum(Xgroup.boot == 1)
    n1
    groupPos.boot <- which(Xgroup.boot == 1)
    groupPos.boot
    groupNeg.boot <- which(Xgroup.boot == 0)
    groupNeg.boot
    
    Dmat.boot <- t(matrix(Data_resample_new$Disease_status, nrow = K, ncol = n))
    
    XgroupPos.boot <- Data_resample_new[groupPos.boot,]
    #XgroupPos.boot
    XgroupNeg.boot <- Data_resample_new[groupNeg.boot, ]
    #XgroupNeg.boot
    #Data_resample_final.boot <- rbind(XgroupPos.boot, XgroupNeg.boot)
    
    p.boot = mean(Data_resample_new$Disease_status == 1)
    p.boot
    #sum(Data_resample_final.boot$Disease_status[which(Data_resample_final.boot $Disease_status == 1)])
    
    Test_Biomarker_1_0_new.boot <- Data_resample_new$Biomarker_1[which(Xgroup.boot == 0)]
    Test_Biomarker_1_1_new.boot <- Data_resample_new$Biomarker_1[which(Xgroup.boot == 1)]
    Test_Biomarker_2_0_new.boot <- Data_resample_new$Biomarker_2[which(Xgroup.boot == 0)]
    Test_Biomarker_2_1_new.boot <- Data_resample_new$Biomarker_2[which(Xgroup.boot == 1)]
    Test_Biomarker_3_0_new.boot <- Data_resample_new$Biomarker_3[which(Xgroup.boot == 0)]
    Test_Biomarker_3_1_new.boot <- Data_resample_new$Biomarker_3[which(Xgroup.boot == 1)]
    
    Bio0.boot <- cbind(Test_Biomarker_1_0_new.boot, Test_Biomarker_2_0_new.boot, Test_Biomarker_3_0_new.boot)
    Bio1.boot <- cbind(Test_Biomarker_1_1_new.boot, Test_Biomarker_2_1_new.boot, Test_Biomarker_3_1_new.boot)
    
    
    # Mean of Non-diseased
    mu01.boot = mean(Test_Biomarker_1_0_new.boot)
    mu02.boot = mean(Test_Biomarker_2_0_new.boot)
    mu03.boot = mean(Test_Biomarker_3_0_new.boot)
    ## Mean of Diseased
    mu11.boot = mean(Test_Biomarker_1_1_new.boot)
    mu12.boot = mean(Test_Biomarker_2_1_new.boot)
    mu13.boot = mean(Test_Biomarker_3_1_new.boot)
    
    ## Sd of Non-diseased
    sd01.boot = sd(Test_Biomarker_1_0_new.boot)
    sd02.boot = sd(Test_Biomarker_2_0_new.boot)
    sd03.boot = sd(Test_Biomarker_3_0_new.boot)
    ## Sd of Diseased
    sd11.boot = sd(Test_Biomarker_1_1_new.boot)
    sd12.boot = sd(Test_Biomarker_2_1_new.boot)
    sd13.boot = sd(Test_Biomarker_3_1_new.boot)
    
    
    
    ## Correlation Coefficients in Non-diseased
    rho0_12.boot = cor(Test_Biomarker_1_0_new.boot, Test_Biomarker_2_0_new.boot)
    rho0_13.boot = cor(Test_Biomarker_1_0_new.boot, Test_Biomarker_3_0_new.boot)
    rho0_23.boot = cor(Test_Biomarker_2_0_new.boot, Test_Biomarker_3_0_new.boot)
    
    
    ## Correlation Coefficients in Diseased
    rho1_12.boot = cor(Test_Biomarker_1_1_new.boot, Test_Biomarker_2_1_new.boot)
    rho1_13.boot = cor(Test_Biomarker_1_1_new.boot, Test_Biomarker_3_1_new.boot)
    rho1_23.boot = cor(Test_Biomarker_2_1_new.boot, Test_Biomarker_3_1_new.boot)
    
    
    #####################################
    mu0.boot = c(mu01.boot, mu02.boot, mu03.boot)
    mu1.boot = c(mu11.boot, mu12.boot, mu13.boot)
    
    sd0.boot = c(sd01.boot, sd02.boot, sd03.boot)
    sd1.boot = c(sd11.boot, sd12.boot, sd13.boot)
    
    rho0.boot = c(rho0_12.boot, rho0_13.boot, rho0_23.boot)
    rho1.boot = c(rho1_12.boot, rho1_13.boot, rho1_23.boot)
    #########################################################################################################################
    
    
    mu0_comb.boot = comb_n(mu0.boot, 2)
    mu0.boot_local <- comb_n(length(mu0.boot), 2)
    mu0.boot_local
    
    mu1_comb.boot = comb_n(mu1.boot, 2)
    mu1.boot_local <- comb_n(length(mu1.boot), 2)
    mu1.boot_local
    
    sd0_comb.boot = comb_n(sd0.boot, 2)
    sd0.boot_local <- comb_n(length(sd0.boot), 2)
    sd0.boot_local
    
    sd1_comb.boot = comb_n(sd1.boot, 2)
    sd1.boot_local <- comb_n(length(sd1.boot), 2)
    sd1.boot_local
    
    
    
    lambda = 0.02
    
    
    
    mu0.mat.boot <- matrix(0, rho_leg, length(mu0))
    mu0.mat.boot
    mu1.mat.boot <- matrix(0, rho_leg, length(mu1))
    mu1.mat.boot
    sd0.mat.boot <- matrix(0, rho_leg, length(sd0))
    sd0.mat.boot
    sd1.mat.boot <- matrix(0, rho_leg, length(sd1))
    sd1.mat.boot
    rho0.mat.boot <- rep(0, rho_leg)
    rho0.mat.boot
    rho1.mat.boot <- rep(0, rho_leg)
    rho1.mat.boot
    p.vec.boot    <- rep(0, rho_leg)
    
    
    
    dVal <- rowSums(Dmat.boot) ###GroupM_true !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ## The case
    groupD <- ifelse(dVal == 0, yes = 0, no = 1)
    print(sum(groupD == 1))
    
    ### For the negative group
    D0Loc <- which(groupD == 0)
    D0M1.num <- sum(rbinom(sum(groupD == 0), size = 1, prob = 1 - pi0) == 1)
    #D0M1.num
    groupM.boot <- groupD
    D0M1Loc <- sample(D0Loc, D0M1.num, replace = F)
    groupM.boot[D0M1Loc] <- 1
    groupM.boot
    
    SenFun <- function(pi1, K, dnum, lambda)
    {
      res = pi1 * dnum /(dnum + lambda * (K - dnum))
      return(res)
    }
    
    ## For the positive group
    pi1.vec <- apply(matrix(1:K, ncol = 1), 1, SenFun, pi1 = pi1, K = K, lambda = lambda)
    for(dd in 1:K)
    {
      tempLoc <- which(dVal == dd)
      if(length(tempLoc) > 0)
      {
        tempNum <- sum(dVal == dd)
        D1M0.num <- sum(rbinom(tempNum, size = 1, prob = pi1.vec[dd]) == 0)
        if(D1M0.num > 0)
        {
          D1M0Loc <- sample(tempLoc, D1M0.num, replace = F)
          groupM.boot[D1M0Loc] <- 0
        }
      }
    }
    groupM.boot
    print(sum(groupM.boot == 1))
    
    pi1_lik = apply(matrix(1:K, ncol = 1), 1, SenFun, pi1 = pi1, K = K, lambda = lambda)
    pi1_lik.mat = as.matrix(pi1_lik)
    pi1_lik.mat
    
    ############ binary matrix ############
    a.mat <- rbind()
    pi1_lik.matrix <- rbind()
    for(d in 1:K)
    {
      for(i in 1 : choose(K, d))
      {
        a <- rep(0, K)
        a[as.matrix(comb_n(K, d))[ ,i]] <- 1
        pi1_lik.matrix <- rbind(pi1_lik.matrix, apply(matrix(d, ncol = 1), 1, SenFun, pi1 = pi1, K = K, lambda = lambda))
        a.mat <- rbind(a.mat,a)
      }
    }
    aa.mat <- matrix(a.mat, dim(a.mat)[1], dim(a.mat)[2])
    pi1_lik.matrix
    ############ Key Matrix #############
    aa.mat
    
    for(l in 1 : rho_leg)
    {
      
      XYvec.boot <- matrix(NA, N, 2)
      
      XYvec.boot[Dmat.boot==1, ] <- Bio1.boot[, Bio1_local[,l]]
      XYvec.boot[Dmat.boot==0, ] <- Bio0.boot[, Bio0_local[,l]]
      
      ####################################
      
      LoglikPMeanSD_est.boot <- function(theta, pi0, pi1, N, n, K, XYvec, lambda, groupM, aa.mat, pi1_lik.matrix)
      {
        mu01_est <- theta[1]
        mu02_est <- theta[2]
        mu11_est <- theta[3]
        mu12_est <- theta[4]
        sd01_est <- theta[5]
        sd02_est <- theta[6]
        sd11_est <- theta[7]
        sd12_est <- theta[8]
        rho0_est <- theta[9]
        rho1_est <- theta[10]
        p_est    <- theta[11]
        
        LoglikPMeanSD = rep(NA, n)
        
        f0.vec  <- rep(NA, N)
        f0.vec  <- (1/(2 * pi * sd01_est * sd02_est * sqrt(1 - rho0_est^2))) * exp(-1/(2 * (1 - rho0_est^2)) * (((XYvec.boot[ ,1] - mu01_est)/sd01_est)^2 - 2 * rho0_est * ((XYvec.boot[ ,1] - mu01_est)/sd01_est) * ((XYvec.boot[ ,2] - mu02_est)/sd02_est) + ((XYvec.boot[ ,2] - mu02_est)/sd02_est)^2))
        f0.mat  <- t(matrix(f0.vec, nrow = K, ncol = n))
        
        f1.vec  <- rep(NA, N)
        f1.vec  <- (1/(2 * pi * sd11_est * sd12_est * sqrt(1 - rho1_est^2))) * exp(-1/(2 * (1 - rho1_est^2)) * (((XYvec.boot[ ,1] - mu11_est)/sd11_est)^2 - 2 * rho1_est * ((XYvec.boot[ ,1] - mu11_est)/sd11_est) * ((XYvec.boot[ ,2] - mu12_est)/sd12_est) + ((XYvec.boot[ ,2] - mu12_est)/sd12_est)^2))
        f1.mat  <- t(matrix(f1.vec, nrow = K, ncol = n))
        
        
        mat.time <- dim(aa.mat)[1]
        f0.mat_final = matrix(rep(f0.mat, mat.time), n, K * mat.time)
        f1.mat_final = matrix(rep(f1.mat, mat.time), n, K * mat.time)
        
        
        aa.mat_new = t(replicate(n, as.vector(t(aa.mat))))
        
        
        if(K==1)
        {
          f0.matrix = f0.mat_final * t(1 - aa.mat_new)
          f1.matrix = f1.mat_final * t(aa.mat_new)
          
          f.matrix = f0.matrix + f1.matrix
        }else
        {
          f0.matrix = f0.mat_final * (1 - aa.mat_new)
          f1.matrix = f1.mat_final * aa.mat_new
          
          f.matrix = f0.matrix + f1.matrix
        }
        
        
        f.matrix_final <- matrix(NA, n, mat.time)
        for(j in 1 : mat.time)
        {
          for(i in 1 : n)
            f.matrix_final[i, j] <- prod(f.matrix[i , (1 + K * (j - 1)) : (K + K * (j - 1))])
        }
        
        
        ############### pi1 * p^d * (1 - p)^(K - d) ###############
        pd <- p_est^(apply(aa.mat, 1, sum)) * (1 - p_est)^(K-apply(aa.mat, 1, sum))
        pd.mat <- as.matrix(pd)
        pd.mat
        
        LoglikPMeanSD_final <- 0
        for(j in 1 : n)
        {
          
          if(groupM[j] == 0)
          {
            tempA <- 0
            
            tempA = tempA + f.matrix_final[j, ] %*% ((1 - pi1_lik.matrix) * pd.mat)
            
            LoglikPMeanSD[j] <- log(pi0 * (1 - p_est)^K * prod(f0.mat[j, ]) + tempA)
            
          }
          else
          {
            tempE <- 0
            
            tempE = tempE + f.matrix_final[j, ] %*% (pi1_lik.matrix * pd.mat)
            
            LoglikPMeanSD[j] <- log((1 - pi0) * (1 - p_est)^K * prod(f0.mat[j, ]) + tempE)
            
          }
          
          LoglikPMeanSD_final =  LoglikPMeanSD_final + LoglikPMeanSD[j]
        }
        
        res = -LoglikPMeanSD_final
        return(res)
      }
      
      print(j)
      EstPMeanSD.boot[j, ] <- as.vector(optim(c(mu0_comb.boot[1,l], mu0_comb.boot[2,l], mu1_comb.boot[1,l], mu1_comb.boot[2,l], sd0_comb.boot[1,l], sd0_comb.boot[2,l], sd1_comb.boot[1,l], sd1_comb.boot[2,l], rho0.boot[l], rho1.boot[l], p.boot),
                                              LoglikPMeanSD_est.boot, pi0 = pi0, pi1 = pi1, N = N, n = n, K = K, XYvec = XYvec.boot, lambda = lambda,
                                              groupM = groupM.boot, aa.mat = aa.mat, pi1_lik.matrix = pi1_lik.matrix)$par)

      
      mu01_final.boot <- EstPMeanSD.boot[j, 1]
      mu02_final.boot <- EstPMeanSD.boot[j, 2]
      mu11_final.boot <- EstPMeanSD.boot[j, 3]
      mu12_final.boot <- EstPMeanSD.boot[j, 4]
      sd01_final.boot <- EstPMeanSD.boot[j, 5]
      sd02_final.boot <- EstPMeanSD.boot[j, 6]
      sd11_final.boot <- EstPMeanSD.boot[j, 7]
      sd12_final.boot <- EstPMeanSD.boot[j, 8]
      rho0_final.boot <- EstPMeanSD.boot[j, 9]
      rho1_final.boot <- EstPMeanSD.boot[j, 10]
      p_est_final.boot <- EstPMeanSD.boot[j, 11]
      
      print(l)
      
      mu0.mat.boot[l, (mu0_local[ ,l])] <- mu0.mat.boot[l, (mu0_local[ ,l])] + c(mu01_final.boot, mu02_final.boot)
      mu1.mat.boot[l, (mu1_local[ ,l])] <- mu1.mat.boot[l, (mu1_local[ ,l])] + c(mu11_final.boot, mu12_final.boot)
      sd0.mat.boot[l, (sd0_local[ ,l])] <- sd0.mat.boot[l, (sd0_local[ ,l])] + c(sd01_final.boot, sd02_final.boot)
      sd1.mat.boot[l, (sd1_local[ ,l])] <- sd1.mat.boot[l, (sd1_local[ ,l])] + c(sd11_final.boot, sd12_final.boot)
      rho0.mat.boot[l] <- rho0.mat.boot[l] + rho0_final.boot
      rho1.mat.boot[l] <- rho1.mat.boot[l] + rho1_final.boot
      p.vec.boot[l]    <- p.vec.boot[l] + p_est_final.boot
      
    }
    
    mu0.final.boot[j, ] <- apply(mu0.mat.boot, 2, sum)/(length(mu0) - 1)
    mu1.final.boot[j, ] <- apply(mu1.mat.boot, 2, sum)/(length(mu1) - 1)
    sd0.final.boot[j, ] <- apply(sd0.mat.boot, 2, sum)/(length(sd0) - 1)
    sd1.final.boot[j, ] <- apply(sd1.mat.boot, 2, sum)/(length(sd1) - 1)
    rho0.final.boot[j, ] <- rho0.mat.boot
    rho1.final.boot[j, ] <- rho1.mat.boot
    p.final.boot[j] <- mean(na.omit(p.vec.boot))
    
    
    sigma_case_final.boot <- sigma_control_final.boot <- matrix(NA, 3, 3)
    sigma_case_final.boot <- matrix(c(sd1.final.boot[j, 1]^2, rho1.final.boot[j, 1] * sd1.final.boot[j, 1] * sd1.final.boot[j, 2], rho1.final.boot[j, 2] * sd1.final.boot[j, 1] * sd1.final.boot[j, 3],
                                      rho1.final.boot[j, 1] * sd1.final.boot[j, 2] * sd1.final.boot[j, 1], sd1.final.boot[j, 2]^2, rho1.final.boot[j, 3] * sd1.final.boot[j, 2] * sd1.final.boot[j, 3],
                                      rho1.final.boot[j, 2] * sd1.final.boot[j, 3] * sd1.final.boot[j, 1], rho1.final.boot[j, 3] * sd1.final.boot[j, 3] * sd1.final.boot[j, 2], sd1.final.boot[j, 3]^2), 3, 3)
    sigma_case_final.boot
    sigma_control_final.boot <- matrix(c(sd0.final.boot[j, 1]^2, rho0.final.boot[j, 1] * sd0.final.boot[j, 1] * sd0.final.boot[j, 2], rho0.final.boot[j, 2] * sd0.final.boot[j, 1] * sd0.final.boot[j, 3],
                                         rho0.final.boot[j, 1] * sd0.final.boot[j, 2] * sd0.final.boot[j, 1], sd0.final.boot[j, 2]^2, rho0.final.boot[j, 3] * sd0.final.boot[j, 2] * sd0.final.boot[j, 3],
                                         rho0.final.boot[j, 2] * sd0.final.boot[j, 3] * sd0.final.boot[j, 1], rho0.final.boot[j, 3] * sd0.final.boot[j, 3] * sd0.final.boot[j, 2], sd0.final.boot[j, 3]^2), 3, 3)
    sigma_control_final.boot
    AUC.value.boot[j] <- pnorm(sqrt(t(mu1.final.boot[j, ] - mu0.final.boot[j, ]) %*% solve(sigma_control_final.boot + sigma_case_final.boot) %*% (mu1.final.boot[j, ] - mu0.final.boot[j, ])))
    print(AUC.value.boot[j])
    # write.table(AUC.value.boot[j], file = paste0("./results/Boot/AUC_boot/TTB_scaleSep_Seed.boot_", j , "_", "K=", K, "_", CurrentData, ".txt"),
    #             sep="\t", row.names=F)
    
    ####################################### Individual ###########################################
    ################################## Bio_1
    XYvec_Bio1 <- rep(NA, N)
    XYvec_Bio1[Dmat.boot==1] <- Test_Biomarker_1_1_new.boot
    XYvec_Bio1[Dmat.boot==0] <- Test_Biomarker_1_0_new.boot
    
    Loglik_Bio1.boot <- function(theta, pi0, pi1, N, n, K, XYvec, lambda, groupM, aa.mat, pi1_lik.matrix)
    {
      mu0_Bio1 <- theta[1]
      mu1_Bio1 <- theta[2]
      sd0_Bio1 <- theta[3]
      sd1_Bio1 <- theta[4]
      p_est    <- theta[5]
      
      LoglikPMeanSD = rep(NA, n)
      
      f0.vec  <- rep(NA, N)
      f0.vec  <- dnorm(XYvec, mean = mu0_Bio1, sd = sd0_Bio1)
      f0.mat  <- t(matrix(f0.vec, nrow = K, ncol = n))
      
      f1.vec  <- rep(NA, N)
      f1.vec  <- dnorm(XYvec, mean = mu1_Bio1, sd = sd1_Bio1)
      f1.mat  <- t(matrix(f1.vec, nrow = K, ncol = n))
      
      
      mat.time <- dim(aa.mat)[1]
      f0.mat_final = matrix(rep(f0.mat, mat.time), n, K * mat.time)
      f1.mat_final = matrix(rep(f1.mat, mat.time), n, K * mat.time)
      
      
      aa.mat_new = t(replicate(n, as.vector(t(aa.mat))))
      
      
      if(K==1)
      {
        f0.matrix = f0.mat_final * t(1 - aa.mat_new)
        f1.matrix = f1.mat_final * t(aa.mat_new)
        
        f.matrix = f0.matrix + f1.matrix
      }else
      {
        f0.matrix = f0.mat_final * (1 - aa.mat_new)
        f1.matrix = f1.mat_final * aa.mat_new
        
        f.matrix = f0.matrix + f1.matrix
      }
      
      
      f.matrix_final <- matrix(NA, n, mat.time)
      for(j in 1 : mat.time)
      {
        for(i in 1 : n)
          f.matrix_final[i, j] <- prod(f.matrix[i , (1 + K * (j - 1)) : (K + K * (j - 1))])
      }
      
      
      ############### pi1 * p^d * (1 - p)^(K - d) ###############
      pd <- p_est^(apply(aa.mat, 1, sum)) * (1 - p_est)^(K-apply(aa.mat, 1, sum))
      pd.mat <- as.matrix(pd)
      pd.mat
      
      LoglikPMeanSD_final <- 0
      for(j in 1 : n)
      {
        
        if(groupM[j] == 0)
        {
          tempA <- 0
          
          tempA = tempA + f.matrix_final[j, ] %*% ((1 - pi1_lik.matrix) * pd.mat)
          
          LoglikPMeanSD[j] <- log(pi0 * (1 - p_est)^K * prod(f0.mat[j, ]) + tempA)
          
        }
        else
        {
          tempE <- 0
          
          tempE = tempE + f.matrix_final[j, ] %*% (pi1_lik.matrix * pd.mat)
          
          LoglikPMeanSD[j] <- log((1 - pi0) * (1 - p_est)^K * prod(f0.mat[j, ]) + tempE)
          
        }
        
        LoglikPMeanSD_final =  LoglikPMeanSD_final + LoglikPMeanSD[j]
      }
      
      res = -LoglikPMeanSD_final
      return(res)
    }
    
    print(j)
    Est_Bio1.boot[j, ] <- as.vector(optim(c(mu01.boot, mu11.boot, sd01.boot, sd11.boot, p.boot), Loglik_Bio1.boot, pi0 = pi0, pi1 = pi1, N = N, n = n, K = K, XYvec = XYvec_Bio1,
                                          lambda = lambda, groupM = groupM.boot, aa.mat = aa.mat, pi1_lik.matrix = pi1_lik.matrix)$par)
    
    mu0_Bio1.boot[j] <- Est_Bio1.boot[j, 1]
    mu1_Bio1.boot[j] <- Est_Bio1.boot[j, 2]
    sd0_Bio1.boot[j] <- Est_Bio1.boot[j, 3]
    sd1_Bio1.boot[j] <- Est_Bio1.boot[j, 4]
    p_Bio1.boot[j]   <- Est_Bio1.boot[j, 5]
    
    
    AUC_1.boot[j] <- pnorm(sqrt(t(mu1_Bio1.boot[j] - mu0_Bio1.boot[j]) %*% solve(sd0_Bio1.boot[j]^2 + sd1_Bio1.boot[j]^2) %*% (mu1_Bio1.boot[j] - mu0_Bio1.boot[j])))
    AUC_1.boot[j]
    
    
    ################################## Bio_2
    XYvec_Bio2 <- rep(NA, N)
    XYvec_Bio2[Dmat.boot==1] <- Test_Biomarker_2_1_new.boot
    XYvec_Bio2[Dmat.boot==0] <- Test_Biomarker_2_0_new.boot
    
    Loglik_Bio2.boot <- function(theta, pi0, pi1, N, n, K, XYvec, lambda, groupM, aa.mat, pi1_lik.matrix)
    {
      mu0_Bio2 <- theta[1]
      mu1_Bio2 <- theta[2]
      sd0_Bio2 <- theta[3]
      sd1_Bio2 <- theta[4]
      p_est    <- theta[5]
      
      LoglikPMeanSD = rep(NA, n)
      
      f0.vec  <- rep(NA, N)
      f0.vec  <- dnorm(XYvec, mean = mu0_Bio2, sd = sd0_Bio2)
      f0.mat  <- t(matrix(f0.vec, nrow = K, ncol = n))
      
      f1.vec  <- rep(NA, N)
      f1.vec  <- dnorm(XYvec, mean = mu1_Bio2, sd = sd1_Bio2)
      f1.mat  <- t(matrix(f1.vec, nrow = K, ncol = n))
      
      
      mat.time <- dim(aa.mat)[1]
      f0.mat_final = matrix(rep(f0.mat, mat.time), n, K * mat.time)
      f1.mat_final = matrix(rep(f1.mat, mat.time), n, K * mat.time)
      
      
      aa.mat_new = t(replicate(n, as.vector(t(aa.mat))))
      
      
      if(K==1)
      {
        f0.matrix = f0.mat_final * t(1 - aa.mat_new)
        f1.matrix = f1.mat_final * t(aa.mat_new)
        
        f.matrix = f0.matrix + f1.matrix
      }else
      {
        f0.matrix = f0.mat_final * (1 - aa.mat_new)
        f1.matrix = f1.mat_final * aa.mat_new
        
        f.matrix = f0.matrix + f1.matrix
      }
      
      
      f.matrix_final <- matrix(NA, n, mat.time)
      for(j in 1 : mat.time)
      {
        for(i in 1 : n)
          f.matrix_final[i, j] <- prod(f.matrix[i , (1 + K * (j - 1)) : (K + K * (j - 1))])
      }
      
      
      ############### pi1 * p^d * (1 - p)^(K - d) ###############
      pd <- p_est^(apply(aa.mat, 1, sum)) * (1 - p_est)^(K-apply(aa.mat, 1, sum))
      pd.mat <- as.matrix(pd)
      pd.mat
      
      LoglikPMeanSD_final <- 0
      for(j in 1 : n)
      {
        
        if(groupM[j] == 0)
        {
          tempA <- 0
          
          tempA = tempA + f.matrix_final[j, ] %*% ((1 - pi1_lik.matrix) * pd.mat)
          
          LoglikPMeanSD[j] <- log(pi0 * (1 - p_est)^K * prod(f0.mat[j, ]) + tempA)
          
        }
        else
        {
          tempE <- 0
          
          tempE = tempE + f.matrix_final[j, ] %*% (pi1_lik.matrix * pd.mat)
          
          LoglikPMeanSD[j] <- log((1 - pi0) * (1 - p_est)^K * prod(f0.mat[j, ]) + tempE)
          
        }
        
        LoglikPMeanSD_final =  LoglikPMeanSD_final + LoglikPMeanSD[j]
      }
      
      res = -LoglikPMeanSD_final
      return(res)
    }
    
    print(j)
    Est_Bio2.boot[j, ] <- as.vector(optim(c(mu02.boot, mu12.boot, sd02.boot, sd12.boot, p.boot), Loglik_Bio2.boot, pi0 = pi0, pi1 = pi1, N = N, n = n, K = K, XYvec = XYvec_Bio2,
                                          lambda = lambda, groupM = groupM.boot, aa.mat = aa.mat, pi1_lik.matrix = pi1_lik.matrix)$par)
    
    mu0_Bio2.boot[j] <- Est_Bio2.boot[j, 1]
    mu1_Bio2.boot[j] <- Est_Bio2.boot[j, 2]
    sd0_Bio2.boot[j] <- Est_Bio2.boot[j, 3]
    sd1_Bio2.boot[j] <- Est_Bio2.boot[j, 4]
    p_Bio2.boot[j]   <- Est_Bio2.boot[j, 5]
    
    
    AUC_2.boot[j] <- pnorm(sqrt(t(mu1_Bio2.boot[j] - mu0_Bio2.boot[j]) %*% solve(sd0_Bio2.boot[j]^2 + sd1_Bio2.boot[j]^2) %*% (mu1_Bio2.boot[j] - mu0_Bio2.boot[j])))
    AUC_2.boot[j]
    
    
    ################################## Bio_3
    XYvec_Bio3 <- rep(NA, N)
    XYvec_Bio3[Dmat.boot==1] <- Test_Biomarker_3_1_new.boot
    XYvec_Bio3[Dmat.boot==0] <- Test_Biomarker_3_0_new.boot
    
    Loglik_Bio3.boot <- function(theta, pi0, pi1, N, n, K, XYvec, lambda, groupM, aa.mat, pi1_lik.matrix)
    {
      mu0_Bio3 <- theta[1]
      mu1_Bio3 <- theta[2]
      sd0_Bio3 <- theta[3]
      sd1_Bio3 <- theta[4]
      p_est    <- theta[5]
      
      LoglikPMeanSD = rep(NA, n)
      
      f0.vec  <- rep(NA, N)
      f0.vec  <- dnorm(XYvec, mean = mu0_Bio3, sd = sd0_Bio3)
      f0.mat  <- t(matrix(f0.vec, nrow = K, ncol = n))
      
      f1.vec  <- rep(NA, N)
      f1.vec  <- dnorm(XYvec, mean = mu1_Bio3, sd = sd1_Bio3)
      f1.mat  <- t(matrix(f1.vec, nrow = K, ncol = n))
      
      
      mat.time <- dim(aa.mat)[1]
      f0.mat_final = matrix(rep(f0.mat, mat.time), n, K * mat.time)
      f1.mat_final = matrix(rep(f1.mat, mat.time), n, K * mat.time)
      
      
      aa.mat_new = t(replicate(n, as.vector(t(aa.mat))))
      
      
      if(K==1)
      {
        f0.matrix = f0.mat_final * t(1 - aa.mat_new)
        f1.matrix = f1.mat_final * t(aa.mat_new)
        
        f.matrix = f0.matrix + f1.matrix
      }else
      {
        f0.matrix = f0.mat_final * (1 - aa.mat_new)
        f1.matrix = f1.mat_final * aa.mat_new
        
        f.matrix = f0.matrix + f1.matrix
      }
      
      
      f.matrix_final <- matrix(NA, n, mat.time)
      for(j in 1 : mat.time)
      {
        for(i in 1 : n)
          f.matrix_final[i, j] <- prod(f.matrix[i , (1 + K * (j - 1)) : (K + K * (j - 1))])
      }
      
      
      ############### pi1 * p^d * (1 - p)^(K - d) ###############
      pd <- p_est^(apply(aa.mat, 1, sum)) * (1 - p_est)^(K-apply(aa.mat, 1, sum))
      pd.mat <- as.matrix(pd)
      pd.mat
      
      LoglikPMeanSD_final <- 0
      for(j in 1 : n)
      {
        
        if(groupM[j] == 0)
        {
          tempA <- 0
          
          tempA = tempA + f.matrix_final[j, ] %*% ((1 - pi1_lik.matrix) * pd.mat)
          
          LoglikPMeanSD[j] <- log(pi0 * (1 - p_est)^K * prod(f0.mat[j, ]) + tempA)
          
        }
        else
        {
          tempE <- 0
          
          tempE = tempE + f.matrix_final[j, ] %*% (pi1_lik.matrix * pd.mat)
          
          LoglikPMeanSD[j] <- log((1 - pi0) * (1 - p_est)^K * prod(f0.mat[j, ]) + tempE)
          
        }
        
        LoglikPMeanSD_final =  LoglikPMeanSD_final + LoglikPMeanSD[j]
      }
      
      res = -LoglikPMeanSD_final
      return(res)
    }
    
    print(j)
    Est_Bio3.boot[j, ] <- as.vector(optim(c(mu03.boot, mu13.boot, sd03.boot, sd13.boot, p.boot), Loglik_Bio3.boot, pi0 = pi0, pi1 = pi1, N = N, n = n, K = K, XYvec = XYvec_Bio3,
                                          lambda = lambda, groupM = groupM.boot, aa.mat = aa.mat, pi1_lik.matrix = pi1_lik.matrix)$par)
    
    mu0_Bio3.boot[j] <- Est_Bio3.boot[j, 1]
    mu1_Bio3.boot[j] <- Est_Bio3.boot[j, 2]
    sd0_Bio3.boot[j] <- Est_Bio3.boot[j, 3]
    sd1_Bio3.boot[j] <- Est_Bio3.boot[j, 4]
    p_Bio3.boot[j]   <- Est_Bio3.boot[j, 5]
    
    
    AUC_3.boot[j] <- pnorm(sqrt(t(mu1_Bio3.boot[j] - mu0_Bio3.boot[j]) %*% solve(sd0_Bio3.boot[j]^2 + sd1_Bio3.boot[j]^2) %*% (mu1_Bio3.boot[j] - mu0_Bio3.boot[j])))
    AUC_3.boot[j]
    
    
    return(list("mu0.final.boot"  = mu0.final.boot[j, ],  "mu1.final.boot" = mu1.final.boot[j, ],
                "sd0.final.boot"  = sd0.final.boot[j, ],  "sd1.final.boot" = sd1.final.boot[j, ],
                "rho0.final.boot" = rho0.final.boot[j, ], "rho1.final.boot" = rho1.final.boot[j, ],
                "AUC.value.boot"  = AUC.value.boot[j], "p.final.boot" = p.final.boot[j], ## Combine
                "mu0_Bio1.boot" = mu0_Bio1.boot[j], "mu1_Bio1.boot" = mu1_Bio1.boot[j], "sd0_Bio1.boot" = sd0_Bio1.boot[j], "sd1_Bio1.boot" = sd1_Bio1.boot[j],
                "p_Bio1.boot" = p_Bio1.boot[j], "AUC_1.boot" = AUC_1.boot[j], ## Bio_1
                "mu0_Bio2.boot" = mu0_Bio2.boot[j], "mu1_Bio2.boot" = mu1_Bio2.boot[j], "sd0_Bio2.boot" = sd0_Bio2.boot[j], "sd1_Bio2.boot" = sd1_Bio2.boot[j],
                "p_Bio2.boot" = p_Bio2.boot[j], "AUC_2.boot" = AUC_2.boot[j], ## Bio_1
                "mu0_Bio3.boot" = mu0_Bio3.boot[j], "mu1_Bio3.boot" = mu1_Bio3.boot[j], "sd0_Bio3.boot" = sd0_Bio3.boot[j], "sd1_Bio3.boot" = sd1_Bio3.boot[j],
                "p_Bio3.boot" = p_Bio3.boot[j], "AUC_3.boot" = AUC_3.boot[j] ## Bio_1
    )
    
    )
  }

# EstPMeanSD.vec.boot <- matrix(NA, seed.num.boot, 11)
p_est_final.vec.boot <- rep(NA, seed.num.boot)
AUC.value_vec.boot <- rep(NA, seed.num.boot)
mu0.final_mat.boot  <- matrix(NA, seed.num.boot, length(mu0))
mu1.final_mat.boot  <- matrix(NA, seed.num.boot, length(mu1))
sd0.final_mat.boot  <- matrix(NA, seed.num.boot, length(sd0))
sd1.final_mat.boot  <- matrix(NA, seed.num.boot, length(sd1))
rho0.final_mat.boot <- matrix(NA, seed.num.boot, length(rho0))
rho1.final_mat.boot <- matrix(NA, seed.num.boot, length(rho1))
############# Individual ###################
## Bio_1
mu0_Bio1_vec.boot <- rep(NA, seed.num.boot)
mu1_Bio1_vec.boot <- rep(NA, seed.num.boot)
sd0_Bio1_vec.boot <- rep(NA, seed.num.boot)
sd1_Bio1_vec.boot <- rep(NA, seed.num.boot)
p_Bio1_vec.boot   <- rep(NA, seed.num.boot)
AUC_1_vec.boot    <- rep(NA, seed.num.boot)
## Bio_2
mu0_Bio2_vec.boot <- rep(NA, seed.num.boot)
mu1_Bio2_vec.boot <- rep(NA, seed.num.boot)
sd0_Bio2_vec.boot <- rep(NA, seed.num.boot)
sd1_Bio2_vec.boot <- rep(NA, seed.num.boot)
p_Bio2_vec.boot   <- rep(NA, seed.num.boot)
AUC_2_vec.boot    <- rep(NA, seed.num.boot)
## Bio_3
mu0_Bio3_vec.boot <- rep(NA, seed.num.boot)
mu1_Bio3_vec.boot <- rep(NA, seed.num.boot)
sd0_Bio3_vec.boot <- rep(NA, seed.num.boot)
sd1_Bio3_vec.boot <- rep(NA, seed.num.boot)
p_Bio3_vec.boot   <- rep(NA, seed.num.boot)
AUC_3_vec.boot    <- rep(NA, seed.num.boot)

for(i in 1 : seed.num.boot)
{
  #EstPMeanSD.vec.boot[i, ] <- Main_Fun[[i]]$EstPMeanSD.boot
  p_est_final.vec.boot[i]  <- Main_Fun[[i]]$p.final.boot
  AUC.value_vec.boot[i]    <- Main_Fun[[i]]$AUC.value.boot
  mu0.final_mat.boot[i, ]  <- Main_Fun[[i]]$mu0.final.boot
  mu1.final_mat.boot[i, ]  <- Main_Fun[[i]]$mu1.final.boot
  sd0.final_mat.boot[i, ]  <- Main_Fun[[i]]$sd0.final.boot
  sd1.final_mat.boot[i, ]  <- Main_Fun[[i]]$sd1.final.boot
  rho0.final_mat.boot[i, ] <- Main_Fun[[i]]$rho0.final.boot
  rho1.final_mat.boot[i, ] <- Main_Fun[[i]]$rho1.final.boot
  ############# Individual ###################
  ## Bio_1
  mu0_Bio1_vec.boot[i] <- Main_Fun[[i]]$mu0_Bio1.boot
  mu1_Bio1_vec.boot[i] <- Main_Fun[[i]]$mu1_Bio1.boot
  sd0_Bio1_vec.boot[i] <- Main_Fun[[i]]$sd0_Bio1.boot
  sd1_Bio1_vec.boot[i] <- Main_Fun[[i]]$sd1_Bio1.boot
  p_Bio1_vec.boot[i]   <- Main_Fun[[i]]$p_Bio1.boot
  AUC_1_vec.boot[i]    <- Main_Fun[[i]]$AUC_1.boot
  ## Bio_2
  mu0_Bio2_vec.boot[i] <- Main_Fun[[i]]$mu0_Bio2.boot
  mu1_Bio2_vec.boot[i] <- Main_Fun[[i]]$mu1_Bio2.boot
  sd0_Bio2_vec.boot[i] <- Main_Fun[[i]]$sd0_Bio2.boot
  sd1_Bio2_vec.boot[i] <- Main_Fun[[i]]$sd1_Bio2.boot
  p_Bio2_vec.boot[i]   <- Main_Fun[[i]]$p_Bio2.boot
  AUC_2_vec.boot[i]    <- Main_Fun[[i]]$AUC_2.boot
  ## Bio_3
  mu0_Bio3_vec.boot[i] <- Main_Fun[[i]]$mu0_Bio3.boot
  mu1_Bio3_vec.boot[i] <- Main_Fun[[i]]$mu1_Bio3.boot
  sd0_Bio3_vec.boot[i] <- Main_Fun[[i]]$sd0_Bio3.boot
  sd1_Bio3_vec.boot[i] <- Main_Fun[[i]]$sd1_Bio3.boot
  p_Bio3_vec.boot[i]   <- Main_Fun[[i]]$p_Bio3.boot
  AUC_3_vec.boot[i]    <- Main_Fun[[i]]$AUC_3.boot
}


p.final_boot <- mean(na.omit(p_est_final.vec.boot))
AUC.est.boot <- mean(na.omit(AUC.value_vec.boot))
AUC.var.boot <- var(na.omit(AUC.value_vec.boot))
AUC.len.boot  <- length(na.omit(AUC.value_vec.boot))
AUC.sd.boot <- sd(na.omit(AUC.value_vec.boot))
AUC.CI.left.boot  <- AUC.est.boot - qnorm(1 - alpha/2) * AUC.sd.boot
AUC.CI.right.boot <- AUC.est.boot + qnorm(1 - alpha/2) * AUC.sd.boot

p.bias.boot <- mean(na.omit(p_est_final.vec.boot)) - p
p.MSE.boot  <- mean((na.omit(p_est_final.vec.boot) - p)^2)
p.var.boot <- var(na.omit(p_est_final.vec.boot))
p.sd.boot <- sd(na.omit(p_est_final.vec.boot))
p.CI.left.boot  <- p.final_boot - qnorm(1 - alpha/2) * p.sd.boot
p.CI.right.boot <- p.final_boot + qnorm(1 - alpha/2) * p.sd.boot
########################## Individual #######################################
## Bio_1
mu0_Bio1_final.boot <- mean(na.omit(mu0_Bio1_vec.boot))
mu1_Bio1_final.boot <- mean(na.omit(mu1_Bio1_vec.boot))
sd0_Bio1_final.boot <- mean(na.omit(sd0_Bio1_vec.boot))
sd1_Bio1_final.boot <- mean(na.omit(sd1_Bio1_vec.boot))
p_Bio1.final_est   <- mean(na.omit(p_Bio1_vec.boot))
p_Bio1.final.var   <- var(na.omit(p_Bio1_vec.boot))
AUC_1.est <- mean(na.omit(AUC_1_vec.boot))
AUC_1.var <- var(na.omit(AUC_1_vec.boot))
## Bio_2
mu0_Bio2_final.boot <- mean(na.omit(mu0_Bio2_vec.boot))
mu1_Bio2_final.boot <- mean(na.omit(mu1_Bio2_vec.boot))
sd0_Bio2_final.boot <- mean(na.omit(sd0_Bio2_vec.boot))
sd1_Bio2_final.boot <- mean(na.omit(sd1_Bio2_vec.boot))
p_Bio2.final_est   <- mean(na.omit(p_Bio2_vec.boot))
p_Bio2.final.var   <- var(na.omit(p_Bio2_vec.boot))
AUC_2.est <- mean(na.omit(AUC_2_vec.boot))
AUC_2.var <- var(na.omit(AUC_2_vec.boot))
## Bio_3
mu0_Bio3_final.boot <- mean(na.omit(mu0_Bio3_vec.boot))
mu1_Bio3_final.boot <- mean(na.omit(mu1_Bio3_vec.boot))
sd0_Bio3_final.boot <- mean(na.omit(sd0_Bio3_vec.boot))
sd1_Bio3_final.boot <- mean(na.omit(sd1_Bio3_vec.boot))
p_Bio3.final_est   <- mean(na.omit(p_Bio3_vec.boot))
p_Bio3.final.var   <- var(na.omit(p_Bio3_vec.boot))
AUC_3.est <- mean(na.omit(AUC_3_vec.boot))
AUC_3.var <- var(na.omit(AUC_3_vec.boot))

sim_stats =
  as.vector(
    list("mu0.final_est.boot" = apply(mu0.final_mat.boot, 2, mean), "mu1.final_est.boot" = apply(mu1.final_mat.boot, 2, mean),
         "sd0.final_est.boot" = apply(sd0.final_mat.boot, 2, mean),  "sd1.final_est.boot" = apply(sd1.final_mat.boot, 2, mean),
         "rho0.final_est.boot" = apply(rho0.final_mat.boot, 2, mean), "rho1.final_est.boot" = apply(rho1.final_mat.boot, 2, mean),
         "p.final.boot" = p.final_boot, "p.bias.boot" = p.bias.boot, "p.MSE.boot" = p.MSE.boot, "p.var.boot" = p.var.boot, "p.CI.left.boot" = p.CI.left.boot, "p.CI.right.boot" = p.CI.right.boot,
         "AUC_true" = AUC_true, "AUC.est.boot" = AUC.est.boot, "AUC.var.boot" = AUC.var.boot, "AUC.len.boot" = AUC.len.boot, "AUC.sd.boot" = AUC.sd.boot,
         "AUC.CI.left.boot" = AUC.CI.left.boot, "AUC.CI.right.boot" = AUC.CI.right.boot, ### Combine
         "mu0_Bio1_final.boot" = mu0_Bio1_final.boot, "mu1_Bio1_final.boot" = mu1_Bio1_final.boot, "sd0_Bio1_final.boot" = sd0_Bio1_final.boot, "sd1_Bio1_final.boot" = sd1_Bio1_final.boot,
         "p_Bio1.final_est" = p_Bio1.final_est, "p_Bio1.final.var" = p_Bio1.final.var, "AUC_1.est" = AUC_1.est, "AUC_1.var" = AUC_1.var, ## Bio1
         "mu0_Bio2_final.boot" = mu0_Bio2_final.boot, "mu1_Bio2_final.boot" = mu1_Bio2_final.boot, "sd0_Bio2_final.boot" = sd0_Bio2_final.boot, "sd1_Bio2_final.boot" = sd1_Bio2_final.boot,
         "p_Bio2.final_est" = p_Bio2.final_est, "p_Bio2.final.var" = p_Bio2.final.var, "AUC_2.est" = AUC_2.est, "AUC_2.var" = AUC_2.var, ## Bio2
         "mu0_Bio3_final.boot" = mu0_Bio3_final.boot, "mu1_Bio3_final.boot" = mu1_Bio3_final.boot, "sd0_Bio3_final.boot" = sd0_Bio3_final.boot, "sd1_Bio3_final.boot" = sd1_Bio3_final.boot,
         "p_Bio3.final_est" = p_Bio3.final_est, "p_Bio3.final.var" = p_Bio3.final.var, "AUC_3.est" = AUC_3.est, "AUC_3.var" = AUC_3.var ## Bio1
    )
  )


write.csv(sim_stats, file = paste0("./results/Boot/TTB_scaleSep_Seed.boot=", seed.num.boot, "_", "N=", N, "_", "K=", K, "_", Pi0, "_", Pi1, "_", CurrentData, ".csv"))


stopCluster(cl)

time2 <- Sys.time()
print(time2-time1)



