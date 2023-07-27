rm(list = ls())
library(foreach)
library(doParallel)
library(combinat)
library(bbmle)
library(Rfast)
library(MASS) 

CurrentData <- Sys.Date()
# SEED.NUM = as.numeric(commandArgs(TRUE)[1])
# pp       = as.numeric(commandArgs(TRUE)[2])
# KK       = as.numeric(commandArgs(TRUE)[3])
# Pi0      = as.numeric(commandArgs(TRUE)[4])
# Pi1      = as.numeric(commandArgs(TRUE)[5])
SEED.NUM = 1000
#pp       = 0.02
#KK       = 1
Pi0      = 0.95
Pi1      = 0.95
time1 <- Sys.time()


seed.num = SEED.NUM
# SIM_DIR = "Results"
p = 0.01
N = 20000
K = 1
n = N/ K

lambda = 0.02
pi0 = Pi0
pi1 = Pi1
alpha = 0.05

#SIM_DIR <- paste0(SIM_DIR, "_", pp, "_", KK, "_", Pi0, "_", Pi1)
mu01 = 0
mu02 = 0
mu03 = 0

mu11 = 1
mu12 = 1.5
mu13 = 1.2

sd01 = 1
sd02 = 1
sd03 = 1

sd11 = 1.1
sd12 = 1.2
sd13 = 1.3

rho01 = 1/3
rho02 = 1/3
rho03 = 1/3

rho11 = 0.5
rho12 = 0.6
rho13 = 0.3

mu0 = c(mu01, mu02, mu03)
mu0_comb = comb_n(mu0, 2)
mu0_local = comb_n(length(mu0), 2)

mu1 = c(mu11, mu12, mu13)
mu1_comb = comb_n(mu1, 2)
mu1_local <- comb_n(length(mu1), 2)

sd0 = c(sd01, sd02, sd03)
sd0_comb = comb_n(sd0, 2)
sd0_local <- comb_n(length(sd0), 2) 

sd1 = c(sd11, sd12, sd13)
sd1_comb = comb_n(sd1, 2)
sd1_local <- comb_n(length(sd1), 2)  


rho0 = c(rho01, rho02, rho03)
rho1 = c(rho11, rho12, rho13)

rho_leg = length(rho0)

sigma_0 = matrix(c(sd01^2, rho01 * sd01 * sd02, rho02 * sd01 * sd03, 
                   rho01 * sd02 * sd01, sd02^2, rho03 * sd02 * sd03,
                   rho02 * sd03 * sd01, rho03 * sd03 * sd02, sd03^2), 3, 3)
sigma_0
sigma_1 = matrix(c(sd11^2, rho11 * sd11 * sd12, rho12 * sd11 * sd13, 
                   rho11 * sd12 * sd11, sd12^2, rho13 * sd12 * sd13,
                   rho12 * sd13 * sd11, rho13 * sd13 * sd12, sd13^2), 3, 3)
sigma_1

mu_true <- as.matrix(mu1 - mu0)
AUC_true <- pnorm(sqrt(t(mu_true) %*% solve(sigma_1 + sigma_0) %*% mu_true))
AUC_true


##########################################################################################################
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
#CDF_true
knotNum <- length(CDF_true)
aa <- CDF_true

ROC.Fun_1 <- function(u)
{
  return(pnorm((sd01 * qnorm(u) + mu11 - mu01)/sd11))
}
#pnorm((CDF_true-mu01)/sd01)
ROC.true_1 <- ROC.Fun_1(aa)
#ROC.Fun_1 <- Vectorize(ROC.Fun_1, "u")

plot(aa,ROC.true_1, type = "l")

abline(0, 1, col = "red", lwd = 0.2, lty = 1)


ROC.Fun_2 <- function(u)
{
  return(pnorm((sd02 * qnorm(u) + mu12 - mu02)/sd12))
}
ROC.true_2 <- ROC.Fun_2(aa)
#ROC.true_2 <- Vectorize(ROC.true_2)

plot(aa,ROC.true_2, type = "l")

abline(0, 1, col = "red", lwd = 0.2, lty = 1)


ROC.Fun_3 <- function(u)
{
  return(pnorm((sd03 * qnorm(u) + mu13 - mu03)/sd13))
}
ROC.true_3 <- ROC.Fun_3(aa)

plot(aa,ROC.true_3, type = "l")

abline(0, 1, col = "red", lwd = 0.2, lty = 1)

# ###### AUC Combination
c <- t(mu_true) %*% solve(sigma_1 + sigma_0)
c

linear_alpha <- (1/c[,1]) * t(mu_true) %*% solve(sigma_1 + sigma_0)
linear_alpha

new_mu <- c(linear_alpha %*% mu_true)
new_mu
new_sigma0 <- c(linear_alpha %*% sigma_0 %*% t(linear_alpha))
new_sigma0
new_sigma1 <- c(linear_alpha %*% sigma_1 %*% t(linear_alpha))
new_sigma1

AUC_try <- pnorm(sqrt(t(new_mu) %*% solve(new_sigma0 + new_sigma1) %*% new_mu))
AUC_try
ROC.Fun <- function(u)
{
  return(pnorm((new_mu + sqrt(new_sigma0) * qnorm(u))/sqrt(new_sigma1)))
}

ROC_true <- ROC.Fun(aa)
plot(aa, ROC_true, type = "l")
abline(0, 1, col = "red", lwd = 1, lty = 1)

##########################################################################################################



EstPMeanSD <- matrix(NA, seed.num, 11)
Est_Bio1 <- matrix(NA, seed.num, 5)
Est_Bio2 <- matrix(NA, seed.num, 5)
Est_Bio3 <- matrix(NA, seed.num, 5)

#p_est_final <- rep(NA, seed.num)
mu0.final  <- matrix(NA, seed.num, length(mu0))
mu1.final  <- matrix(NA, seed.num, length(mu1))
sd0.final  <- matrix(NA, seed.num, length(sd0))
sd1.final  <- matrix(NA, seed.num, length(sd1))
rho0.final <- matrix(NA, seed.num, length(rho0))
rho1.final <- matrix(NA, seed.num, length(rho1))
p.vec.final <- rep(NA, seed.num)
AUC.value  <- rep(NA, seed.num)
############# Individual ###########

mu0_Bio1_final <- rep(NA, seed.num)
mu1_Bio1_final <- rep(NA, seed.num)
sd0_Bio1_final <- rep(NA, seed.num)
sd1_Bio1_final <- rep(NA, seed.num)
p_Bio1_final   <- rep(NA, seed.num)
AUC_1.value  <- rep(NA, seed.num)

mu0_Bio2_final <- rep(NA, seed.num)
mu1_Bio2_final <- rep(NA, seed.num)
sd0_Bio2_final <- rep(NA, seed.num)
sd1_Bio2_final <- rep(NA, seed.num)
p_Bio2_final   <- rep(NA, seed.num)
AUC_2.value  <- rep(NA, seed.num)

mu0_Bio3_final <- rep(NA, seed.num)
mu1_Bio3_final <- rep(NA, seed.num)
sd0_Bio3_final <- rep(NA, seed.num)
sd1_Bio3_final <- rep(NA, seed.num)
p_Bio3_final   <- rep(NA, seed.num)
AUC_3.value  <- rep(NA, seed.num)
###################################


UseCores <- 100
#UseCores
cl <- makeCluster(UseCores)
registerDoParallel(cl)

Main_Fun <- foreach(j = 1 : seed.num, .packages= c('MASS', 'Rfast', 'bbmle', 'combinat')) %dopar%
  {
    set.seed(j)
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
    
    D.status <- rbinom(N, size = 1, prob = p)
    Dmat <- t(matrix(D.status, nrow = K, ncol = n))
    Dmat
    
    Num.case <- sum(D.status==1)
    Num.control <- sum(D.status==0)
    
    ################################# Combine ##################################
    XYvec_ini <- matrix(NA, N, 3)
    
    XYvec_ini[D.status==1, ] <- mvrnorm(Num.case, mu = mu1, Sigma = sigma_1)
    XYvec_ini[D.status==0, ] <- mvrnorm(Num.control, mu = mu0, Sigma = sigma_0)
    
    XYvec_ini_D  <- XYvec_ini[D.status==1, ]
    XYvec_ini_ND <- XYvec_ini[D.status==0, ]
    
    
    ##############################################################################
    
    
    ########### GroupM ############
    dVal <- rowSums(Dmat)
    ## The case
    groupD <- ifelse(dVal == 0, yes = 0, no = 1)
    print(sum(groupD == 1))
    
    ### For the negative group
    D0Loc <- which(groupD == 0)
    D0M1.num <- sum(rbinom(sum(groupD == 0), size = 1, prob = 1 - pi0) == 1)
    D0M1.num
    groupM <- groupD
    D0M1Loc <- sample(D0Loc, D0M1.num, replace = F)
    groupM[D0M1Loc] <- 1
    groupM
    
    SenFun <- function(pi1, K, dnum, lambda)
    {
      res = pi1 * dnum /(dnum + lambda * (K - dnum))
      return(res)
    }
    
    ### For the positive group
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
    groupM
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
    for(l in 1 : rho_leg)
    {
      XYvec <- matrix(NA, N, 2)
      
      XYvec[D.status==1, ] <- XYvec_ini_D[, mu1_local[,l]]
      XYvec[D.status==0, ] <- XYvec_ini_ND[, mu0_local[,l]]
      
      #LoglikPMeanSD_est <- function(mu0_est, mu1_est, sd0_est, sd1_est, p_est, theta, pi0, pi1, N, n, K, X_0.vec, X_1.vec, lambda, groupM)
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
        #X_0.mat <- t(matrix(X_0.vec, nrow = K, ncol = n))
        
        f1.vec  <- rep(NA, N)
        f1.vec  <- (1/(2 * pi * sd11_est * sd12_est * sqrt(1 - rho1_est^2))) * exp(-1/(2 * (1 - rho1_est^2)) * (((XYvec[ ,1] - mu11_est)/sd11_est)^2 - 2 * rho1_est * ((XYvec[ ,1] - mu11_est)/sd11_est) * ((XYvec[ ,2] - mu12_est)/sd12_est) + ((XYvec[ ,2] - mu12_est)/sd12_est)^2))
        f1.mat  <- t(matrix(f1.vec, nrow = K, ncol = n))
        #X_1.mat <- t(matrix(X_1.vec, nrow = K, ncol = n))
        
        mat.time <- dim(aa.mat)[1]
        f0.mat_final = matrix(rep(f0.mat, mat.time), n, K * mat.time)
        #print(f0.mat_final)
        f1.mat_final = matrix(rep(f1.mat, mat.time), n, K * mat.time)
        #print(f1.mat_final)
        
        aa.mat_new = t(replicate(n, as.vector(t(aa.mat))))
        # aa.mat_new
        #print(aa.mat_new)
        
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
      EstPMeanSD[j, ] <- as.vector(optim(c(mu0_comb[1,l], mu0_comb[2,l], mu1_comb[1,l], mu1_comb[2,l], sd0_comb[1,l], sd0_comb[2,l], sd1_comb[1,l], sd1_comb[2,l], rho0[l], rho1[l], p), LoglikPMeanSD_est, pi0 = pi0, pi1 = pi1, N = N, n = n, K = K, XYvec = XYvec, lambda = lambda,
                                         groupM = groupM, aa.mat = aa.mat, pi1_lik.matrix = pi1_lik.matrix)$par) ### issue!!!
      # EstPMeanSD[j, ] <- as.vector(optim(c(mu0_comb[1,l], mu0_comb[2,l], mu1_comb[1,l], mu1_comb[2,l], sd0_comb[1,l], sd0_comb[2,l], sd1_comb[1,l], sd1_comb[2,l], rho0[l], rho1[l], p), LoglikPMeanSD_est, pi0 = pi0, pi1 = pi1, N = N, n = n, K = K, XYvec = XYvec, lambda = lambda,
      #                                    groupM = groupM, aa.mat = aa.mat, pi1_lik.matrix = pi1_lik.matrix, method="L-BFGS-B",
      #                                    lower = c(mu0_comb[1,l]-0.1, mu0_comb[2,l]-0.1, mu1_comb[1,l]-0.1, mu1_comb[2,l]-0.1, sd0_comb[1,l]-0.1, sd0_comb[2,l]-0.1, sd1_comb[1,l]-0.1, sd1_comb[2,l]-0.1, rho0[l]-0.05, rho1[l]-0.05, 1e-5),
      #                                    upper = c(mu0_comb[1,l]+0.1, mu0_comb[2,l]+0.1, mu1_comb[1,l]+0.1, mu1_comb[2,l]+0.1, sd0_comb[1,l]+0.1, sd0_comb[2,l]+0.1, sd1_comb[1,l]+0.1, sd1_comb[2,l]+0.1, rho0[l]+0.05, rho1[l]+0.05, 1 - 1e-5))$par)
      p_est_final <- EstPMeanSD[j, 11]
      mu01_final <- EstPMeanSD[j, 1]
      mu02_final <- EstPMeanSD[j, 2]
      mu11_final <- EstPMeanSD[j, 3]
      mu12_final <- EstPMeanSD[j, 4]
      sd01_final <- EstPMeanSD[j, 5]
      sd02_final <- EstPMeanSD[j, 6]
      sd11_final <- EstPMeanSD[j, 7]
      sd12_final <- EstPMeanSD[j, 8]
      rho0_final <- EstPMeanSD[j, 9]
      rho1_final <- EstPMeanSD[j, 10]
      
      print(l)
      
      mu0.mat[l, (mu0_local[ ,l])] <- mu0.mat[l, (mu0_local[ ,l])] + c(mu01_final, mu02_final)
      mu1.mat[l, (mu1_local[ ,l])] <- mu1.mat[l, (mu1_local[ ,l])] + c(mu11_final, mu12_final)
      sd0.mat[l, (sd0_local[ ,l])] <- sd0.mat[l, (sd0_local[ ,l])] + c(sd01_final, sd02_final)
      sd1.mat[l, (sd1_local[ ,l])] <- sd1.mat[l, (sd1_local[ ,l])] + c(sd11_final, sd12_final)
      rho0.mat[l] <- rho0.mat[l] + rho0_final
      rho1.mat[l] <- rho1.mat[l] + rho1_final
      p.vec[l] <- p.vec[l] + p_est_final
    }
    
    #print(list("mu0.mat" = mu0.mat, "mu1.mat" = mu1.mat, "sd0.mat" = sd0.mat, "sd1.mat" = sd1.mat, "rho0.mat" = rho0.mat, "rho1.mat" = rho1.mat))
    
    mu0.final[j, ] <- apply(mu0.mat, 2, sum)/(length(mu0) - 1)
    mu1.final[j, ] <- apply(mu1.mat, 2, sum)/(length(mu1) - 1)
    sd0.final[j, ] <- apply(sd0.mat, 2, sum)/(length(sd0) - 1)
    sd1.final[j, ] <- apply(sd1.mat, 2, sum)/(length(sd1) - 1)
    rho0.final[j, ] <- rho0.mat
    rho1.final[j, ] <- rho1.mat
    p.vec.final[j] <- mean(na.omit(p.vec))
    
    sigma0.final <- sigma1.final <- matrix(NA, 3, 3)
    sigma0.final <- matrix(c(sd0.final[j,1]^2, rho0.final[j,1] * sd0.final[j,1] * sd0.final[j,2], rho0.final[j,2] * sd0.final[j,1] * sd0.final[j,3],
                             rho0.final[j,1] * sd0.final[j,1] * sd0.final[j,2], sd0.final[j,2]^2, rho0.final[j,3] * sd0.final[j,2] * sd0.final[j,3],
                             rho0.final[j,2] * sd0.final[j,1] * sd0.final[j,3], rho0.final[j,3] * sd0.final[j,2] * sd0.final[j,3], sd0.final[j,3]^2), 3, 3)
    
    sigma1.final <- matrix(c(sd1.final[j,1]^2, rho1.final[j,1] * sd1.final[j,1] * sd1.final[j,2], rho1.final[j,2] * sd1.final[j,1] * sd1.final[j,3],
                             rho1.final[j,1] * sd1.final[j,1] * sd1.final[j,2], sd1.final[j,2]^2, rho1.final[j,3] * sd1.final[j,2] * sd1.final[j,3],
                             rho1.final[j,2] * sd1.final[j,1] * sd1.final[j,3], rho1.final[j,3] * sd1.final[j,2] * sd1.final[j,3], sd1.final[j,3]^2), 3, 3)
    
    AUC.value[j] <- pnorm(sqrt(t(mu1.final[j, ] - mu0.final[j, ]) %*% solve(sigma1.final + sigma0.final) %*% (mu1.final[j, ] - mu0.final[j, ])))
    # AUC_1.value[j] <- pnorm(sqrt(t(mu1.final[j, 1] - mu0.final[j, 1]) %*% solve(sd0.final[j, 1]^2 + sd1.final[j, 1]^2) %*% (mu1.final[j, 1] - mu0.final[j, 1])))
    # AUC_1.value[j]
    # AUC_2.value[j] <- pnorm(sqrt(t(mu1.final[j, 2] - mu0.final[j, 2]) %*% solve(sd0.final[j, 2]^2 + sd1.final[j, 2]^2) %*% (mu1.final[j, 2] - mu0.final[j, 2])))
    # AUC_2.value[j]
    # AUC_3.value[j] <- pnorm(sqrt(t(mu1.final[j, 3] - mu0.final[j, 3]) %*% solve(sd0.final[j, 3]^2 + sd1.final[j, 3]^2) %*% (mu1.final[j, 3] - mu0.final[j, 3])))
    # AUC_3.value[j]
    
    
    ######################################################################## Individual ###################################################################################
    ########################## Bio_1
    XYvec_Bio1 <- rep(NA, N)
    XYvec_Bio1[D.status==1] <- rnorm(Num.case, mean = mu11, sd = sd11)
    XYvec_Bio1[D.status==0] <- rnorm(Num.control, mean = mu01, sd = sd01)
    
    Loglik_Bio1_est <- function(theta, pi0, pi1, N, n, K, XYvec, lambda, groupM, aa.mat, pi1_lik.matrix)
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
    Est_Bio1[j, ] <- as.vector(optim(c(mu01, mu11, sd01, sd11, p), Loglik_Bio1_est, pi0 = pi0, pi1 = pi1, N = N, n = n, K = K, XYvec = XYvec_Bio1, lambda = lambda,
                                     groupM = groupM, aa.mat = aa.mat, pi1_lik.matrix = pi1_lik.matrix)$par)
    
    mu0_Bio1_final[j] = Est_Bio1[j, 1]
    mu1_Bio1_final[j] = Est_Bio1[j, 2]
    sd0_Bio1_final[j] = Est_Bio1[j, 3]
    sd1_Bio1_final[j] = Est_Bio1[j, 4]
    p_Bio1_final[j]   = Est_Bio1[j, 5]
    
    AUC_1.value[j] <- pnorm(sqrt(t(mu1_Bio1_final[j] - mu0_Bio1_final[j]) %*% solve(sd0_Bio1_final[j]^2 + sd1_Bio1_final[j]^2) %*% (mu1_Bio1_final[j] - mu0_Bio1_final[j])))
    AUC_1.value[j]
    
    
    ########################## Bio_2
    XYvec_Bio2 <- rep(NA, N)
    XYvec_Bio2[D.status==1] <- rnorm(Num.case, mean = mu12, sd = sd12)
    XYvec_Bio2[D.status==0] <- rnorm(Num.control, mean = mu02, sd = sd02)
    
    Loglik_Bio2_est <- function(theta, pi0, pi1, N, n, K, XYvec, lambda, groupM, aa.mat, pi1_lik.matrix)
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
    Est_Bio2[j, ] <- as.vector(optim(c(mu02, mu12, sd02, sd12, p), Loglik_Bio2_est, pi0 = pi0, pi1 = pi1, N = N, n = n, K = K, XYvec = XYvec_Bio2, lambda = lambda,
                                     groupM = groupM, aa.mat = aa.mat, pi1_lik.matrix = pi1_lik.matrix)$par)
    
    mu0_Bio2_final[j] = Est_Bio2[j, 1]
    mu1_Bio2_final[j] = Est_Bio2[j, 2]
    sd0_Bio2_final[j] = Est_Bio2[j, 3]
    sd1_Bio2_final[j] = Est_Bio2[j, 4]
    p_Bio2_final[j]   = Est_Bio2[j, 5]
    
    AUC_2.value[j] <- pnorm(sqrt(t(mu1_Bio2_final[j] - mu0_Bio2_final[j]) %*% solve(sd0_Bio2_final[j]^2 + sd1_Bio2_final[j]^2) %*% (mu1_Bio2_final[j] - mu0_Bio2_final[j])))
    AUC_2.value[j]
    
    
    ########################## Bio_3
    XYvec_Bio3 <- rep(NA, N)
    XYvec_Bio3[D.status==1] <- rnorm(Num.case, mean = mu13, sd = sd13)
    XYvec_Bio3[D.status==0] <- rnorm(Num.control, mean = mu03, sd = sd03)
    
    Loglik_Bio3_est <- function(theta, pi0, pi1, N, n, K, XYvec, lambda, groupM, aa.mat, pi1_lik.matrix)
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
    Est_Bio3[j, ] <- as.vector(optim(c(mu03, mu13, sd03, sd13, p), Loglik_Bio3_est, pi0 = pi0, pi1 = pi1, N = N, n = n, K = K, XYvec = XYvec_Bio3, lambda = lambda,
                                     groupM = groupM, aa.mat = aa.mat, pi1_lik.matrix = pi1_lik.matrix)$par)
    
    mu0_Bio3_final[j] = Est_Bio3[j, 1]
    mu1_Bio3_final[j] = Est_Bio3[j, 2]
    sd0_Bio3_final[j] = Est_Bio3[j, 3]
    sd1_Bio3_final[j] = Est_Bio3[j, 4]
    p_Bio3_final[j]   = Est_Bio3[j, 5]
    
    AUC_3.value[j] <- pnorm(sqrt(t(mu1_Bio3_final[j] - mu0_Bio3_final[j]) %*% solve(sd0_Bio3_final[j]^2 + sd1_Bio3_final[j]^2) %*% (mu1_Bio3_final[j] - mu0_Bio3_final[j])))
    AUC_3.value[j]
    
    return(list("mu0.final"  = mu0.final[j, ],  "mu1.final" = mu1.final[j, ], 
                "sd0.final"  = sd0.final[j, ],  "sd1.final" = sd1.final[j, ], 
                "rho0.final" = rho0.final[j, ], "rho1.final" = rho1.final[j, ],
                "AUC.value"  = AUC.value[j],   
                "p.vec.final" = p.vec.final[j], ### Combine
                "mu0_Bio1_final" = mu0_Bio1_final[j], "mu1_Bio1_final" = mu1_Bio1_final[j], "sd0_Bio1_final" = sd0_Bio1_final[j], "sd1_Bio1_final" = sd1_Bio1_final[j], 
                "AUC_1.value" = AUC_1.value[j], "p_Bio1_final" = p_Bio1_final[j], ## Bio_1
                "mu0_Bio2_final" = mu0_Bio2_final[j], "mu1_Bio2_final" = mu1_Bio2_final[j], "sd0_Bio2_final" = sd0_Bio2_final[j], "sd1_Bio2_final" = sd1_Bio2_final[j], 
                "AUC_2.value" = AUC_2.value[j], "p_Bio2_final" = p_Bio2_final[j], ## Bio_2
                "mu0_Bio3_final" = mu0_Bio3_final[j], "mu1_Bio3_final" = mu1_Bio3_final[j], "sd0_Bio3_final" = sd0_Bio3_final[j], "sd1_Bio3_final" = sd1_Bio3_final[j], 
                "AUC_3.value" = AUC_3.value[j], "p_Bio3_final" = p_Bio3_final[j]  ## Bio_3
    )
    )
  }



p_est_final.vec <- rep(NA, seed.num)
AUC.value_vec <- rep(NA, seed.num)
mu0.final_mat  <- matrix(NA, seed.num, length(mu0))
mu1.final_mat  <- matrix(NA, seed.num, length(mu1))
sd0.final_mat  <- matrix(NA, seed.num, length(sd0))
sd1.final_mat  <- matrix(NA, seed.num, length(sd1))
rho0.final_mat <- matrix(NA, seed.num, length(rho0))
rho1.final_mat <- matrix(NA, seed.num, length(rho1))
############# Individual ###################
mu0_Bio1_final.vec <- rep(NA, seed.num)
mu1_Bio1_final.vec <- rep(NA, seed.num)
sd0_Bio1_final.vec <- rep(NA, seed.num)
sd1_Bio1_final.vec <- rep(NA, seed.num)
p_Bio1_final.vec <- rep(NA, seed.num)
AUC_1.value_vec <- rep(NA, seed.num)

mu0_Bio2_final.vec <- rep(NA, seed.num)
mu1_Bio2_final.vec <- rep(NA, seed.num)
sd0_Bio2_final.vec <- rep(NA, seed.num)
sd1_Bio2_final.vec <- rep(NA, seed.num)
p_Bio2_final.vec <- rep(NA, seed.num)
AUC_2.value_vec <- rep(NA, seed.num)

mu0_Bio3_final.vec <- rep(NA, seed.num)
mu1_Bio3_final.vec <- rep(NA, seed.num)
sd0_Bio3_final.vec <- rep(NA, seed.num)
sd1_Bio3_final.vec <- rep(NA, seed.num)
p_Bio3_final.vec <- rep(NA, seed.num)
AUC_3.value_vec <- rep(NA, seed.num)


for(i in 1:seed.num)
{
  p_est_final.vec[i]  <- Main_Fun[[i]]$p.vec.final
  AUC.value_vec[i]    <- Main_Fun[[i]]$AUC.value
  mu0.final_mat[i, ]  <- Main_Fun[[i]]$mu0.final
  mu1.final_mat[i, ]  <- Main_Fun[[i]]$mu1.final
  sd0.final_mat[i, ]  <- Main_Fun[[i]]$sd0.final
  sd1.final_mat[i, ]  <- Main_Fun[[i]]$sd1.final
  rho0.final_mat[i, ] <- Main_Fun[[i]]$rho0.final
  rho1.final_mat[i, ] <- Main_Fun[[i]]$rho1.final
  ## Bio_1
  mu0_Bio1_final.vec[i] <- Main_Fun[[i]]$mu0_Bio1_final
  mu1_Bio1_final.vec[i] <- Main_Fun[[i]]$mu1_Bio1_final
  sd0_Bio1_final.vec[i] <- Main_Fun[[i]]$sd0_Bio1_final
  sd1_Bio1_final.vec[i] <- Main_Fun[[i]]$sd1_Bio1_final
  p_Bio1_final.vec[i]   <- Main_Fun[[i]]$p_Bio1_final
  AUC_1.value_vec[i]    <- Main_Fun[[i]]$AUC_1.value
  ## Bio_2
  mu0_Bio2_final.vec[i] <- Main_Fun[[i]]$mu0_Bio2_final
  mu1_Bio2_final.vec[i] <- Main_Fun[[i]]$mu1_Bio2_final
  sd0_Bio2_final.vec[i] <- Main_Fun[[i]]$sd0_Bio2_final
  sd1_Bio2_final.vec[i] <- Main_Fun[[i]]$sd1_Bio2_final
  p_Bio2_final.vec[i]   <- Main_Fun[[i]]$p_Bio2_final
  AUC_2.value_vec[i]    <- Main_Fun[[i]]$AUC_2.value
  ## Bio_3
  mu0_Bio3_final.vec[i] <- Main_Fun[[i]]$mu0_Bio3_final
  mu1_Bio3_final.vec[i] <- Main_Fun[[i]]$mu1_Bio3_final
  sd0_Bio3_final.vec[i] <- Main_Fun[[i]]$sd0_Bio3_final
  sd1_Bio3_final.vec[i] <- Main_Fun[[i]]$sd1_Bio3_final
  p_Bio3_final.vec[i]   <- Main_Fun[[i]]$p_Bio3_final
  AUC_3.value_vec[i]    <- Main_Fun[[i]]$AUC_3.value
}

p.final <- mean(na.omit(p_est_final.vec))

AUC.est  <- mean(na.omit(AUC.value_vec))
AUC.bias <- AUC.est - AUC_true
AUC.var  <- var(na.omit(AUC.value_vec))
AUC.len <- length(na.omit(AUC.value_vec))
AUC.SE   <- sqrt(var(na.omit(AUC.value_vec)))
AUC.CI.left  <- AUC.est - qnorm(1 - alpha/2) * AUC.SE
AUC.CI.right <- AUC.est + qnorm(1 - alpha/2) * AUC.SE

p.bias <- mean(na.omit(p_est_final.vec)) - p
p.var  <- var(na.omit(p_est_final.vec))
p.MSE  <- mean((na.omit(p_est_final.vec) - p)^2)
p.SE   <- sqrt(var(na.omit(p_est_final.vec)))
p.CI.left  <- mean(na.omit(p_est_final.vec)) - qnorm(1 - alpha/2) * p.SE
p.CI.right <- mean(na.omit(p_est_final.vec)) + qnorm(1 - alpha/2) * p.SE
########################## Individual #######################################
## Bio_1
mu0_Bio1.final_est <- mean(na.omit(mu0_Bio1_final.vec))
mu1_Bio1.final_est <- mean(na.omit(mu1_Bio1_final.vec))
sd0_Bio1.final_est <- mean(na.omit(sd0_Bio1_final.vec))
sd1_Bio1.final_est <- mean(na.omit(sd1_Bio1_final.vec))
p_Bio1.final_est   <- mean(na.omit(p_Bio1_final.vec))
p_Bio1.final.var   <- var(na.omit(p_Bio1_final.vec))
AUC_1.est <- mean(na.omit(AUC_1.value_vec))
AUC_1.var <- var(na.omit(AUC_1.value_vec))
## Bio_2
mu0_Bio2.final_est <- mean(na.omit(mu0_Bio2_final.vec))
mu1_Bio2.final_est <- mean(na.omit(mu1_Bio2_final.vec))
sd0_Bio2.final_est <- mean(na.omit(sd0_Bio2_final.vec))
sd1_Bio2.final_est <- mean(na.omit(sd1_Bio2_final.vec))
p_Bio2.final_est   <- mean(na.omit(p_Bio2_final.vec))
p_Bio2.final.var   <- var(na.omit(p_Bio2_final.vec))
AUC_2.est <- mean(na.omit(AUC_2.value_vec))
AUC_2.var <- var(na.omit(AUC_2.value_vec))
## Bio_3
mu0_Bio3.final_est <- mean(na.omit(mu0_Bio3_final.vec))
mu1_Bio3.final_est <- mean(na.omit(mu1_Bio3_final.vec))
sd0_Bio3.final_est <- mean(na.omit(sd0_Bio3_final.vec))
sd1_Bio3.final_est <- mean(na.omit(sd1_Bio3_final.vec))
p_Bio3.final_est   <- mean(na.omit(p_Bio3_final.vec))
p_Bio3.final.var   <- var(na.omit(p_Bio3_final.vec))
AUC_3.est <- mean(na.omit(AUC_3.value_vec))
AUC_3.var <- var(na.omit(AUC_3.value_vec))

sim_stats =
  as.matrix(
    list("mu0.final_est" = apply(mu0.final_mat, 2, mean), "mu1.final_est" = apply(mu1.final_mat, 2, mean),
         "sd0.final_est" = apply(sd0.final_mat, 2, mean),  "sd1.final_est" = apply(sd1.final_mat, 2, mean),
         "rho0.final_est" = apply(rho0.final_mat, 2, mean), "rho1.final_est" = apply(rho1.final_mat, 2, mean),
         "p.final" = p.final, "p.bias" = p.bias, "p.var" = p.var, "p.MSE" = p.MSE, "p.SE" = p.SE, "p.CI.left" = p.CI.left, "p.CI.right" = p.CI.right, 
         "AUC_true" = AUC_true, "AUC.est" = AUC.est, "AUC.bias" = AUC.bias, "AUC.var" = AUC.var, 
         "AUC.len" = AUC.len, "AUC.SE" = AUC.SE, "AUC.CI.left" = AUC.CI.left, "AUC.CI.right" = AUC.CI.right, ### Combine
         ############## Individual ##############
         ## Bio_1
         "mu0_Bio1.final_est" = mu0_Bio1.final_est, "mu1_Bio1.final_est" = mu1_Bio1.final_est,
         "sd0_Bio1.final_est" = sd0_Bio1.final_est, "sd1_Bio1.final_est" = sd1_Bio1.final_est,
         "p_Bio1.final_est" = p_Bio1.final_est, "p_Bio1.final.var" = p_Bio1.final.var,
         "AUC_1.est" = AUC_1.est, "AUC_1.var" = AUC_1.var,
         ## Bio_2
         "mu0_Bio2.final_est" = mu0_Bio2.final_est, "mu1_Bio2.final_est" = mu1_Bio2.final_est,
         "sd0_Bio2.final_est" = sd0_Bio2.final_est, "sd1_Bio2.final_est" = sd1_Bio2.final_est,
         "p_Bio2.final_est" = p_Bio2.final_est, "p_Bio2.final.var" = p_Bio2.final.var,
         "AUC_2.est" = AUC_2.est, "AUC_2.var" = AUC_2.var,
         ## Bio_3
         "mu0_Bio3.final_est" = mu0_Bio3.final_est, "mu1_Bio3.final_est" = mu1_Bio3.final_est,
         "sd0_Bio3.final_est" = sd0_Bio3.final_est, "sd1_Bio3.final_est" = sd1_Bio3.final_est,
         "p_Bio3.final_est" = p_Bio3.final_est, "p_Bio3.final.var" = p_Bio3.final.var,
         "AUC_3.est" = AUC_3.est, "AUC_3.var" = AUC_3.var
    )  
  )

write.csv(sim_stats, file = paste0("simulation/Fix_sample/20000/Revised_K=", K, "_SEED=", SEED.NUM, "_p=", p, "_", pi0, "_", pi1, "_", CurrentData, ".csv"))




#write.csv(data.frame(sim_stats), "result.csv")
stopCluster(cl)

time2 <- Sys.time()
print(time2-time1)
