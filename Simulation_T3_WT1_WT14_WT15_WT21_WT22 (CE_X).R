##########################################################
# This R program help to reproduces 
# Table 3,  Web Table 1, Web Table 14, Web Table 15,
# Web Table 21 and Web Table 22 in the simulation results
##########################################################

setwd("...") #Need to change the directory here
source("gendata") #Should use the respective data generating function
source("CE_X.R")

require(nlme)

#Change the effect size
#delta <- 0.2 # 0.4

#########################################
# Generate the left part of these tables
#########################################

m_bar <- c(rep(20,12), rep(50,12), rep(100,12))
rho <- rep(c(rep(0.02,4), rep(0.05,4), rep(0.1,4)), 3)
CV <- rep(c(0,0.3,0.6,0.9), 9)
table <- cbind(m_bar, rho, CV)

n <- numeric(36)
pred.power <- numeric(36)
for (i in 1:nrow(table)){
  m_bar <- as.numeric(table[i,][1])
  rho <- as.numeric(table[i,][2])
  CV <- as.numeric(table[i,][3])
  
  n[i] <- n_calc_beta2(eff=delta, m_bar=m_bar, rho=rho, CV=CV)
  pred.power[i] <- power_calc_beta2(n=n[i], eff=delta, m_bar=m_bar, rho=rho, CV=CV)
}

table <- cbind(table, n)


#function to compute empirical power or empirical type I error
empirical_beta2 <- function(nullcase=F, parameter, nsims=5000){
  
  #parameter=table[1,]
  m_bar <- as.numeric(parameter[1])
  rho <- as.numeric(parameter[2])
  CV <- as.numeric(parameter[3])
  n <- as.numeric(parameter[4])
  
  beta1 <- 1
  beta2 <- delta
  beta3 <- 0.05
  beta4 <- 0.05
  
  if (nullcase==T){
    beta2 <- 0
  }
  
  pvalue <- NULL
  count <- NULL
  for (i in 1:nsims){
    set.seed(510+2021*i)
    simdata <- gendata(beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4,
                       m_bar=m_bar, rho=rho, CV=CV, n=n)
    
    fit <- try(lme(Y ~ ctrt + itrt + ctrt:itrt, data=simdata, random= ~ 1 | cluster.id), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    R <- c(0,1,0,0)
    beta <- fit$coef$fixed
    #Use Wald test, which follows the chi-squared distribution with 1 degree of freedom
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    pvalue[i] <- 1-pchisq(test.stat,1)
  }
  
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  return(c(empirical, error.rate))
}

#Compute empirical power (and error rate)
empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power <- rbind(empirical.power, empirical_beta2(parameter=table[i,]))
}

#Compute empirical type I error (and error rate)
empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_beta2(nullcase=T, parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
names(result)[5:9] <- c("emp.tIe", "nconv.rate", "emp.power", "nconv.rate.p", "pred.power") 

mainDir = '...'
write.table(result, paste0(mainDir, "/CE_X_", delta, ".csv"), sep =",", row.names=F)


##########################################
# Generate the right part of these tables
##########################################

n <- numeric(36)
pred.power <- numeric(36)
for (i in 1:nrow(table)){
  m_bar <- as.numeric(table[i,][1])
  rho <- as.numeric(table[i,][2])
  CV <- as.numeric(table[i,][3])
  
  n[i] <- calc_beta2_corrected(eff=delta, m_bar=m_bar, rho=rho, CV=CV)[1]
  pred.power[i] <- calc_beta2_corrected(eff=delta, m_bar=m_bar, rho=rho, CV=CV)[2]
}

table <- cbind(table, n)


#function to compute empirical power or empirical type I error
empirical_beta2_corrected <- function(nullcase=F, parameter, nsims=5000){
  
  #parameter=table[1,]
  m_bar <- as.numeric(parameter[1])
  rho <- as.numeric(parameter[2])
  CV <- as.numeric(parameter[3])
  n <- as.numeric(parameter[4])
  
  beta1 <- 1
  beta2 <- delta
  beta3 <- 0.05
  beta4 <- 0.05
  
  if (nullcase==T){
    beta2 <- 0
  }
  
  pvalue <- NULL
  count <- NULL
  for (i in 1:nsims){
    set.seed(510+2021*i)
    simdata <- gendata(beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4,
                       m_bar=m_bar, rho=rho, CV=CV, n=n)
    
    fit <- try(lme(Y ~ ctrt + itrt + ctrt:itrt, data=simdata, random= ~ 1 | cluster.id), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    R <- c(0,1,0,0)
    beta <- fit$coef$fixed
    #Use Wald test, which follows the chi-squared distribution with 1 degree of freedom
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    pvalue[i] <- 1-pf(test.stat, 1, n-2)
  }
  
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  return(c(empirical, error.rate))
}

#Compute empirical power (and error rate)
empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power <- rbind(empirical.power, empirical_beta2_corrected(parameter=table[i,]))
}

#Compute empirical type I error (and error rate)
empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_beta2_corrected(nullcase=T, parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
names(result)[5:9] <- c("emp.tIe", "nconv.rate", "emp.power", "nconv.rate.p", "pred.power") 

mainDir = '...'
write.table(result, paste0(mainDir, "/CE_X_corrected_", delta, ".csv"), sep =",", row.names=F)
