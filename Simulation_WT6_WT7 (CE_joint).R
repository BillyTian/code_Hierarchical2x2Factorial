##########################################################
# This R program help to reproduces 
# Web Table 6 and Web Table 7
# in the simulation results
##########################################################

setwd("...") #Need to change the directory here
source("gendata") #Should use the respective data generating function
source("CE_joint.R")

require(nlme)

#Change the effect size
#delta_x <- 0.2 # 0.5
#delta_z <- 0.1 # 0.25

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
  
  n[i] <- calc_mjoint(eff_x=delta_x, eff_z=delta_z, m_bar=m_bar, rho=rho, CV=CV)[1]
  pred.power[i] <- calc_mjoint(eff_x=delta_x, eff_z=delta_z, m_bar=m_bar, rho=rho, CV=CV)[2]
}

table <- cbind(table, n)


#function to compute empirical power or empirical type I error
empirical_mjoint <- function(nullcase=F, parameter, nsims=5000){
  
  #parameter=table[1,]
  m_bar <- as.numeric(parameter[1])
  rho <- as.numeric(parameter[2])
  CV <- as.numeric(parameter[3])
  n <- as.numeric(parameter[4])
  
  beta1 <- 1
  beta2 <- delta_x
  beta3 <- delta_z
  beta4 <- 0.05
  
  if (nullcase==T){
    beta2 <- 0
    beta3 <- 0
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
    R <- matrix(c(0,1,0,0,0,0,1,0), nrow=2, byrow=T)
    beta <- fit$coef$fixed
    #Use Wald test, which follows the chi-squared distribution with 1 degree of freedom
    test.stat <- as.numeric(t(R%*%beta) %*% solve(R%*%vcov(fit)%*%t(R)) %*% (R%*%beta))
    pvalue[i] <- 1-pchisq(test.stat,2)
    
  }
  
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  return(c(empirical, error.rate))
}

#Compute empirical power (and error rate)
empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power <- rbind(empirical.power, empirical_mjoint(parameter=table[i,]))
}

#Compute empirical type I error (and error rate)
empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_mjoint(nullcase=T, parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
names(result)[5:9] <- c("emp.tIe", "nconv.rate", "emp.power", "nconv.rate.p", "pred.power") 

mainDir = '...'
write.table(result, paste0(mainDir, "/CEjoint_x_", delta_x, "_z_", delta_z, ".csv"), sep =",", row.names=F)


##########################################
# Generate the right part of these tables
##########################################

n <- numeric(36)
pred.power <- numeric(36)
for (i in 1:nrow(table)){
  m_bar <- as.numeric(table[i,][1])
  rho <- as.numeric(table[i,][2])
  CV <- as.numeric(table[i,][3])
  
  n[i] <- calc_mjoint(eff_x=delta_x, eff_z=delta_z, m_bar=m_bar, rho=rho, CV=CV)[1]
  pred.power[i] <- calc_mjoint(eff_x=delta_x, eff_z=delta_z, m_bar=m_bar, rho=rho, CV=CV)[2]
}

table <- cbind(table, n)


#function to compute empirical power or empirical type I error
empirical_mjoint <- function(nullcase=F, parameter, nsims=5000){
  
  #parameter=table[1,]
  m_bar <- as.numeric(parameter[1])
  rho <- as.numeric(parameter[2])
  CV <- as.numeric(parameter[3])
  n <- as.numeric(parameter[4])
  
  beta1 <- 1
  beta2 <- delta_x
  beta3 <- delta_z
  beta4 <- 0.05
  
  if (nullcase==T){
    beta2 <- 0
    beta3 <- 0
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
    R <- matrix(c(0,1,0,0,0,0,1,0), nrow=2, byrow=T)
    beta <- fit$coef$fixed
    #Use Wald test, which follows the chi-squared distribution with 1 degree of freedom
    test.stat <- as.numeric(t(R%*%beta) %*% solve(R%*%vcov(fit)%*%t(R)) %*% (R%*%beta))/2
    pvalue[i] <- 1-pf(test.stat, 2, n-2)
    
  }
  
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  return(c(empirical, error.rate))
}

#Compute empirical power (and error rate)
empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power <- rbind(empirical.power, empirical_mjoint(parameter=table[i,]))
}

#Compute empirical type I error (and error rate)
empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_mjoint(nullcase=T, parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
names(result)[5:9] <- c("emp.tIe", "nconv.rate", "emp.power", "nconv.rate.p", "pred.power") 

mainDir = '...'
write.table(result, paste0(mainDir, "/CEjoint_Ftry_x_", delta_x, "_z_", delta_z, ".csv"), sep =",", row.names=F)
