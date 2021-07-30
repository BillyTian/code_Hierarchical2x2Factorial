##########################################################
# This R program help to reproduces 
# Web Table 12 and Web Table 13
# in the simulation results
##########################################################

setwd("...") #Need to change the directory here
source("gendata") #Should use the respective data generating function
source("NE_IU.R")

require(nlme)

#Change the effect size
#delta_x <- 0.2 # 0.4
#delta_z <- 0.1 # 0.2

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
  
  n[i] <- calc_IU(eff_x=delta_x, eff_z=delta_z, m_bar=m_bar, rho=rho, CV=CV)[1]
  pred.power[i] <- calc_IU(eff_x=delta_x, eff_z=delta_z, m_bar=m_bar, rho=rho, CV=CV)[2]
}

table <- cbind(table, n)


#function to compute empirical power or empirical type I error
empirical_IU <- function(nullcase=F, parameter, nsims=5000){
  
  #parameter=table[1,]
  m_bar <- as.numeric(parameter[1])
  rho <- as.numeric(parameter[2])
  CV <- as.numeric(parameter[3])
  n <- as.numeric(parameter[4])
  
  if (nullcase==T){
    delta_x <- 0
  }
  
  #set true regression parameters
  beta1 <- 1
  beta2 <- 0.15
  beta4 <- 2*(delta_x-beta2)
  beta3 <- delta_z-0.5*beta4
  
  p.c <- c(NA, nsims)
  p.i <- c(NA, nsims)
  count <- NULL
  for (i in 1:nsims){
    set.seed(510+2021*i)
    simdata <- gendata(beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4,
                       m_bar=m_bar, rho=rho, CV=CV, n=n)
    
    fit <- try(lme(Y ~ ctrt + itrt + ctrt:itrt, data=simdata, random= ~ 1 | cluster.id), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    A.c <- c(0,1,0,0.5)
    A.i <- c(0,0,1,0.5)
    
    #extract parameter estimates
    beta <- fit$coef$fixed
    
    #compute 2 test statistics independently
    test.stat.c <- as.numeric((t(A.c)%*%beta)^2/(t(A.c)%*%vcov(fit)%*%A.c))
    test.stat.i <- as.numeric((t(A.i)%*%beta)^2/(t(A.i)%*%vcov(fit)%*%A.i))
    
    #compute 2 p-values independently, require both significance simultaneously
    p.c[i] <- 1-pchisq(test.stat.c,1)
    p.i[i] <- 1-pchisq(test.stat.i,1)
  }
  
  empirical <- mean(p.c<0.05 & p.i<0.05, na.rm = T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  return(c(empirical, error.rate))
}

#Compute empirical power (and error rate)
empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power <- rbind(empirical.power, empirical_IU(parameter=table[i,]))
}

#Compute empirical type I error (and error rate)
empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_IU(nullcase=T, parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
names(result)[5:9] <- c("emp.tIe", "nconv.rate", "emp.power", "nconv.rate.p", "pred.power") 

mainDir = '...'
write.table(result, paste0(mainDir, "/NEIU_x_",delta_x,"_z_",delta_z,".csv"), sep =",", row.names=F)


#########################################
# Generate the left part of these tables
#########################################

n <- numeric(36)
pred.power <- numeric(36)
for (i in 1:nrow(table)){
  m_bar <- as.numeric(table[i,][1])
  rho <- as.numeric(table[i,][2])
  CV <- as.numeric(table[i,][3])
  
  n[i] <- calc_IU_corrected(eff_x=delta_x, eff_z=delta_z, m_bar=m_bar, rho=rho, CV=CV)[1]
  pred.power[i] <- calc_IU_corrected(eff_x=delta_x, eff_z=delta_z, m_bar=m_bar, rho=rho, CV=CV)[2]
}

table <- cbind(table, n)


#function to compute empirical power or empirical type I error
empirical_IU_corrected <- function(nullcase=F, parameter, nsims=5000){
  
  #parameter=table[1,]
  m_bar <- as.numeric(parameter[1])
  rho <- as.numeric(parameter[2])
  CV <- as.numeric(parameter[3])
  n <- as.numeric(parameter[4])
  
  if (nullcase==T){
    delta_x <- 0
  }
  
  #set true regression parameters
  beta1 <- 1
  beta2 <- 0.15
  beta4 <- 2*(delta_x-beta2)
  beta3 <- delta_z-0.5*beta4
  
  p.c <- c(NA, nsims)
  p.i <- c(NA, nsims)
  count <- NULL
  for (i in 1:nsims){
    set.seed(510+2021*i)
    simdata <- gendata(beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4,
                       m_bar=m_bar, rho=rho, CV=CV, n=n)
    
    fit <- try(lme(Y ~ ctrt + itrt + ctrt:itrt, data=simdata, random= ~ 1 | cluster.id), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    A.c <- c(0,1,0,0.5)
    A.i <- c(0,0,1,0.5)
    
    #extract parameter estimates
    beta <- fit$coef$fixed
    
    #compute 2 test statistics independently
    test.stat.c <- as.numeric((t(A.c)%*%beta)^2/(t(A.c)%*%vcov(fit)%*%A.c))
    test.stat.i <- as.numeric((t(A.i)%*%beta)^2/(t(A.i)%*%vcov(fit)%*%A.i))
    
    #compute 2 p-values independently, require both significance simultaneously
    p.c[i] <- 1-pf(test.stat.c,1,n-2)
    p.i[i] <- 1-pchisq(test.stat.i,1)
  }
  
  empirical <- mean(p.c<0.05 & p.i<0.05, na.rm = T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  return(c(empirical, error.rate))
}

#Compute empirical power (and error rate)
empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power <- rbind(empirical.power, empirical_IU_corrected(parameter=table[i,]))
}

#Compute empirical type I error (and error rate)
empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_IU_corrected(nullcase=T, parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
names(result)[5:9] <- c("emp.tIe", "nconv.rate", "emp.power", "nconv.rate.p", "pred.power") 

mainDir = '...'
write.table(result, paste0(mainDir, "/NEIU_corrected_x_",delta_x,"_z_",delta_z,".csv"), sep =",", row.names=F)

