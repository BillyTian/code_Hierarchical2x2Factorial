##########################################################
# This R program help to reproduces 
# Web Table 4,  Web Table 19, and Web Table 26
# in the simulation results
##########################################################

setwd("...") #Need to change the directory here
source("gendata") #Should use the respective data generating function
source("NE_Z.R")

require(nlme)

#Change the effect size
#delta_z <- 0.1 # 0.15

########################
# Generate these tables
########################

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
  
  n[i] <- calcn_mind(eff=delta_z, m_bar=m_bar, rho=rho, CV=CV)
  pred.power[i] <- calcpower_mind(n=n[i], eff=delta_z, m_bar=m_bar, rho=rho, CV=CV)
}

table <- cbind(table, n)


#function to compute empirical power or empirical type I error
empirical_mind <- function(nullcase=F, parameter, nsims=5000){
  
  #parameter=table[1,]
  m_bar <- as.numeric(parameter[1])
  rho <- as.numeric(parameter[2])
  CV <- as.numeric(parameter[3])
  n <- as.numeric(parameter[4])
  
  if (nullcase==F){
    beta1 <- 1
    beta2 <- 0.15
    beta3 <- 0.05
    beta4 <- 2*(delta_z-beta3)
  } else if (nullcase==T){ 
    beta1 <- 1
    beta2 <- 0.15
    beta3 <- 0.05
    beta4 <- -2*beta3
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
    R <- c(0,0,1,0.5)
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
  empirical.power <- rbind(empirical.power, empirical_mind(parameter=table[i,]))
}

#Compute empirical type I error (and error rate)
empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_mind(nullcase=T, parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
names(result)[5:8] <- c("emp.tIe", "nconv.rate", "emp.power", "nconv.rate.p", "pred.power") 

mainDir = '...'
write.table(result, paste0(mainDir, "/NE_Z_",delta_z,".csv"), sep =",", row.names=F)

