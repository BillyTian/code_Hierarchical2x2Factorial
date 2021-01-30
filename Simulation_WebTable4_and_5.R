#########################################################################
# This R program reproduces Web Table 4 and 5 in the simulation results
#########################################################################

setwd("...") #need to change the directory here

source("joint.R")
source("joint_corrected.R")

require(nlme)
require(Matrix)

##########################################
# Left part of Web Table 4
##########################################

delta_x <- 0.2
delta_z <- 0.1

m_bar <- c(rep(50,12), rep(100,12))
ICC <- rep(c(rep(0.02,4), rep(0.05,4), rep(0.1,4)), 2)
CV <- rep(c(0,0.3,0.6,0.9), 6)
joint.data <- data.frame(cbind(m_bar, ICC, CV))

joint.n <- function(parameter, desired.power=0.8){
  m_bar <- parameter[1]
  ICC <- parameter[2]
  CV <- parameter[3]
  a <- 0.05
  omega_x <- 4*(1+(m_bar-1)*ICC)/(m_bar*(1-CV^2*m_bar*ICC*(1-ICC)/(1+(m_bar-1)*ICC)^2))
  omega_z <- 4*(1-ICC)*(1+(m_bar-1)*ICC)^3/(m_bar*((1+(m_bar-2)*ICC)*(1+(m_bar-1)*ICC)^2+CV^2*m_bar*ICC^2*(1-ICC)))
  
  n <- 2
  for (i in 1:500){
    #non-centrality parameter (include the unknown of interest, n)
    theta <- n*(delta_x^2/omega_x + delta_z^2/omega_z)
    power <- pchisq(qchisq(1-a, 2), 2, ncp = theta, lower.tail = F)
    
    if (power>=desired.power) {break}
    n <- n+2
  }
  return(c(n, power))
}

joint.pred <- NULL
for (i in 1:nrow(joint.data)){
  joint.pred <- rbind(joint.pred, joint.n(parameter=unlist(joint.data[i,])))
}
n <- joint.pred[,1]
predicted.power <- joint.pred[,2]
joint_table <- data.frame(cbind(joint.data, n))

empirical.power <- NULL
for (l in 1:nrow(joint_table)){
  empirical.power <- rbind(empirical.power, joint(choice=1, parameter=unlist(joint_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(joint_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, joint(choice=2, parameter=unlist(joint_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

joint_table <- data.frame(cbind(joint_table, empirical.typeIerror, empirical.power, predicted.power, diff))
mainDir = '...' #need to change the directory here
write.table(joint_table, paste0(mainDir, "/joint_x_",delta_x,"_z_",delta_z,".csv"), sep =",", row.names=F)


##########################################
# Right part of Web Table 4
##########################################

j.data <- data.frame(cbind(m_bar, rho, CV))

joint.n2 <- function(parameter, desired.power=0.8){
  m_bar <- parameter[1]
  rho <- parameter[2]
  CV <- parameter[3]
  a <- 0.05
  omega_x <- 4*(1+(m_bar-1)*rho)/(m_bar*(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2))
  omega_z <- 4*(1-rho)*(1+(m_bar-1)*rho)^3/(m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2+CV^2*m_bar*rho^2*(1-rho)))
  
  n <- 4
  for (i in 1:500){
    set.seed(20201231)
    
    #Simulate the mixed distribution (CENTRAL) to identify rejection region bound
    f.distn <- rf(10000, 1, n-2)
    chisq.distn <- rchisq(10000, 1)
    mix.distn <- f.distn + chisq.distn
    crt.value <- quantile(mix.distn, 1-a)
    
    #Simulate the mixed distribution (NONCENTRAL VERSION) for power calculation
    nc.f.distn <- rf(10000, 1, n-2, ncp = n*delta_x^2/omega_x)
    nc.chisq.distn <- rchisq(10000, 1, ncp = n*delta_z^2/omega_z)
    nc.mix.distn <- nc.f.distn + nc.chisq.distn
    
    power <- mean(nc.mix.distn>crt.value)
    
    if (power>=desired.power) {break}
    n <- n+2
  }
  return(c(n, power))
}

joint.pred <- NULL
for (i in 1:nrow(j.data)){
  joint.pred <- rbind(joint.pred, joint.n2(parameter=unlist(j.data[i,])))
}
n <- joint.pred[,1]
predicted.power <- joint.pred[,2]
joint_table <- data.frame(cbind(j.data, n))

empirical.power <- NULL
for (l in 1:nrow(joint_table)){
  empirical.power <- rbind(empirical.power, joint_corrected(choice=1, parameter=unlist(joint_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(joint_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, joint_corrected(choice=2, parameter=unlist(joint_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

joint_table <- data.frame(cbind(joint_table, empirical.typeIerror, empirical.power, predicted.power, diff))
write.table(joint_table, paste0(mainDir, "/joint_x_",delta_x,"_z_",delta_z,"_corrected.csv"), sep =",", row.names=F)


##########################################
# Left part of Web Table 5
##########################################

delta_x <- 0.25
delta_z <- 0.15

m_bar <- c(rep(50,12), rep(100,12))
ICC <- rep(c(rep(0.02,4), rep(0.05,4), rep(0.1,4)), 2)
CV <- rep(c(0,0.3,0.6,0.9), 6)
joint.data <- data.frame(cbind(m_bar, ICC, CV))

joint.pred <- NULL
for (i in 1:nrow(joint.data)){
  joint.pred <- rbind(joint.pred, joint.n(parameter=unlist(joint.data[i,])))
}
n <- joint.pred[,1]
predicted.power <- joint.pred[,2]
joint_table <- data.frame(cbind(joint.data, n))

empirical.power <- NULL
for (l in 1:nrow(joint_table)){
  empirical.power <- rbind(empirical.power, joint(choice=1, parameter=unlist(joint_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(joint_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, joint(choice=2, parameter=unlist(joint_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

joint_table <- data.frame(cbind(joint_table, empirical.typeIerror, empirical.power, predicted.power, diff))
write.table(joint_table, paste0(mainDir, "/joint_x_",delta_x,"_z_",delta_z,".csv"), sep =",", row.names=F)


##########################################
# Right part of Web Table 5
##########################################

j.data <- data.frame(cbind(m_bar, rho, CV))

joint.pred <- NULL
for (i in 1:nrow(j.data)){
  joint.pred <- rbind(joint.pred, joint.n2(parameter=unlist(j.data[i,])))
}
n <- joint.pred[,1]
predicted.power <- joint.pred[,2]
joint_table <- data.frame(cbind(j.data, n))

empirical.power <- NULL
for (l in 1:nrow(joint_table)){
  empirical.power <- rbind(empirical.power, joint_corrected(choice=1, parameter=unlist(joint_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(joint_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, joint_corrected(choice=2, parameter=unlist(joint_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

joint_table <- data.frame(cbind(joint_table, empirical.typeIerror, empirical.power, predicted.power, diff))
write.table(joint_table, paste0(mainDir, "/joint_x_",delta_x,"_z_",delta_z,"_corrected.csv"), sep =",", row.names=F)