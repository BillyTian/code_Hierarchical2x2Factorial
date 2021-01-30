###############################################################################
# This R program reproduces Web Table 6 and Table 4 in the simulation results
###############################################################################

setwd("...") #need to change the directory here

source("IU.R")
source("IU_corrected.R")

require(nlme)

##########################################
# Left part of Web Table 6
##########################################

delta_x <- 0.2
delta_z <- 0.1

#list all combinations of parameters in a factorial sense 
m_bar <- c(rep(50,12), rep(100,12))
ICC <- rep(c(rep(0.02,4), rep(0.05,4), rep(0.1,4)), 2)
CV <- rep(c(0,0.3,0.6,0.9), 6)
IU.data <- data.frame(cbind(m_bar, ICC, CV))

#Create a function to compute predicted power through attempting series of cluster numbers, which also gives the predicted cluster number
IU.n <- function(parameter, desired.power=0.8){
  m_bar <- parameter[1]
  ICC <- parameter[2]
  CV <- parameter[3]
  a <- 0.05
  omega_x <- 4*(1+(m_bar-1)*ICC)/(m_bar*(1-CV^2*m_bar*ICC*(1-ICC)/(1+(m_bar-1)*ICC)^2))
  omega_z <- 4*(1-ICC)*(1+(m_bar-1)*ICC)^3/(m_bar*((1+(m_bar-2)*ICC)*(1+(m_bar-1)*ICC)^2+CV^2*m_bar*ICC^2*(1-ICC)))
  
  n <- 2
  for (i in 1:500){
    wmean.c <- sqrt(n)*delta_x/sqrt(omega_x)
    wmean.i <- sqrt(n)*delta_z/sqrt(omega_z)
    power <- pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F) + pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(a/2), mean = wmean.i) + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F) + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(a/2), mean = wmean.i)
    if (power>=desired.power) {break}
    n <- n+2
  }
  return(c(n, power))
}

IU.pred <- NULL
for (i in 1:nrow(IU.data)){
  IU.pred <- rbind(IU.pred, IU.n(parameter=unlist(IU.data[i,])))
}
n <- IU.pred[,1]
predicted.power <- IU.pred[,2]
IU_table <- data.frame(cbind(IU.data, n))

empirical.power <- NULL
for (l in 1:nrow(IU_table)){
  empirical.power <- rbind(empirical.power, IU(choice=1, parameter=unlist(IU_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(IU_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, IU(choice=2, parameter=unlist(IU_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

IU_table <- data.frame(cbind(IU_table, empirical.typeIerror, empirical.power, predicted.power, diff))
mainDir = '...' #need to change the directory here
write.table(IU_table, paste0(mainDir, "/IU_x_",delta_x,"_z_",delta_z,".csv"), sep =",", row.names=F)


##########################################
# Right part of Web Table 6
##########################################

#list all combinations of parameters in a factorial sense 
m_bar <- c(rep(50,12), rep(100,12))
rho <- rep(c(rep(0.02,4), rep(0.05,4), rep(0.1,4)), 2)
CV <- rep(c(0,0.3,0.6,0.9), 6)
IU.data <- data.frame(cbind(m_bar, rho, CV))

#Create a function to compute predicted power through attempting series of cluster numbers, which also gives the predicted cluster number
IU.n2 <- function(parameter, desired.power=0.8){
  m_bar <- parameter[1]
  rho <- parameter[2]
  CV <- parameter[3]
  a <- 0.05
  omega_x <- 4*(1+(m_bar-1)*rho)/(m_bar*(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2))
  omega_z <- 4*(1-rho)*(1+(m_bar-1)*rho)^3/(m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2+CV^2*m_bar*rho^2*(1-rho)))
  
  n <- 4
  for (i in 1:500){
    c.ncp <- sqrt(n)*delta_x/sqrt(omega_x)
    i.mean <- sqrt(n)*delta_z/sqrt(omega_z)
    power <- pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F) + pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(a/2), mean = i.mean) + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F) + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(a/2), mean = i.mean)
    if (power>=desired.power) {break}
    n <- n+2
  }
  return(c(n, power))
}

IU.pred <- NULL
for (i in 1:nrow(IU.data)){
  IU.pred <- rbind(IU.pred, IU.n2(parameter=unlist(IU.data[i,])))
}
n <- IU.pred[,1]
predicted.power <- IU.pred[,2]
IU_table <- data.frame(cbind(IU.data, n))

empirical.power <- NULL
for (l in 1:nrow(IU_table)){
  empirical.power <- rbind(empirical.power, IU_corrected(choice=1, parameter=unlist(IU_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(IU_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, IU_corrected(choice=2, parameter=unlist(IU_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

IU_table <- data.frame(cbind(IU_table, empirical.typeIerror, empirical.power, predicted.power, diff))
write.table(IU_table, paste0(mainDir, "/IU_x_",delta_x,"_z_",delta_z,"_corrected.csv"), sep =",", row.names=F)




##########################################
# Left part of Table 4
##########################################

delta_x <- 0.4
delta_z <- 0.2

IU.data <- data.frame(cbind(m_bar, ICC, CV))

IU.pred <- NULL
for (i in 1:nrow(IU.data)){
  IU.pred <- rbind(IU.pred, IU.n(parameter=unlist(IU.data[i,])))
}
n <- IU.pred[,1]
predicted.power <- IU.pred[,2]
IU_table <- data.frame(cbind(IU.data, n))

empirical.power <- NULL
for (l in 1:nrow(IU_table)){
  empirical.power <- rbind(empirical.power, IU(choice=1, parameter=unlist(IU_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(IU_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, IU(choice=2, parameter=unlist(IU_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

IU_table <- data.frame(cbind(IU_table, empirical.typeIerror, empirical.power, predicted.power, diff))
write.table(IU_table, paste0(mainDir, "/IU_x_",delta_x,"_z_",delta_z,".csv"), sep =",", row.names=F)


##########################################
# Right part of Table 4
##########################################

IU.data <- data.frame(cbind(m_bar, rho, CV))

IU.pred <- NULL
for (i in 1:nrow(IU.data)){
  IU.pred <- rbind(IU.pred, IU.n2(parameter=unlist(IU.data[i,])))
}
n <- IU.pred[,1]
predicted.power <- IU.pred[,2]
IU_table <- data.frame(cbind(IU.data, n))

empirical.power <- NULL
for (l in 1:nrow(IU_table)){
  empirical.power <- rbind(empirical.power, IU_corrected(choice=1, parameter=unlist(IU_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(IU_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, IU_corrected(choice=2, parameter=unlist(IU_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

IU_table <- data.frame(cbind(IU_table, empirical.typeIerror, empirical.power, predicted.power, diff))
write.table(IU_table, paste0(mainDir, "/IU_x_",delta_x,"_z_",delta_z,"_corrected.csv"), sep =",", row.names=F)
