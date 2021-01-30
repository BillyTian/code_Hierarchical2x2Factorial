#########################################################################
# This R program reproduces Web Table 1 and 2 in the simulation results
#########################################################################

setwd("...") #need to change the directory here

source("marginal_cluster.R")
source("marginal_cluster_corrected.R")

require(nlme)

##########################################
# Left part of Web Table 1
##########################################

delta_x <- 0.2

a <- 0.05
b <- 0.2
#z_{1-alpha/2}
z_a <- qnorm(1-a/2)
#z_{1-beta}
z_b <- qnorm(1-b)

#create data set for the information combinations
m_bar <- c(rep(50,12), rep(100,12))
rho <- rep(c(rep(0.02,4), rep(0.05,4), rep(0.1,4)), 2)
CV <- rep(c(0,0.3,0.6,0.9), 6)

omega_x <- 4*1*(1+(m_bar-1)*rho)/m_bar/(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2)

#cluster number for detecting cluster-level treatment effect
n <- (z_a+z_b)^2*omega_x/delta_x^2
#rounded to the nearest even integer above
n <- ceiling(ceiling(n)/2)*2
#theoretical/predicted power
predicted.power <- pnorm(sqrt(n*delta_x^2/omega_x)-z_a)
marginal.cluster_table <- data.frame(cbind(m_bar, rho, CV, n))

empirical.power <- NULL
for (l in 1:nrow(marginal.cluster_table)){
  empirical.power <- rbind(empirical.power, marginal.cluster(choice=1, parameter=unlist(marginal.cluster_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(marginal.cluster_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, marginal.cluster(choice=2, parameter=unlist(marginal.cluster_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

marginal.cluster_table <- data.frame(cbind(marginal.cluster_table, empirical.typeIerror, empirical.power, predicted.power, diff))
mainDir = '...' #need to change directory here
write.table(marginal.cluster_table, paste0(mainDir, "/marginal.cluster_",delta_x,".csv"), sep =",", row.names=F)

###########################################
# Right part of Web Table 1
###########################################

c.data <- data.frame(cbind(m_bar, rho, CV))

cluster.n <- function(parameter, desired.power=0.8){
  m_bar <- parameter[1]
  rho <- parameter[2]
  CV <- parameter[3]
  a <- 0.05
  omega_x <- 4*(1+(m_bar-1)*rho)/m_bar/(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2)
  
  n <- 4
  for (i in 1:20000){
    power <- pt(qt(1-a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n), lower.tail = F) + pt(qt(a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n))
    if (power>=desired.power) {break}
    n <- n+2
  }
  return(c(n, power))
}

cluster.pred <- NULL
for (i in 1:nrow(c.data)){
  cluster.pred <- rbind(cluster.pred, cluster.n(parameter=unlist(c.data[i,])))
}
n <- cluster.pred[,1]
predicted.power <- cluster.pred[,2]
marginal.cluster_table <- data.frame(cbind(c.data, n))

empirical.power <- NULL
for (l in 1:nrow(marginal.cluster_table)){
  empirical.power <- rbind(empirical.power, marginal.cluster_corrected(choice=1, parameter=unlist(marginal.cluster_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(marginal.cluster_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, marginal.cluster_corrected(choice=2, parameter=unlist(marginal.cluster_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

marginal.cluster_table <- data.frame(cbind(marginal.cluster_table, empirical.typeIerror, empirical.power, predicted.power, diff))
write.table(marginal.cluster_table, paste0(mainDir, "/marginal.cluster_",delta_x,"_corrected.csv"), sep =",", row.names=F)


##########################################
# Left part of Web Table 2
##########################################

delta_x <- 0.4

omega_x <- 4*1*(1+(m_bar-1)*rho)/m_bar/(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2)

#cluster number for detecting cluster-level treatment effect
n <- (z_a+z_b)^2*omega_x/delta_x^2
#rounded to the nearest even integer above
n <- ceiling(ceiling(n)/2)*2
#theoretical/predicted power
predicted.power <- pnorm(sqrt(n*delta_x^2/omega_x)-z_a)
marginal.cluster_table <- data.frame(cbind(m_bar, rho, CV, n))

empirical.power <- NULL
for (l in 1:nrow(marginal.cluster_table)){
  empirical.power <- rbind(empirical.power, marginal.cluster(choice=1, parameter=unlist(marginal.cluster_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(marginal.cluster_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, marginal.cluster(choice=2, parameter=unlist(marginal.cluster_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

marginal.cluster_table <- data.frame(cbind(marginal.cluster_table, empirical.typeIerror, empirical.power, predicted.power, diff))
write.table(marginal.cluster_table, paste0(mainDir, "/marginal.cluster_",delta_x,".csv"), sep =",", row.names=F)

###########################################
# Right part of Web Table 2
###########################################

cluster.pred <- NULL
for (i in 1:nrow(c.data)){
  cluster.pred <- rbind(cluster.pred, cluster.n(parameter=unlist(c.data[i,])))
}
n <- cluster.pred[,1]
predicted.power <- cluster.pred[,2]
marginal.cluster_table <- data.frame(cbind(c.data, n))

empirical.power <- NULL
for (l in 1:nrow(marginal.cluster_table)){
  empirical.power <- rbind(empirical.power, marginal.cluster_corrected(choice=1, parameter=unlist(marginal.cluster_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(marginal.cluster_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, marginal.cluster_corrected(choice=2, parameter=unlist(marginal.cluster_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

marginal.cluster_table <- data.frame(cbind(marginal.cluster_table, empirical.typeIerror, empirical.power, predicted.power, diff))
write.table(marginal.cluster_table, paste0(mainDir, "/marginal.cluster_",delta_x,"_corrected.csv"), sep =",", row.names=F)
