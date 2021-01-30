#########################################################################
# This R program reproduces Table 3 in the simulation results
#########################################################################

setwd("...") #need to change the directory here

source("interaction.R")

require(nlme)

##########################################
# Left part of Table 3
##########################################

delta_xz <- 0.2

a <- 0.05
b <- 0.2
#z_{1-alpha/2}
z_a <- qnorm(1-a/2)
#z_{1-beta}
z_b <- qnorm(1-b)
#assume in the cluster number formula, sigma^2=1
const <- (z_a+z_b)^2*16*1

m_bar <- c(rep(50,12), rep(100,12))
ICC <- rep(c(rep(0.02,4), rep(0.05,4), rep(0.1,4)), 2)
CV <- rep(c(0,0.3,0.6,0.9), 6)

#cluster number for detecting cluster-level treatment effect
n <- const*(1-ICC)*(1+(m_bar-1)*ICC)^3/((m_bar*delta_xz^2)*((1+(m_bar-2)*ICC)*(1+(m_bar-1)*ICC)^2+CV^2*m_bar*ICC^2*(1-ICC)))
#rounded to the nearest even integer above
n <- ceiling(ceiling(n)/2)*2
#theoretical/predicted power
predicted.power <- pnorm(sqrt(n*m_bar*delta_xz^2*((1+(m_bar-2)*ICC)*(1+(m_bar-1)*ICC)^2+CV^2*m_bar*ICC^2*(1-ICC))/(16*(1-ICC)*(1+(m_bar-1)*ICC)^3))-z_a)
interaction_table <- data.frame(cbind(m_bar, ICC, CV, n))

empirical.power <- NULL
for (l in 1:nrow(interaction_table)){
  empirical.power <- rbind(empirical.power, interact(choice=1, parameter=unlist(interaction_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(interaction_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, interact(choice=2, parameter=unlist(interaction_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

interaction_table <- data.frame(cbind(interaction_table, empirical.typeIerror, empirical.power, predicted.power, diff))
mainDir = '...' #need to change the directory here
write.table(interaction_table, paste0(mainDir, "/interaction_",delta_xz,".csv"), sep =",", row.names=F)

##########################################
# Left part of Table 3
##########################################

delta_xz <- 0.3

#cluster number for detecting cluster-level treatment effect
n <- const*(1-ICC)*(1+(m_bar-1)*ICC)^3/((m_bar*delta_xz^2)*((1+(m_bar-2)*ICC)*(1+(m_bar-1)*ICC)^2+CV^2*m_bar*ICC^2*(1-ICC)))
#rounded to the nearest even integer above
n <- ceiling(ceiling(n)/2)*2
#theoretical/predicted power
predicted.power <- pnorm(sqrt(n*m_bar*delta_xz^2*((1+(m_bar-2)*ICC)*(1+(m_bar-1)*ICC)^2+CV^2*m_bar*ICC^2*(1-ICC))/(16*(1-ICC)*(1+(m_bar-1)*ICC)^3))-z_a)
interaction_table <- data.frame(cbind(m_bar, ICC, CV, n))

empirical.power <- NULL
for (l in 1:nrow(interaction_table)){
  empirical.power <- rbind(empirical.power, interact(choice=1, parameter=unlist(interaction_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(interaction_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, interact(choice=2, parameter=unlist(interaction_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

interaction_table <- data.frame(cbind(interaction_table, empirical.typeIerror, empirical.power, predicted.power, diff))
write.table(interaction_table, paste0(mainDir, "/interaction_",delta_xz,".csv"), sep =",", row.names=F)