#########################################################################
# This R program reproduces Web Table 3 in the simulation results
#########################################################################

setwd("...") #need to change the directory here

source("marginal_ind.R")

require(nlme)

##########################################
# Left part of Web Table 3
##########################################

delta_z <- 0.1

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

lambda <- m_bar/(m_bar+(1-rho)/rho)
omega_z <- 4*(1-rho)/(m_bar-lambda*(1-(1-lambda)*lambda*CV^2))

#cluster number for detecting cluster-level treatment effect
n <- (z_a+z_b)^2*omega_z/delta_z^2
#rounded to the nearest even integer above
n <- ceiling(ceiling(n)/2)*2
#theoretical/predicted power
predicted.power <- pnorm(sqrt(n*delta_z^2/omega_z)-z_a)

marginal.ind_table <- data.frame(cbind(m_bar, rho, CV, n))

empirical.power <- NULL
for (l in 1:nrow(marginal.ind_table)){
  empirical.power <- rbind(empirical.power, marginal.ind(choice=1, parameter=unlist(marginal.ind_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(marginal.ind_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, marginal.ind(choice=2, parameter=unlist(marginal.ind_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

marginal.ind_table <- data.frame(cbind(marginal.ind_table, empirical.typeIerror, empirical.power, predicted.power, diff))
mainDir = '...' #need to change directory here
write.table(marginal.ind_table, paste0(mainDir, "/marginal.ind_",delta_z,".csv"), sep =",", row.names=F)

##########################################
# Right part of Web Table 3
##########################################

delta_z <- 0.15

lambda <- m_bar/(m_bar+(1-rho)/rho)
omega_z <- 4*(1-rho)/(m_bar-lambda*(1-(1-lambda)*lambda*CV^2))

#cluster number for detecting cluster-level treatment effect
n <- (z_a+z_b)^2*omega_z/delta_z^2
#rounded to the nearest even integer above
n <- ceiling(ceiling(n)/2)*2

#theoretical/predicted power
predicted.power <- pnorm(sqrt(n*delta_z^2/omega_z)-z_a)

marginal.ind_table <- data.frame(cbind(m_bar, rho, CV, n))

empirical.power <- NULL
for (l in 1:nrow(marginal.ind_table)){
  empirical.power <- rbind(empirical.power, marginal.ind(choice=1, parameter=unlist(marginal.ind_table[l,])))
}
diff <- predicted.power-empirical.power

empirical.typeIerror <- NULL
for (l in 1:nrow(marginal.ind_table)){
  empirical.typeIerror <- rbind(empirical.typeIerror, marginal.ind(choice=2, parameter=unlist(marginal.ind_table[l,])))
}
empirical.typeIerror <- round(empirical.typeIerror, 2)
empirical.power <- round(empirical.power, 2)
predicted.power <- round(predicted.power, 2)

marginal.ind_table <- data.frame(cbind(marginal.ind_table, empirical.typeIerror, empirical.power, predicted.power, diff))
write.table(marginal.ind_table, paste0(mainDir, "/marginal.ind_",delta_z,".csv"), sep =",", row.names=F)