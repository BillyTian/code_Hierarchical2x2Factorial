##############################################################
# Functions regarding the cluster-level controlled effect test 
# without finite-sample correction
##############################################################

n_calc_beta2 <- function(alpha=0.05, beta=0.2, sigma2y=1, eff, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  nvar2 <- sigma2y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) * 
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*nvar2/eff^2
  n <- 2*(ceiling(n/2))
  return (n)
  
}

power_calc_beta2 <- function(n, alpha=0.05, beta=0.2, sigma2y=1, eff, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  nvar2 <- sigma2y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) * 
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
  power <- pnorm( sqrt(n*eff^2/nvar2)-qnorm(1-alpha/2) )
  return(power)
  
}


##############################################################
# Functions regarding the cluster-level controlled effect test 
# with finite-sample correction
##############################################################

calc_beta2_corrected <- function(alpha=0.05, beta=0.2, sigma2y=1, eff, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  nvar2 <- sigma2y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) * 
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
  
  n <- 2
  power <- 0
  while (power < 1-beta){
    n <- n+2
    power <- pt(qt(1-alpha/2, n-2), n-2, ncp=eff/sqrt(nvar2/n), lower.tail = F) + pt(qt(alpha/2, n-2), n-2, ncp=eff/sqrt(nvar2/n))
  }
  return(c(n, power))
}