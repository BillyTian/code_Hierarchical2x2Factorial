##############################################################
# Functions regarding the joint test of the controlled effects
# without finite-sample correction
##############################################################

calc_mjoint <- function(alpha=0.05, beta=0.2, sigma2y=1, eff_x, eff_z, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  nvar2 <- sigma2y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) * 
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
  nvar3 <- sigma2y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
  ncov23 <- pi_z*nvar3
  
  power <- 0
  n <- 0
  while (power < 1-beta){
    n <- n+2
    #non-centrality parameter (include the unknown of interest, n)
    beta23 <- c(eff_x, eff_z)
    omega <- matrix(c(nvar2, ncov23, ncov23, nvar3), nrow=2, byrow=T)
    
    theta <- n*t(beta23) %*% solve(omega) %*% beta23
    power <- pchisq(qchisq(1-alpha, 2), 2, ncp = theta, lower.tail = F)

  }
  return(c(n, power))
}


##############################################################
# Functions regarding the joint test of the controlled effects 
# with finite-sample correction
##############################################################

calc_mjoint <- function(alpha=0.05, beta=0.2, sigma2y=1, eff_x, eff_z, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  nvar2 <- sigma2y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) * 
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
  nvar3 <- sigma2y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
  ncov23 <- pi_z*nvar3
  
  power <- 0
  n <- 2
  while (power < 1-beta){
    n <- n+2
    #non-centrality parameter (include the unknown of interest, n)
    beta23 <- c(eff_x, eff_z)
    omega <- matrix(c(nvar2, ncov23, ncov23, nvar3), nrow=2, byrow=T)
    
    theta <- n*t(beta23) %*% solve(omega) %*% beta23
    power <- pf(qf(1-alpha, 2, n-2), 2, n-2, ncp = theta, lower.tail = F)
    
  }
  return(c(n, power))
}
