##################################################################
# Functions regarding the individual-level controlled effect test 
##################################################################

n_calc_beta3 <- function(alpha=0.05, beta=0.2, sigma2y=1, eff, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  nvar3 <- sigma2y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*nvar3/eff^2
  n <- 2*(ceiling(n/2))
  return (n)
  
}

power_calc_beta3 <- function(n, alpha=0.05, beta=0.2, sigma2y=1, eff, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  nvar3 <- sigma2y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
  power <- pnorm( sqrt(n*eff^2/nvar3)-qnorm(1-alpha/2) )
  return(power)
  
}
