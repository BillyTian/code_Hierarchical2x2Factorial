##################################################################
# Functions regarding the individual-level natural effect test 
##################################################################

calcn_mind <- function(alpha=0.05, beta=0.2, sigma2y=1, eff, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  mbar_plus_eta1bar <- ( m_bar*(1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar^2*rho^2*(1-rho) )/(1+(m_bar-1)*rho)^3
  omega_z <- sigma2y*(1-rho)/( mbar_plus_eta1bar*pi_z*(1-pi_z) )
  
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*omega_z/eff^2
  n <- 2*(ceiling(n/2))
  return (n)
}


calcpower_mind <- function(n, alpha=0.05, beta=0.2, sigma2y=1, eff, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  mbar_plus_eta1bar <- ( m_bar*(1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar^2*rho^2*(1-rho) )/(1+(m_bar-1)*rho)^3
  omega_z <- sigma2y*(1-rho)/( mbar_plus_eta1bar*pi_z*(1-pi_z) )
  
  power <- pnorm( sqrt(n*eff^2/omega_z)-qnorm(1-alpha/2) )
  return(power)
  
}


