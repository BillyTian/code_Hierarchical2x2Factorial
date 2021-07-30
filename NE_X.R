##############################################################
# Functions regarding the cluster-level natural effect test 
# without finite-sample correction
##############################################################

calcn_mcluster <- function(alpha=0.05, beta=0.2, sigma2y=1, eff, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  mbar_plus_eta2bar <- m_bar*(1-rho)/( 1+(m_bar-1)*rho )*(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2)
  omega_x <- sigma2y*(1-rho)/( mbar_plus_eta2bar*pi_x*(1-pi_x) )
  
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*omega_x/eff^2
  n <- 2*(ceiling(n/2))
  return (n)
}


calcpower_mcluster <- function(n, alpha=0.05, beta=0.2, sigma2y=1, eff, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  mbar_plus_eta2bar <- m_bar*(1-rho)/( 1+(m_bar-1)*rho )*(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2)
  omega_x <- sigma2y*(1-rho)/( mbar_plus_eta2bar*pi_x*(1-pi_x) )  
  
  power <- pnorm( sqrt(n*eff^2/omega_x)-qnorm(1-alpha/2) )
  return(power)
  
}


##############################################################
# Functions regarding the cluster-level natural effect test 
# with finite-sample correction
##############################################################


calc_mcluster_corrected <- function(alpha=0.05, beta=0.2, sigma2y=1, eff, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  mbar_plus_eta2bar <- m_bar*(1-rho)/( 1+(m_bar-1)*rho )*(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2)
  omega_x <- sigma2y*(1-rho)/( mbar_plus_eta2bar*pi_x*(1-pi_x) )
  
  n <- 2
  power <- 0
  while (power < 1-beta){
    n <- n+2
    power <- pt(qt(1-alpha/2, n-2), n-2, ncp=eff/sqrt(omega_x/n), lower.tail = F) + pt(qt(alpha/2, n-2), n-2, ncp=eff/sqrt(omega_x/n))
  }
  return(c(n, power))
}