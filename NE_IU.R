##############################################################
# Functions regarding the I-U test of the natural effects
# without finite-sample correction
##############################################################


calc_IU <- function(alpha=0.05, beta=0.2, sigma2y=1, eff_x, eff_z, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  mbar_plus_eta2bar <- m_bar*(1-rho)/( 1+(m_bar-1)*rho )*(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2)
  omega_x <- sigma2y*(1-rho)/( mbar_plus_eta2bar*pi_x*(1-pi_x) )
  mbar_plus_eta1bar <- ( m_bar*(1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar^2*rho^2*(1-rho) )/(1+(m_bar-1)*rho)^3
  omega_z <- sigma2y*(1-rho)/( mbar_plus_eta1bar*pi_z*(1-pi_z) )
  
  n <- 0
  power <- 0
  while (power < 1-beta){
    n <- n+2
    wmean.c <- sqrt(n)*eff_x/sqrt(omega_x)
    wmean.i <- sqrt(n)*eff_z/sqrt(omega_z)
    power <- pnorm(qnorm(1-alpha/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(1-alpha/2), mean = wmean.i, lower.tail = F) + 
      pnorm(qnorm(1-alpha/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(alpha/2), mean = wmean.i) + 
      pnorm(qnorm(alpha/2), mean = wmean.c)*pnorm(qnorm(1-alpha/2), mean = wmean.i, lower.tail = F) + 
      pnorm(qnorm(alpha/2), mean = wmean.c)*pnorm(qnorm(alpha/2), mean = wmean.i)
  }
  return(c(n, power))
}


##############################################################
# Functions regarding the I-U test of the natural effects
# with finite-sample correction
##############################################################


calc_IU_corrected <- function(alpha=0.05, beta=0.2, sigma2y=1, eff_x, eff_z, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  mbar_plus_eta2bar <- m_bar*(1-rho)/( 1+(m_bar-1)*rho )*(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2)
  omega_x <- sigma2y*(1-rho)/( mbar_plus_eta2bar*pi_x*(1-pi_x) )
  mbar_plus_eta1bar <- ( m_bar*(1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar^2*rho^2*(1-rho) )/(1+(m_bar-1)*rho)^3
  omega_z <- sigma2y*(1-rho)/( mbar_plus_eta1bar*pi_z*(1-pi_z) )
  
  n <- 2
  power <- 0
  while (power < 1-beta){
    n <- n+2
    c.ncp <- sqrt(n)*delta_x/sqrt(omega_x)
    i.mean <- sqrt(n)*delta_z/sqrt(omega_z)
    power <- pt(qt(1-alpha/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(1-alpha/2), mean = i.mean, lower.tail = F) + 
      pt(qt(1-alpha/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(alpha/2), mean = i.mean) + 
      pt(qt(alpha/2, n-2), n-2, c.ncp)*pnorm(qnorm(1-alpha/2), mean = i.mean, lower.tail = F) + 
      pt(qt(alpha/2, n-2), n-2, c.ncp)*pnorm(qnorm(alpha/2), mean = i.mean)
  }
  return(c(n, power))
}