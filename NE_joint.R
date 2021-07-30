
##############################################################
# Functions regarding the joint test of the natural effects
# without finite-sample correction
##############################################################

calc_joint <- function(alpha=0.05, beta=0.2, sigma2y=1, eff_x, eff_z, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  mbar_plus_eta2bar <- m_bar*(1-rho)/( 1+(m_bar-1)*rho )*(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2)
  omega_x <- sigma2y*(1-rho)/( mbar_plus_eta2bar*pi_x*(1-pi_x) )
  mbar_plus_eta1bar <- ( m_bar*(1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar^2*rho^2*(1-rho) )/(1+(m_bar-1)*rho)^3
  omega_z <- sigma2y*(1-rho)/( mbar_plus_eta1bar*pi_z*(1-pi_z) )
  
  n <- 0
  power <- 0
  while (power < 1-beta){
    n <- n+2
    theta <- n*(eff_x^2/omega_x + eff_z^2/omega_z)
    power <- pchisq(qchisq(1-alpha, 2), 2, ncp = theta, lower.tail = F)
  }
  return(c(n, power))
}


##############################################################
# Functions regarding the joint test of the natural effects
# with finite-sample correction
##############################################################

calc_joint_corrected <- function(alpha=0.05, beta=0.2, sigma2y=1, eff_x, eff_z, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  mbar_plus_eta2bar <- m_bar*(1-rho)/( 1+(m_bar-1)*rho )*(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2)
  omega_x <- sigma2y*(1-rho)/( mbar_plus_eta2bar*pi_x*(1-pi_x) )
  mbar_plus_eta1bar <- ( m_bar*(1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar^2*rho^2*(1-rho) )/(1+(m_bar-1)*rho)^3
  omega_z <- sigma2y*(1-rho)/( mbar_plus_eta1bar*pi_z*(1-pi_z) )
  
  n <- 2
  power <- 0
  while (power < 1-beta){
    n <- n+2
    set.seed(2021)
    
    #Simulate the mixed distribution (CENTRAL) to identify rejection region bound
    f.distn <- rf(10000, 1, n-2)
    chisq.distn <- rchisq(10000, 1)
    mix.distn <- f.distn + chisq.distn
    crt.value <- quantile(mix.distn, 1-alpha)
    
    #Simulate the mixed distribution (NONCENTRAL VERSION) for power calculation
    nc.f.distn <- rf(10000, 1, n-2, ncp = n*eff_x^2/omega_x)
    nc.chisq.distn <- rchisq(10000, 1, ncp = n*eff_z^2/omega_z)
    nc.mix.distn <- nc.f.distn + nc.chisq.distn
    
    power <- mean(nc.mix.distn > crt.value)  
  }
  return(c(n, power))
}

