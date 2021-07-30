##############################################################
# Functions regarding the I-U test of the controlled effects
# without finite-sample correction
##############################################################


calc_mIU <- function(alpha=0.05, beta=0.2, sigma2y=1, eff_x, eff_z, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  nvar2 <- sigma2y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) * 
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
  nvar3 <- sigma2y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
  ncov23 <- pi_z*nvar3
  
  power <- 0
  n <- 0
  while (power < 1-beta){
    n <- n+2
    #non-centrality parameter (include the unknown of interest, n)
    mean_W <- c(eff_x/sqrt(nvar2/n), eff_z/sqrt(nvar3/n))
    cov_W <- (ncov23/n) / (sqrt(nvar2/n)*sqrt(nvar3/n))
    sigma_W <- matrix(c(1, cov_W, cov_W, 1), nrow=2, byrow=T)

    power <- pmvnorm(lower=rep(qnorm(1-alpha/2),2), upper=rep(Inf,2), mean=mean_W, sigma=sigma_W) + 
      pmvnorm(lower=c(qnorm(1-alpha/2),-Inf), upper=c(Inf,qnorm(alpha/2)), mean=mean_W, sigma=sigma_W) +
      pmvnorm(lower=rep(-Inf,2), upper=rep(qnorm(alpha/2),2), mean=mean_W, sigma=sigma_W) +
      pmvnorm(lower=c(-Inf,qnorm(alpha/2)), upper=c(qnorm(alpha/2),Inf), mean=mean_W, sigma=sigma_W)
  }
  return(c(n, power))
}


##############################################################
# Functions regarding the I-U test of the controlled effects
# with finite-sample correction
##############################################################


calc_mIU <- function(alpha=0.05, beta=0.2, sigma2y=1, eff_x, eff_z, m_bar, rho, CV, pi_x=0.5, pi_z=0.5){
  
  nvar2 <- sigma2y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) * 
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
  nvar3 <- sigma2y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
  ncov23 <- pi_z*nvar3
  
  power <- 0
  n <- 2
  while (power < 1-beta){
    n <- n+2
    #non-centrality parameter (include the unknown of interest, n)
    mean_W <- c(eff_x/sqrt(nvar2/n), eff_z/sqrt(nvar3/n))
    cov_W <- (ncov23/n) / (sqrt(nvar2/n)*sqrt(nvar3/n))
    sigma_W <- matrix(c(1, cov_W, cov_W, 1), nrow=2, byrow=T)
    
    power <- pmvt(df=n-2, lower=c(qt(1-alpha/2, n-2), qt(1-alpha/2, n-2)), upper=rep(Inf,2), delta=mean_W, sigma=sigma_W) + 
      pmvt(df=n-2, lower=c(qt(1-alpha/2, n-2),-Inf), upper=c(Inf,qt(alpha/2, n-2)), delta=mean_W, sigma=sigma_W) +
      pmvt(df=n-2, lower=rep(-Inf,2), upper=c(qt(alpha/2, n-2), qt(alpha/2, n-2)), delta=mean_W, sigma=sigma_W) +
      pmvt(df=n-2, lower=c(-Inf,qt(alpha/2, n-2)), upper=c(qt(alpha/2, n-2),Inf), delta=mean_W, sigma=sigma_W)
  }
  return(c(n, power))
}