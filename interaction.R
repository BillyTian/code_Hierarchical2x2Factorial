##########################################################################################
# Function for calculating empirical power or empirical type I error via 5000 simulations
# 
# for interaction test
##########################################################################################

interact <- function(choice, parameter, nsims=5000){
  
  #number of iterations
  nsims <- nsims
  #cluster size
  m_bar <- parameter[1]
  #outcome ICC
  rho <- parameter[2]
  #Coefficient of variation
  cv <- parameter[3]
  #number of clusters
  n <- parameter[4]
  
  if (choice==1){
    #set assumed parameter estimates for generating data 
    beta1 <- 1
    beta2 <- 0.15
    beta3 <- 0.05
    beta4 <- delta_xz #effect size
  } else if (choice==2){
    beta1 <- 1
    beta2 <- 0.15
    beta3 <- 0.05
    beta4 <- 0
  }

  #set within- and between-cluster variability's standard deviations (assume total variance to be 1)
  sd.b <- sqrt(rho)
  sd.w <- sqrt(1-rho)
  
  #initialize the p-value recorder for computing empirical power
  p <- c(NA, nsims)
  
  for (i in 1:nsims){
    #set seed for random generator
    set.seed(2020+1231*i)
    
    #randomly generate a vector of varying cluster sizes, cluster sizes should be positive and total sample size should be even
    if (cv==0){
      m_vector <- rep(m_bar, n)
    } else{
      m_vector <- round(rgamma(n, shape=cv^(-2), rate=m_bar^(-1)*cv^(-2)))
      m_vector[m_vector<2] <- 2
    }
    
    #create individual identifications
    id <- seq(1,sum(m_vector),1)
    #assign cluster id, 1 for the first m_1 individuals, to n for the last m_n individuals
    cluster.id <- rep(1:n, m_vector)
    
    simulated.data <- data.frame(cbind(id, cluster.id))
    
    #equally and randomly assign 0 or 1 to n clusters, representing the cluster-level treatment
    ctrt.assign <- sample(1:n, n/2)
    simulated.data$ctrt <- ifelse(simulated.data$cluster.id %in% ctrt.assign, 1, 0)
    
    #randomize the individual-level treatment WITHIN each cluster
    itrt <- NULL
    for (k in 1:n){
      if (m_vector[k]%%2==0){ #if the current cluster size is even
        itrt_by_cluster <- sample(c(rep(0, floor(m_vector[k]/2)), rep(1, floor(m_vector[k]/2))))
      } else {
        itrt_by_cluster <- sample(c(rep(0, floor(m_vector[k]/2)), rep(1, floor(m_vector[k]/2)), rbinom(1,1,0.5)))
      }
      
      itrt <- c(itrt, itrt_by_cluster)
    }
    
    simulated.data$itrt <- itrt
    
    #randomly generate within-cluster variability for each individual
    simulated.data$epsilon <- rnorm(sum(m_vector), 0, sd.w)
    
    #randomly generate between-cluster variability for each cluster
    gamma <- NULL
    #assign cluster id, 1 for the first m_1 individuals, to n for the last m_n individuals
    for (l in 1:n){
      gamma <- c(gamma, rep(rnorm(1, 0, sd.b), m_vector[l]))
    }
    simulated.data$gamma <- gamma
      
    #create outcome Y_ij's values
    simulated.data$Y <- beta1+beta2*simulated.data$ctrt+beta3*simulated.data$itrt+beta4*simulated.data$itrt*simulated.data$ctrt+simulated.data$gamma+simulated.data$epsilon
      
    fit = try(lme(Y ~ ctrt + itrt + ctrt:itrt, data=simulated.data, random= ~ 1 | cluster.id), silent=T)
    #cope with potentially problematic simulated data set
    if(class(fit)=="try-error"){next}
      
    R <- c(0,0,0,1)
      
    beta <- fit$coef$fixed
    #Use Wald test, which follows the chi-squared distribution with 1 degree of freedom
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    p[i] <- 1-pchisq(test.stat,1)
  }
  empirical <- mean(p<0.05, na.rm = T)
  return(empirical)
}

