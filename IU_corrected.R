##########################################################################################
# Function for calculating empirical power or empirical type I error via 5000 simulations
# 
# for intersection-union test WITH finite-sample consideration
##########################################################################################

IU_corrected <- function(choice, parameter, nsims=5000){
  
  #number of iterations
  nsims <- nsims
  
  if (choice==2){ #null case
    delta_x <- 0
  }
  
  #mean cluster size
  m_bar <- parameter[1]
  #outcome ICC
  rho <- parameter[2]
  #coefficient of variation
  cv <- parameter[3]
  #number of clusters
  n <- parameter[4]

  #set true regression parameters
  beta1 <- 1
  #fix beta2 at 0.15 reltative to other test scenarios
  beta2 <- 0.15
  #two constraints among beta2, 3, and 4, decide beta4 firstly, then decide beta3
  beta4 <- 2*(delta_x-beta2)
  beta3 <- delta_z-0.5*beta4
  
  #set within- and between-cluster variability's standard deviations (assume total variance to be 1)
  sd.b <- sqrt(rho)
  sd.w <- sqrt(1-rho)
  
  #initialize the p-value recorder for computing empirical power
  p.c <- c(NA, nsims)
  p.i <- c(NA, nsims)
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
    
    #randomize the individual-level treatment WITHIN each cluster, using completely randomization
    itrt <- NULL
    for (j in 1:n){
      if (m_vector[j]%%2==0){ 
        #if the current cluster size is even, randomly assign treatment to exactly half of the individuals within the cluster 
        itrt_by_cluster <- sample(c(rep(0, floor(m_vector[j]/2)), rep(1, floor(m_vector[j]/2))))
      } else { 
        #if the current cluster size is odd: minus 1 first to get even cluster size -> do complete randomization -> do a Bernoulli trial for the rest 1 individual
        itrt_by_cluster <- sample(c(rep(0, floor(m_vector[j]/2)), rep(1, floor(m_vector[j]/2)), rbinom(1,1,0.5)))
      }
      itrt <- c(itrt, itrt_by_cluster)
    }
    simulated.data$itrt <- itrt
    
    #randomly generate within-cluster variability for each individual
    simulated.data$epsilon <- rnorm(sum(m_vector), 0, sd.w)
    
    #randomly generate between-cluster variability for each cluster
    gamma <- NULL
    #assign gamma based on cluster id
    for (k in 1:n){
      gamma <- c(gamma, rep(rnorm(1, 0, sd.b), m_vector[k]))
    }
    simulated.data$gamma <- gamma
      
    #create outcome Y_ij's values
    simulated.data$Y <- beta1+beta2*simulated.data$ctrt+beta3*simulated.data$itrt+beta4*simulated.data$itrt*simulated.data$ctrt+simulated.data$gamma+simulated.data$epsilon
      
    #cope with potentially problematic simulated data set with try()
    fit = try(lme(Y ~ ctrt + itrt + ctrt:itrt, data=simulated.data, random= ~ 1 | cluster.id), silent=T)
    #if non-convergent, jump to next iteration 
    if(class(fit)=="try-error"){next}
      
    A.c <- c(0,1,0,0.5)
    A.i <- c(0,0,1,0.5)
    
    #extract parameter estimates
    beta <- fit$coef$fixed
    
    #compute 2 test statistics independently
    test.stat.c <- as.numeric((t(A.c)%*%beta)^2/(t(A.c)%*%vcov(fit)%*%A.c))
    test.stat.i <- as.numeric((t(A.i)%*%beta)^2/(t(A.i)%*%vcov(fit)%*%A.i))
    
    #compute 2 p-values independently, require both significance simultaneously
    p.c[i] <- 1-pf(test.stat.c,1,n-2)
    p.i[i] <- 1-pchisq(test.stat.i,1)
  }
  
  empirical <- mean(p.c<0.05 & p.i<0.05, na.rm = T)
  return(empirical)
}
