############################################################################################################
# Function to generate the simulated data sets for linear mixed modelling (CV follows a Gamma distribution)
############################################################################################################

gendata <- function(beta1, beta2, beta3, beta4, rho, sigma2y=1, n, m_bar, CV){
  
  #randomly generate a vector of varying cluster sizes, cluster sizes should be positive and total sample size should be even
  if (CV==0){
    m_vector <- rep(m_bar, n)
  } else{
    m_vector <- round(rgamma(n, shape=CV^(-2), rate=m_bar^(-1)*CV^(-2)))
    m_vector[m_vector<2] <- 2
  }
  
  #create individual identifications
  id <- seq(1, sum(m_vector), 1)
  #assign cluster id, 1 for the first m_1 individuals, to n for the last m_n individuals
  cluster.id <- rep(1:n, m_vector)
  data <- data.frame(cbind(id, cluster.id))
  
  #equally and randomly assign 0 or 1 to n clusters, representing the cluster-level treatment
  ctrt.assign <- sample(1:n, n/2)
  data$ctrt <- ifelse(data$cluster.id %in% ctrt.assign, 1, 0)
  
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
  data$itrt <- itrt
  
  #randomly generate within-cluster variability for each individual
  data$epsilon <- rnorm(sum(m_vector), 0, sqrt(1-rho))
  
  #randomly generate between-cluster variability for each cluster
  data$gamma <- rep(rnorm(n, 0, sqrt(rho)), m_vector)
  
  #create outcome Y_ij's values
  data$Y <- beta1 + beta2*data$ctrt + beta3*data$itrt + beta4*data$itrt*data$ctrt + data$gamma + data$epsilon
  return(data)
}


###################################################################
# Data generation function when CV follows a Normal distribution
###################################################################


gendata_normal <- function(beta1, beta2, beta3, beta4, rho, sigma2y=1, n, m_bar, CV){
  
  #randomly generate a vector of varying cluster sizes, cluster sizes should be positive and total sample size should be even
  if (CV==0){
    m_vector <- rep(m_bar, n)
  } else{
    m_vector <- round(rnorm(n, mean=m_bar, sd=m_bar*CV))
    m_vector[m_vector<2] <- 2
  }
  
  #create individual identifications
  id <- seq(1, sum(m_vector), 1)
  #assign cluster id, 1 for the first m_1 individuals, to n for the last m_n individuals
  cluster.id <- rep(1:n, m_vector)
  data <- data.frame(cbind(id, cluster.id))
  
  #equally and randomly assign 0 or 1 to n clusters, representing the cluster-level treatment
  ctrt.assign <- sample(1:n, n/2)
  data$ctrt <- ifelse(data$cluster.id %in% ctrt.assign, 1, 0)
  
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
  data$itrt <- itrt
  
  #randomly generate within-cluster variability for each individual
  data$epsilon <- rnorm(sum(m_vector), 0, sqrt(1-rho))
  
  #randomly generate between-cluster variability for each cluster
  data$gamma <- rep(rnorm(n, 0, sqrt(rho)), m_vector)
  
  #create outcome Y_ij's values
  data$Y <- beta1 + beta2*data$ctrt + beta3*data$itrt + beta4*data$itrt*data$ctrt + data$gamma + data$epsilon
  return(data)
}


###################################################################
# Data generation function when CV follows a Uniform distribution
###################################################################


gendata_uniform <- function(beta1, beta2, beta3, beta4, rho, sigma2y=1, n, m_bar, CV){
  
  #randomly generate a vector of varying cluster sizes, cluster sizes should be positive and total sample size should be even
  if (CV==0){
    m_vector <- rep(m_bar, n)
  } else{
    m_vector <- round(runif(n, min=m_bar-sqrt(3)*m_bar*CV, max=m_bar+sqrt(3)*m_bar*CV))
    m_vector[m_vector<2] <- 2
  }
  
  #create individual identifications
  id <- seq(1, sum(m_vector), 1)
  #assign cluster id, 1 for the first m_1 individuals, to n for the last m_n individuals
  cluster.id <- rep(1:n, m_vector)
  data <- data.frame(cbind(id, cluster.id))
  
  #equally and randomly assign 0 or 1 to n clusters, representing the cluster-level treatment
  ctrt.assign <- sample(1:n, n/2)
  data$ctrt <- ifelse(data$cluster.id %in% ctrt.assign, 1, 0)
  
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
  data$itrt <- itrt
  
  #randomly generate within-cluster variability for each individual
  data$epsilon <- rnorm(sum(m_vector), 0, sqrt(1-rho))
  
  #randomly generate between-cluster variability for each cluster
  data$gamma <- rep(rnorm(n, 0, sqrt(rho)), m_vector)
  
  #create outcome Y_ij's values
  data$Y <- beta1 + beta2*data$ctrt + beta3*data$itrt + beta4*data$itrt*data$ctrt + data$gamma + data$epsilon
  return(data)
}