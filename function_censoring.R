# Find the values of the parameter under the censoring rates
n_large <- 1000000  # Very large samples

find_censoring <- function(target_rate) {
  set.seed(2024)
  lambda <- 1  # initial value
  found <- FALSE
  
  while (!found) {
    T_large <- rexp(n_large, rate = rho)^(1 / gamma)  # Generate very large Ts
    C_large <- rexp(n_large, rate = lambda)           # Generate very large Cs
    
    censoring_rate <- mean(C_large < T_large)  # Calculate the censoring rate
    
    if (abs(censoring_rate - target_rate) < 0.01) {
      found <- TRUE
    } else {
      lambda <- lambda * (target_rate / censoring_rate)  
    }
  }
  
  return(lambda)
}

lambda_c <- find_censoring(censoring_rate)