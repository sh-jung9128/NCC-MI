library(survival)
library(MASS)
library(Epi)
library(openxlsx)

# fixed setting
n <- 30000
rho <- 0.5; gamma <- 1.5
beta1 <- log(2); beta2 <- log(2); beta3 <- -1

# changeable setting
iter <- 3000  # number of simulations
n_imputations <- 10  # Multiple imputations
control_size <- 3

# Find the values of the parameter under the censoring rates
censoring_rate <- 0.99 # target rate
lambda_c <- find_censoring(censoring_rate)


# Initialize estimates
means_full <- matrix(NA, nrow = iter, ncol = 3); ses_full <- matrix(NA, nrow = iter, ncol = 3)
means_thomas <- matrix(NA, nrow = iter, ncol = 3); ses_thomas <- matrix(NA, nrow = iter, ncol = 3)
means_ipw <- matrix(NA, nrow = iter, ncol = 3); ses_ipw <- matrix(NA, nrow = iter, ncol = 3)
pooled_means <- matrix(NA, nrow = iter, ncol = 3); pooled_ses <- matrix(NA, nrow = iter, ncol = 3)


# Simulations
set.seed(2024) 
for (i in 1:iter) {
  Z <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, 0.75, 0.75, 1), nrow = 2))
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  W1 <- rnorm(n, mean = 0, sd = 1)
  W2 <- rbinom(n, 1, 0.5)
  u <- runif(n)
  T <- (-log(1-u) / (rho*exp(beta1*Z1+beta2*W1+beta3*W2)))^(1/gamma)
  C <- rexp(n, rate = lambda_c)
  X <- pmin(T, C)
  Delta <- as.numeric(T < C)
  
  data <- data.frame(X = X, Delta = Delta, Z1 = Z1, W1 = W1, W2 = W2)
  ncc_data <- ccwc(exit = X, fail = Delta, data = data, include = list(Z1, Z2, W1, W2), controls = control_size, silent = TRUE) # ncc for thomas
  ncc_data2 <- ncc_data_f(data = data, control_size = control_size) # ncc for ipw
  
  # Fitting
  cox_model_full <- coxph(Surv(X, Delta) ~ Z1 + W1 + factor(W2), data = data) # full cohort
  cox_model_thomas <- clogit(Fail ~ Z1 + W1 + factor(W2) + strata(Set), data=ncc_data)  # NCC-thomas
  cox_model_ipw <- coxph(Surv(X, Delta) ~ Z1 + W1 + factor(W2), data = ncc_data2, weights = weight)  # NCC-ipw
  
  # Estimate coefficients and standard error
  means_full[i, ] <- store_results(cox_model_full)$coefficients; ses_full[i, ] <- store_results(cox_model_full)$se
  means_thomas[i, ] <- store_results(cox_model_thomas)$coefficients; ses_thomas[i, ] <- store_results(cox_model_thomas)$se
  means_ipw[i, ] <- store_results(cox_model_ipw)$coefficients; ses_ipw[i, ] <- store_results(cox_model_ipw)$se
  
  # Imputation Step
  imputed_datasets <- impute_data(data, baseline_hazard, n_imputations, ncc_data2)
  analysis_results <- perform_analysis(imputed_datasets)
  pooled_means[i, ] <- analysis_results$mean_beta; pooled_ses[i, ] <- analysis_results$se_beta
}


true_beta = c(beta1, beta2, beta3)

result_full <- calculate_statistics(means_full, ses_full, true_beta)
result_thomas <- calculate_statistics(means_thomas, ses_thomas, true_beta)
result_ipw <- calculate_statistics(means_ipw, ses_ipw, true_beta)
result_mi <- calculate_statistics(pooled_means, pooled_ses, true_beta)

# Export
#write.xlsx(list(Full = result_full, Thomas = result_thomas, IPW = result_ipw, MI = result_mi), file = "Results.xlsx", overwrite = TRUE)
