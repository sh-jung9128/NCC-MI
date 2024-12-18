library(survival)
library(MASS)

# Function for multiple imputation
# Generate baseline hazard (Nelson-Aalen estimate approximation)
baseline_hazard <- function(data) {
  order.time <- order(data$X)  # Sort the entire dataset by survival time
  nafit <- coxph(Surv(X, Delta) ~ 1, data = data) # A Cox model without covariates
  cum.haz <- basehaz(nafit)$hazard  # Calculate the cumulative hazard function (Nelson-Aalen estimator)
  
  # Assign cumulative hazard values to H0
  n.s <- dim(data)[1]
  H0<-rep(0, n.s)
  for (k in 1:n.s) {
    H0[order.time[k]] <- c(cum.haz, rep(cum.haz[length(cum.haz)], n.s - length(cum.haz)))[k]
  }
  return(H0)
}

# Multiple imputation function with posterior sampling
impute_data <- function(data, baseline_hazard, n_imputations, ncc_data) {
  imputed_datasets <- vector("list", n_imputations)
  
  for (m in 1:n_imputations) {
    # Calculate Baseline Hazard
    H0 <- baseline_hazard(data)
    
    # Add interaction terms to data
    data$W1_H0 <- data$W1 * H0
    data$W2_H0 <- data$W2 * H0
    
    # Fit the imputation model
    fit <- lm(Z1 ~ Z2 + W1 + factor(W2) + factor(Delta) + H0 + W1_H0 + W2_H0, data = data)
    
    # Extract posterior draws of parameters
    posterior_draws <- MASS::mvrnorm(
      n = 1,  # Single draw for each imputation
      mu = coef(fit),  # Mean from estimated coefficients
      Sigma = vcov(fit)  # Variance-covariance matrix of estimates
    )
    theta0 <- posterior_draws["(Intercept)"]
    thetaZ2 <- posterior_draws["Z2"]
    thetaW1 <- posterior_draws["W1"]
    thetaW2 <- posterior_draws["factor(W2)1"]
    thetaDelta <- posterior_draws["factor(Delta)1"]
    thetaT <- posterior_draws["H0"]
    thetaW1_H0 <- posterior_draws["W1_H0"]
    thetaW2_H0 <- posterior_draws["W2_H0"]
    
    # Impute missing values for Z1
    data$Z1_imp <- with(data, ifelse(
      rownames(data) %in% rownames(ncc_data), Z1,  # If Z1 is observed, retain the existing Z1
      theta0 + thetaZ2 * Z2 + thetaW1 * W1 + thetaW2 * W2 +
        thetaDelta * Delta + thetaT * H0 + thetaW1_H0 * W1_H0 + thetaW2_H0 * W2_H0 +
        rnorm(nrow(data), mean = 0, sd = summary(fit)$sigma)  # Generate replacement values.
    ))
    imputed_datasets[[m]] <- data
  }
  
  return(imputed_datasets)
}


# Analysis function for imputed datasets
perform_analysis <- function(imputed_datasets) {
  beta_estimates <- matrix(NA, nrow = n_imputations, ncol = 3)
  variances <- matrix(NA, nrow = n_imputations, ncol = 3)  
  
  for (m in seq_along(imputed_datasets)) {
    data <- imputed_datasets[[m]]
    
    # Cox regression on nested case-control sample
    cox_model <- coxph(Surv(X, Delta) ~ Z1_imp + W1 + factor(W2), data = data)
    beta_estimates[m, ] <- c(cox_model$coefficients["Z1_imp"],
                             cox_model$coefficients["W1"],
                             cox_model$coefficients["factor(W2)1"])
    variances[m, ] <- diag(vcov(cox_model))  # Variance of the coefficients
  }
  # Rubin's Rule for pooling
  pooled_means <- colMeans(beta_estimates, na.rm = TRUE)
  within_variance <- colMeans(variances, na.rm = TRUE)
  between_variance <- apply(beta_estimates, 2, var, na.rm = TRUE)  # Variance of estimates
  
  total_variance <- within_variance + (1 + 1 / n_imputations) * between_variance
  pooled_se <- sqrt(total_variance)  # Standard error using total variance
  
  return(list(mean_beta = pooled_means, se_beta = pooled_se))
}
