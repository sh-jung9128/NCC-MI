# Functions for summary
store_results <- function(model) {
  coefficients <- c(model$coefficients["Z1"],
                    model$coefficients["W1"],
                    model$coefficients["factor(W2)1"])
  se <- summary(model)$coef[, "se(coef)"][c("Z1", "W1", "factor(W2)1")]
  
  return(list(coefficients = coefficients, se = se))
}

calculate_statistics <- function(means, ses, true_beta) {
  mean_beta <- colMeans(means, na.rm = TRUE)
  sd_beta <- apply(means, 2, sd, na.rm = TRUE)
  mean_se <- sqrt(colMeans(ses^2, na.rm = TRUE))
  
  result <- data.frame(
    True_beta = true_beta,
    Mean_beta = mean_beta,
    Bias = mean_beta - true_beta,
    Mean_SE_beta = mean_se,
    SD_beta = sd_beta
  )
  return(result)
}
